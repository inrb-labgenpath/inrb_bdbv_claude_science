#!/usr/bin/env python3
"""build_alignments.py - dated, masked alignments for Delphy / BEAST X.

Recovered from the code that wrote cleaned_20260912.fasta, all_dated_20260912.fasta
and minus6aug_20260912.fasta (12 Sep 2026 round) and parameterised.

Sets written (suffix from --out-suffix):
  cleaned_<suffix>.fasta      dated, completeness >= --min-completeness, screen class
                              not in --exclude-classes (default E,E-,D,M,H)
  all_dated_<suffix>.fasta    dated, completeness >= --min-completeness (flagged genomes kept)
  <name>_<suffix>.fasta       cleaned set without the release batch(es) given by
                              --drop-release (default name 'minus6aug' for 2026-08-06)
  <name>_<suffix>.fasta       optional: set built under an external exclusion list
                              (--exclusion-list FILE --exclusion-set-name NAME): dated,
                              completeness filter, genomes on the list removed; the screen
                              classes are NOT applied unless --exclusion-list-adds-to-classes
  alignment_sets_<suffix>.csv per-genome membership, exclusion reasons, masked positions
  taxa_dates_<suffix>.csv     taxon label, date, decimal date (BEAST X tip dates)

Rules applied to every sequence written
  * ADAR-type clusters (>= --adar-min derived T>C, or A>G, changes, single linkage within
    --adar-window nt): all positions of a cluster except the first are set to N
  * positions after --max-pos (1-based 18,901..end) are set to N
  * header: >accessionVersion|YYYY-MM-DD ; sequence on one line
  * 'dated' = collection date given to the day (YYYY-MM-DD); month-precision dates and
    ranges are not dated
  * completeness = the screen table's called_frac (fraction A/C/G/T over positions
    1..max-pos, rounded to 3 decimals); the threshold acts on the rounded value
  * version rule: only genomes present in the screen table are eligible (latest,
    non-revoked versions)
  * duplicate rule (--duplicate-rule): 'none' (12 Sep behaviour: both copies of a sample
    released under two accessions are kept), 'keep-best' or 'drop-discordant'

Set-specification mode (24 Sep 2026 round): --set-spec FILE.json
  The JSON file holds {"primary": NAME, "sets": {NAME: {filters}}}.  Only the sets of the file are
  written.  Filters (all optional; defaults in brackets):
    dated [true]                 collection date given to the day
    min_completeness [0.9]       null = no completeness filter; acts on the EXACT called fraction
                                 with --completeness-precision exact
    exclude_classes [E,E-,D,M,H] list of screen classes removed; [] = flagged genomes kept
    duplicate_rule [none]        none | keep-best | drop-discordant
                                 keep-best: one copy per pair (unflagged class, then completeness,
                                 then earlier release, then accession)
                                 drop-discordant: copies identical at all sites called in both ->
                                 the more complete copy is kept (ties: earlier release, accession);
                                 copies that differ -> both removed
    data_use [exclude]           exclude | include | only-add : genomes of --data-use-exclude
    drop_batches []              list of batch labels, or "anomalous" (column 'anomalous' of --batch-table)
    drop_hold_loo [false]        remove genomes with hold_loo in the screen table
    country [null]               keep genomes of this geoLocCountry only
    groups [null]                keep genomes of these submitting groups only (groupName)
    collected_from / collected_to / released_to [null]   inclusive dates (release = --release-column)
    exclusion_list [false]       remove the accessions of --exclusion-list
    class_column ["class"]       column of the screen table that holds the class used by this set
                                 (e.g. the class of an earlier rule set attached with screen.py --attach-class)
    exclude_accessions_file [null]   file with accession(Version)s removed from this set, one per line
  Key "external_data_set" (with "neutral_wording"): {"country", "main_group_only", "released_to", "collected_to",
  "include_undated"} describes the data set of the analysis that published --exclusion-list; the set table then
  holds the column state_in_external_analysis (excluded / retained / not in that data set) instead of a yes/no column.
  Keys beside "sets": "primary" (name of the primary set) and "neutral_wording" [false]: records that
  share a sample identifier are named as such in the tables (no statement that one sample was
  sequenced twice), and the reasons for absence are written out.
  Additional outputs: alignment_set_summary_<suffix>.csv, alignment_sets_by_zone_<suffix>.csv,
  alignment_sets_by_month_<suffix>.csv, alignment_sets_last_weeks_<suffix>.csv,
  alignment_duplicate_pairs_<suffix>.csv
"""
import argparse
import gzip
import hashlib
import json
import os
import re
import sys

import numpy as np
import pandas as pd

BASES = [b"A", b"C", b"G", b"T"]


def read_fasta(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    seqs, name, buf = {}, None, []
    with opener(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\n").rstrip("\r")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(buf).upper()
                name = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if name is not None:
        seqs[name] = "".join(buf).upper()
    return seqs


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def decimal_year(ts):
    """Decimal date used for BEAST tip dates and the Skygrid geometry (12 Sep definition):
    year + (days since 1 Jan) / (days in that year)."""
    ts = pd.Timestamp(ts)
    start = pd.Timestamp(ts.year, 1, 1)
    end = pd.Timestamp(ts.year + 1, 1, 1)
    return ts.year + (ts - start).days / (end - start).days


def adar_clusters(positions_tc, positions_ag, window, minn):
    """Single-linkage clusters; returns list of sorted 0-based position arrays."""
    out = []
    for cand in (positions_tc, positions_ag):
        if len(cand) < minn:
            continue
        cand = np.sort(cand)
        groups = np.split(cand, np.where(np.diff(cand) > window)[0] + 1)
        out += [g for g in groups if len(g) >= minn]
    return out


def read_exclusion_list(path):
    """Accessions (with or without version suffix), one per line or a CSV with a column
    named accession / accessionVersion.  Returns set of accessions WITHOUT version."""
    txt = open(path).read().strip().splitlines()
    if not txt:
        return set()
    if "," in txt[0] or "\t" in txt[0]:
        df = pd.read_csv(path, sep=None, engine="python")
        col = next((c for c in df.columns if c.lower() in ("accession", "accessionversion")), None)
        if col is None:
            raise SystemExit(f"{path}: no column named accession/accessionVersion")
        vals = df[col].dropna().astype(str)
    else:
        vals = pd.Series([t.strip() for t in txt if t.strip() and not t.startswith("#")])
    return set(re.sub(r"\.\d+$", "", v.strip()) for v in vals)


def site_counts(rows):
    """Number of variable and parsimony-informative columns among A/C/G/T of the rows (2-D S1 array)."""
    cnt = np.stack([(rows == b).sum(0) for b in BASES])
    n_states = (cnt > 0).sum(0)
    n_states2 = (cnt >= 2).sum(0)
    return int((n_states >= 2).sum()), int((n_states2 >= 2).sum())


def build_from_spec(args, ctx):
    """Set-specification mode (24 Sep 2026 round)."""
    ids, arr, sub, valid = ctx["ids"], ctx["arr"], ctx["sub"], ctx["valid"]
    meta, screen, masked = ctx["meta"], ctx["screen"], ctx["masked"]
    dated, parsed, date_str = ctx["dated"], ctx["parsed"], ctx["date_str"]
    cf_round, n_called, cls, rel = ctx["cf"], ctx["n_called"], ctx["cls"], ctx["rel"]
    n_clusters = ctx["n_clusters"]
    MAXPOS, sfx = args.max_pos, args.out_suffix
    spec = json.load(open(args.set_spec))
    primary = spec["primary"]
    NEUTRAL = bool(spec.get("neutral_wording", False))
    LAB = {"undated": "undated", "complete": "below completeness", "class": "class ",
           "disc": "duplicate discordant", "low": "duplicate lower-quality copy", "du": "data-use D7",
           "other": "other filter of the primary set"}
    if NEUTRAL:
        LAB = {"undated": "collection date not given to the day", "complete": "called fraction below threshold",
               "class": "class ", "disc": "shares a sample identifier with a discordant record",
               "low": "shares a sample identifier with a concordant, more complete record",
               "du": "restricted use, agreement of the submitters pending",
               "other": "other filter of the primary set"}

    def class_of(f):
        col = f.get("class_column", "class")
        assert col in screen.columns, f"screen table has no column {col}"
        return screen.loc[ids, col].astype(str).values
    assert primary in spec["sets"], "primary set missing from the specification"
    idx = {a: k for k, a in enumerate(ids)}
    acc_nov = np.array([re.sub(r"\.\d+$", "", a) for a in ids])
    EXACT = args.completeness_precision == "exact"

    def complete_at(thr):
        if thr is None:
            return np.ones(len(ids), bool)
        if EXACT:
            return n_called >= int(np.ceil(round(thr * MAXPOS, 6)))
        return cf_round >= thr

    # ---- duplicate pairs ---------------------------------------------------------------------
    if args.duplicate_pairs:
        dp = pd.read_csv(args.duplicate_pairs)
        pairs = [(a, b) for a, b in zip(dp["accessionVersion_1"], dp["accessionVersion_2"])]
        missing = [a for p in pairs for a in p if a not in idx]
        assert not missing, f"duplicate pairs name accessions absent from the screen table: {missing[:5]}"
    else:
        sid = screen.loc[ids, "sid"].fillna("").astype(str).values
        grp = pd.Series(np.arange(len(ids))).groupby(sid).apply(list)
        pairs = []
        for s, mem in grp.items():
            if s == "" or len(mem) < 2:
                continue
            assert len(mem) == 2, "sample id with more than two accessions: extend the duplicate rule"
            pairs.append((ids[mem[0]], ids[mem[1]]))
    in_pair = pd.Series(pd.Series([a for p in pairs for a in p]).value_counts())
    assert (in_pair == 1).all(), "an accession occurs in more than one duplicate pair"

    def rule_masks(rule, bad, cls=cls):
        """Returns (dropped, reason) for one duplicate rule."""
        drop = np.zeros(len(ids), bool)
        why = np.array([""] * len(ids), dtype=object)
        if rule == "none":
            return drop, why
        for a, b in pairs:
            i, j = idx[a], idx[b]
            nd = int(((sub[i] != sub[j]) & valid[i] & valid[j]).sum())
            if rule == "drop-discordant":
                if nd > 0:
                    drop[[i, j]] = True
                    why[i] = why[j] = LAB["disc"]
                    continue
                order = sorted([i, j], key=lambda k: (-int(n_called[k]), rel[k], ids[k]))
            else:                                   # keep-best
                order = sorted([i, j], key=lambda k: (cls[k] in bad, -int(n_called[k]), rel[k], ids[k]))
            drop[order[1]] = True
            why[order[1]] = LAB["low"]
        return drop, why

    # ---- other per-genome attributes ---------------------------------------------------------
    du = set()
    if args.data_use_exclude:
        txt = args.data_use_exclude
        vals = open(txt).read().replace("\n", ",").split(",") if os.path.exists(txt) else txt.split(",")
        du = set(re.sub(r"\.\d+$", "", v.strip()) for v in vals if v.strip())
    is_du = np.isin(acc_nov, list(du)) if du else np.zeros(len(ids), bool)
    assert int(is_du.sum()) == len(du), "data-use accessions absent from the snapshot"
    batch = screen.loc[ids, "batch"].astype(str).values
    anomalous = []
    if args.batch_table:
        bt = pd.read_csv(args.batch_table, index_col=0)
        anomalous = sorted(bt.index[bt["anomalous"].astype(bool)])
    hold_loo = screen.loc[ids, "hold_loo"].fillna(False).astype(bool).values if "hold_loo" in screen.columns \
        else np.zeros(len(ids), bool)
    country = meta.loc[ids, "geoLocCountry"].astype("string").fillna("").values
    group = meta.loc[ids, "groupName"].astype("string").fillna("").values
    reldate = pd.to_datetime(pd.Series(rel), format="ISO8601", errors="coerce").values
    on_list = np.zeros(len(ids), bool)
    absent = []
    if args.exclusion_list:
        excl = read_exclusion_list(args.exclusion_list)
        on_list = np.isin(acc_nov, list(excl))
        absent = sorted(excl - set(acc_nov))
        print(f"[build] exclusion list: {len(excl)} accessions; {int(on_list.sum())} present in snapshot; "
              f"{len(absent)} not in snapshot")
    pdate = parsed.values

    def accession_file(path):
        vals = [t.strip() for t in open(path).read().replace(",", "\n").splitlines() if t.strip() and not t.startswith("#")]
        nov = set(re.sub(r"\.\d+$", "", v) for v in vals)
        hit = np.isin(acc_nov, list(nov))
        assert int(hit.sum()) == len(nov), f"{path}: {len(nov) - int(hit.sum())} accessions absent from the screen table"
        return hit

    sets, used = {}, {}
    for name, f in spec["sets"].items():
        bad = set(f.get("exclude_classes", ["E", "E-", "D", "M", "H"]))
        cls_f = class_of(f)
        keep = np.ones(len(ids), bool)
        if f.get("dated", True):
            keep &= dated
        keep &= complete_at(f.get("min_completeness", 0.9))
        keep &= ~np.isin(cls_f, list(bad))
        drop, _ = rule_masks(f.get("duplicate_rule", "none"), bad, cls_f)
        keep &= ~drop
        if f.get("exclude_accessions_file"):
            keep &= ~accession_file(f["exclude_accessions_file"])
        du_mode = f.get("data_use", "exclude")
        assert du_mode in ("exclude", "include")
        if du_mode == "exclude":
            keep &= ~is_du
        db = f.get("drop_batches", [])
        db = anomalous if db == "anomalous" else list(db)
        if f.get("drop_batches") == "anomalous":
            assert args.batch_table, "drop_batches = anomalous needs --batch-table"
        keep &= ~np.isin(batch, db)
        if f.get("drop_hold_loo", False):
            assert "hold_loo" in screen.columns, "screen table without hold_loo"
            keep &= ~hold_loo
        if f.get("country"):
            keep &= (country == f["country"])
        if f.get("groups"):
            keep &= np.isin(group, f["groups"])
        if f.get("collected_from"):
            keep &= dated & (pdate >= np.datetime64(f["collected_from"]))
        if f.get("collected_to"):
            keep &= dated & (pdate <= np.datetime64(f["collected_to"]))
        if f.get("released_to"):
            keep &= (reldate <= np.datetime64(f["released_to"]))
        if f.get("exclusion_list", False):
            assert args.exclusion_list, "exclusion_list = true needs --exclusion-list"
            keep &= ~on_list
        sets[name] = keep
        used[name] = dict(f, drop_batches_resolved=db)

    # ---- reason for absence from the primary set (first applicable) -----------------------------
    pf = spec["sets"][primary]
    pbad = set(pf.get("exclude_classes", ["E", "E-", "D", "M", "H"]))
    cls = class_of(pf)                                 # class of the primary set in the tables below
    pdrop, pwhy = rule_masks(pf.get("duplicate_rule", "none"), pbad, cls)
    pcomplete = complete_at(pf.get("min_completeness", 0.9))
    flagged = np.isin(cls, list(pbad))
    reason = np.array([""] * len(ids), dtype=object)
    for k in range(len(ids)):
        if sets[primary][k]:
            continue
        if not dated[k]:
            reason[k] = LAB["undated"]
        elif not pcomplete[k]:
            reason[k] = LAB["complete"]
        elif flagged[k]:
            reason[k] = LAB["class"] + cls[k]
        elif pdrop[k]:
            reason[k] = pwhy[k]
        elif is_du[k] and pf.get("data_use", "exclude") == "exclude":
            reason[k] = LAB["du"]
        else:
            reason[k] = LAB["other"]

    # ---- write ------------------------------------------------------------------------------------
    def masked_row(k):
        row = arr[k].copy()
        if not args.no_adar_mask and k in masked:
            row[np.array(masked[k])] = b"N"
        row[MAXPOS:] = b"N"
        return row

    zone = None
    if args.validation:
        V = pd.read_csv(args.validation, low_memory=False).set_index("accessionVersion")
        zone = V["zone"].reindex(ids).astype("string").fillna("(not recorded)").values
    else:
        zone = meta.loc[ids, "geoLocAdmin2"].astype("string").fillna("(not recorded)").values
    month = np.where(dated, pd.Series(date_str.values).str[:7].values, "(undated)")

    summary, rows_zone, rows_month, rows_tail = [], [], [], []
    for name, keep in sets.items():
        fn = os.path.join(args.outdir, f"{name}_{sfx}.fasta")
        kk = np.where(keep)[0]
        undated_in = [ids[k] for k in kk if not dated[k]]
        assert not undated_in, f"set {name} holds undated genomes; FASTA headers need a day-precision date"
        M = np.stack([masked_row(k) for k in kk])
        with open(fn, "w") as fh:
            for r, k in zip(M, kk):
                fh.write(f">{ids[k]}|{date_str.iloc[k]}\n{r.tobytes().decode()}\n")
        nvar, npi = site_counts(M[:, :MAXPOS])
        d = parsed[keep]
        last = d.max()
        summary.append(dict(set=name, file=os.path.basename(fn), n=int(keep.sum()), first_tip=str(d.min().date()),
                            last_tip=str(last.date()), last_tip_decimal=round(decimal_year(last), 6),
                            n_tips_on_last_date=int((d == last).sum()),
                            second_last_tip=str(d[d < last].max().date()),
                            n_variable_sites=nvar, n_parsimony_informative_sites=npi,
                            n_with_adar_mask=int(sum(1 for k in masked if keep[k])),
                            n_masked_adar_positions=int(sum(len(masked[k]) for k in masked if keep[k])),
                            sha256=sha256(fn)))
        for z, c in pd.Series(zone[keep]).value_counts().items():
            rows_zone.append(dict(set=name, health_zone=z, n=int(c)))
        for m_, c in pd.Series(month[keep]).value_counts().sort_index().items():
            rows_month.append(dict(set=name, collection_month=m_, n=int(c)))
        cut = last - pd.Timedelta(days=args.tail_days - 1)
        tail = d[d >= cut]
        for day, c in tail.dt.strftime("%Y-%m-%d").value_counts().sort_index().items():
            rows_tail.append(dict(set=name, collection_date=day, n=int(c),
                                  days_before_last_tip=int((last - pd.Timestamp(day)).days)))
        print(f"[build] {name}: n={int(keep.sum())} first tip {d.min().date()} last tip {last.date()} "
              f"variable sites {nvar} -> {fn}")
    S = pd.DataFrame(summary).set_index("set")
    S.to_csv(os.path.join(args.outdir, f"alignment_set_summary_{sfx}.csv"))
    Z = pd.DataFrame(rows_zone).pivot(index="health_zone", columns="set", values="n").fillna(0).astype(int)
    Z = Z[list(sets)].sort_values(primary, ascending=False)
    Z.to_csv(os.path.join(args.outdir, f"alignment_sets_by_zone_{sfx}.csv"))
    Mo = pd.DataFrame(rows_month).pivot(index="collection_month", columns="set", values="n").fillna(0).astype(int)
    Mo[list(sets)].to_csv(os.path.join(args.outdir, f"alignment_sets_by_month_{sfx}.csv"))
    T = pd.DataFrame(rows_tail).pivot(index="collection_date", columns="set", values="n").fillna(0).astype(int)
    T[[s for s in sets if s in T.columns]].to_csv(os.path.join(args.outdir, f"alignment_sets_last_weeks_{sfx}.csv"))

    tab = pd.DataFrame(dict(
        accessionVersion=ids, collection_date=date_str.values, dated=dated,
        called_fraction_exact=n_called / MAXPOS, screen_class=cls, release=rel,
        in_duplicate_pair=np.isin(ids, [a for p in pairs for a in p]),
        data_use_D7=is_du, on_exclusion_list=on_list, hold_loo=hold_loo,
        in_anomalous_batch=np.isin(batch, anomalous), adar_clusters=n_clusters,
        masked_positions_1based=[",".join(str(q + 1) for q in masked.get(k, [])) for k in range(len(ids))]))
    if NEUTRAL:
        tab = tab.rename(columns={"in_duplicate_pair": "shares_sample_identifier",
                                  "data_use_D7": "restricted_use_agreement_pending",
                                  "on_exclusion_list": "on_external_exclusion_table",
                                  "in_anomalous_batch": "in_release_with_excess_of_genomes_without_derived_alleles"})
        tab = tab.drop(columns=["hold_loo"])
        eds = spec.get("external_data_set")
        if eds:
            # three states with respect to the analysis that published the exclusion table:
            # excluded (on the table) / retained (in its data set, not on the table) / not in that data set
            in_ds = np.ones(len(ids), bool)
            if eds.get("country"):
                in_ds &= country == eds["country"]
            if eds.get("main_group_only", False):
                in_ds &= group == pd.Series(group).value_counts().index[0]
            if eds.get("released_to"):
                in_ds &= reldate <= np.datetime64(eds["released_to"])
            if eds.get("collected_to"):
                late = dated & (pdate > np.datetime64(eds["collected_to"]))
                in_ds &= ~late
            if not eds.get("include_undated", True):
                in_ds &= dated
            assert not (on_list & ~in_ds).any(), "accessions of the exclusion table outside the external data set"
            state = np.where(on_list, "excluded", np.where(in_ds, "retained", "not in that data set"))
            pos = list(tab.columns).index("on_external_exclusion_table")
            tab = tab.drop(columns=["on_external_exclusion_table"])
            tab.insert(pos, "state_in_external_analysis", state)
    for name, keep in sets.items():
        tab[f"in_{name}"] = keep
    tab[f"reason_not_in_{primary}"] = reason
    tab.to_csv(os.path.join(args.outdir, f"alignment_sets_{sfx}.csv"), index=False)

    prow = []
    for a, b in pairs:
        i, j = idx[a], idx[b]
        nd = int(((sub[i] != sub[j]) & valid[i] & valid[j]).sum())
        prow.append(dict(accessionVersion_1=a, accessionVersion_2=b, collection_date_1=date_str.iloc[i],
                         collection_date_2=date_str.iloc[j], release_1=rel[i], release_2=rel[j],
                         called_fraction_1=n_called[i] / MAXPOS, called_fraction_2=n_called[j] / MAXPOS,
                         class_1=cls[i], class_2=cls[j], n_sites_both_called=int((valid[i] & valid[j]).sum()),
                         n_differences=nd, concordant=nd == 0,
                         **{f"in_{primary}_1": bool(sets[primary][i]), f"in_{primary}_2": bool(sets[primary][j]),
                            "reason_1": reason[i], "reason_2": reason[j]}))
    if NEUTRAL:
        pd.DataFrame(prow).rename(columns={"concordant": "records_concordant"}).to_csv(
            os.path.join(args.outdir, f"records_sharing_sample_identifier_{sfx}.csv"), index=False)
    else:
        pd.DataFrame(prow).to_csv(os.path.join(args.outdir, f"alignment_duplicate_pairs_{sfx}.csv"), index=False)

    anyset = np.zeros(len(ids), bool)
    for keep in sets.values():
        anyset |= keep
    td = pd.DataFrame(dict(taxon_acc=np.array(ids)[anyset], date=date_str.values[anyset]))
    td["decimal"] = [f"{decimal_year(x):.6f}" for x in td.date]
    td["taxon"] = td.taxon_acc + "|" + td.date
    td.to_csv(os.path.join(args.outdir, f"taxa_dates_{sfx}.csv"), index=False)
    with open(os.path.join(args.outdir, f"alignment_parameters_{sfx}.json"), "w") as fh:
        json.dump(dict(parameters=vars(args), primary=primary, sets=used,
                       summary=S.reset_index().to_dict(orient="records"), n_screened=len(ids),
                       n_duplicate_pairs=len(pairs), anomalous_batches=anomalous,
                       data_use_excluded=sorted(du), exclusion_list_absent_from_snapshot=absent,
                       inputs={k: dict(file=os.path.basename(getattr(args, k)), sha256=sha256(getattr(args, k)))
                               for k in ("metadata", "alignment", "screen")}), fh, indent=1, default=str)
    return 0


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--screen", required=True, help="per-genome table written by screen.py")
    ap.add_argument("--out-suffix", required=True)
    ap.add_argument("--outdir", default=".")
    ap.add_argument("--exclude-classes", default="E,E-,D,M,H")
    ap.add_argument("--min-completeness", type=float, default=0.9)
    ap.add_argument("--polarise-completeness", type=float, default=0.95,
                    help="majority base computed among genomes with completeness >= this "
                         "(0 = all genomes, as in the 12 Sep alignment code; identical result on 12 Sep data)")
    ap.add_argument("--genome-length", type=int, default=18940)
    ap.add_argument("--max-pos", type=int, default=18900, help="positions 1..max-pos are kept")
    ap.add_argument("--adar-window", type=int, default=300)
    ap.add_argument("--adar-min", type=int, default=3)
    ap.add_argument("--no-adar-mask", action="store_true")
    ap.add_argument("--release-column", default="earliestReleaseDate",
                    help="12 Sep code used releasedDate (same sets on 12 Sep data)")
    ap.add_argument("--drop-release", default="2026-08-06",
                    help="comma-separated release dates removed from the cleaned set for the "
                         "leave-release-out set; empty string = do not write it")
    ap.add_argument("--drop-release-name", default="minus6aug")
    ap.add_argument("--exclusion-list", default=None)
    ap.add_argument("--exclusion-set-name", default="external_exclusions")
    ap.add_argument("--exclusion-list-adds-to-classes", action="store_true",
                    help="apply the list on top of the class exclusions instead of replacing them")
    ap.add_argument("--duplicate-rule", choices=["none", "keep-best", "drop-discordant"], default="none")
    ap.add_argument("--expected-n-cleaned", type=int, default=None)
    n = ap.add_argument_group("24 Sep 2026 round (defaults = 12 Sep behaviour)")
    n.add_argument("--set-spec", default=None, help="JSON file with the sets to write (see module docstring)")
    n.add_argument("--completeness-precision", choices=["rounded", "exact"], default="rounded",
                   help="exact = the completeness threshold acts on the unrounded called fraction")
    n.add_argument("--duplicate-pairs", default=None,
                   help="table of same-sample pairs (accessionVersion_1, accessionVersion_2); "
                        "default: groups of equal 'sid' in the screen table")
    n.add_argument("--data-use-exclude", default=None,
                   help="comma-separated accession(Version)s, or a file, of genomes left out for data-use reasons")
    n.add_argument("--batch-table", default=None, help="screen_batches_<suffix>.csv (column 'anomalous')")
    n.add_argument("--validation", default=None, help="per-genome table with the harmonised 'zone'")
    n.add_argument("--tail-days", type=int, default=42, help="window of the tip-date table before the last tip")
    args = ap.parse_args(argv)
    os.makedirs(args.outdir, exist_ok=True)
    sfx = args.out_suffix
    MAXPOS = args.max_pos
    BAD = set(c.strip() for c in args.exclude_classes.split(",") if c.strip())

    meta = pd.read_csv(args.metadata, low_memory=False).set_index("accessionVersion")
    screen = pd.read_csv(args.screen, index_col=0)
    seqs = read_fasta(args.alignment)
    ids = [a for a in seqs if a in screen.index]
    not_screened = [a for a in seqs if a not in screen.index]
    assert len(ids) > 0, "no alignment record is present in the screen table"
    if not_screened:
        print(f"[build] WARNING {len(not_screened)} alignment records absent from the screen table are ignored")
    L = len(seqs[ids[0]])
    assert L == args.genome_length and all(len(seqs[a]) == L for a in ids)
    arr = np.frombuffer("".join(seqs[a] for a in ids).encode(), dtype="S1").reshape(len(ids), L).copy()
    sub = arr[:, :MAXPOS]
    valid = np.isin(sub, BASES)

    date_str = meta.loc[ids, "sampleCollectionDate"].astype("string")
    dated = date_str.str.fullmatch(r"\d{4}-\d{2}-\d{2}").fillna(False).astype(bool).values
    parsed = pd.to_datetime(date_str, format="ISO8601", errors="coerce")
    assert parsed[dated].notna().all()
    cf = screen.loc[ids, "called_frac"].astype(float).values          # rounded to 3 decimals by screen.py
    cf_raw = valid.mean(1)
    assert np.allclose(np.round(cf_raw, 3), cf, atol=1e-9), "screen called_frac does not match alignment"
    cls = screen.loc[ids, "class"].astype(str).values
    n_called = valid.sum(1)
    if args.completeness_precision == "exact":
        complete = n_called >= int(np.ceil(round(args.min_completeness * MAXPOS, 6)))
        pol_ok = n_called >= int(np.ceil(round(args.polarise_completeness * MAXPOS, 6)))
    else:
        complete = cf >= args.min_completeness
        pol_ok = cf >= args.polarise_completeness

    # ---- polarisation and ADAR masks --------------------------------------------
    pol = pol_ok
    counts = np.stack([((sub[pol] == b)).sum(0) for b in BASES])
    maj = np.array(BASES)[counts.argmax(0)]
    nocov = counts.sum(0) == 0
    masked = {}
    n_clusters = np.zeros(len(ids), int)
    for k in range(len(ids)):
        d = np.where(valid[k] & (sub[k] != maj) & ~nocov)[0]
        if len(d) == 0:
            continue
        tc = d[(maj[d] == b"T") & (sub[k, d] == b"C")]
        ag = d[(maj[d] == b"A") & (sub[k, d] == b"G")]
        cl = adar_clusters(tc, ag, args.adar_window, args.adar_min)
        n_clusters[k] = len(cl)
        m = [int(q) for g in cl for q in g[1:]]
        if m:
            masked[k] = m
    if "adar_clusters" in screen.columns:
        n_scr = screen.loc[ids, "adar_clusters"].fillna(0).astype(int).values
        n_dis = int((n_scr != n_clusters).sum())
        print(f"[build] ADAR cluster counts vs screen table: {n_dis} disagreements")
        assert n_dis == 0, "ADAR clusters disagree with the screen table"

    if args.set_spec:
        return build_from_spec(args, dict(ids=ids, arr=arr, sub=sub, valid=valid, meta=meta, screen=screen,
                                          masked=masked, dated=dated, parsed=parsed, date_str=date_str, cf=cf,
                                          n_called=n_called, cls=cls, n_clusters=n_clusters,
                                          rel=meta.loc[ids, args.release_column].astype(str).values))

    # ---- duplicate rule ------------------------------------------------------------
    sid = screen.loc[ids, "sid"].fillna("").astype(str).values if "sid" in screen.columns else np.array([""] * len(ids))
    rel = meta.loc[ids, args.release_column].astype(str).values
    dup_drop = np.zeros(len(ids), bool)
    dup_note = [""] * len(ids)
    groups = pd.Series(np.arange(len(ids))).groupby(sid).apply(list)
    n_dup_groups = 0
    for s, members in groups.items():
        if s == "" or len(members) < 2:
            continue
        n_dup_groups += 1
        for k in members:
            dup_note[k] = "same-sample-id group of %d" % len(members)
        if args.duplicate_rule == "none":
            continue
        # rank: unflagged class first, then completeness (high), release (early), accession
        ranked = sorted(members, key=lambda k: (cls[k] in BAD, -cf[k], rel[k], ids[k]))
        best = ranked[0]
        discordant = any(int(((sub[best] != sub[k]) & valid[best] & valid[k]).sum()) > 0 for k in ranked[1:])
        for k in ranked[1:]:
            dup_drop[k] = True
            dup_note[k] += "; dropped (duplicate of %s)" % ids[best]
        if args.duplicate_rule == "drop-discordant" and discordant:
            dup_drop[best] = True
            dup_note[best] += "; dropped (copies of the sample disagree)"

    flagged = np.isin(cls, list(BAD))
    base = dated & complete & ~dup_drop
    sets = {"all_dated": base.copy(), "cleaned": base & ~flagged}
    drop_rel = [r.strip() for r in args.drop_release.split(",") if r.strip()]
    if drop_rel:
        sets[args.drop_release_name] = sets["cleaned"] & ~np.isin(rel, drop_rel)
    excl = None
    if args.exclusion_list:
        excl = read_exclusion_list(args.exclusion_list)
        acc_nov = np.array([re.sub(r"\.\d+$", "", a) for a in ids])
        on_list = np.isin(acc_nov, list(excl))
        absent = sorted(excl - set(acc_nov))
        print(f"[build] exclusion list: {len(excl)} accessions; {int(on_list.sum())} present in snapshot; "
              f"{len(absent)} not in snapshot")
        sets[args.exclusion_set_name] = (base & ~on_list & ~flagged) if args.exclusion_list_adds_to_classes \
            else (base & ~on_list)
    else:
        on_list = np.zeros(len(ids), bool)

    def write(name, keep):
        fn = os.path.join(args.outdir, f"{name}_{sfx}.fasta")
        with open(fn, "w") as fh:
            for k in np.where(keep)[0]:
                row = arr[k].copy()
                if not args.no_adar_mask and k in masked:
                    row[np.array(masked[k])] = b"N"
                row[MAXPOS:] = b"N"
                fh.write(f">{ids[k]}|{date_str.iloc[k]}\n{row.tobytes().decode()}\n")
        return fn

    summary = {}
    for name, keep in sets.items():
        fn = write(name, keep)
        d = parsed[keep]
        summary[name] = dict(file=os.path.basename(fn), n=int(keep.sum()), sha256=sha256(fn),
                             first_tip=str(d.min().date()), last_tip=str(d.max().date()),
                             last_tip_decimal=round(decimal_year(d.max()), 10),
                             n_with_adar_mask=int(sum(1 for k in masked if keep[k])),
                             n_masked_positions=int(sum(len(masked[k]) for k in masked if keep[k])),
                             n_tips_in_same_sample_groups=int(sum(1 for k in np.where(keep)[0] if dup_note[k])))
        print(f"[build] {name}: n={summary[name]['n']} first tip {summary[name]['first_tip']} "
              f"last tip {summary[name]['last_tip']} -> {fn}")
    if args.expected_n_cleaned is not None:
        assert summary["cleaned"]["n"] == args.expected_n_cleaned

    tab = pd.DataFrame(dict(
        accessionVersion=ids, collection_date=date_str.values, dated=dated, completeness=cf,
        screen_class=cls, release=rel, flagged_class=flagged, on_exclusion_list=on_list,
        duplicate_note=dup_note, dropped_by_duplicate_rule=dup_drop, adar_clusters=n_clusters,
        masked_positions_1based=[",".join(str(q + 1) for q in masked.get(k, [])) for k in range(len(ids))]))
    for name, keep in sets.items():
        tab[f"in_{name}"] = keep
    reason = np.where(~dated, "not dated to the day",
                      np.where(~complete, f"completeness < {args.min_completeness}",
                               np.where(dup_drop, "duplicate rule",
                                        np.where(flagged, "screen class " + pd.Series(cls).astype(str), ""))))
    tab["reason_not_in_cleaned"] = reason
    tab.to_csv(os.path.join(args.outdir, f"alignment_sets_{sfx}.csv"), index=False)

    td = pd.DataFrame(dict(taxon_acc=np.array(ids)[base], date=date_str.values[base]))
    td["decimal"] = [f"{decimal_year(x):.6f}" for x in td.date]
    td["taxon"] = td.taxon_acc + "|" + td.date
    td.to_csv(os.path.join(args.outdir, f"taxa_dates_{sfx}.csv"), index=False)

    with open(os.path.join(args.outdir, f"alignment_parameters_{sfx}.json"), "w") as fh:
        json.dump(dict(parameters=vars(args), sets=summary, n_screened=len(ids),
                       n_same_sample_groups=int(n_dup_groups),
                       counts=dict(not_dated=int((~dated).sum()),
                                   dated_low_completeness=int((dated & ~complete).sum()),
                                   dated_complete_flagged=int((dated & complete & flagged).sum())),
                       inputs={k: dict(file=os.path.basename(getattr(args, k)), sha256=sha256(getattr(args, k)))
                               for k in ("metadata", "alignment", "screen")}), fh, indent=1, default=str)
    return 0


if __name__ == "__main__":
    sys.exit(main())
