#!/usr/bin/env python3
"""make_set_tables.py (round of October 2026)

Compares the genome sets that the pipeline wrote (build_completeness_sets.py -> build_alignments.py) with the second
computation of sets_independent.py, tests the alignments, and writes the tables of step 2:
  genome_sets_by_called_fraction_<sfx>.csv     (a) one row per genome: membership in each set, first reason for absence
                                                (order of the pipeline) and first step failed (order of the article)
  genome_set_alignments_<sfx>.csv              (b) per set: file, SHA-256, number of genomes, earliest and most recent date
  set_steps_by_set_<sfx>.csv                   (c) per set: steps from 810 genomes to the set in the order of the article, each
                                                genome counted once at the first step that it fails; the same in the order of
                                                the pipeline; number of genomes that fail each rule whatever the order
  set_composition_<sfx>.csv                    (e) per set: genomes by month of collection, by health zone, by range of called fraction
  checks_step2_<sfx>.csv                       (d) tests
All tables hold public accession versions only.
usage: python make_set_tables.py --sets-dir RUN/sets --independent RUN/independent --metadata META --alignment ALN --validation VAL
                                 --suffix SFX --outdir OUT
"""
import argparse
import gzip
import hashlib
import os

import numpy as np
import pandas as pd

SETS = ["C95", "C90", "C80", "C00"]
THR = {"C95": 17955, "C90": 17010, "C80": 15120, "C00": 0}
RANGES4 = ["0.95 or more", "0.90 to below 0.95", "0.80 to below 0.90", "below 0.80"]
LEFT_OUT = ["E", "E-", "D", "M", "H"]
CHECKS = []


def check(name, expected, observed, ok):
    CHECKS.append(dict(step="2 sets", check=name, expected=str(expected), observed=str(observed),
                       passed="yes" if ok is True else "no" if ok is False else ok))


def sha(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for c in iter(lambda: fh.read(1 << 20), b""):
            h.update(c)
    return h.hexdigest()


def read_fasta_records(path):
    out, name, buf = [], None, []
    op = gzip.open if str(path).endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:], []
            else:
                buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sets-dir", required=True)
    ap.add_argument("--independent", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--validation", required=True)
    ap.add_argument("--suffix", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    sfx = a.suffix
    os.makedirs(a.outdir, exist_ok=True)
    P = pd.read_csv(os.path.join(a.sets_dir, f"genome_sets_by_called_fraction_{sfx}.csv"), dtype=str, keep_default_na=False).set_index("accessionVersion")
    I = pd.read_csv(os.path.join(a.independent, f"independent_membership_{sfx}.csv"), dtype=str, keep_default_na=False).set_index("accessionVersion")
    meta = pd.read_csv(a.metadata, dtype=str, keep_default_na=False, low_memory=False).set_index("accessionVersion")
    val = pd.read_csv(a.validation, dtype=str, keep_default_na=False, low_memory=False).set_index("accessionVersion")
    raw = dict(read_fasta_records(a.alignment))
    assert len(P) == 810 and set(P.index) == set(I.index) == set(raw)
    I = I.loc[P.index]
    called = I.called_positions.astype(int)

    # ---- (d) membership and reasons: pipeline against the second computation ----
    for s in SETS:
        eq_m = int((P[f"in_{s}"] == I[f"in_{s}"]).sum())
        eq_r = int((P[f"reason_not_in_{s}"] == I[f"reason_pipeline_order_{s}"]).sum())
        check(f"{s}: membership written by the pipeline against the logical AND of the five rules evaluated by code that does not call the pipeline",
              "equal for 810 of 810 genomes", f"equal for {eq_m} of 810 genomes; in the set: pipeline {int((P[f'in_{s}'] == 'True').sum())}, second computation {int((I[f'in_{s}'] == 'True').sum())}", eq_m == 810)
        check(f"{s}: first reason for absence written by the pipeline against the second computation (order of the pipeline: date, called fraction, class, shared sample identifier, restricted use)",
              "equal for 810 of 810 genomes", f"equal for {eq_r} of 810 genomes", eq_r == 810)
        mem = I[f"in_{s}"] == "True"
        one = int(((I[f"step_article_order_{s}"] != "") == ~mem).sum())
        check(f"{s}: every absent genome has one first step in the order of the article and every genome of the set has none",
              "810 of 810", f"{one} of 810", one == 810)
    M = {s: (P[f"in_{s}"] == "True").values for s in SETS}
    for hi, lo in [("C95", "C90"), ("C90", "C80"), ("C80", "C00"), ("C95", "C80"), ("C95", "C00"), ("C90", "C00")]:
        n_out = int((M[hi] & ~M[lo]).sum())
        check(f"nested sets: genomes of {hi} that are not in {lo}", "0", n_out, n_out == 0)
    r_date = (I.rule_date == "True").values
    r_id = (I.rule_identifier == "True").values
    r_use = (I.rule_restricted_use == "True").values
    r_cls = (I.rule_class == "True").values
    for s in SETS:
        r_cf = (called >= THR[s]).values
        bad = dict(undated=int((M[s] & ~r_date).sum()), below_threshold=int((M[s] & ~r_cf).sum()), left_out_by_the_rule_on_shared_identifiers=int((M[s] & ~r_id).sum()),
                   restricted_use=int((M[s] & ~r_use).sum()), class_E_Eminus_D_M_H=int((M[s] & ~r_cls).sum()))
        check(f"{s}: genomes of the set that fail a rule", "0 for each of the five rules", bad, sum(bad.values()) == 0)
        miss = int((~M[s] & r_date & r_cf & r_id & r_use & r_cls).sum())
        check(f"{s}: genomes that meet all five rules and are absent from the set", "0", miss, miss == 0)
    cls_p = P.screen_class
    check("class in the table of the pipeline against the class of the screen table (read by the second computation)", "equal for 810 of 810",
          f"equal for {int((cls_p == I['class']).sum())} of 810", bool((cls_p == I['class']).all()))
    check("called fraction in the table of the pipeline against called positions / 18,900 of the second computation", "equal for 810 of 810 (absolute difference below 1e-12)",
          f"equal for {int(np.isclose(P.called_fraction_exact.astype(float).values, called.values / 18900, atol=1e-12).sum())} of 810",
          bool(np.isclose(P.called_fraction_exact.astype(float).values, called.values / 18900, atol=1e-12).all()))
    check("positions masked as ADAR-type clusters: table of the pipeline against the second computation", "equal for 810 of 810 genomes",
          f"equal for {int((P.masked_positions_1based == I.masked_positions_1based).sum())} of 810; genomes with a mask: {int((I.masked_positions_1based != '').sum())}",
          bool((P.masked_positions_1based == I.masked_positions_1based).all()))
    lv = meta.loc[P.index, "versionStatus"].eq("LATEST_VERSION") & meta.loc[P.index, "isRevocation"].ne("true")
    check("every screened genome is a latest version without revocation in the metadata of the snapshot", "810 of 810", f"{int(lv.sum())} of 810", bool(lv.all()))

    # ---- (b) alignments ----
    rowsA = []
    for s in SETS:
        f1 = os.path.join(a.sets_dir, s, f"{s}_{sfx}.fasta")
        f2 = os.path.join(a.sets_dir, "all", f"{s}_{sfx}.fasta")
        recs = read_fasta_records(f1)
        names = [h.split("|")[0] for h, _ in recs]
        dates = [h.split("|")[1] for h, _ in recs]
        members = [x for x in P.index if M[s][list(P.index).index(x)]] if False else list(P.index[M[s]])
        check(f"{s}: records of the alignment against the genomes of the set", "the same accession versions, each once",
              f"{len(recs)} records, {len(set(names))} distinct; genomes of the set {len(members)}; in both {len(set(names) & set(members))}",
              len(recs) == len(set(names)) == len(members) and set(names) == set(members))
        ok_date = sum(d == meta.loc[n, "sampleCollectionDate"] for n, d in zip(names, dates))
        check(f"{s}: collection date in the header of every record against the metadata of the new snapshot", f"{len(recs)} of {len(recs)} equal", f"{ok_date} of {len(recs)} equal", ok_date == len(recs))
        n_len = sum(len(q) == 18940 for _, q in recs)
        n_seq_ok, n_mask_ok, n_end_ok, chars = 0, 0, 0, set()
        for n, (_, q) in zip(names, recs):
            src = raw[n].upper()
            mp = [int(x) for x in I.loc[n, "masked_positions_1based"].split(",") if x]
            exp = list(src)
            for p_ in mp:
                exp[p_ - 1] = "N"
            exp[18900:] = "N" * 40
            n_seq_ok += ("".join(exp) == q)
            n_mask_ok += all(q[p_ - 1] == "N" for p_ in mp)
            n_end_ok += (q[18900:] == "N" * 40)
            chars |= set(q)
        check(f"{s}: sequence of every record against the sequence of the snapshot with the masks of the second computation (ADAR-type clusters except their first position; positions 18,901 to 18,940)",
              f"{len(recs)} of {len(recs)} identical, length 18,940", f"{n_seq_ok} of {len(recs)} identical; length 18,940: {n_len}; characters: {''.join(sorted(chars))}", n_seq_ok == len(recs) == n_len)
        h1, h2 = sha(f1), sha(f2)
        summ = pd.read_csv(os.path.join(a.sets_dir, s, f"alignment_set_summary_{sfx}.csv")).set_index("set").loc[s]
        check(f"{s}: SHA-256 of the alignment: recomputed, written by the pipeline (alignment_set_summary), and of the run with the four sets in one specification",
              "three equal values", f"{h1}; summary equal: {summ['sha256'] == h1}; run 'all' equal: {h2 == h1}", bool(summ['sha256'] == h1 and h2 == h1))
        d = pd.Series(dates)
        dm = meta.loc[members, "sampleCollectionDate"]
        check(f"{s}: earliest and most recent collection date: headers of the alignment, metadata of the members, summary of the pipeline",
              "three equal pairs", f"headers {d.min()} to {d.max()}; metadata {dm.min()} to {dm.max()}; pipeline {summ['first_tip']} to {summ['last_tip']}",
              (d.min(), d.max()) == (dm.min(), dm.max()) == (summ['first_tip'], summ['last_tip']))
        rowsA.append(dict(set=s, called_fraction="no threshold" if s == "C00" else f"{int(s[1:]) / 100:.2f} or more", called_positions_needed_of_18900=THR[s] if s != "C00" else "",
                          alignment=f"{s}_{sfx}.fasta", sha256=h1, bytes=os.path.getsize(f1), genomes=len(recs), alignment_length=18940,
                          earliest_collection_date=d.min(), most_recent_collection_date=d.max(), genomes_on_the_most_recent_date=int((d == d.max()).sum()),
                          second_most_recent_collection_date=d[d < d.max()].max(),
                          variable_sites_counted_by_the_pipeline=int(summ["n_variable_sites"]), parsimony_informative_sites_counted_by_the_pipeline=int(summ["n_parsimony_informative_sites"]),
                          genomes_with_masked_adar_positions=int(summ["n_with_adar_mask"]), masked_adar_positions=int(summ["n_masked_adar_positions"])))
    A = pd.DataFrame(rowsA)
    A.to_csv(os.path.join(a.outdir, f"genome_set_alignments_{sfx}.csv"), index=False)

    # ---- (a) table of all genomes ----
    T = pd.DataFrame(index=pd.Index(P.index, name="accession_version"))
    T["collection_date_as_recorded"] = meta.loc[P.index, "sampleCollectionDate"].values
    T["collection_date_given_to_the_day"] = r_date
    T["collection_month"] = [d[:7] if ok else "(no date to the day)" for d, ok in zip(T.collection_date_as_recorded, r_date)]
    T["country"] = val.loc[P.index, "geoLocCountry"].replace("", "(not recorded)").values
    T["admin1_as_recorded"] = val.loc[P.index, "geoLocAdmin1"].replace("", "(not recorded)").values
    T["health_zone"] = val.loc[P.index, "zone"].replace("", "(not recorded)").values
    T["earliest_release_date"] = meta.loc[P.index, "earliestReleaseDate"].values
    T["called_positions_1_to_18900"] = called.values
    T["called_fraction"] = called.values / 18900
    T["range_of_called_fraction"] = np.where(called >= 17955, RANGES4[0], np.where(called >= 17010, RANGES4[1], np.where(called >= 15120, RANGES4[2], RANGES4[3])))
    T["class_of_the_quality_checks"] = I["class"].values
    T["class_is_left_out"] = ~r_cls
    T["shares_sample_identifier"] = (I.partner != "").values
    T["record_with_the_same_sample_identifier"] = I.partner.values
    T["positions_that_differ_from_that_record"] = I.differences_with_partner.values
    T["left_out_by_the_rule_on_shared_sample_identifiers"] = ~r_id
    T["restricted_use_agreement_pending"] = ~r_use
    T["adar_type_clusters"] = I.adar_clusters.values
    T["masked_positions_1based"] = I.masked_positions_1based.values
    for s in SETS:
        T[f"in_{s}"] = M[s]
    for s in SETS:
        T[f"reason_not_in_{s}_order_of_the_pipeline"] = P[f"reason_not_in_{s}"].values
    for s in SETS:
        T[f"first_step_failed_{s}_order_of_the_article"] = I[f"step_article_order_{s}"].values
    T["in_no_set"] = ~M["C00"]
    T.to_csv(os.path.join(a.outdir, f"genome_sets_by_called_fraction_{sfx}.csv"))

    # ---- (c) steps ----
    rows = []
    id_disc = (I.partner != "").values & ~r_id & (I.differences_with_partner.replace("", "0").astype(int) > 0).values
    id_conc = ~r_id & ~id_disc
    for s in SETS:
        r_cf = (called >= THR[s]).values
        n = 810
        rows.append(dict(set=s, order="article", step=0, rule="genomes of the outbreak set", genomes_that_fail_at_this_step_first="", genomes_left_after_the_step=n, genomes_that_fail_the_rule_whatever_the_order=""))
        left = np.ones(810, bool)
        steps_article = [("no collection date to the day", ~r_date), ("called fraction below the threshold", ~r_cf),
                         ("shares a sample identifier with a record that differs, or with a concordant record that is more complete", ~r_id),
                         ("restricted use, agreement of the submitters pending", ~r_use), ("did not pass the quality checks (class E, E-, D, M or H)", ~r_cls)]
        for k, (nm, fail) in enumerate(steps_article, 1):
            f = left & fail
            left &= ~fail
            rows.append(dict(set=s, order="article", step=k, rule=nm, genomes_that_fail_at_this_step_first=int(f.sum()), genomes_left_after_the_step=int(left.sum()),
                             genomes_that_fail_the_rule_whatever_the_order=int(fail.sum())))
            if k == 3:
                for nm2, f2 in [("of these: with a record that differs (both records are left out)", id_disc), ("of these: with a concordant record that is more complete (the less complete record is left out)", id_conc)]:
                    rows.append(dict(set=s, order="article", step="3 (part)", rule=nm2, genomes_that_fail_at_this_step_first=int((f & f2).sum()), genomes_left_after_the_step="",
                                     genomes_that_fail_the_rule_whatever_the_order=int(f2.sum())))
            if k == 5:
                for c in LEFT_OUT:
                    cm = (I["class"] == c).values
                    rows.append(dict(set=s, order="article", step="5 (part)", rule=f"of these: class {c}", genomes_that_fail_at_this_step_first=int((f & cm).sum()), genomes_left_after_the_step="",
                                     genomes_that_fail_the_rule_whatever_the_order=int(cm.sum())))
        assert int(left.sum()) == int(M[s].sum()), s
        rows.append(dict(set=s, order="article", step="set", rule=f"genomes of the set {s}", genomes_that_fail_at_this_step_first="", genomes_left_after_the_step=int(M[s].sum()), genomes_that_fail_the_rule_whatever_the_order=""))
        left = np.ones(810, bool)
        rows.append(dict(set=s, order="pipeline", step=0, rule="genomes of the outbreak set", genomes_that_fail_at_this_step_first="", genomes_left_after_the_step=n, genomes_that_fail_the_rule_whatever_the_order=""))
        steps_pipe = [("collection date not given to the day", ~r_date), ("called fraction below threshold", ~r_cf), ("class E, E-, D, M or H", ~r_cls),
                      ("shares a sample identifier (discordant record, or concordant and more complete record)", ~r_id), ("restricted use, agreement of the submitters pending", ~r_use)]
        for k, (nm, fail) in enumerate(steps_pipe, 1):
            f = left & fail
            left &= ~fail
            rows.append(dict(set=s, order="pipeline", step=k, rule=nm, genomes_that_fail_at_this_step_first=int(f.sum()), genomes_left_after_the_step=int(left.sum()),
                             genomes_that_fail_the_rule_whatever_the_order=int(fail.sum())))
        rows.append(dict(set=s, order="pipeline", step="set", rule=f"genomes of the set {s}", genomes_that_fail_at_this_step_first="", genomes_left_after_the_step=int(M[s].sum()), genomes_that_fail_the_rule_whatever_the_order=""))
    ST = pd.DataFrame(rows)
    ST.to_csv(os.path.join(a.outdir, f"set_steps_by_set_{sfx}.csv"), index=False)
    for s in SETS:
        for order in ("article", "pipeline"):
            x = ST[(ST.set == s) & (ST.order == order) & ST.step.astype(str).isin(list("12345"))]
            tot = int(x.genomes_that_fail_at_this_step_first.astype(int).sum()) + int(M[s].sum())
            check(f"{s}: genomes counted at their first step ({order} order) plus genomes of the set", "810", tot, tot == 810)
        # the table of steps against the reasons of the pipeline table
        vc = P[f"reason_not_in_{s}"].value_counts()
        x = ST[(ST.set == s) & (ST.order == "pipeline")].set_index("step")
        obs = dict(date=int(x.loc[1, "genomes_that_fail_at_this_step_first"]), called=int(x.loc[2, "genomes_that_fail_at_this_step_first"]), cls=int(x.loc[3, "genomes_that_fail_at_this_step_first"]),
                   ident=int(x.loc[4, "genomes_that_fail_at_this_step_first"]), use=int(x.loc[5, "genomes_that_fail_at_this_step_first"]))
        exp = dict(date=int(vc.get("collection date not given to the day", 0)), called=int(vc.get("called fraction below threshold", 0)),
                   cls=int(sum(v for k_, v in vc.items() if k_.startswith("class "))), ident=int(sum(v for k_, v in vc.items() if k_.startswith("shares a sample identifier"))),
                   use=int(vc.get("restricted use, agreement of the submitters pending", 0)))
        check(f"{s}: counts of the steps in the order of the pipeline against the counts of the reasons in the table that the pipeline wrote", exp, obs, exp == obs)

    # ---- (e) composition ----
    rows = []
    for s in SETS + ["all 810 genomes"]:
        m = M[s] if s in M else np.ones(810, bool)
        sub = T[m]
        for k, v in sub.collection_month.value_counts().sort_index().items():
            rows.append(dict(set=s, kind="month of collection", country="", admin1_as_recorded="", value=k, genomes=int(v)))
        z = sub.groupby(["country", "admin1_as_recorded", "health_zone"]).size().sort_values(ascending=False)
        for (c, a1, hz), v in z.items():
            rows.append(dict(set=s, kind="health zone", country=c, admin1_as_recorded=a1, value=hz, genomes=int(v)))
        for g in RANGES4:
            rows.append(dict(set=s, kind="range of called fraction", country="", admin1_as_recorded="", value=g, genomes=int((sub.range_of_called_fraction == g).sum())))
        for c in ["A", "A*", "0", "C", "H", "E", "E-", "D", "M"]:
            rows.append(dict(set=s, kind="class of the quality checks", country="", admin1_as_recorded="", value=c, genomes=int((sub.class_of_the_quality_checks == c).sum())))
    C = pd.DataFrame(rows)
    C.to_csv(os.path.join(a.outdir, f"set_composition_{sfx}.csv"), index=False)
    for s in SETS:
        for kind in ["month of collection", "health zone", "range of called fraction", "class of the quality checks"]:
            tot = int(C[(C.set == s) & (C.kind == kind)].genomes.sum())
            check(f"{s}: composition by {kind} adds up to the size of the set", int(M[s].sum()), tot, tot == int(M[s].sum()))
        # against the tables by month and by zone that the pipeline wrote
        pm = pd.read_csv(os.path.join(a.sets_dir, s, f"alignment_sets_by_month_{sfx}.csv"), dtype=str).set_index("collection_month")[s].astype(int)
        om = C[(C.set == s) & (C.kind == "month of collection")].set_index("value").genomes
        check(f"{s}: genomes by month of collection against the table of the pipeline", pm.to_dict(), om.to_dict(), pm.to_dict() == om.to_dict())
        pz = pd.read_csv(os.path.join(a.sets_dir, s, f"alignment_sets_by_zone_{sfx}.csv"), dtype=str, keep_default_na=False).set_index("health_zone")[s].astype(int)
        oz = T[M[s]].groupby("health_zone").size()
        check(f"{s}: genomes by health zone (zones of the same name in two provinces counted together, as the pipeline does) against the table of the pipeline",
              f"{len(pz)} zones, equal counts", f"{len(oz)} zones; equal counts: {pz.sort_index().to_dict() == oz.sort_index().to_dict()}", pz.sort_index().to_dict() == oz.sort_index().to_dict())
    # genomes next to a threshold
    for s in ["C95", "C90", "C80"]:
        r_all = r_date & r_id & r_use & r_cls
        above = called[(called >= THR[s]).values].min(); below = called[(called < THR[s]).values].max()
        check(f"{s}: called positions of the genomes next to the threshold of {THR[s]} positions", f"smallest value at or above the threshold >= {THR[s]}; largest value below it < {THR[s]}",
              f"smallest at or above: {int(above)}; largest below: {int(below)}", bool(above >= THR[s] > below))
    pd.DataFrame(CHECKS).to_csv(os.path.join(a.outdir, f"checks_step2_{sfx}.csv"), index=False)
    print(pd.DataFrame(CHECKS).passed.value_counts().to_dict())
    bad = [c for c in CHECKS if c["passed"] != "yes"]
    for c in bad:
        print("NOT PASSED:", c)
    print(A.to_string(index=False))


if __name__ == "__main__":
    main()
