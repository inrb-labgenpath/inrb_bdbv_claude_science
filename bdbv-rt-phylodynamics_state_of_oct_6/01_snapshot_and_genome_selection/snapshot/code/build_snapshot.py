"""Step 2 - the derivations for ONE download, with the names of the files as parameters.

This is the build script of 2 October (03_build_20261002.py) turned into a function, so that
the SAME code runs (a) on the download of 2 October, as the control, and (b) on the download of 7 October.
The derivations themselves are in snapshot_lib.py, which is byte for byte the file of 2 October.
What differs from the script of 2 October:
  * the paths, the tag of the file names, the earlier zone map and the baseline are parameters;
  * the baseline may hold fewer database fields than the download (the comparison columns use the shared fields);
  * names of facts that carried a date or a count in their name are neutral ("aligned_alphabet", "n_stored_map", ...);
  * the list of spellings that the earlier map does not hold is also written as a file of its own.
"""
import json
import pathlib
from collections import Counter

import numpy as np
import pandas as pd

import snapshot_lib as sl

UNMAPPED_COLS = ["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "geoLocAdmin2_as_python_string", "n_genomes", "accessionVersions",
                 "zone_assigned", "flag", "information_only_zone_that_the_rule_of_the_code_of_20260924_would_give_NOT_APPLIED",
                 "information_only_rule_NOT_APPLIED"]


def run(tag, raw_dir, out_dir, freeze, zone_map_path, baseline_meta=None, baseline_fasta=None, baseline_tag=None, verbose=True):
    RAW = pathlib.Path(raw_dir)
    OUT = pathlib.Path(out_dir)
    OUT.mkdir(parents=True, exist_ok=True)
    TAG = tag
    assert freeze["all_sha256_match"] and freeze["all_sizes_match"]
    T_NEW = freeze["details_allversions_query_start_utc"]
    DV_NEW = freeze["data_version"]
    F = {"tag": TAG, "snapshotDatetimeUtc": T_NEW, "lapisDataVersion": DV_NEW, "baseline_tag": baseline_tag or ""}

    full = sl.read_csv_str(RAW / "details_allversions.csv")
    filt = sl.read_csv_str(RAW / "details_latest_nonrevoked.csv")
    al, al_info = sl.read_fasta(RAW / "aligned_latest_nonrevoked.fasta")
    un, un_info = sl.read_fasta(RAW / "unaligned_latest_nonrevoked.fasta")
    F["n_records_all_versions"] = len(full)
    F["n_records_latest_nonrevoked"] = len(filt)
    F["n_fields"] = len(full.columns)
    F["fields"] = list(full.columns)
    F["fields_same_in_both_details_files"] = bool(list(full.columns) == list(filt.columns))
    F["n_accessions"] = int(full.accession.nunique())
    F["accessionVersion_unique_all_versions"] = bool(full.accessionVersion.is_unique)
    F["accessionVersion_equals_accession_dot_version"] = bool((full.accessionVersion == full.accession + "." + full.version).all())
    F["version_status_x_revocation"] = {f"{a}|isRevocation={b}": int(n) for (a, b), n in full.groupby(["versionStatus", "isRevocation"]).size().items()}

    # ---- the filtered download is the subset LATEST_VERSION & not revocation of the unfiltered one
    sub = full[(full.versionStatus == "LATEST_VERSION") & (full.isRevocation.str.lower() == "false")]
    a = sub.set_index("accessionVersion").sort_index()
    b = filt.set_index("accessionVersion").sort_index()
    F["filtered_is_subset_same_accessionVersions"] = bool(list(a.index) == list(b.index))
    F["filtered_cells_different_from_unfiltered"] = int((a.values != b.values).sum()) if a.shape == b.shape else None
    F["one_latest_version_per_accession"] = bool(full[full.versionStatus == "LATEST_VERSION"].accession.is_unique)

    # ---- records and sequences (hard requirement: nothing dropped silently)
    al_heads = [h for h, _ in al]
    un_heads = [h for h, _ in un]
    seqs = dict(al)
    useqs = dict(un)
    F["aligned"] = dict(al_info, n_records=len(al), n_unique_headers=len(set(al_heads)),
                        lengths=dict(Counter(len(s) for _, s in al)))
    F["unaligned"] = dict(un_info, n_records=len(un), n_unique_headers=len(set(un_heads)),
                          n_distinct_lengths=len(set(len(s) for _, s in un)),
                          min_length=min(len(s) for _, s in un), max_length=max(len(s) for _, s in un))
    F["aligned_headers_repeated"] = sorted(h for h, n in Counter(al_heads).items() if n > 1)
    F["unaligned_headers_repeated"] = sorted(h for h, n in Counter(un_heads).items() if n > 1)
    F["latest_nonrevoked_records_without_aligned_sequence"] = sorted(set(filt.accessionVersion) - set(seqs))
    F["aligned_sequences_without_record"] = sorted(set(seqs) - set(filt.accessionVersion))
    F["aligned_sequences_of_other_length"] = sorted(h for h, s in al if len(s) != sl.L_ALN)
    F["aligned_sequences_empty"] = sorted(h for h, s in al if len(s) == 0)
    F["aligned_sequences_without_any_called_base"] = sorted(h for h, s in al if sl.count_acgt(s, sl.L_ALN) == 0)
    F["latest_nonrevoked_records_without_unaligned_sequence"] = sorted(set(filt.accessionVersion) - set(useqs))
    F["unaligned_sequences_without_record"] = sorted(set(useqs) - set(filt.accessionVersion))
    F["unaligned_sequences_empty"] = sorted(h for h, s in un if len(s) == 0)
    alpha = Counter()
    for _, s in al:
        alpha.update(s)
    F["aligned_alphabet"] = dict(sorted(alpha.items()))
    F["aligned_lowercase_characters"] = int(sum(n for c, n in alpha.items() if c.islower()))
    ualpha = Counter()
    for _, s in un:
        ualpha.update(s)
    F["unaligned_alphabet"] = dict(sorted(ualpha.items()))

    # ---- the unaligned file: what it holds (nothing below uses it)
    len_field = filt.set_index("accessionVersion")["length"].astype(int)
    ins_field = filt.set_index("accessionVersion")["totalInsertedNucs"].astype(int)
    ob_field = set(filt.accessionVersion[filt.outbreak == "Bdbv-2026"])
    core = lambda x: x.strip("Nn")  # noqa: E731
    U = pd.DataFrame({"len_un": {k: len(v) for k, v in useqs.items()}})
    U["length_field"] = len_field
    U["acgt_un"] = pd.Series({k: sl.count_acgt(v, len(v)) for k, v in useqs.items()})
    U["acgt_al"] = pd.Series({k: sl.count_acgt(v, sl.L_ALN) for k, v in seqs.items()})
    U["equal_to_aligned_without_gaps_and_outer_N"] = pd.Series({k: core(useqs[k]) == core(seqs[k].replace("-", "")) for k in useqs if k in seqs})
    U["totalInsertedNucs"] = ins_field
    U["field_outbreak_Bdbv2026"] = U.index.isin(ob_field)
    uf = {}
    for nm, usub in (("all latest non-revoked records", U), ("records with field outbreak = Bdbv-2026", U[U.field_outbreak_Bdbv2026]),
                    ("other records", U[~U.field_outbreak_Bdbv2026])):
        extra = usub.acgt_un - usub.acgt_al
        eq = usub.equal_to_aligned_without_gaps_and_outer_N.fillna(False).astype(bool)
        uf[nm] = dict(n=len(usub), length_min=int(usub.len_un.min()), length_median=float(usub.len_un.median()), length_max=int(usub.len_un.max()),
                      n_length_equals_field_length=int((usub.len_un == usub.length_field).sum()),
                      n_equal_to_aligned_without_gaps_and_outer_N=int(eq.sum()),
                      n_with_more_called_bases_than_aligned=int((extra > 0).sum()), n_with_fewer_called_bases_than_aligned=int((extra < 0).sum()),
                      extra_called_bases_distribution={int(k): int(v) for k, v in extra.value_counts().sort_index().items()},
                      n_with_insertions_by_field=int((usub.totalInsertedNucs > 0).sum()),
                      n_not_equal_and_no_insertion_by_field=int(((~eq) & (usub.totalInsertedNucs == 0)).sum()),
                      n_equal_but_insertion_by_field=int((eq & (usub.totalInsertedNucs > 0)).sum()),
                      n_shorter_than_15000=int((usub.len_un < 15000).sum()),
                      accessionVersions_shorter_than_15000=sorted(usub.index[usub.len_un < 15000]))
    F["unaligned_file"] = uf
    F["unaligned_gap_characters"] = int(sum(v.count("-") for v in useqs.values()))

    # ---- unreadable dates, all records of all versions
    date_fields = ["sampleCollectionDate", "sampleCollectionDateRangeLower", "sampleCollectionDateRangeUpper", "earliestReleaseDate",
                   "releasedDate", "submittedDate", "sampleReceivedDate", "sequencingDate", "dataUseTermsRestrictedUntil",
                   "dataBecameOpenAt", "ncbiReleaseDate", "ncbiUpdateDate"]
    unread = {}
    for c in date_fields:
        v = full[c]
        bad = full.accessionVersion[(v != "") & sl.iso(v).isna()].tolist()
        unread[c] = dict(n_filled=int((v != "").sum()), n_not_read_by_iso8601_parser=len(bad), accessionVersions=bad)
    F["dates_not_read_by_iso8601_parser"] = unread
    prec_all = full.sampleCollectionDate.map(sl.precision)
    F["collection_date_precision_all_records"] = prec_all.value_counts().to_dict()
    F["collection_date_precision_latest_nonrevoked"] = prec_all[full.accessionVersion.isin(set(filt.accessionVersion))].value_counts().to_dict()
    F["collection_date_precision_other"] = full.accessionVersion[prec_all == "other"].tolist()

    # ---- derivations (the baseline feeds the comparison columns only)
    zm_stored, look = sl.load_zone_map(zone_map_path)
    derived = sl.DERIVED16 + sl.NEW_COLS
    if baseline_meta is not None:
        om = sl.read_csv_str(baseline_meta)
        old_fields = [c for c in om.columns if c not in derived]
        F["baseline_n_fields"] = len(old_fields)
        F["fields_added_since_baseline"] = [c for c in full.columns if c not in old_fields]
        F["fields_removed_since_baseline"] = [c for c in old_fields if c not in full.columns]
        F["shared_fields_in_same_order"] = bool([c for c in full.columns if c in old_fields] == [c for c in old_fields if c in full.columns])
        old_seqs = dict(sl.read_fasta(baseline_fasta)[0])
        res = sl.build(full, seqs, look, T_NEW, DV_NEW, baseline_full=om[old_fields], baseline_seqs=old_seqs, add_new_columns=True)
    else:
        res = sl.build(full, seqs, look, T_NEW, DV_NEW, baseline_full=None, baseline_seqs=None, add_new_columns=True)
    o, M, V = res["o"], res["M"], res["V"]

    # ---- the outbreak set and its cross-checks
    set_def = set(o.index)
    set_iso = set(res["ob_iso"].accessionVersion)
    F["n_outbreak"] = len(set_def)
    F["outbreak_iso8601_parser_same_set"] = bool(set_def == set_iso)
    F["outbreak_definition_arm_first10"] = o.definition_arm_first10.value_counts().to_dict()
    F["outbreak_definition_arm_iso8601"] = o.definition_arm_iso8601.value_counts().to_dict()
    F["outbreak_through_release_arm_iso8601"] = sorted(o.index[o.definition_arm_iso8601 == "release"])
    F["outbreak_through_release_arm_first10"] = sorted(o.index[o.definition_arm_first10 == "release"])
    F["outbreak_field_in_set"] = o.outbreak.replace("", "(unassigned)").value_counts().to_dict()
    rest = sub[~sub.accessionVersion.isin(set_def)]
    F["outbreak_field_latest_nonrevoked_outside_set"] = rest.outbreak.replace("", "(unassigned)").value_counts().to_dict()
    F["latest_nonrevoked_outside_set_with_field_Bdbv2026"] = sorted(rest.accessionVersion[rest.outbreak == "Bdbv-2026"])
    F["in_set_without_field_Bdbv2026"] = sorted(o.index[o.outbreak != "Bdbv-2026"])
    F["outbreak_field_all_records"] = {f"{a}|{b}|isRevocation={c}": int(n) for (a, b, c), n in
                                       full.assign(outbreak=full.outbreak.replace("", "(unassigned)")).groupby(["outbreak", "versionStatus", "isRevocation"]).size().items()}
    F["outbreak_without_aligned_sequence"] = sorted(k for k in set_def if k not in seqs)
    F["outbreak_aligned_length_not_18940"] = sorted(k for k in set_def if k in seqs and len(seqs[k]) != sl.L_ALN)
    # records of the latest versions that the definition leaves out, with the reason (every record accounted for)
    lat = full[full.versionStatus == "LATEST_VERSION"].copy()
    lat["_coll"] = lat.sampleCollectionDate.map(sl.norm_date)
    lat["_rel"] = lat.earliestReleaseDate.map(sl.norm_date)

    def reason(av, rev, coll):
        if av in set_def:
            return "in the outbreak set"
        if rev.lower() == "true":
            return "revocation entry"
        if pd.isna(coll):
            return "collection date missing or not readable, and earliest release before 2026-05-01"
        return "collection date before 2026-01-01"

    lat["reason"] = [reason(av, rev, coll) for av, rev, coll in zip(lat.accessionVersion, lat.isRevocation, lat["_coll"])]
    F["latest_versions_by_reason"] = lat.reason.value_counts().to_dict()
    F["n_latest_versions"] = len(lat)
    F["records_by_status_and_place_in_set"] = {f"{a}|isRevocation={b}|in set={c}": int(n) for (a, b, c), n in
                                               full.assign(_in=full.accessionVersion.isin(set_def)).groupby(["versionStatus", "isRevocation", "_in"]).size().items()}

    # ---- rule on blanks: every record (all versions) and every place field on which it acts
    OB = sl.outer_blanks(full)
    OB["latest_nonrevoked"] = OB.accessionVersion.isin(set(sub.accessionVersion))
    OB["in_outbreak_set"] = OB.accessionVersion.isin(set(o.index))
    OB["field_used_by_zone_map"] = OB.field.isin(["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2"])
    OB["zone_assigned"] = OB.accessionVersion.map(o.zone).fillna("")
    OB["zone_rule"] = OB.accessionVersion.map(o.zone_rule).fillna("")
    F["place_fields_checked"] = [f for f in sl.PLACE_FIELDS if f in full.columns]
    F["place_fields_filled"] = {f: int((full[f] != "").sum()) for f in F["place_fields_checked"]}
    F["rule_on_blanks_records"] = OB.to_dict("records")
    F["rule_on_blanks_n_records"] = int(OB.accessionVersion[OB.rule_acts].nunique())
    F["rule_on_blanks_n_fields"] = int(OB.rule_acts.sum())
    F["other_white_space_at_ends_n_fields"] = int(OB.other_white_space_left.sum())
    F["rule_on_blanks_acted_in_outbreak_set"] = sorted(o.index[o.rule_on_blanks_acted])
    assert set(F["rule_on_blanks_acted_in_outbreak_set"]) == set(OB.accessionVersion[OB.rule_acts & OB.in_outbreak_set & OB.field_used_by_zone_map])
    inner = {f: int(full[f].str.contains(r"[\t\r\n\u00a0\u2000-\u200b\u202f\u3000]|  ", regex=True).sum()) for f in F["place_fields_checked"]}
    F["place_fields_with_tab_nonbreaking_or_double_blank_inside"] = {k: v for k, v in inner.items() if v}

    # ---- zone map: raw spellings that the earlier map does not hold (flagged and listed; no zone assigned)
    nm = o[o.flag_zone_not_in_map]
    unm_rows = []
    for (c, a1, a2), g in nm.groupby(["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2"]):
        info = sl.rule_of_20260924_code(c, a1, a2)
        assert (g.zone == "").all(), "a spelling that the map does not hold was given a zone"
        unm_rows.append(dict(geoLocCountry=c, geoLocAdmin1=a1, geoLocAdmin2=a2, geoLocAdmin2_as_python_string=repr(a2), n_genomes=int(len(g)),
                             accessionVersions=";".join(sorted(g.index)), zone_assigned="", flag="zone_spelling_not_in_map",
                             information_only_zone_that_the_rule_of_the_code_of_20260924_would_give_NOT_APPLIED=info["zone"],
                             information_only_rule_NOT_APPLIED=info["zone_rule"]))
    UNM = pd.DataFrame(unm_rows, columns=UNMAPPED_COLS)
    F["zone_keys_not_in_stored_map"] = UNM.to_dict("records")
    F["n_zone_keys_not_in_stored_map"] = len(UNM)
    F["n_genomes_zone_not_in_stored_map"] = int(o.flag_zone_not_in_map.sum())
    assert F["n_genomes_zone_not_in_stored_map"] == (int(UNM.n_genomes.sum()) if len(UNM) else 0)
    zmap = res["zmap"].copy()
    used = set(zip(zmap.geoLocCountry, zmap.geoLocAdmin1, zmap.geoLocAdmin2))
    unused = zm_stored[[k not in used for k in zip(zm_stored.geoLocCountry, zm_stored.geoLocAdmin1, zm_stored.geoLocAdmin2)]].copy()
    F["stored_map_rows"] = len(zm_stored)
    F["stored_map_rows_without_genome_in_this_set"] = unused[["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "zone"]].to_dict("records")
    zmap["n_genomes"] = zmap["n_genomes"].astype(int)
    info_rule = [sl.rule_of_20260924_code(c, a1, a2)["zone"] for c, a1, a2 in zip(zmap.geoLocCountry, zmap.geoLocAdmin1, zmap.geoLocAdmin2)]
    mapped = (zmap.flag_zone_not_in_map.astype(str) == "False").values
    assert all(z == i for z, i, m in zip(zmap.zone, info_rule, mapped) if m), "rule of the earlier code and stored map disagree on a mapped key"
    # every row of the new map whose key the earlier map holds carries the values of the earlier map
    st_rows = {k: r for k, r in zip(zip(zm_stored.geoLocCountry, zm_stored.geoLocAdmin1, zm_stored.geoLocAdmin2), zm_stored.to_dict("records"))}
    cmp_cols = ["zone", "zone_rule", "zone_ref_province", "flag_province_conflict", "flag_zone_not_in_reference", "flag_zone_missing", "flag_zone_interpreted"]
    n_same = n_other = 0
    for r in zmap.to_dict("records"):
        k = (r["geoLocCountry"], r["geoLocAdmin1"], r["geoLocAdmin2"])
        if k in st_rows:
            same = all(str(r[c]) == str(st_rows[k][c]) for c in cmp_cols)
            n_same += same
            n_other += (not same)
    F["zone_map_rows_with_key_in_stored_map_same_values"] = int(n_same)
    F["zone_map_rows_with_key_in_stored_map_other_values"] = int(n_other)
    old_counts = {k: int(n) for k, n in zip(zip(zm_stored.geoLocCountry, zm_stored.geoLocAdmin1, zm_stored.geoLocAdmin2), zm_stored.n_genomes)}
    new_counts = {k: int(n) for k, n in zip(zip(zmap.geoLocCountry, zmap.geoLocAdmin1, zmap.geoLocAdmin2), zmap.n_genomes)}
    F["zone_map_counts_changed_vs_stored_map"] = [dict(geoLocCountry=k[0], geoLocAdmin1=k[1], geoLocAdmin2=k[2], geoLocAdmin2_repr=repr(k[2]),
                                                       zone=(st_rows[k]["zone"] if k in st_rows else ""),
                                                       n_stored_map=old_counts.get(k, 0), n_this_snapshot=new_counts.get(k, 0))
                                                  for k in sorted(set(old_counts) | set(new_counts)) if old_counts.get(k, 0) != new_counts.get(k, 0)]
    F["zone_map_rows"] = len(zmap)

    # ---- sample identifiers: exact matching vs matching without case and punctuation; identifiers occurring more than twice
    sid = o.specimenCollectorSampleId
    exact_groups = {k: sorted(g.index) for k, g in o[sid != ""].groupby("specimenCollectorSampleId") if len(g) > 1}
    norm_groups = {k: sorted(g.index) for k, g in o[o.sid_norm != ""].groupby("sid_norm") if len(g) > 1}
    F["dup_n_identifiers_exact"] = len(exact_groups)
    F["dup_n_identifiers_normalised"] = len(norm_groups)
    F["dup_same_groups_exact_and_normalised"] = bool(sorted(map(tuple, exact_groups.values())) == sorted(map(tuple, norm_groups.values())))
    F["dup_max_records_per_identifier"] = max([len(v) for v in exact_groups.values()], default=0)
    F["dup_records_per_identifier"] = {int(k): int(v) for k, v in sorted(Counter(len(v) for v in exact_groups.values()).items())}
    F["dup_n_genomes"] = int(o.duplicate_sample_id.sum())
    F["dup_n_pairs"] = len(res["dups"])
    F["sample_id_missing_in_set"] = sorted(o.index[sid == ""])

    # ---- write the output tables
    M.to_csv(OUT / f"pp_snapshot_{TAG}_metadata.csv", index=False)
    sl.write_fasta_gz(OUT / f"pp_snapshot_{TAG}_outbreak_aligned.fasta.gz", {k: seqs[k] for k in o.index})
    V.to_csv(OUT / f"snapshot_validation_{TAG}.csv")
    res["dups"].to_csv(OUT / f"duplicate_sample_pairs_{TAG}.csv", index=False)
    zmap.to_csv(OUT / f"zone_harmonisation_map_{TAG}.csv", index=False)
    UNM.to_csv(OUT / f"unmapped_health_zone_spellings_{TAG}.csv", index=False)
    res["byzone"].to_csv(OUT / f"genomes_by_zone_{TAG}.csv", index=False)
    res["zm"].to_csv(OUT / f"genomes_by_zone_month_{TAG}.csv")
    # working table for the later scripts
    o.drop(columns=["collection_date", "release_date", "submitted_date", "_dayd"]).to_csv(OUT / "work_outbreak_table_INTERNAL.csv")

    # ---- read back and check what was written
    Mb = sl.read_csv_str(OUT / f"pp_snapshot_{TAG}_metadata.csv")
    assert Mb.shape == (len(full), len(full.columns) + 16 + 2), Mb.shape
    assert list(Mb.columns) == sl.DERIVED16[:2] + list(full.columns) + sl.DERIVED16[2:] + sl.NEW_COLS
    F["metadata_database_fields_unchanged_cells_different"] = int((Mb[list(full.columns)].values != full.values).sum())
    F["metadata_rows"] = len(Mb)
    F["metadata_columns"] = len(Mb.columns)
    F["metadata_outbreak2026_true"] = int((Mb.outbreak2026 == "True").sum())
    F["metadata_latest_nonrevoked_true"] = int((Mb.latest_nonrevoked == "True").sum())
    F["metadata_called_frac_filled"] = int((Mb.called_frac != "").sum())
    F["metadata_n_acgt_filled"] = int((Mb.n_acgt_1_18900 != "").sum())
    cf_chk = Mb[Mb.called_frac != ""]
    F["called_frac_equals_n_acgt_1_18900_over_18900"] = bool(np.all(cf_chk.called_frac.astype(float).values == cf_chk.n_acgt_1_18900.astype(int).values / 18900))
    fb, _ = sl.read_fasta(OUT / f"pp_snapshot_{TAG}_outbreak_aligned.fasta.gz")
    F["fasta_records"] = len(fb)
    F["fasta_headers_sorted"] = bool([h for h, _ in fb] == sorted(h for h, _ in fb))
    F["fasta_headers_equal_outbreak_set"] = bool(set(h for h, _ in fb) == set_def)
    F["fasta_sequences_equal_raw_download"] = bool(all(seqs[h] == s for h, s in fb))
    F["fasta_lengths"] = dict(Counter(len(s) for _, s in fb))
    Vb = sl.read_csv_str(OUT / f"snapshot_validation_{TAG}.csv")
    F["validation_rows"] = len(Vb)
    F["validation_columns"] = len(Vb.columns)
    F["status_vs_baseline"] = Vb.status_vs_baseline.replace("", "(no baseline given)").value_counts().to_dict()
    F["metadata_changed_vs_baseline"] = Vb.metadata_changed_vs_baseline.replace("", "(not in baseline)").value_counts().to_dict()
    F["metadata_fields_changed_vs_baseline"] = Vb.metadata_fields_changed_vs_baseline.replace("", "(none or not in baseline)").value_counts().to_dict()
    dcf = Vb.called_frac.astype(float) - Vb.completeness_pathoplexus.astype(float)
    F["called_frac_minus_completeness_field"] = dict(min=float(dcf.min()), max=float(dcf.max()), n=int(dcf.notna().sum()))
    F["sha256"] = {p.name: sl.sha256_file(p) for p in sorted(OUT.glob(f"*{TAG}*")) if p.is_file()}
    json.dump(F, open(OUT / "build_facts.json", "w"), indent=1)
    if verbose:
        hide = ("aligned_alphabet", "unaligned_alphabet", "sha256", "dates_not_read_by_iso8601_parser", "fields", "unaligned_file",
                "stored_map_rows_without_genome_in_this_set")
        print(json.dumps({k: v for k, v in F.items() if k not in hide}, indent=1)[:6000])
        print({k: (v["n_filled"], v["n_not_read_by_iso8601_parser"]) for k, v in unread.items()})
    return F
