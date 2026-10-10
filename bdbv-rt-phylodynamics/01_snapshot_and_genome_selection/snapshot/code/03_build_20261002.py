"""Step 2b - apply the validated derivations to the download of 2 October 2026.
Input:  raw/pathoplexus/ (verified by 01_freeze.py), zone map of 24 September,
        metadata and aligned genomes of 24 September (baseline of the comparison columns only).
Output: out/pp_snapshot_20261002_metadata.csv, out/pp_snapshot_20261002_outbreak_aligned.fasta.gz,
        out/snapshot_validation_20261002.csv, out/duplicate_sample_pairs_20261002.csv,
        out/zone_harmonisation_map_20261002.csv, out/genomes_by_zone_20261002.csv,
        out/genomes_by_zone_month_20261002.csv, out/build_facts.json
Run from snapshot_20261002/ :  python code/03_build_20261002.py
"""
import json
import pathlib
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, "code")
import snapshot_lib as sl  # noqa: E402

CFG = json.load(open("code/inputs.json"))
P = CFG["paths"]
RAW = pathlib.Path("raw/pathoplexus")
OUT = pathlib.Path("out")
TAG = "20261002"
freeze = json.load(open(OUT / "freeze_facts.json"))
assert freeze["all_sha256_match"] and freeze["all_sizes_match"]
T_NEW = freeze["details_allversions_query_start_utc"]
DV_NEW = freeze["data_version"]
F = {"snapshotDatetimeUtc": T_NEW, "lapisDataVersion": DV_NEW}

full = sl.read_csv_str(RAW / "details_allversions.csv")
filt = sl.read_csv_str(RAW / "details_latest_nonrevoked.csv")
al, al_info = sl.read_fasta(RAW / "aligned_latest_nonrevoked.fasta")
un, un_info = sl.read_fasta(RAW / "unaligned_latest_nonrevoked.fasta")
F["n_records_all_versions"] = len(full)
F["n_records_latest_nonrevoked"] = len(filt)
F["n_fields"] = len(full.columns)
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
F["latest_nonrevoked_records_without_aligned_sequence"] = sorted(set(filt.accessionVersion) - set(seqs))
F["aligned_sequences_without_record"] = sorted(set(seqs) - set(filt.accessionVersion))
F["aligned_sequences_of_other_length"] = sorted(h for h, s in al if len(s) != sl.L_ALN)
F["aligned_sequences_empty"] = sorted(h for h, s in al if len(s) == 0)
F["aligned_sequences_without_any_called_base"] = sorted(h for h, s in al if sl.count_acgt(s, sl.L_ALN) == 0)
F["latest_nonrevoked_records_without_unaligned_sequence"] = sorted(set(filt.accessionVersion) - set(useqs))
F["unaligned_sequences_without_record"] = sorted(set(useqs) - set(filt.accessionVersion))
alpha = Counter()
for _, s in al:
    alpha.update(s)
F["aligned_alphabet_854"] = dict(sorted(alpha.items()))
F["aligned_lowercase_characters"] = int(sum(n for c, n in alpha.items() if c.islower()))
ualpha = Counter()
for _, s in un:
    ualpha.update(s)
F["unaligned_alphabet_854"] = dict(sorted(ualpha.items()))

# ---- the unaligned file: what it holds (nothing below uses it)
len_field = filt.set_index("accessionVersion")["length"].astype(int)
ins_field = filt.set_index("accessionVersion")["totalInsertedNucs"].astype(int)
ob_field = set(filt.accessionVersion[filt.outbreak == "Bdbv-2026"])
core = lambda x: x.strip("Nn")
U = pd.DataFrame({"len_un": {k: len(v) for k, v in useqs.items()}})
U["length_field"] = len_field
U["acgt_un"] = pd.Series({k: sl.count_acgt(v, len(v)) for k, v in useqs.items()})
U["acgt_al"] = pd.Series({k: sl.count_acgt(v, sl.L_ALN) for k, v in seqs.items()})
U["equal_to_aligned_without_gaps_and_outer_N"] = pd.Series({k: core(useqs[k]) == core(seqs[k].replace("-", "")) for k in useqs})
U["totalInsertedNucs"] = ins_field
U["field_outbreak_Bdbv2026"] = U.index.isin(ob_field)
uf = {}
for nm, usub in (("all 854 latest non-revoked records", U), ("records with field outbreak = Bdbv-2026", U[U.field_outbreak_Bdbv2026]),
                ("other records", U[~U.field_outbreak_Bdbv2026])):
    extra = usub.acgt_un - usub.acgt_al
    uf[nm] = dict(n=len(usub), length_min=int(usub.len_un.min()), length_median=float(usub.len_un.median()), length_max=int(usub.len_un.max()),
                  n_length_equals_field_length=int((usub.len_un == usub.length_field).sum()),
                  n_equal_to_aligned_without_gaps_and_outer_N=int(usub.equal_to_aligned_without_gaps_and_outer_N.sum()),
                  n_with_more_called_bases_than_aligned=int((extra > 0).sum()), n_with_fewer_called_bases_than_aligned=int((extra < 0).sum()),
                  extra_called_bases_distribution={int(k): int(v) for k, v in extra.value_counts().sort_index().items()},
                  n_with_insertions_by_field=int((usub.totalInsertedNucs > 0).sum()),
                  n_not_equal_and_no_insertion_by_field=int(((~usub.equal_to_aligned_without_gaps_and_outer_N) & (usub.totalInsertedNucs == 0)).sum()),
                  n_equal_but_insertion_by_field=int((usub.equal_to_aligned_without_gaps_and_outer_N & (usub.totalInsertedNucs > 0)).sum()),
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

# ---- derivations (baseline = snapshot of 24 September, comparison columns only)
om = sl.read_csv_str(P["old_metadata"])
old_fields = [c for c in om.columns if c not in sl.DERIVED16]
assert old_fields == list(full.columns), "database fields differ from those of 24 September"
old_seqs = dict(sl.read_fasta(P["old_fasta"])[0])
zm_stored, look = sl.load_zone_map(P["old_zone_map"])
res = sl.build(full, seqs, look, T_NEW, DV_NEW, baseline_full=om[old_fields], baseline_seqs=old_seqs, add_new_columns=True)
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

# ---- zone map: raw spellings that the map of 24 September does not hold
nm = o[o.flag_zone_not_in_map]
F["zone_keys_not_in_stored_map"] = [dict(geoLocCountry=c, geoLocAdmin1=a1, geoLocAdmin2=a2, geoLocAdmin2_repr=repr(a2), n_genomes=int(n),
                                         accessionVersions=";".join(sorted(g.index)),
                                         under_rule_of_earlier_code_NOT_APPLIED=sl.rule_of_20260924_code(c, a1, a2))
                                    for (c, a1, a2), n, g in ((k, len(g), g) for k, g in nm.groupby(["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2"]))]
F["n_genomes_zone_not_in_stored_map"] = int(o.flag_zone_not_in_map.sum())
zmap = res["zmap"].copy()
used = set(zip(zmap.geoLocCountry, zmap.geoLocAdmin1, zmap.geoLocAdmin2))
unused = zm_stored[[k not in used for k in zip(zm_stored.geoLocCountry, zm_stored.geoLocAdmin1, zm_stored.geoLocAdmin2)]].copy()
F["stored_map_rows"] = len(zm_stored)
F["stored_map_rows_without_genome_in_new_set"] = unused[["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "zone"]].to_dict("records")
zmap["n_genomes"] = zmap["n_genomes"].astype(int)
info_rule = [sl.rule_of_20260924_code(c, a1, a2)["zone"] for c, a1, a2 in zip(zmap.geoLocCountry, zmap.geoLocAdmin1, zmap.geoLocAdmin2)]
mapped = (zmap.flag_zone_not_in_map.astype(str) == "False").values
assert all(z == i for z, i, m in zip(zmap.zone, info_rule, mapped) if m), "rule of the earlier code and stored map disagree on a mapped key"
old_counts = {k: int(n) for k, n in zip(zip(zm_stored.geoLocCountry, zm_stored.geoLocAdmin1, zm_stored.geoLocAdmin2), zm_stored.n_genomes)}
new_counts = {k: int(n) for k, n in zip(zip(zmap.geoLocCountry, zmap.geoLocAdmin1, zmap.geoLocAdmin2), zmap.n_genomes)}
F["zone_map_counts_changed_vs_stored_map"] = [dict(geoLocCountry=k[0], geoLocAdmin1=k[1], geoLocAdmin2=k[2], geoLocAdmin2_repr=repr(k[2]),
                                                   n_20260924=old_counts.get(k, 0), n_20261002=new_counts.get(k, 0))
                                              for k in sorted(set(old_counts) | set(new_counts)) if old_counts.get(k, 0) != new_counts.get(k, 0)]
F["zone_map_rows_20261002"] = len(zmap)

# ---- sample identifiers: exact matching vs matching without case and punctuation; identifiers occurring more than twice
sid = o.specimenCollectorSampleId
exact_groups = {k: sorted(g.index) for k, g in o[sid != ""].groupby("specimenCollectorSampleId") if len(g) > 1}
norm_groups = {k: sorted(g.index) for k, g in o[o.sid_norm != ""].groupby("sid_norm") if len(g) > 1}
F["dup_n_identifiers_exact"] = len(exact_groups)
F["dup_n_identifiers_normalised"] = len(norm_groups)
F["dup_same_groups_exact_and_normalised"] = bool(sorted(map(tuple, exact_groups.values())) == sorted(map(tuple, norm_groups.values())))
F["dup_max_records_per_identifier"] = max([len(v) for v in exact_groups.values()], default=0)
F["dup_n_genomes"] = int(o.duplicate_sample_id.sum())
F["dup_n_pairs"] = len(res["dups"])

# ---- write the output tables
M.to_csv(OUT / f"pp_snapshot_{TAG}_metadata.csv", index=False)
sl.write_fasta_gz(OUT / f"pp_snapshot_{TAG}_outbreak_aligned.fasta.gz", {k: seqs[k] for k in o.index})
V.to_csv(OUT / f"snapshot_validation_{TAG}.csv")
res["dups"].to_csv(OUT / f"duplicate_sample_pairs_{TAG}.csv", index=False)
zmap.to_csv(OUT / f"zone_harmonisation_map_{TAG}.csv", index=False)
res["byzone"].to_csv(OUT / f"genomes_by_zone_{TAG}.csv", index=False)
res["zm"].to_csv(OUT / f"genomes_by_zone_month_{TAG}.csv")
# working table for the later scripts
o.drop(columns=["collection_date", "release_date", "submitted_date", "_dayd"]).to_csv(OUT / "work_outbreak_table_INTERNAL.csv")

# ---- read back and check what was written
Mb = sl.read_csv_str(OUT / f"pp_snapshot_{TAG}_metadata.csv")
assert Mb.shape == (len(full), 132 + 16 + 2), Mb.shape
assert list(Mb.columns) == sl.DERIVED16[:2] + list(full.columns) + sl.DERIVED16[2:] + sl.NEW_COLS
F["metadata_database_fields_unchanged_cells_different"] = int((Mb[list(full.columns)].values != full.values).sum())
F["metadata_rows"] = len(Mb); F["metadata_columns"] = len(Mb.columns)
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
F["validation_rows"] = len(Vb); F["validation_columns"] = len(Vb.columns)
F["status_vs_baseline"] = Vb.status_vs_baseline.value_counts().to_dict()
F["metadata_changed_vs_baseline"] = Vb.metadata_changed_vs_baseline.replace("", "(not in baseline)").value_counts().to_dict()
dcf = Vb.called_frac.astype(float) - Vb.completeness_pathoplexus.astype(float)
F["called_frac_minus_completeness_field"] = dict(min=float(dcf.min()), max=float(dcf.max()), n=int(dcf.notna().sum()))
F["sha256"] = {p.name: sl.sha256_file(p) for p in sorted(OUT.glob(f"*{TAG}*")) if p.is_file()}
json.dump(F, open(OUT / "build_facts.json", "w"), indent=1)
show = {k: v for k, v in F.items() if k not in ("aligned_alphabet_854", "unaligned_alphabet_854", "sha256", "dates_not_read_by_iso8601_parser")}
print(json.dumps(show, indent=1)[:7000])
print({k: (v["n_filled"], v["n_not_read_by_iso8601_parser"]) for k, v in unread.items()})
