#!/usr/bin/env python3
"""build_completeness_sets.py (round of October 2026) - runs the frozen build_alignments.py with the specifications of the sets
C95, C90, C80 and C00 and joins the four tables into one table of genome sets with the first reason for every absence.

The script changes no rule: it calls pipeline_20260924/build_alignments.py once per specification, with the command-line
options of the round of 24 Sep 2026 (--completeness-precision exact, --duplicate-pairs, --data-use-exclude, --batch-table,
--validation, --tail-days 42; --exclusion-list if given), and reads the tables that the pipeline writes.

  run 'C95', 'C90', 'C80', 'C00'  one specification each (alignment_set_spec_<set>_20261002.json); the pipeline writes
                                  <set>_<suffix>.fasta and alignment_sets_<suffix>.csv with reason_not_in_<set>
  run 'all'                       the four sets in one specification (alignment_set_spec_completeness_all_20261002.json);
                                  used here as a check: its alignments must be the same files as those of the single runs

Output in --outdir:
  <set>/...                        files of the pipeline, one directory per run
  genome_sets_by_called_fraction_<suffix>.csv   one row per screened genome: descriptive columns of the pipeline, in_C95 ...
                                  in_C00, reason_not_in_C95 ... reason_not_in_C00, in_any_set, first_reason_in_no_set
  genome_sets_summary_<suffix>.csv  per set: number of genomes, first and last collection date, SHA-256 of the alignment,
                                  and whether the alignment of the run 'all' is the same file
The accession versions in these tables are public accession numbers of the database; the tables hold no laboratory sample
identifier and no name of a submitter.

usage (from anywhere):
  python build_completeness_sets.py --pipeline TREE/pipeline_20260924 --specs TREE/set_specs --metadata META.csv
         --alignment ALN.fasta.gz --screen public_quality_screen_<sfx>.csv --batch-table screen_batches_<sfx>.csv
         --duplicate-pairs PAIRS.csv --data-use-exclude ACC1,ACC2 --validation VALIDATION.csv [--exclusion-list LIST.csv]
         --suffix <sfx> --outdir OUT
"""
import argparse
import hashlib
import os
import subprocess
import sys

import pandas as pd

SETS = ["C95", "C90", "C80", "C00"]
SPEC_SUFFIX = "20261002"


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--pipeline", required=True, help="directory with the frozen build_alignments.py")
    ap.add_argument("--specs", required=True, help="directory with the set specifications")
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--screen", required=True)
    ap.add_argument("--batch-table", required=True)
    ap.add_argument("--duplicate-pairs", required=True)
    ap.add_argument("--data-use-exclude", required=True, help="comma-separated accession(Version)s or a file; '' = none")
    ap.add_argument("--validation", required=True)
    ap.add_argument("--exclusion-list", default="")
    ap.add_argument("--suffix", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--python", default=sys.executable)
    a = ap.parse_args()
    ab = os.path.abspath
    env = dict(os.environ, OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1")
    runs = {s: os.path.join(a.specs, f"alignment_set_spec_{s}_{SPEC_SUFFIX}.json") for s in SETS}
    runs["all"] = os.path.join(a.specs, f"alignment_set_spec_completeness_all_{SPEC_SUFFIX}.json")
    for name, spec in runs.items():
        od = os.path.join(a.outdir, name)
        os.makedirs(od, exist_ok=True)
        cmd = [a.python, ab(os.path.join(a.pipeline, "build_alignments.py")), "--metadata", ab(a.metadata), "--alignment", ab(a.alignment),
               "--screen", ab(a.screen), "--out-suffix", a.suffix, "--outdir", ".", "--set-spec", ab(spec),
               "--completeness-precision", "exact", "--duplicate-pairs", ab(a.duplicate_pairs),
               "--batch-table", ab(a.batch_table), "--validation", ab(a.validation), "--tail-days", "42"]
        if a.data_use_exclude:
            cmd += ["--data-use-exclude", ab(a.data_use_exclude) if os.path.exists(a.data_use_exclude) else a.data_use_exclude]
        if a.exclusion_list:
            cmd += ["--exclusion-list", ab(a.exclusion_list)]
        with open(os.path.join(od, "build_alignments.log"), "w") as log:
            log.write(" ".join(cmd) + "\n")
            log.flush()
            rc = subprocess.run(cmd, cwd=od, stdout=log, stderr=subprocess.STDOUT, env=env).returncode
        assert rc == 0, f"build_alignments.py failed for the specification {spec} (see {od}/build_alignments.log)"
    # ---- one table ----
    tabs = {n: pd.read_csv(os.path.join(a.outdir, n, f"alignment_sets_{a.suffix}.csv"), dtype=str, keep_default_na=False) for n in runs}
    base_cols = [c for c in tabs["C90"].columns if not c.startswith("in_C") and not c.startswith("reason_not_in_")]
    T = tabs["C90"][base_cols].copy()
    for s in SETS:
        t = tabs[s]
        assert list(t["accessionVersion"]) == list(T["accessionVersion"]), f"row order differs in the table of the run {s}"
        for c in base_cols:
            assert (t[c].values == T[c].values).all(), f"descriptive column {c} differs between the runs C90 and {s}"
        T[f"in_{s}"] = t[f"in_{s}"].values
    for s in SETS:
        T[f"reason_not_in_{s}"] = tabs[s][f"reason_not_in_{s}"].values
        assert ((T[f"in_{s}"] == "True") == (T[f"reason_not_in_{s}"] == "")).all(), f"set {s}: a genome of the set has a reason, or an absent genome has none"
        assert (tabs["all"][f"in_{s}"].values == T[f"in_{s}"].values).all(), f"set {s}: membership differs between the single run and the run 'all'"
    any_set = (T[[f"in_{s}" for s in SETS]] == "True").any(axis=1)
    T["in_any_set"] = any_set.map({True: "True", False: "False"})
    # a genome that is in no set is absent from every set; the reason given is that for the set without threshold (C00)
    T["first_reason_in_no_set"] = [r if not k else "" for r, k in zip(T["reason_not_in_C00"], any_set)]
    T.to_csv(os.path.join(a.outdir, f"genome_sets_by_called_fraction_{a.suffix}.csv"), index=False)
    rows = []
    for s in SETS:
        f1 = os.path.join(a.outdir, s, f"{s}_{a.suffix}.fasta")
        f2 = os.path.join(a.outdir, "all", f"{s}_{a.suffix}.fasta")
        d = T.loc[T[f"in_{s}"] == "True", "collection_date"]
        rows.append(dict(set=s, n_genomes=int((T[f"in_{s}"] == "True").sum()), first_collection_date=d.min(), last_collection_date=d.max(),
                         alignment=os.path.basename(f1), sha256=sha256(f1), same_file_in_run_all=sha256(f1) == sha256(f2)))
    S = pd.DataFrame(rows)
    S.to_csv(os.path.join(a.outdir, f"genome_sets_summary_{a.suffix}.csv"), index=False)
    print(S[["set", "n_genomes", "first_collection_date", "last_collection_date", "same_file_in_run_all"]].to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
