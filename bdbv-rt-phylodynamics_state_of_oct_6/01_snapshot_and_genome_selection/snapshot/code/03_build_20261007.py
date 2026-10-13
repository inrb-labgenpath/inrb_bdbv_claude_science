"""Step 2b - the derivations for the download of 7 October 2026, by the code that passed the control (build_snapshot.run).
Input:  raw/pathoplexus/ (verified by 01_freeze.py); the zone map of 2 October; metadata and aligned genomes of
        2 October (baseline of the comparison columns of the validation table only).
Output: out/pp_snapshot_20261007_metadata.csv, out/pp_snapshot_20261007_outbreak_aligned.fasta.gz,
        out/snapshot_validation_20261007.csv, out/duplicate_sample_pairs_20261007.csv,
        out/zone_harmonisation_map_20261007.csv, out/unmapped_health_zone_spellings_20261007.csv,
        out/genomes_by_zone_20261007.csv, out/genomes_by_zone_month_20261007.csv, out/build_facts.json
Run from snapshot_20261007/ :  python code/03_build_20261007.py     (after 01_freeze.py and 02_control_on_20261002.py)
"""
import json
import pathlib
import sys

sys.path.insert(0, "code")
import build_snapshot as bs  # noqa: E402

CFG = json.load(open("code/inputs.json"))
P = CFG["paths"]
control = json.load(open("out/control_facts.json"))
assert control["control_passed"], "the control on the download of 2 October did not pass: the new tables are not built"
freeze = json.load(open("out/freeze_facts.json"))
bs.run("20261007", pathlib.Path("raw/pathoplexus"), pathlib.Path("out"), freeze, P["zone_map_20261002"],
       baseline_meta=P["metadata_20261002"], baseline_fasta=P["fasta_20261002"], baseline_tag="20261002")
