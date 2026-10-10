#!/usr/bin/env python3
"""make_completeness_specs.py (round of October 2026) - writes the specifications of the genome sets C95, C90, C80 and C00
(analysis plan of 2 October 2026, section 2) in the format read by build_alignments.py --set-spec.

Every set takes the filters of the set 'primary' of the frozen specification of 24 Sep 2026
(pipeline_20260924/alignment_set_spec_20260924.json) and changes ONE key, min_completeness:
    C95  0.95      C90  0.9 (the rule of 24 Sep)      C80  0.8      C00  null (no threshold)
Nothing else is changed, and no rule of the pipeline is changed.  The key min_completeness is always written: a
specification WITHOUT the key means 0.9 in build_alignments.py (default of the code), not 'no threshold'.

Files written (--outdir):
  alignment_set_spec_C95_20261002.json ... alignment_set_spec_C00_20261002.json
        one set per file, named as primary set; the table alignment_sets_<suffix>.csv of a run then holds the column
        reason_not_in_<set> (first applicable reason for the absence of a genome from THAT set)
  alignment_set_spec_completeness_all_20261002.json
        the same four sets in one file (primary: C90), for a run that writes the four alignments and one table with the
        four membership columns; its column reason_not_in_C90 refers to C90 only

Top-level keys: 'neutral_wording' is kept (true) so that column names and reasons are worded as in the table of 24 Sep.
'external_data_set' of the frozen specification is NOT copied: it describes the data set of another analysis and
enters no rule of membership; without it the table holds the yes/no column on_external_exclusion_table.

usage: python round_20261002/make_completeness_specs.py --frozen pipeline_20260924/alignment_set_spec_20260924.json --outdir set_specs
"""
import argparse
import json
import os
import sys

SETS = [("C95", 0.95), ("C90", 0.9), ("C80", 0.8), ("C00", None)]
SUFFIX = "20261002"


def dump(spec):
    """one set per line, as the frozen specification is written"""
    lines = ["{", f' "primary": {json.dumps(spec["primary"])},', f' "neutral_wording": {json.dumps(spec["neutral_wording"])},', ' "sets": {']
    names = list(spec["sets"])
    for i, n in enumerate(names):
        lines.append(f'  {json.dumps(n)}: {json.dumps(spec["sets"][n])}' + ("," if i < len(names) - 1 else ""))
    lines += [" }", "}"]
    return "\n".join(lines) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--frozen", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    frozen = json.load(open(a.frozen))
    base = frozen["sets"][frozen["primary"]]
    assert base.get("min_completeness") == 0.9, "the primary set of the frozen specification is expected to have min_completeness 0.9"
    os.makedirs(a.outdir, exist_ok=True)
    sets = {}
    for name, thr in SETS:
        f = dict(base)
        f["min_completeness"] = thr
        assert list(f) == list(base) and all(f[k] == base[k] for k in base if k != "min_completeness")
        sets[name] = f
    for name, _ in SETS:
        spec = dict(primary=name, neutral_wording=bool(frozen.get("neutral_wording", False)), sets={name: sets[name]})
        txt = dump(spec)
        assert json.loads(txt) == spec
        open(os.path.join(a.outdir, f"alignment_set_spec_{name}_{SUFFIX}.json"), "w").write(txt)
    spec = dict(primary="C90", neutral_wording=bool(frozen.get("neutral_wording", False)), sets=sets)
    txt = dump(spec)
    assert json.loads(txt) == spec
    open(os.path.join(a.outdir, f"alignment_set_spec_completeness_all_{SUFFIX}.json"), "w").write(txt)
    print("filters of the frozen primary set:", json.dumps(base))
    for name, _ in SETS:
        print(name, json.dumps(sets[name]))
    return 0


if __name__ == "__main__":
    sys.exit(main())
