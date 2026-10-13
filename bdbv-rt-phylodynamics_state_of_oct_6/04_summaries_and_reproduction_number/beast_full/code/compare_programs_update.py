"""compare_programs_update.py - written for this repository; it is not one of the programs that ran in the analysis.

The two tables of the comparison of the two programs (results/program_comparison_values_20261007.csv and
results/program_comparison_trajectories_20261007.csv) were written by compare_programs_v3.py, which stands beside
this file. In the analysis its two input files were written by a module that is not in the repository. This program
writes the same two input files from the tables of the repository and runs compare_programs_v3.py on them:

  periods.csv   six rows, P1 to P6, from settings/periods_of_the_update_20261007.csv (for P3 the end 'END', for P6
                the end that the table gives for the four sets), each with the status below;
  pairs.json    the three pairs below (A, B and B_estimated, each with the summary of the analysis): for each side
                the analysis, its table of estimates and its table of intervals in results/, and the labels under
                which its table of estimates names the quantities.
The words of the pairs (pair_in_the_form, differs_within_the_pair_in, what_differs_from_the_primary_analysis_of_Delphy,
the labels of the quantities), the note of the pair B_estimated and the statuses of the periods are those of the input
files that the module wrote. Two cells of the pair B_estimated (differs_within_the_pair_in, note_on_the_pair) stand
here in the shortened wording of the tables of results/.

usage, from the top folder of the repository:

    python 04_summaries_and_reproduction_number/beast_full/code/compare_programs_update.py OUTDIR

It writes OUTDIR/periods.csv, OUTDIR/pairs.json, copies of the tables read (OUTDIR/comparison/inputs/) and the two
tables, and prints whether each table is byte for byte the table of results/. The tables of results/ are those of
the analysis in the form of this repository: the cells that named a version of a table are empty, and cells of text
columns are reworded (README, 'What differs from the code as it ran').
exit status: 0 if both tables are those of results/, 3 if one is another, 2 if a step failed.
"""
import csv
import json
import os
import re
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(os.path.dirname(os.path.dirname(HERE)))
RESULTS = os.path.join(TOP, "results")
PERIODS_TABLE = "periods_of_the_update_20261007.csv"
ENDS_TABLE = "ends_by_set_20261007.csv"
VALUES = "program_comparison_values_20261007.csv"
TRAJECTORIES = "program_comparison_trajectories_20261007.csv"
DAY = re.compile(r"^\d{4}-\d\d-\d\d$")
ROWS_OF_THE_PERIODS = (("P1", "P1"), ("P2", "P2"), ("P3", "P3 (own end)"), ("P4", "P4"), ("P5", "P5"), ("P6", "P6 (own end)"))
STATUS = {
    "P1": "planned; carried over unchanged from the round of 2 October (plan of the update, section 1)",
    "P2": "planned; carried over unchanged from the round of 2 October (plan of the update, section 1)",
    "P3": "planned; carried over unchanged from the round of 2 October (plan of the update, section 1)",
    "P4": "planned; carried over unchanged from the round of 2 October (plan of the update, section 1)",
    "P5": "exploratory; carried over unchanged from the round of 2 October (plan of the update, section 1: 'P5, 15 June - 1 August (exploratory)')",
    "P6": "defined on 2 October 2026 as '1 August - 1 September', before the snapshot of that round was seen and with knowledge of the earlier estimates; in this update 1 September does not lie beyond the end of any of the four sets, and the period ends on 1 September (plan of the update, section 4 and its entries of section 8)"
}
PAIRS = [
    {
        "pair": "A",
        "pair_in_the_form": "A",
        "model": "skygrid",
        "summary": "the summary of the analysis",
        "differs_within_the_pair_in": "priors that each program sets for the other parameters (evolutionary rate: CTMC scale prior in BEAST X, uniform prior in Delphy)",
        "note_on_the_pair": "",
        "Delphy": {
            "analysis": "delphy_fixed_C00",
            "estimates": "estimates_delphy_fixed_C00_20261007.csv",
            "intervals": "intervals_delphy_fixed_C00_20261007.csv",
            "quantity_labels": {},
            "what_differs_from_the_primary_analysis_of_Delphy": "reference (primary analysis)"
        },
        "BEAST X": {
            "analysis": "beast_A_C00",
            "estimates": "estimates_beast_update_20261007.csv",
            "intervals": "intervals_beast_update_20261007.csv",
            "quantity_labels": {
                "reproduction number, P1": "reproduction number P1",
                "reproduction number, P2": "reproduction number P2",
                "reproduction number, P3": "reproduction number P3 (own end)",
                "reproduction number, P4": "reproduction number P4",
                "reproduction number, P5": "reproduction number P5",
                "reproduction number, P6": "reproduction number P6 (own end)"
            },
            "what_differs_from_the_primary_analysis_of_Delphy": "the program"
        }
    },
    {
        "pair": "B",
        "pair_in_the_form": "B",
        "model": "skygrid",
        "summary": "the summary of the analysis",
        "differs_within_the_pair_in": "smoothing: fixed in Delphy, estimated in BEAST X; priors of the other parameters as in A",
        "note_on_the_pair": "",
        "Delphy": {
            "analysis": "delphy_fixed_C00",
            "estimates": "estimates_delphy_fixed_C00_20261007.csv",
            "intervals": "intervals_delphy_fixed_C00_20261007.csv",
            "quantity_labels": {},
            "what_differs_from_the_primary_analysis_of_Delphy": "reference (primary analysis)"
        },
        "BEAST X": {
            "analysis": "beast_B_C00",
            "estimates": "estimates_beast_update_20261007.csv",
            "intervals": "intervals_beast_update_20261007.csv",
            "quantity_labels": {
                "reproduction number, P1": "reproduction number P1",
                "reproduction number, P2": "reproduction number P2",
                "reproduction number, P3": "reproduction number P3 (own end)",
                "reproduction number, P4": "reproduction number P4",
                "reproduction number, P5": "reproduction number P5",
                "reproduction number, P6": "reproduction number P6 (own end)"
            },
            "what_differs_from_the_primary_analysis_of_Delphy": "the program; smoothing estimated"
        }
    },
    {
        "pair": "B_estimated",
        "pair_in_the_form": "not in the stored form; added in this round",
        "model": "skygrid",
        "summary": "the summary of the analysis",
        "differs_within_the_pair_in": "the program; the smoothing is estimated in both analyses (settings of the two tables of estimates, whole cells: Delphy 'Skygrid, smoothing estimated (hyperprior at the program defaults)'; BEAST X 'B: Skygrid, precision estimated (Gamma(0.001, 1000) prior)'); prior of the smoothing in Delphy: not stated beyond that cell; priors of the other parameters in the two programs: not stated (the input file of the analysis of BEAST X and the command lines of the analysis of Delphy were not available to the comparison)",
        "note_on_the_pair": "Delphy with the smoothing estimated beside BEAST X with the smoothing estimated (setting B). The form of 26 September does not hold this pair; it was added in the round of 2 October and is run for the update beside pair B. The cell pair_in_the_form of this pair is that of the pairs file of 2 October: 'this round' in it is the round of 2 October. The two cells 'setting' that the cell differs_within_the_pair_in quotes are, character for character, those of the two tables of estimates of the update.",
        "Delphy": {
            "analysis": "delphy_estimated_C00",
            "estimates": "estimates_delphy_estimated_C00_20261007.csv",
            "intervals": "intervals_delphy_estimated_C00_20261007.csv",
            "quantity_labels": {},
            "what_differs_from_the_primary_analysis_of_Delphy": "smoothing estimated (setting of the table of estimates of the analysis, whole cell: 'Skygrid, smoothing estimated (hyperprior at the program defaults)')"
        },
        "BEAST X": {
            "analysis": "beast_B_C00",
            "estimates": "estimates_beast_update_20261007.csv",
            "intervals": "intervals_beast_update_20261007.csv",
            "quantity_labels": {
                "reproduction number, P1": "reproduction number P1",
                "reproduction number, P2": "reproduction number P2",
                "reproduction number, P3": "reproduction number P3 (own end)",
                "reproduction number, P4": "reproduction number P4",
                "reproduction number, P5": "reproduction number P5",
                "reproduction number, P6": "reproduction number P6 (own end)"
            },
            "what_differs_from_the_primary_analysis_of_Delphy": "the program; smoothing estimated"
        }
    }
]


def table(path):
    with open(path, encoding="utf-8", newline="") as fh:
        return list(csv.DictReader(fh))


def periods_file(path):
    """the periods in the form that compare_programs_v3.py reads (period, start, end, status)"""
    periods = {r["period"]: r for r in table(os.path.join(TOP, "settings", PERIODS_TABLE))}
    ends = table(os.path.join(TOP, "settings", ENDS_TABLE))
    if len(ends) != 4 or not all(DAY.match(r["own_end_mid_point_of_the_most_recent_interval"])
                                 and r["own_end_mid_point_of_the_most_recent_interval"] >= "2026-09-01" for r in ends):
        raise SystemExit("an own end of a set lies before 1 September, or the table of the ends has not four rows")
    rows = []
    for k, key in ROWS_OF_THE_PERIODS:
        r = periods[key]
        start, end = r["start"], r["end"]
        if k == "P3":
            if not end.startswith("END of the set and grid of the analysis"):
                raise SystemExit("P3: the end is not END: " + end)
            end = "END"
        elif k == "P6":
            if not end.startswith("2026-09-01 if 1 September 2026 does not lie beyond END") \
                    or r["in_the_four_sets_of_this_round"] != "2026-08-01 to 2026-09-01":
                raise SystemExit("P6: the table of the periods does not give 1 August to 1 September for the four sets")
            end = "2026-09-01"
        elif r["rule"] != "fixed dates":
            raise SystemExit(k + ": the rule is not 'fixed dates'")
        if not DAY.match(start) or not (end == "END" or DAY.match(end)):
            raise SystemExit("%s: start %s, end %s" % (k, start, end))
        rows.append(dict(period=k, start=start, end=end, status=STATUS[k], source="%s, row '%s'" % (PERIODS_TABLE, key)))
    with open(path, "w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["period", "start", "end", "status", "source"], lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def one(rows, column):
    values = sorted(set(r[column] for r in rows))
    if len(values) != 1:
        raise SystemExit("the column %s holds %d values in the rows of the analysis" % (column, len(values)))
    return values[0]


def pairs_file(outdir):
    """pairs.json for the three pairs; the tables are copied to OUTDIR/comparison/inputs/"""
    os.makedirs(os.path.join(outdir, "comparison", "inputs"), exist_ok=True)
    out = []
    for p in PAIRS:
        q = {k: p[k] for k in ("pair", "pair_in_the_form", "model", "summary", "differs_within_the_pair_in", "note_on_the_pair")}
        for program in ("Delphy", "BEAST X"):
            s0 = p[program]
            s = dict(analysis=s0["analysis"], analysis_in_the_table=s0["analysis"])
            for kind in ("estimates", "intervals"):
                rel = "comparison/inputs/" + s0[kind]
                shutil.copyfile(os.path.join(RESULTS, s0[kind]), os.path.join(outdir, rel))
                s[kind], s[kind + "_name"], s[kind + "_version"] = rel, s0[kind], ""
            estimates = [r for r in table(os.path.join(outdir, s["estimates"])) if r.get("analysis", s["analysis"]) == s["analysis"]]
            intervals = [r for r in table(os.path.join(outdir, s["intervals"])) if r.get("analysis", s["analysis"]) == s["analysis"]]
            if not estimates or not intervals:
                raise SystemExit("no row of the analysis %s in %s or %s" % (s["analysis"], s0["estimates"], s0["intervals"]))
            s["quantity_labels"] = dict(s0["quantity_labels"])
            s["genomes"], s["most_recent_collection_date"] = one(estimates, "genomes"), one(estimates, "most_recent_collection_date")
            if "setting" in estimates[0]:
                s["setting"] = one(estimates, "setting")
            s["n_intervals"] = str(len(intervals))
            s["what_differs_from_the_primary_analysis_of_Delphy"] = s0["what_differs_from_the_primary_analysis_of_Delphy"]
            q[program] = s
        out.append(q)
    with open(os.path.join(outdir, "pairs.json"), "w", encoding="utf-8") as fh:
        json.dump(dict(pairs=out, written_by="compare_programs_update.py"), fh, indent=1, ensure_ascii=False)


def main(argv):
    if len(argv) != 2:
        print("usage: python 04_summaries_and_reproduction_number/beast_full/code/compare_programs_update.py OUTDIR", file=sys.stderr)
        return 2
    outdir = os.path.abspath(argv[1])
    os.makedirs(outdir, exist_ok=True)
    periods_file(os.path.join(outdir, "periods.csv"))
    pairs_file(outdir)
    step = subprocess.run([sys.executable, os.path.join(HERE, "compare_programs_v3.py"), "--pairs", "pairs.json", "--periods", "periods.csv",
                           "--out-values", VALUES, "--out-trajectories", TRAJECTORIES], cwd=outdir)
    if step.returncode != 0:
        print("compare_programs_v3.py ended with status", step.returncode, file=sys.stderr)
        return 2
    same = True
    for name in (VALUES, TRAJECTORIES):
        ok = open(os.path.join(outdir, name), "rb").read() == open(os.path.join(RESULTS, name), "rb").read()
        same = same and ok
        print("%s: the table of results/: %s" % (name, "the same" if ok else "NOT the same"))
    return 0 if same else 3


if __name__ == "__main__":
    sys.exit(main(sys.argv))
