"""regenerate_case_tables.py - written for this repository; it is not one of the programs that ran in the analysis.

The repository holds no table of case counts. This program regenerates them from the public source, the
transcription of the situation reports of the Institut National de Sante Publique (INSP) in the repository
INRB-UMIE/BDBV2026-Data at the commit of the analysis:

    git clone https://github.com/INRB-UMIE/BDBV2026-Data
    git -C BDBV2026-Data checkout 8eb57154cf967e7dfdb7e690230f61f2a298d2b0

usage, from the top folder of this repository:

    python 06_case_counts_and_response_events/regenerate_case_tables.py BDBV2026-Data

What it does:
  1. runs build_case_series_20261002.py (the program of the analysis, unchanged) on the folder given, with the
     arguments of the analysis; its tables go to results/case_counts_regenerated/;
  2. copies the weekly table of confirmed cases by province to results/;
  3. runs 07_figures/add_case_counts_to_figure_1_values.py, which writes the file of values of Figure 1
     (results/fig_ne_R_cases_20261006_values.csv), from which draw_figure_1.py draws.
It compares the SHA-256 of the weekly table and of the file of values with those of the files of the analysis and
prints whether they are the same. The files it writes are not listed in MANIFEST.csv.

The terms of the source for the situation reports ask for attribution to the INSP, citation of the number and date
of the report, and confirmation of the terms of distribution with the INSP before the counts are published elsewhere
(see README_case_series_20261002.md).

exit status: 0 if both files are those of the analysis, 3 if one is another, 2 if a step failed.
"""
import hashlib
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
SUFFIX = "20261002"
WEEKLY = f"weekly_confirmed_cases_by_province_{SUFFIX}.csv"
SHA256_OF_THE_WEEKLY_TABLE_OF_THE_ANALYSIS = "7dad27e7551c8a6cf4f324a9d75c8356fdff6a4f18568cad1cf6452ccf4582ff"
NEEDED = ("data/insp_sitrep/processed/insp_sitrep__cumulative_confirmed_cases__daily.csv",
          "data/insp_sitrep/processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv")


def main(argv):
    if len(argv) != 2:
        print("usage: python 06_case_counts_and_response_events/regenerate_case_tables.py SOURCE_FOLDER", file=sys.stderr)
        return 2
    source = argv[1]
    missing = [n for n in NEEDED if not os.path.exists(os.path.join(source, n))]
    if missing:
        print("not found in %s: %s" % (source, missing), file=sys.stderr)
        return 2
    outdir = os.path.join(TOP, "results", "case_counts_regenerated")
    os.makedirs(outdir, exist_ok=True)
    step = subprocess.run([sys.executable, os.path.join(HERE, "build_case_series_20261002.py"), "--indir", source,
                           "--outdir", outdir, "--suffix", SUFFIX,
                           "--zones-file", os.path.join(HERE, "zones_of_analysed_set_20261002.txt")],
                          stdout=subprocess.DEVNULL)
    if step.returncode != 0:
        print("build_case_series_20261002.py ended with status", step.returncode, file=sys.stderr)
        return 2
    weekly = os.path.join(TOP, "results", WEEKLY)
    shutil.copyfile(os.path.join(outdir, WEEKLY), weekly)
    got = hashlib.sha256(open(weekly, "rb").read()).hexdigest()
    same = got == SHA256_OF_THE_WEEKLY_TABLE_OF_THE_ANALYSIS
    print("written:", os.path.relpath(outdir, TOP) + os.sep, "and", os.path.relpath(weekly, TOP))
    print("SHA-256 of the weekly table:", got)
    print("the weekly table of the analysis:",
          "the same" if same else "NOT the same (" + SHA256_OF_THE_WEEKLY_TABLE_OF_THE_ANALYSIS + ")")
    sys.stdout.flush()
    step = subprocess.run([sys.executable, os.path.join(TOP, "07_figures", "add_case_counts_to_figure_1_values.py"), weekly])
    if step.returncode not in (0, 3):
        return 2
    return 0 if same and step.returncode == 0 else 3


if __name__ == "__main__":
    sys.exit(main(sys.argv))
