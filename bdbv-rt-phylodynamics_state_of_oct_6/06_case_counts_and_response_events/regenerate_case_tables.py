"""regenerate_case_tables.py - written for this repository; it is not one of the programs that ran in the analysis.

The repository holds no table of case counts. This program regenerates them from the public source, the
transcription of the situation reports of the Institut National de Sante Publique (INSP) in the repository
INRB-UMIE/BDBV2026-Data at the commit of the analysis (commit 71e06bda, state of 7 October 2026):

    git clone https://github.com/INRB-UMIE/BDBV2026-Data
    git -C BDBV2026-Data checkout 71e06bdaeda49c068ba499cab41359179cc4b9d0

usage, from the top folder of this repository:

    python 06_case_counts_and_response_events/regenerate_case_tables.py BDBV2026-Data

What it does:
  1. builds the case series as the program of the analysis did (case_series_update_20261007.py, function main, the
     parts 'the script of 2 October, checked by its SHA-256' and '2. the series from the new files'), by the same calls:
       sha256(...) of case_series_update_20261007.py   build_case_series_20261002.py is the script whose SHA-256 that
                                                        program holds; the constant is the SHA-256 of the script as it
                                                        stands in this repository, where it differs from the script
                                                        that ran in comments and documentation strings only
       read_zone_list(...) of build_case_series_20261002.py   the list of health zones
       build(...) of build_case_series_20261002.py     the series, from the folder given
       write_step1(...) of case_series_update_20261007.py     the nine tables of step 1, written as that program
                                                        writes them
     The tables go to results/case_counts_regenerated/. The cumulative counts by province that the text of the
     article prints stand in the last row of cumulative_confirmed_cases_by_province_20261007.csv;
  2. copies the weekly table of confirmed cases by province to results/;
  3. runs 07_figures/add_case_counts_to_figure_1_values.py, which writes the file of values of Figure 1
     (results/fig_ne_R_cases_20261007_values.csv), from which draw_figure_1.py draws.
It compares the SHA-256 of every table with that of the table of the analysis and prints whether they are the same.
07_figures/add_case_counts_to_figure_1_values.py does the same for the file of values; there the file of the analysis
is meant in the form in which this repository holds such files (the cells that named a version of a table empty,
the tables named as the members of results/ are named, three cells of the column note shortened), as the head of
that program says.

case_series_update_20261007.py is not run as a whole (README says how the analysis called it). Beside the source it
takes the fetch log of the download, a folder with tables against which it checks its output and the sentence of
the article on the case counts, none of which is in this repository.
Left out here are therefore: the check of the files against the fetch log and against the list of files of the
commit; the lists of the report dates, of the situation reports and of the falls of a cumulative count; the
comparison with those tables; and the table of the numbers that the text of the article prints, which that program
compares with that sentence.

The licence of the source repository is the MIT licence (its LICENSE.md). For the situation reports the source
asks for attribution to the INSP, citation of the number and date of the report, and confirmation of the terms of
distribution with the INSP before the counts are published elsewhere (its data/insp_sitrep/metadata.yaml).

exit status: 0 if every file is that of the analysis, 3 if one is another, 2 if a step failed.
"""
import hashlib
import os
import shutil
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
SUFFIX = "20261007"
WEEKLY = f"weekly_confirmed_cases_by_province_{SUFFIX}.csv"
ZONES = "zones_of_analysed_set_20261002.txt"
SHA256_OF_THE_TABLES_OF_THE_ANALYSIS = {
    "weekly_confirmed_cases_by_province_20261007.csv": "d6fa06fde049120851efaa2cf607ed78b094640bf19d7cbdc02c3725da513082",
    "cumulative_confirmed_cases_by_province_20261007.csv": "5a24664d16595c254f34b67bae618cf8a94804e915cf37533c7b62449cbaf3d9",
    "zone_province_map_20261007.csv": "358df2061c23c715dd458658567d51c4b08c33966f9b90b88e5183f59dacd70c",
    "mean_weekly_cases_by_zone_group_20261007.csv": "eda4095f9c524c9ff8bee514f65fb09e985cb9f276f88500bca0b0171e8c8975",
    "mean_weekly_cases_seven_zones_20261007.csv": "e8b47bd9061fd684bcdc67fea5af4af33e06731a02f550325566ac7671f71e83",
    "mean_weekly_cases_seven_zones_by_treatment_20261007.csv": "2339f1c0e74257c7c1cf00a613cfbc987e9533c4d9ec862960767866c1aa12d4",
    "weekly_confirmed_cases_seven_zones_20261007.csv": "025ae482c72086dee8609ef513b0f1ca84bd414964f823ab6ff05cb25f179c31",
    "case_series_crosscheck_20261007.csv": "4204a8d548962f2ec837d4204c59a09caa2cc1e04ae601c99f0a4d1883379245",
    "case_series_values_set_aside_20261007.csv": "dea1e4d619ab1aa450036d5a96379a83fd3a94ee0370087a835fdd7d9d8afccb",
}
NEEDED = ("data/insp_sitrep/processed/insp_sitrep__cumulative_confirmed_cases__daily.csv",
          "data/insp_sitrep/processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv",
          "data/shapefiles/DRC_Health_zones.dbf")


def main(argv):
    if len(argv) != 2:
        print("usage: python 06_case_counts_and_response_events/regenerate_case_tables.py SOURCE_FOLDER", file=sys.stderr)
        return 2
    source = argv[1]
    missing = [n for n in NEEDED if not os.path.exists(os.path.join(source, n))]
    if missing:
        print("not found in %s: %s" % (source, missing), file=sys.stderr)
        return 2
    sys.path.insert(0, HERE)
    import case_series_update_20261007 as update          # the program of the analysis; nothing runs on import
    if update.sha256(os.path.join(HERE, update.STORED_SCRIPT)) != update.STORED_SCRIPT_SHA256:
        print(update.STORED_SCRIPT + " is not the script whose SHA-256 case_series_update_20261007.py holds", file=sys.stderr)
        return 2
    import build_case_series_20261002 as stored           # the stored script
    outdir = os.path.join(TOP, "results", "case_counts_regenerated")
    try:
        zones = stored.read_zone_list(os.path.join(HERE, ZONES))
        tables = stored.build(source, zone_list=zones)
        files = update.write_step1(tables, outdir, SUFFIX)
    except SystemExit as stop:
        print("the build of the case series stopped:", stop, file=sys.stderr)
        return 2
    written = sorted(files.values())
    if written != sorted(SHA256_OF_THE_TABLES_OF_THE_ANALYSIS):
        print("the tables written are not the nine tables of step 1: %s" % written, file=sys.stderr)
        return 2
    same = True
    print("written to", os.path.relpath(outdir, TOP) + os.sep)
    for name in written:
        got = hashlib.sha256(open(os.path.join(outdir, name), "rb").read()).hexdigest()
        ok = got == SHA256_OF_THE_TABLES_OF_THE_ANALYSIS[name]
        same = same and ok
        print("  %s  SHA-256 %s  the table of the analysis: %s" % (
            name, got, "the same" if ok else "NOT the same (" + SHA256_OF_THE_TABLES_OF_THE_ANALYSIS[name] + ")"))
    weekly = os.path.join(TOP, "results", WEEKLY)
    shutil.copyfile(os.path.join(outdir, WEEKLY), weekly)
    print("copied:", os.path.relpath(weekly, TOP))
    print("every table is that of the analysis:", "yes" if same else "NO")
    sys.stdout.flush()
    step = subprocess.run([sys.executable, os.path.join(TOP, "07_figures", "add_case_counts_to_figure_1_values.py"), weekly])
    if step.returncode not in (0, 3):
        return 2
    return 0 if same and step.returncode == 0 else 3


if __name__ == "__main__":
    sys.exit(main(sys.argv))
