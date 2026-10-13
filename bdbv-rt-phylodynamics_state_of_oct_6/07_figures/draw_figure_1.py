"""draw_figure_1.py - written for this repository; it is not one of the programs that ran in the analysis.

Draws Figure 1 (panels a to d) of the article. It does nothing but call run(...) of
fig_ne_R_cases_short_version_zoom_update_20261007.py with the arguments with which the picture of the article was
drawn: the three further programs of this folder that it loads, the folder code_main_figures_20261007, the two
files of values, the analyses and labels below, and the four tables from which panel d is drawn. As in the run of
the article, the environment variable MAIN_FIGURES_ROUND names the file of the names of the tables
(names_of_the_round_main_figures_20261007.json) before any program is loaded.

The file of values results/fig_ne_R_cases_20261007_values.csv is not in the repository, because it holds the
weekly case counts. Write it first:

    python 06_case_counts_and_response_events/regenerate_case_tables.py SOURCE_FOLDER

usage, from the top folder of the repository:

    MPLBACKEND=Agg PYTHONHASHSEED=0 OMP_NUM_THREADS=1 python 07_figures/draw_figure_1.py OUTDIR

It writes OUTDIR/fig_ne_R_cases_short_version_zoom_update_20261007.png and .pdf and, as the program does on every
run, the values of panel d (OUTDIR/fig_ne_R_cases_short_version_zoom_update_20261007_values_panel_d.csv), prints
their paths, and says whether the values of panel d are those of results/. The program is handed each table under
the name of the member and with no version of a file, which is the form in which results/ holds the values
of panel d of the analysis.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
RESULTS = os.path.join(TOP, "results")
STEM = "fig_ne_R_cases_short_version_zoom_update_20261007"
ANALYSES = [("beast_A_C00", "a"), ("beast_B_C00", "b")]
LABELS = ["Delphy, primary analysis", "Delphy, before the median tMRCA", "BEAST X, same model",
          "BEAST X, smoothing estimated"]
TABLES_OF_PANEL_D = {
    "primary_pairs": "primary_pairs_sensitivity_20261007.csv",
    "weekly_series": "trajectories_delphy_sens_gridWeekly32_20261007.csv",
    "weekly_computed": "pairs_of_adjacent_intervals_delphy_sens_gridWeekly32_update_20261007_computed_20261008.csv",
    "events": "figure1b_response_events_20261006.csv",
    "values_abc": "fig_ne_R_cases_20261007_values.csv",
}
ENVIRONMENT_OF_THE_ARTICLE = {"MPLBACKEND": "Agg", "PYTHONHASHSEED": "0", "OMP_NUM_THREADS": "1"}


def main(outdir):
    values = os.path.join(RESULTS, "fig_ne_R_cases_20261007_values.csv")
    if not os.path.exists(values):
        raise SystemExit("results/fig_ne_R_cases_20261007_values.csv does not exist: run "
                         "06_case_counts_and_response_events/regenerate_case_tables.py first (see README.md)")
    other = {k: os.environ.get(k) for k, v in ENVIRONMENT_OF_THE_ARTICLE.items() if os.environ.get(k) != v}
    if other:
        print("note: the picture of the article was drawn with %s; here: %s" % (ENVIRONMENT_OF_THE_ARTICLE, other),
              file=sys.stderr)
    os.environ["MAIN_FIGURES_ROUND"] = os.path.join(HERE, "names_of_the_round_main_figures_20261007.json")
    inp = {key: {"path": os.path.join(RESULTS, name), "file": name, "version": ""}
           for key, name in TABLES_OF_PANEL_D.items()}
    inp["pipeline_dir"] = os.path.join(TOP, "04_summaries_and_reproduction_number", "pipeline_round_20261002",
                                       "work", "pipeline_ext")
    spec = importlib.util.spec_from_file_location(STEM, os.path.join(HERE, STEM + ".py"))
    program = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(program)
    os.makedirs(outdir, exist_ok=True)
    result = program.run(os.path.join(HERE, "fig_ne_R_cases_short_version_event_legend_20261007.py"),
                         os.path.join(HERE, "fig_ne_R_cases_short_version_20261007.py"),
                         os.path.join(HERE, "code_main_figures_20261007"),
                         os.path.join(HERE, "fig_ne_R_cases_lanes_update_20261007.py"),
                         values,
                         os.path.join(RESULTS, "fig_program_trajectories_pairs_A_and_B_estimated_20261007_values.csv"),
                         outdir, STEM, ANALYSES, LABELS, inp)
    for kind in ("png", "pdf"):
        print(kind, result["saved"][kind])
    written = result["values_of_panel_d"]
    held = os.path.join(RESULTS, STEM + "_values_panel_d.csv")
    same = open(written, "rb").read() == open(held, "rb").read()
    print("values of panel d", written)
    print("the values of panel d of results/:", "the same" if same else "NOT the same")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python 07_figures/draw_figure_1.py OUTDIR")
    main(sys.argv[1])
