"""draw_figure_1.py - written for this repository; it is not one of the programs that ran in the analysis.

Draws Figure 1 of the short version of the article from the two files of values in results/. It does nothing
but call run(...) of fig_ne_R_cases_short_version_event_legend_20261006.py with the arguments with which the
picture of the article was drawn (analyses and labels below).

The file of values results/fig_ne_R_cases_20261006_values.csv is not in the repository, because it holds the
weekly case counts. Write it first:

    python 06_case_counts_and_response_events/regenerate_case_tables.py SOURCE_FOLDER

usage, from the top folder of the repository:

    MPLBACKEND=Agg PYTHONHASHSEED=0 python 07_figures/draw_figure_1.py OUTDIR

It writes OUTDIR/fig_ne_R_cases_short_version_event_legend_20261006.png and .pdf and prints their paths.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
STEM = "fig_ne_R_cases_short_version_event_legend_20261006"
ANALYSES = [["beast_A_C00", "a"], ["beast_B_C00", "b"]]
LABELS = ["Delphy, primary analysis", "Delphy, before the median tMRCA", "BEAST X, same model",
          "BEAST X, smoothing estimated"]


def main(outdir):
    values = os.path.join(TOP, "results", "fig_ne_R_cases_20261006_values.csv")
    if not os.path.exists(values):
        raise SystemExit("results/fig_ne_R_cases_20261006_values.csv does not exist: run "
                         "06_case_counts_and_response_events/regenerate_case_tables.py first (see README.md)")
    path = os.path.join(HERE, STEM + ".py")
    spec = importlib.util.spec_from_file_location(STEM, path)
    program = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(program)
    os.makedirs(outdir, exist_ok=True)
    result = program.run(os.path.join(HERE, "fig_ne_R_cases_short_version_20261006.py"),
                         os.path.join(HERE, "code_main_figures_20261006"),
                         os.path.join(HERE, "fig_ne_R_cases_lanes_20261006.py"),
                         values,
                         os.path.join(TOP, "results", "fig_program_trajectories_20261006_values.csv"),
                         outdir, STEM, ANALYSES, LABELS)
    for kind in ("png", "pdf"):
        print(kind, result["saved"][kind])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python 07_figures/draw_figure_1.py OUTDIR")
    main(sys.argv[1])
