"""draw_figure_2.py - written for this repository; it is not one of the programs that ran in the analysis.

Draws Figure 2 of the article from its file of values in results/. It does nothing but call
run(...) of fig_R_robustness_focus_short_version_periods_20261007.py with the arguments with which the picture of
the article was drawn (the folder code_main_figures_20261007 and the file of values). In the run of the article no
file of the names of a round was named for this figure; the environment variable MAIN_FIGURES_ROUND is therefore
taken out of the environment before the program is loaded.

usage, from the top folder of the repository:

    MPLBACKEND=Agg PYTHONHASHSEED=0 OMP_NUM_THREADS=1 python 07_figures/draw_figure_2.py OUTDIR

It writes OUTDIR/fig_R_robustness_focus_short_version_periods_20261007.png and .pdf and prints their paths.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
STEM = "fig_R_robustness_focus_short_version_periods_20261007"


def main(outdir):
    os.environ.pop("MAIN_FIGURES_ROUND", None)
    spec = importlib.util.spec_from_file_location(STEM, os.path.join(HERE, STEM + ".py"))
    program = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(program)
    os.makedirs(outdir, exist_ok=True)
    result = program.run(os.path.join(HERE, "code_main_figures_20261007"),
                         os.path.join(TOP, "results", "fig_R_robustness_focus_20261007_values.csv"), outdir)
    for kind in ("png", "pdf"):
        print(kind, result["saved"][kind])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python 07_figures/draw_figure_2.py OUTDIR")
    main(sys.argv[1])
