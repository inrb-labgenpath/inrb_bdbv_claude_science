"""draw_figure_2.py - written for this repository; it is not one of the programs that ran in the analysis.

Draws Figure 2 of the short version of the article from its file of values in results/. It does nothing but call
run(...) of fig_R_robustness_focus_short_version_periods_20261006.py.

usage, from the top folder of the repository:

    MPLBACKEND=Agg PYTHONHASHSEED=0 python 07_figures/draw_figure_2.py OUTDIR

It writes OUTDIR/fig_R_robustness_focus_short_version_periods_20261006.png and .pdf and prints their paths.
"""
import importlib.util
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
STEM = "fig_R_robustness_focus_short_version_periods_20261006"


def main(outdir):
    path = os.path.join(HERE, STEM + ".py")
    spec = importlib.util.spec_from_file_location(STEM, path)
    program = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(program)
    result = program.run(os.path.join(HERE, "code_main_figures_20261006"),
                         os.path.join(TOP, "results", "fig_R_robustness_focus_20261006_values.csv"), outdir)
    for kind in ("png", "pdf"):
        print(kind, result["saved"][kind])


if __name__ == "__main__":
    if len(sys.argv) != 2:
        raise SystemExit("usage: python 07_figures/draw_figure_2.py OUTDIR")
    main(sys.argv[1])
