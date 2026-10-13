"""
fig_R_robustness_focus_short_version_periods_20261007.py. Figure 2 of the article, with the three periods named P1, P2,
P3 and the labels of the thresholds of the called fraction in the form "95 % or more (n)" (likewise 90 % and 80 %),
whatever the number of genomes n in the brackets is; the program stops unless exactly one label of each of the three
thresholds is renamed. The name of the picture carries the label 20261007.

The figure is drawn by code_main_figures_20261007/fig_R_robustness_focus_20261007.py (which imports
main_fig_common_20261007.py and fig_ne_R_cases_20261007.py) from its file of values, neither of which is changed. Before
the picture is saved, every text of the drawing that holds W1, W2 or W3 as a word gets P1, P2, P3, and four further
texts get the words of the legend of the article (WORDS). The program stops if the texts renamed are not exactly those
expected (the heads of the three columns, one label of the key, and the four of WORDS), or if a text of the drawing
still holds one of the stored names.

usage: python fig_R_robustness_focus_short_version_periods_20261007.py CODE_DIR VALUES.csv OUTDIR
"""
import os
import re
import sys

STEM = "fig_R_robustness_focus_short_version_periods_20261007"
EXPECTED = ("W1", "W2", "W3", "W3 ends at the date printed")
# the words of the picture as in the legend of the article
WORDS = {
            "Coalescent-prior cells": "Numerical precision"}
# the labels of the three thresholds of the called fraction, with whatever number of genomes the file of values gives
THRESHOLD = re.compile(r"^0\.(95|90|80) or more \((\d[\d,]*)\)$")


class Stop(SystemExit):
    pass


def run(code_dir, values_path, outdir):
    from matplotlib.text import Text
    if code_dir not in sys.path:
        sys.path.insert(0, code_dir)
    import fig_R_robustness_focus_20261007 as F
    C = F.C
    fig, extra = F.draw(values_path)
    before = [t.get_text() for t in fig.findobj(Text)]
    renamed = []
    for t in fig.findobj(Text):
        new = THRESHOLD.sub(r"\1 % or more (\2)", WORDS.get(t.get_text(), re.sub(r"\bW([123])\b", r"P\1", t.get_text())))
        if new != t.get_text():
            renamed.append([t.get_text(), new])
            t.set_text(new)
    thresholds = sorted(THRESHOLD.match(a).group(1) for a, _ in renamed if THRESHOLD.match(a))
    others = sorted(a for a, _ in renamed if not THRESHOLD.match(a))
    if others != sorted(EXPECTED + tuple(WORDS)) or thresholds != ["80", "90", "95"]:
        raise Stop("the texts renamed are %s; expected %s and one label of each threshold 0.95, 0.90 and 0.80 or more"
                   % (sorted(a for a, _ in renamed), sorted(EXPECTED + tuple(WORDS))))
    left = [t.get_text() for t in fig.findobj(Text) if re.search(r"\bW\d\b|[Ww]indow", t.get_text())]
    if left:
        raise Stop("texts of the picture still hold a name of a window: %s" % left)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    k = fig.dpi / 72.0
    h = fig.bbox.height
    boxes = []
    for t in fig.findobj(Text):
        if t.get_text() in [b for _, b in renamed] and t.get_visible():
            b = t.get_window_extent(rend)
            boxes.append({"text": t.get_text(), "box_pt_from_the_upper_left": [b.x0 / k, (h - b.y1) / k, b.x1 / k, (h - b.y0) / k],
                          "size": float(t.get_fontsize()), "weight": str(t.get_fontweight())})
    os.makedirs(outdir, exist_ok=True)
    saved = C.save_figure(fig, STEM, outdir)
    after = [t.get_text() for t in fig.findobj(Text)]
    F.plt.close(fig)
    return {"saved": saved, "renamed": renamed, "boxes": boxes, "texts_before": before, "texts_after": after,
            "drawing": extra, "size_in": [float(v) for v in fig.get_size_inches()]}


if __name__ == "__main__":
    if len(sys.argv) != 4:
        raise Stop(__doc__)
    r = run(sys.argv[1], sys.argv[2], sys.argv[3])
    print({"png": r["saved"]["png"], "renamed": r["renamed"]})
