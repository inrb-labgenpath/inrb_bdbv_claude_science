"""fig_ne_R_cases_short_version_event_legend_20261006.py. Figure 1 of the short version of the article as it is
printed, with a legend of the response events above panel b. It is called by draw_figure_1.py.

The figure is drawn by fig_ne_R_cases_short_version_20261006.py, which this program calls unchanged and whose
function _run it follows step by step; that program draws by code_main_figures_20261006/fig_ne_R_cases_20261006.py
and by fig_ne_R_cases_lanes_20261006.py. What this program adds:

  (1) the marks of the events at 9/8 of their size (their numbers at the size of the numbers of the periods), in at
      most 4 rows;
  (2) the names of the periods (P1, P2, P3; W1 to W3 in the program that draws) in bold, and 'period' for 'window'
      in the key;
  (3) the keys inside panels b and c in the first form (labels broken, size) that stands on white;
  (4) the legend of the nine events above panel b, at the right of the strip of marks: a title, and in two columns
      for each event its mark with the properties that the strip gives it, its number, and its short name;
  (5) more room between the panels, the title of panel b on the line of the title of the legend, every tick label
      at the size of the label of the time axis, and the bar of the root date in panel a a little lower.

The size of the legend is the largest of LEGEND_SIZES_PT at which the legend stands free of the strip, of the ticks of
panel a and of the names of the periods; these sizes lie below the smallest size of the conventions of the figures.

Nothing is read but the two files of values that fig_ne_R_cases_short_version_20261006.py reads. The names and the
title of the legend are constants of this program.

It also holds what the checks need of a live drawing: the count of the texts that overlap, the count of the pixels
that are not white under an artist, and the comparison of two listings of the artists of a drawing.
"""
import importlib.util
import math
import os
import re

MAX_LANES = 4                     # rows of marks at most
STRIP_SCALE = 9.0 / 8.0           # the numbers of the marks at the size of the numbers of the periods
MARK_BOX_PT = 10.0                # the box of the mark
MARK_LINE_MOST_PT = 7.5           # a bar or a dotted line in the legend is at most this long (three dots)
PAD_PT = 1.5                      # from the box of the mark to the number
GAP_PT = 3.5                      # from the box of the number to the words
LEGEND_TITLE = "Response events"  # the title of the legend
LEGEND_SIZES_PT = (7.5, 7.0, 6.5, 6.0)                 # the sizes of the legend that are tried, in this order
LEGEND_COLUMNS = ((1, 2, 3, 4, 5), (6, 7, 8, 9))       # the events of the columns of the legend
LEGEND_BROKEN = {8: ("MSF centre,", "Mongbwalu")}      # the names that stand in two lines
LINE_PITCH = 1.2                  # from base line to base line, in sizes of the font
TITLE_GAP_PT = 1.5                # more than that between the title and the first entry
COLUMN_GAP_PT = 8.0               # between the columns
LEAST_FROM_THE_STRIP_PT = 6.0     # free between the ink of the strip and the legend, at least
LEAST_ABOVE_AND_BELOW_PT = 6.0    # free between the legend and the ticks of panel a, the names of the periods, at least
CONTROL_PAD_PT = 20.0             # control J2: the box of the legend widened by this must not be white
LEAST_BETWEEN_PERIOD_LABELS_PT = 4.0
WHITE_PAD_PT = 1.0                # the box of a key or of the legend is widened by this for the count of white
OVERLAP_W_PT, OVERLAP_H_PT = 0.5, 1.0               # two boxes overlap where they intersect by more than this
NAMES = ("P1", "P2", "P3")        # the names of the periods (W1 to W3 in the original figure program)
MARKS_UP_PT = 4.08                # the strip of the marks is moved up by this (17 pixels at 300 dpi)
GROWTH_PT = 12.24                 # the picture grows by this between panels a and b (51 pixels)
GAP_B_C_PT = 11.52                # the picture grows by this between panels b and c (48 pixels)
HEIGHT_LIMIT_IN = 8.7             # the largest height for the saving of this figure (the conventions hold 8.5 in
                                  # for a drawing 7 in wide; in the document the figure stands 6.5 in wide)
PANEL_B_DOWN_PT = 6.0             # panel b with its head, legend and marks stands this much lower (25 pixels)
HEAD_ABOVE_LEGEND_PT = 4.0        # from the upper edge of the legend to the lower edge of the head of b

WORDS = (                         # the short names of the nine events; no dates
    "First alert to WHO",
    "Field investigation",
    "WHO declares emergency",
    "Treatment tents set on fire",
    "Ituri limits gatherings",
    "MSF centre, Goma",
    "MSF centre, Bunia",
    "MSF centre, Mongbwalu",
    "Burial team attacked",
)
KEY_B = (
    ("period: median and 95 % HPD", "period: median\nand 95 % HPD"),
    ("period: 95 % range, smoothing prior alone", "period: 95 % range,\nsmoothing prior alone"),
    ("pair of adjacent intervals, 95 % HPD", "pair of adjacent\nintervals, 95 % HPD"),
    ("the same, posterior SD / prior SD \u2265 0.90", "the same, posterior SD /\nprior SD \u2265 0.90"),
)
KEY_B_ORDER = ((), (1,), (1, 3), (1, 3, 2), (1, 3, 2, 0))        # the labels that are broken, form by form
KEY_C_LAST = ("province split pooled, 2\u20134 weeks", "province split pooled,\n2\u20134 weeks")


class Stop(SystemExit):
    pass


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


# ---------------------------------------------------------------------------------------------------------------
# what the checks need of a live drawing
# ---------------------------------------------------------------------------------------------------------------
def boxes_of_the_texts(fig):
    """every visible text of the figure that holds a character, with its box of the renderer in points
    (left, bottom, right, top; from the lower left corner of the picture)"""
    from matplotlib.text import Text
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    k = fig.dpi / 72.0
    out = []
    for t in fig.findobj(Text):
        if not t.get_visible() or not t.get_text().strip():
            continue
        b = t.get_window_extent(rend)
        out.append({"text": t.get_text(), "size": float(t.get_fontsize()), "weight": str(t.get_fontweight()),
                    "box": [b.x0 / k, b.y0 / k, b.x1 / k, b.y1 / k]})
    return out


def texts_that_overlap(fig):
    """F5 (i): the pairs of visible texts whose boxes intersect by more than OVERLAP_W_PT in width and more than
    OVERLAP_H_PT in height, and the texts whose box does not lie inside the picture"""
    bx = boxes_of_the_texts(fig)
    w, h = [v * 72.0 for v in fig.get_size_inches()]
    pairs, outside = [], []
    for i, a in enumerate(bx):
        x0, y0, x1, y1 = a["box"]
        if x0 < -0.01 or y0 < -0.01 or x1 > w + 0.01 or y1 > h + 0.01:
            outside.append(a["text"])
        for b in bx[i + 1:]:
            dx = min(x1, b["box"][2]) - max(x0, b["box"][0])
            dy = min(y1, b["box"][3]) - max(y0, b["box"][1])
            if dx > OVERLAP_W_PT and dy > OVERLAP_H_PT:
                pairs.append([a["text"], b["text"], round(dx, 2), round(dy, 2)])
    return {"texts": len(bx), "pairs": pairs, "outside": outside}


def not_white_under(fig, artists, white_from, pad_pt=WHITE_PAD_PT):
    """F5 (ii): in the picture drawn WITHOUT the artists, the number of pixels that are not white in the box that
    holds the artists, widened by pad_pt. The artists are shown again afterwards."""
    import numpy as np
    from matplotlib.transforms import Bbox
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    bb = Bbox.union([a.get_window_extent(rend) for a in artists])
    seen = [a.get_visible() for a in artists]
    for a in artists:
        a.set_visible(False)
    fig.canvas.draw()
    buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()
    for a, v in zip(artists, seen):
        a.set_visible(v)
    fig.canvas.draw()
    hpx, wpx = buf.shape[:2]
    m = pad_pt * fig.dpi / 72.0
    r0, r1 = int(math.floor(hpx - bb.y1 - m)), int(math.ceil(hpx - bb.y0 + m))
    c0, c1 = int(math.floor(bb.x0 - m)), int(math.ceil(bb.x1 + m))
    inside = bool(r0 >= 0 and c0 >= 0 and r1 <= hpx and c1 <= wpx)
    r0, c0, r1, c1 = max(r0, 0), max(c0, 0), min(r1, hpx), min(c1, wpx)
    k = fig.dpi / 72.0
    return {"pixels_not_white": int((buf[r0:r1, c0:c1].min(axis=2) < white_from).sum()),
            "box_widened_inside_the_picture": inside,
            "box_pt": [bb.x0 / k, bb.y0 / k, bb.x1 / k, bb.y1 / k]}


def differences(held0, held1):
    """the listing (function artists_of of the original figure program) of one drawing against that of another.
    Returns a list of differences, each with the exception under which it falls or 'other':
      E1 the axes of the strip are not compared
      E2 the set-off of a word of panel a (axes 0)
      E3 the weight of one of the names W1, W2, W3
      E4 in a key: labels that are equal once the line breaks are taken out; their size"""
    out = []
    if len(held0) != len(held1):
        return [{"exception": "other", "what": "number of axes", "was": len(held0), "is": len(held1)}]
    for d0, d1 in zip(held0, held1):
        n = d0["axes"]
        if d0["strip"] != d1["strip"]:
            out.append({"exception": "other", "axes": n, "what": "strip", "was": d0["strip"], "is": d1["strip"]})
            continue
        if d0["strip"]:
            out.append({"exception": "E1", "axes": n, "what": "the axes of the strip are not compared",
                        "was": len(d0["held"]), "is": len(d1["held"])})
            continue
        for f in ("xlim", "ylim", "xscale", "yscale", "xlabel", "ylabel", "xticklabels", "yticklabels"):
            if d0[f] != d1[f]:
                out.append({"exception": "other", "axes": n, "what": f, "was": d0[f], "is": d1[f]})
        if len(d0["held"]) != len(d1["held"]):
            out.append({"exception": "other", "axes": n, "what": "number of artists", "was": len(d0["held"]),
                        "is": len(d1["held"])})
            continue
        for j, (a0, a1) in enumerate(zip(d0["held"], d1["held"])):
            for f in sorted(set(a0) | set(a1)):
                if a0.get(f) == a1.get(f):
                    continue
                ex = "other"
                if f == "set_off" and a0["kind"] == "word" and n == 0:
                    ex = "E2"
                elif (f == "looks" and a0["kind"] == "word" and a0.get("text") in ("W1", "W2", "W3")
                      and a0.get("text") == a1.get("text") and len(a0[f]) == len(a1[f])
                      and [i for i, (p, q) in enumerate(zip(a0[f], a1[f])) if p != q] == [2]
                      and a0[f][2] == "normal" and a1[f][2] == "bold"):
                    ex = "E3"
                elif a0["kind"] == "key" and a1["kind"] == "key" and f == "text":
                    if [t.replace("\n", " ") for t in a0[f]] == [t.replace("\n", " ") for t in a1[f]]:
                        ex = "E4"
                elif a0["kind"] == "key" and a1["kind"] == "key" and f == "looks":
                    if (len(a0[f]) == 1 and len(a1[f]) == 1 and len(a0[f][0]) == len(a1[f][0])
                            and all(q in (8.0, 9.0) for q in a1[f][0])):
                        ex = "E4"
                out.append({"exception": ex, "axes": n, "artist": j, "kind": a0.get("kind"),
                            "text": a0.get("text"), "what": f, "was": a0.get(f), "is": a1.get(f)})
    return out


# ---------------------------------------------------------------------------------------------------------------
# the drawing
# ---------------------------------------------------------------------------------------------------------------
def legend_of_the_events(M, fig, size):
    """draws the legend above panel b at one size and returns what was drawn and whether it stands free: a title in
    bold and the nine events in the columns LEGEND_COLUMNS; the right edge at the right edge of the key above panel a;
    between the upper edge of the names of the periods and the head of panel b, which it moves above itself
"""
    from matplotlib.lines import Line2D
    I = fig.interventions
    axb, plan, strip = I["axb"], I["plan"], I["strip"]
    axa = fig.short_version["axa"]
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    k = fig.dpi / 72.0
    if [e["number"] for e in plan["events"]] != list(range(1, 10)) or len(strip["marks"]) != 9 \
            or [m["number"] for m in strip["marks"]] != list(range(1, 10)):
        raise Stop("the plan or the strip does not hold the marks 1 to 9 in order")
    if [i for col in LEGEND_COLUMNS for i in col] != list(range(1, 10)):
        raise Stop("the columns of the legend do not hold the events 1 to 9 in order")
    lines_of = {}
    for i in range(1, 10):
        lines_of[i] = list(LEGEND_BROKEN.get(i, (WORDS[i - 1],)))
        if " ".join(lines_of[i]) != WORDS[i - 1]:
            raise Stop("the lines of entry %d are not its name" % i)
    t1 = [t.tick1line for t in axa.xaxis.get_major_ticks() + axa.xaxis.get_minor_ticks() if t.tick1line.get_visible()]
    if not t1:
        raise Stop("panel a shows no tick at its lower edge")
    above = axa.get_window_extent(rend).y0 - max(float(l.get_markersize()) for l in t1) * k    # lower end of the ticks
    names = [t for t in axb.texts if re.fullmatch(r"P\d", t.get_text())]
    if sorted(t.get_text() for t in names) != list(NAMES):
        raise Stop("panel b holds the names %s, P1 to P3 expected" % sorted(t.get_text() for t in names))
    heads = [t for t in I["axes"].texts if hasattr(t, "xyann") and tuple(t.xy) == (0, 1) and t.xycoords == "axes fraction"]
    if sorted(len(t.get_text()) > 1 for t in heads) != [False, True]:
        raise Stop("the strip holds %d annotations at its upper left corner; the letter and the title of panel b are expected" % len(heads))
    hb = [t.get_window_extent(rend) for t in heads]
    head_h = max(b.y1 for b in hb) - min(b.y0 for b in hb)
    below = max(t.get_window_extent(rend).y1 for t in names)                 # the upper edge of the names
    ink = [a.get_window_extent(rend) for m in strip["marks"] for _, a in m["artists"]] \
        + [t["artist"].get_window_extent(rend) for t in strip["texts"]]
    strip_box = [min(b.x0 for b in ink), min(b.y0 for b in ink), max(b.x1 for b in ink), max(b.y1 for b in ink)]
    right = fig.short_version["legend"].get_window_extent(rend).x1          # the right edge of the key above panel a
    numbers = {t["numbers"][0]: t["artist"] for t in strip["texts"]}
    tr = fig.dpi_scale_trans
    per_day = plan["per_day_pt"]

    def measure(text, weight="normal"):
        t = fig.text(0, 0, text, fontsize=size, fontweight=weight)
        w = t.get_window_extent(rend).width
        fp = t.get_fontproperties()
        t.remove()
        return w, fp

    numw = max(measure(str(i))[0] for i in range(1, 10))
    _, fp = measure("lp")
    _, h, desc = rend.get_text_width_height_descent("lp", fp, ismath=False)
    asc = h - desc
    _, digit, _ = rend.get_text_width_height_descent("1", fp, ismath=False)
    lead = (MARK_BOX_PT + PAD_PT + GAP_PT) * k + numw                        # from the left of a column to its names
    widths = [max(measure(l)[0] for i in col for l in lines_of[i]) for col in LEGEND_COLUMNS]
    total = max(sum(lead + w for w in widths) + COLUMN_GAP_PT * k * (len(LEGEND_COLUMNS) - 1),
                measure(LEGEND_TITLE, "bold")[0])
    step = LINE_PITCH * size * k
    n_lines = [sum(len(lines_of[i]) for i in col) for col in LEGEND_COLUMNS]
    height = asc + TITLE_GAP_PT * k + max(n_lines) * step + desc
    x0 = right - total
    # from below the names, a free space, the legend, HEAD_ABOVE_LEGEND_PT, the head of panel b, the same
    # free space, the ticks of panel a
    free = ((above - below) - height - head_h - (HEAD_ABOVE_LEGEND_PT + PANEL_B_DOWN_PT) * k) / 2.0
    y1 = below + free + height
    # the title of panel b stands on the baseline of the title of the legend; the legend keeps the place
    # of the line above, and the space above it stays free
    title_b = [t for t in heads if len(t.get_text()) > 1][0]
    d_head = rend.get_text_width_height_descent("lp", title_b.get_fontproperties(), ismath=False)[2]
    head_up = ((y1 - asc) - (title_b.get_window_extent(rend).y0 + d_head)) / k
    for t in heads:
        t.xyann = (t.xyann[0], t.xyann[1] + head_up)
    base_t = y1 - asc
    title = fig.text(x0 / fig.dpi, base_t / fig.dpi, LEGEND_TITLE, transform=tr, ha="left", va="baseline",
                     fontsize=size, fontweight="bold", zorder=5)
    artists, entries = [title], []
    x = x0
    for c, (col, w) in enumerate(zip(LEGEND_COLUMNS, widths)):
        n = 1
        for i in col:
            e, m = plan["events"][i - 1], strip["marks"][i - 1]
            if (e["number"], e["mark"], e["first"], e["last"]) != (m["number"], m["mark"], m["first"], m["last"]) \
                    or e["number"] != i:
                raise Stop("mark %s of the strip is not event %s of the plan" % (m["number"], e["number"]))
            base = base_t - TITLE_GAP_PT * k - n * step
            words = [fig.text((x + lead) / fig.dpi, (base - j * step) / fig.dpi, l, transform=tr, ha="left",
                              va="baseline", fontsize=size, zorder=5) for j, l in enumerate(lines_of[i])]
            src_n = numbers[i]
            num = fig.text((x + (MARK_BOX_PT + PAD_PT) * k) / fig.dpi, base / fig.dpi, str(i), transform=tr,
                           ha="left", va="baseline", fontsize=size, color=src_n.get_color(), zorder=5)
            kind, src = m["artists"][0]
            xc, yc = x + MARK_BOX_PT * k / 2.0, base + digit / 2.0
            if kind in ("triangle", "cross"):
                xs, length_pt = [xc / fig.dpi], 0.0
            elif kind in ("bar", "dotted line"):
                length_pt = min(MARK_LINE_MOST_PT, (m["x1"] - m["x0"]) * per_day)
                xs = [(xc - length_pt * k / 2.0) / fig.dpi, (xc + length_pt * k / 2.0) / fig.dpi]
            else:
                raise Stop("a mark of the kind %r is not known to the legend" % kind)
            a = Line2D(xs, [yc / fig.dpi] * len(xs))
            a.update_from(src)                                               # the looks of the mark in the strip
            a.set_transform(tr)
            a.set_clip_on(False)
            a.set_clip_box(None)
            a.set_clip_path(None)
            a.set_zorder(5)
            fig.add_artist(a)
            artists += words + [num, a]
            bbs = [t.get_window_extent(rend) for t in words]
            entries.append({"number": i, "column": c, "mark": e["mark"], "kind_in_the_strip": kind,
                            "first": e["first"], "last": e["last"], "lines": list(lines_of[i]),
                            "length_of_the_mark_pt": length_pt,
                            "length_in_the_strip_pt": (m["x1"] - m["x0"]) * per_day
                            if kind in ("bar", "dotted line") else 0.0,
                            "size_of_the_number_pt": float(num.get_fontsize()),
                            "size_of_the_words_pt": [float(t.get_fontsize()) for t in words],
                            "colour_of_the_number": str(num.get_color()),
                            "colour_of_the_number_in_the_strip": str(src_n.get_color()),
                            "colour_of_the_mark": str(a.get_color()), "marker": str(a.get_marker()),
                            "markersize_pt": float(a.get_markersize()), "linewidth_pt": float(a.get_linewidth()),
                            "box_of_the_words_pt": [min(b.x0 for b in bbs) / k, min(b.y0 for b in bbs) / k,
                                                    max(b.x1 for b in bbs) / k, max(b.y1 for b in bbs) / k],
                            "base_line_pt": base / k, "centre_of_the_mark_pt": [xc / k, yc / k]})
            n += len(lines_of[i])
        x += lead + w + COLUMN_GAP_PT * k
    from_the_strip = (x0 - strip_box[2]) / k
    above_and_below = free / k
    tb = title.get_window_extent(rend)
    return {"artists": artists, "entries": entries, "size_pt": float(size),
            "title": {"text": title.get_text(), "weight": str(title.get_fontweight()),
                      "size_pt": float(title.get_fontsize()), "box_pt": [tb.x0 / k, tb.y0 / k, tb.x1 / k, tb.y1 / k]},
            "columns": [list(col) for col in LEGEND_COLUMNS], "lines": n_lines,
            "width_pt": total / k, "height_pt": height / k, "room_pt": (above - below) / k,
            "box_pt": [x0 / k, (y1 - height) / k, right / k, y1 / k],
            "from_the_strip_pt": from_the_strip, "above_and_below_pt": above_and_below,
            "fits": bool(from_the_strip >= LEAST_FROM_THE_STRIP_PT - 1e-9
                         and above_and_below >= LEAST_ABOVE_AND_BELOW_PT - 1e-9),
            "lower_end_of_the_ticks_of_panel_a_pt": above / k, "upper_edge_of_the_names_pt": below / k,
            "box_of_the_strip_pt": [v / k for v in strip_box],
            "head_of_panel_b": {"texts": [t.get_text() for t in heads], "height_pt": head_h / k, "moved_up_pt": head_up,
                                "on_the_baseline_of_the_title_of_the_legend": True, "baseline_pt": (y1 - asc) / k,
                                "free_above_and_below_pt": free / k},
            "right_edge_of_the_key_above_panel_a_pt": right / k}


ROOT_DOWN_PT = 3.6                # the bar of the root date in panel a and its dates stand this much lower
ROOT_Y, ROOT_CAP = 0.075, 0.022   # the height of that bar and half its cap in the original figure program (share of the frame)


def root_date_down(fig, down_pt):
    """the bar of the root date (bar, cap, arrow or second cap, diamond: four lines that the original figure program draws
    at a share of the height of the frame) and the printed dates above it move down by down_pt"""
    axa = fig.interventions["axa"]
    d = down_pt / (axa.get_position().height * fig.get_figheight() * 72.0)
    lines = [l for l in axa.lines if l.get_transform() is not axa.transData and len(l.get_ydata())
             and all(ROOT_Y - ROOT_CAP - 1e-9 <= float(y) <= ROOT_Y + ROOT_CAP + 1e-9 for y in l.get_ydata())]
    texts = [t for t in axa.texts if t.get_text().startswith("tMRCA ")]
    if len(lines) != 4 or len(texts) != 1:
        raise Stop("the bar of the root date: %d lines and %d texts found, 4 and 1 expected" % (len(lines), len(texts)))
    for l in lines:
        l.set_ydata([float(y) - d for y in l.get_ydata()])
    x, y = texts[0].get_position()
    if abs(y - (ROOT_Y + 0.035)) > 1e-9:
        raise Stop("the dates of the root stand at %s of the frame, %s expected" % (y, ROOT_Y + 0.035))
    texts[0].set_position((x, y - d))
    return {"down_pt": float(down_pt), "share_of_the_frame": float(d), "lines_moved": len(lines), "texts_moved": [t.get_text() for t in texts],
            "lowest_point_of_a_line_pt_above_the_axis": float((ROOT_Y - ROOT_CAP - d) * down_pt / d)}


def tick_labels(fig, size):
    """every tick label of the figure (the months and the numbers of the vertical axes) at the size of
    the label of the time axis"""
    for ax in fig.axes:
        ax.tick_params(axis="both", which="both", labelsize=size)
    return float(size)


def periods_and_room(fig, dpi):
    """the names W1 to W3 of the periods become P1 to P3 and 'window:' becomes 'period:' in the key
    inside panel b; the strip of the marks moves up by MARKS_UP_PT; panel a (with the axes that share its place)
    moves up by GROWTH_PT, and the picture grows in height by GROWTH_PT. Panels b and c keep place and size."""
    from matplotlib.text import Text, Annotation
    I = fig.interventions
    axb, axs = I["axb"], I["axes"]
    axa = fig.short_version["axa"]
    from matplotlib.patches import Rectangle
    own = list(fig.texts) + list(fig.legends) + list(fig.lines) + list(fig.patches) + list(fig.images) \
        + [a for a in fig.artists if not (isinstance(a, Rectangle) and a.get_data_transform() == fig.transFigure
                                          and [round(float(v), 9) for v in a.get_bbox().bounds] == [0.0, 0.0, 1.0, 1.0])]
    if own:      # a rectangle over the whole figure in coordinates of the figure (its ground) grows with it
        raise Stop("the figure holds artists of its own before the key above panel a is placed; they would not move: %s"
                   % [type(a).__name__ for a in own])
    follows = []                 # a twin axes is placed by a locator on its twin and moves with it
    for ax in fig.axes:
        if ax.get_axes_locator() is not None:
            twins = [o for o in ax._twinned_axes.get_siblings(ax) if o is not ax and o.get_axes_locator() is None]
            if len(twins) != 1:
                raise Stop("an axes of the figure is placed by a locator and is not the twin of one axes")
            follows.append(ax)
        for t in ax.texts:
            if isinstance(t, Annotation):
                for c in (t.xycoords, t.anncoords):
                    for one in (c if isinstance(c, (tuple, list)) else (c,)):
                        if isinstance(one, str) and one.startswith(("figure", "subfigure")):
                            raise Stop("the annotation %r is placed in coordinates of the figure" % t.get_text())
    names = sorted([t for t in axb.texts if re.fullmatch(r"W\d", t.get_text())], key=lambda t: t.get_text())
    if [t.get_text() for t in names] != ["W1", "W2", "W3"]:
        raise Stop("panel b holds the names %s, W1 to W3 expected" % [t.get_text() for t in names])
    for t, n in zip(names, NAMES):
        t.set_text(n)
    changed = []
    for t in axb.get_legend().get_texts():
        if t.get_text().startswith("window: "):
            changed.append([t.get_text(), "period: " + t.get_text()[len("window: "):]])
            t.set_text(changed[-1][1])
    if len(changed) != 2:
        raise Stop("%d labels of the key of panel b begin with 'window: ', 2 expected" % len(changed))
    left = [t.get_text() for t in fig.findobj(Text) if re.search(r"\bW\d\b|[Ww]indow", t.get_text())]
    if left:
        raise Stop("texts of the picture still hold a name of a window: %s" % left)
    w_in, h0 = [float(v) for v in fig.get_size_inches()]
    px = [v * dpi / 72.0 for v in (GROWTH_PT, MARKS_UP_PT, PANEL_B_DOWN_PT, GAP_B_C_PT)]
    if max(abs(v - round(v)) for v in px) > 1e-6:
        raise Stop("GROWTH_PT and MARKS_UP_PT are not whole pixels at %s dpi: %s" % (dpi, px))
    rows = int(round(h0 * dpi)) + int(round(px[0])) + int(round(px[3]))      # the height of the saved picture in pixels
    h1 = rows / float(dpi)
    while int(h1 * dpi) < rows:                               # the renderer cuts the height to whole pixels
        h1 = math.nextafter(h1, math.inf)
    ya, ys, yb = axa.get_position().y0, axs.get_position().y0, axb.get_position().y0
    if not ya > ys > axb.get_position().y0:
        raise Stop("panel a, the strip and panel b do not stand in this order from above")
    before = [(ax, [float(v) for v in ax.get_position().bounds]) for ax in fig.axes]
    fig.set_size_inches(w_in, h1)
    moved = []
    for ax, (x0, y0, w, h) in before:
        if any(ax is f for f in follows):
            moved.append("with its twin")
            continue
        up = GROWTH_PT + GAP_B_C_PT if y0 >= ya - 1e-9 else MARKS_UP_PT - PANEL_B_DOWN_PT + GAP_B_C_PT if y0 >= ys - 1e-9 \
            else GAP_B_C_PT - PANEL_B_DOWN_PT if y0 >= yb - 1e-9 else 0.0      # panel c keeps its place from below
        ax.set_position([x0, (y0 * h0 + up / 72.0) / h1, w, h * h0 / h1])
        moved.append(up)
    return {"names": list(NAMES), "labels_of_the_key_of_panel_b": changed, "height_before_in": h0, "height_in": h1,
            "axes_moved_up_pt": moved, "MARKS_UP_PT": MARKS_UP_PT, "GROWTH_PT": GROWTH_PT,
            "height_px": rows}


def draw(S, code_dir, lanes_path, values_path, s2_values_path, analyses, labels):
    """one drawing at the full width of the frames: the steps of _run of the original figure program, with the marks at
    STRIP_SCALE in at most MAX_LANES rows; the legend is drawn afterwards"""
    W = S.load_lanes(lanes_path)
    M = W.install(W.load(code_dir))
    C = M.C
    number_size = C.SIZES[2] * STRIP_SCALE
    state = {}
    S.install_search(M, W, MAX_LANES, state, STRIP_SCALE)
    sets = S.beast_rows(s2_values_path, analyses, C.Stop)
    fig, extra = M.draw(values_path)
    tick_size = tick_labels(fig, C.SIZES[0])
    root = root_date_down(fig, ROOT_DOWN_PT)
    axa, floor, drawn = S.add_beast(M, fig, sets)
    delphy = [r for r in C.read_csv(values_path) if r["element"] == M.E_STAIR]
    fig.short_version = {"axa": axa, "beast": drawn, "legend": None,
                         "delphy": [{"x": M.D(r["x"]), "x_end": M.D(r["x_end"]), "value": float(r["value"]),
                                     "lower": float(r["lower"]), "upper": float(r["upper"])} for r in delphy]}
    room = periods_and_room(fig, C.DPI)
    moved = S.clear_the_words(M, fig)
    leg, key = S.key_above(M, fig, axa, labels)
    fig.short_version["legend"] = leg
    over = S.touches(fig)
    if state["sizes"]["size_of_the_numbers_pt"] != number_size:
        raise Stop("the numbers of the strip are at %s pt, %s pt expected" % (state["sizes"]["size_of_the_numbers_pt"],
                                                                             number_size))
    return {"S": S, "W": W, "M": M, "fig": fig, "extra": extra, "state": state, "moved": moved, "key": key,
            "touches": over, "floor": floor, "drawn": drawn, "axa": axa, "room": room, "tick_size": tick_size, "root": root,
            "frames_in": float(C.FIG_WIDTH * (M.GRID["right"] - M.GRID["left"]))}


def names_in_bold(fig):
    from matplotlib import font_manager as fm
    axb = fig.interventions["axb"]
    names = [t for t in axb.texts if re.fullmatch(r"P\d", t.get_text())]
    if sorted(t.get_text() for t in names) != list(NAMES):
        raise Stop("panel b holds the names %s, P1 to P3 expected" % sorted(t.get_text() for t in names))
    numbers = [t for t in axb.texts if re.match(r"^\d+\.\d+ \(", t.get_text())]
    if len(numbers) != 3:
        raise Stop("panel b holds %d labels of period numbers, 3 expected" % len(numbers))
    before = [str(t.get_fontweight()) for t in names]
    for t in names:
        t.set_fontweight("bold")
    return {"names": [t.get_text() for t in names], "weight_before": before,
            "weight": [str(t.get_fontweight()) for t in names], "size_pt": [float(t.get_fontsize()) for t in names],
            "font_file_of_the_names": sorted(set(os.path.basename(fm.findfont(t.get_fontproperties())) for t in names)),
            "font_file_of_the_numbers": sorted(set(os.path.basename(fm.findfont(t.get_fontproperties()))
                                                   for t in numbers)),
            "size_of_the_numbers_pt": [float(t.get_fontsize()) for t in numbers]}


def keys_on_white(S, M, fig):
    """rule: the first form of the key of b, and then of the key of c at the size that
    the key of b got, under which the picture drawn without the key is white"""
    C = M.C
    I = fig.interventions
    leg_b, leg_c = I["axb"].get_legend(), I["axc"].get_legend()
    tb, tc = leg_b.get_texts(), leg_c.get_texts()
    if [t.get_text() for t in tb] != [a for a, _ in KEY_B]:
        raise Stop("the key of panel b holds %s" % [t.get_text() for t in tb])
    if tc[-1].get_text() != KEY_C_LAST[0]:
        raise Stop("the last label of the key of panel c is %r" % tc[-1].get_text())
    rec = {"b": {"size_before_pt": [float(t.get_fontsize()) for t in tb], "tried": []},
           "c": {"size_before_pt": [float(t.get_fontsize()) for t in tc], "tried": []}}
    got = None
    for size in (C.SIZES[1], C.SIZES[2]):
        for form, broken in enumerate(KEY_B_ORDER):
            for i, t in enumerate(tb):
                t.set_text(KEY_B[i][1] if i in broken else KEY_B[i][0])
                t.set_fontsize(size)
            w = not_white_under(fig, [leg_b], S.WHITE_FROM)
            rec["b"]["tried"].append({"size_pt": float(size), "form": form, "pixels_not_white": w["pixels_not_white"],
                                      "inside": w["box_widened_inside_the_picture"]})
            if w["pixels_not_white"] == 0 and w["box_widened_inside_the_picture"]:
                got = (size, form)
                break
        if got:
            break
    if not got:
        raise Stop("no form of the key of panel b stands on white")
    rec["b"].update(size_pt=float(got[0]), form=got[1], labels=[t.get_text() for t in tb])
    got_c = None
    for t in tc:
        t.set_fontsize(got[0])
    for form, last in enumerate(KEY_C_LAST):
        tc[-1].set_text(last)
        w = not_white_under(fig, [leg_c], S.WHITE_FROM)
        rec["c"]["tried"].append({"size_pt": float(got[0]), "form": form, "pixels_not_white": w["pixels_not_white"],
                                  "inside": w["box_widened_inside_the_picture"]})
        if w["pixels_not_white"] == 0 and w["box_widened_inside_the_picture"]:
            got_c = form
            break
    if got_c is None:
        raise Stop("no form of the key of panel c stands on white")
    rec["c"].update(size_pt=float(got[0]), form=got_c, labels=[t.get_text() for t in tc])
    return rec


def run(short_program, code_dir, lanes_path, values_path, s2_values_path, outdir, stem, analyses, labels):
    from matplotlib.text import Text
    S = load(short_program, "fig1_short_version_stored")
    d = draw(S, code_dir, lanes_path, values_path, s2_values_path, analyses, labels)
    M, W, fig, C = d["M"], d["W"], d["fig"], d["M"].C
    I = fig.interventions
    axa, axb, axc, axs, plan, strip = d["axa"], I["axb"], I["axc"], I["axes"], I["plan"], I["strip"]
    os.makedirs(outdir, exist_ok=True)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    k = fig.dpi / 72.0
    num = sorted([t for t in axb.texts if re.match(r"^\d+\.\d+ \(", t.get_text())],
                 key=lambda t: t.get_window_extent(rend).x0)
    nb = [t.get_window_extent(rend) for t in num]
    gaps = [(nb[i + 1].x0 - nb[i].x1) / k for i in range(len(nb) - 1)]
    if len(num) != 3 or min(gaps) < LEAST_BETWEEN_PERIOD_LABELS_PT:
        raise Stop("the labels of the period numbers (%d) stand %s pt apart; at least %.1f pt are asked"
                   % (len(num), [round(g, 1) for g in gaps], LEAST_BETWEEN_PERIOD_LABELS_PT))
    names = names_in_bold(fig)
    keys = keys_on_white(S, M, fig)
    sizes, lg = [], None
    for size in LEGEND_SIZES_PT:
        lg = legend_of_the_events(M, fig, size)
        sizes.append({"size_pt": float(size), "width_pt": lg["width_pt"], "height_pt": lg["height_pt"],
                      "from_the_strip_pt": lg["from_the_strip_pt"], "above_and_below_pt": lg["above_and_below_pt"],
                      "fits": lg["fits"]})
        if lg["fits"]:
            break
        for a in lg["artists"]:
            a.remove()
        lg = None
    if lg is None:
        raise Stop("at no size of %s pt does the legend stand free above panel b" % (LEGEND_SIZES_PT,))
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    # control J3: a text over the title of panel b must be counted as an overlap
    title = [t for t in axs.texts if hasattr(t, "xyann") and len(t.get_text()) > 1][0]
    tb = title.get_window_extent(rend)
    extra_text = fig.text((tb.x0 + 2.0) / fig.dpi, (tb.y0 + 1.0) / fig.dpi, "Control", transform=fig.dpi_scale_trans,
                          ha="left", va="bottom", fontsize=C.SIZES[0])
    j3 = texts_that_overlap(fig)
    extra_text.remove()
    overlap = texts_that_overlap(fig)
    white = {"legend of the events": not_white_under(fig, lg["artists"], S.WHITE_FROM),
             "key of panel b": not_white_under(fig, [axb.get_legend()], S.WHITE_FROM),
             "key of panel c": not_white_under(fig, [axc.get_legend()], S.WHITE_FROM),
             "key above panel a": not_white_under(fig, [fig.short_version["legend"]], S.WHITE_FROM)}
    # control J2: the count for the legend with its box widened far enough to reach the strip and the names
    j2 = not_white_under(fig, lg["artists"], S.WHITE_FROM, pad_pt=CONTROL_PAD_PT)
    touches_at_the_end = S.touches(fig)
    texts = boxes_of_the_texts(fig)
    of_the_legend = set(id(a) for a in lg["artists"])
    shown = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    sizes_of_the_texts = {"outside_the_legend": sorted(set(float(t.get_fontsize()) for t in shown
                                                           if id(t) not in of_the_legend)),
                          "of_the_legend": sorted(set(float(t.get_fontsize()) for t in shown
                                                      if id(t) in of_the_legend)),
                          "texts_of_the_legend": sum(1 for t in shown if id(t) in of_the_legend),
                          "texts": len(shown)}
    rend = fig.canvas.get_renderer()
    name_boxes = [[b.x0 / k, b.y0 / k, b.x1 / k, b.y1 / k] for b in
                  (t.get_window_extent(rend) for t in axb.texts if re.fullmatch(r"P\d", t.get_text()))]
    held = S.artists_of(fig, [axs])
    boxes = W.boxes_in_pixels(fig, axs, plan, C.DPI)
    numbers_of_the_strip = [{"number": t["numbers"][0], "size_pt": float(t["artist"].get_fontsize()),
                             "colour": str(t["artist"].get_color())} for t in strip["texts"]]
    pos = {n_: [float(v) for v in a.get_position().bounds] for n_, a in (("a", axa), ("b", axb), ("c", axc), ("strip", axs))}
    limit_of_the_conventions = C.MAX_HEIGHT
    C.MAX_HEIGHT = max(C.MAX_HEIGHT, HEIGHT_LIMIT_IN)        # for the saving of this figure alone
    try:
        saved = C.save_figure(fig, stem, outdir)
    finally:
        C.MAX_HEIGHT = limit_of_the_conventions
    saved["limit_of_the_conventions_in"] = float(limit_of_the_conventions)
    saved["limit_for_this_figure_in"] = float(HEIGHT_LIMIT_IN)
    size_in = [float(v) for v in fig.get_size_inches()]
    M.plt_close(fig)
    return {"saved": dict(saved), "size_in": size_in, "frames_in": d["frames_in"], "sizes_tried": sizes,
            "legend": {k_: v for k_, v in lg.items() if k_ != "artists"}, "names": names, "keys": keys,
            "key_above_a": d["key"], "words_moved": d["moved"], "touches": d["touches"], "room": d["room"],
            "tick_labels_pt": d["tick_size"], "root_date": d["root"],
            "touches_at_the_end": touches_at_the_end,
            "keys_or_words_over_data_by_the_track": d["extra"].get("keys_or_words_over_data"),
            "period_labels": {"texts": [t.get_text() for t in num], "gaps_pt": gaps,
                              "size_pt": [float(t.get_fontsize()) for t in num]},
            "strip": d["state"]["sizes"], "numbers_of_the_strip": numbers_of_the_strip,
            "plan": {"rule": plan["rule"], "lanes": plan["lanes"], "strip_height_pt": plan["strip_height_pt"],
                     "stack_pt": plan["stack_pt"], "layout": plan["layout"], "room_total_pt": plan["room_total_pt"],
                     "per_day_pt": plan["per_day_pt"],
                     "events": [{"number": e["number"], "mark": e["mark"], "first": e["first"], "last": e["last"],
                                 "lane": e["lane"]["lane"], "side": e["lane"]["side"]} for e in plan["events"]]},
            "boxes": boxes, "positions_of_the_axes": pos, "texts": texts, "boxes_of_the_names_pt": name_boxes,
            "sizes_of_the_texts": sizes_of_the_texts,
            "overlap": overlap, "control_J3": {"pairs": j3["pairs"], "outside": j3["outside"]}, "control_J2": j2,
            "white": white, "held": held,
            "constants": {"MAX_LANES": MAX_LANES, "STRIP_SCALE": STRIP_SCALE, "MARK_BOX_PT": MARK_BOX_PT,
                          "MARK_LINE_MOST_PT": MARK_LINE_MOST_PT, "PAD_PT": PAD_PT, "GAP_PT": GAP_PT,
                          "LEGEND_TITLE": LEGEND_TITLE, "LEGEND_SIZES_PT": list(LEGEND_SIZES_PT), "NAMES": list(NAMES),
                          "MARKS_UP_PT": MARKS_UP_PT, "GROWTH_PT": GROWTH_PT, "HEAD_ABOVE_LEGEND_PT": HEAD_ABOVE_LEGEND_PT,
                          "PANEL_B_DOWN_PT": PANEL_B_DOWN_PT,
                          "LEGEND_COLUMNS": [list(col) for col in LEGEND_COLUMNS],
                          "LEGEND_BROKEN": {str(i): list(v) for i, v in LEGEND_BROKEN.items()},
                          "LINE_PITCH": LINE_PITCH, "TITLE_GAP_PT": TITLE_GAP_PT, "COLUMN_GAP_PT": COLUMN_GAP_PT,
                          "LEAST_FROM_THE_STRIP_PT": LEAST_FROM_THE_STRIP_PT,
                          "LEAST_ABOVE_AND_BELOW_PT": LEAST_ABOVE_AND_BELOW_PT, "CONTROL_PAD_PT": CONTROL_PAD_PT,
                          "LEAST_BETWEEN_PERIOD_LABELS_PT": LEAST_BETWEEN_PERIOD_LABELS_PT,
                          "WHITE_PAD_PT": WHITE_PAD_PT, "OVERLAP_W_PT": OVERLAP_W_PT, "OVERLAP_H_PT": OVERLAP_H_PT,
                          "WORDS": list(WORDS), "KEY_B": [list(x) for x in KEY_B],
                          "KEY_B_ORDER": [list(x) for x in KEY_B_ORDER], "KEY_C_LAST": list(KEY_C_LAST)}}
