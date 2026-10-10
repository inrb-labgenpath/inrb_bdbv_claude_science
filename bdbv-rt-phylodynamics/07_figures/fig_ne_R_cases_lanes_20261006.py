"""fig_ne_R_cases_lanes_20261006.py. Figure 1 of the article with the marks of the response events drawn
in LANES: no line runs along the lower edge of the strip, and no two marks overlap, so that every symbol can be
told from the others.

The figure is drawn by the program fig_ne_R_cases_20261006.py (folder code_main_figures_20261006/)
from its file of values, which is not changed. Two functions of that program are replaced while it draws:
  the plan of the strip   every mark gets its own number beside it and a lane (a height) in which it touches no
                          other mark and no other number; the lowest lane that is free is taken, the number to the
                          right of the mark, or to the left where the right is not free;
  the drawing of the strip no line along the lower edge; a mark is drawn in the middle of its lane (triangle and
                          cross at the day; a solid bar from the first to the last day; a dotted
                          line between the two days where the day is not known).
Nothing else of the program is replaced: the three panels, their sizes and what they hold are drawn by its code.
The checks of that program on the strip are NOT run here (they read the marks back on the lower edge of the strip);
the picture itself was checked instead (pixels outside the strip against the picture of that program; the marks read
back from the pixels of the strip).
"""
import importlib
import math
import os
import sys

LANE_GAP_PT = 1.6       # between two lanes
EXTENT_GAP_PT = 2.0     # smallest distance, in one lane, between two marks with their numbers
NUMBER_PAD_PT = 1.3     # between a mark and its number
BAR_LW_PT = 2.6         # thickness of the bar of a period (first run of this program: a thin bar with a tick at
                        # each end; for a period of two days it read like the letter H beside its number)
STATE = {}


def load(code_dir):
    code_dir = os.path.abspath(code_dir)
    if code_dir not in sys.path:
        sys.path.insert(0, code_dir)
    return importlib.reload(importlib.import_module("fig_ne_R_cases_20261006"))     # afresh at every run


def ink_height_of_the_digits(C):
    import matplotlib
    from matplotlib.textpath import TextPath
    from matplotlib.font_manager import FontProperties
    fp = FontProperties(family=matplotlib.rcParams["font.family"], size=C.SIZES[2])
    return float(TextPath((0, 0), "0123456789", prop=fp).get_extents().y1)


def install(M):
    C = M.C
    original_plan = M.plan_of_the_strip

    def plan_in_lanes(mark_rows, number_rows, xmin_t, xmax_t, declared_t):
        plan = original_plan(mark_rows, number_rows, xmin_t, xmax_t, declared_t)
        if plan["rule"] != "numbered marks":
            raise C.Stop("the lanes are defined for the rule 'numbered marks' only; the plan gives %r" % plan["rule"])
        ev, per_day, width = plan["events"], plan["per_day_pt"], plan["geometry"]["axes_width_pt"]
        XMIN = M.D(xmin_t)
        x_decl = (M.D(declared_t) - XMIN) * per_day
        size = M._measure()
        ink = ink_height_of_the_digits(C)
        body = max(M.MARK_PT, ink)
        placed = []
        for e in ev:
            if e["mark"] in ("triangle", "cross"):
                g0, g1 = e["x_pt"] - M.MARK_PT / 2.0, e["x_pt"] + M.MARK_PT / 2.0
            else:
                a, b = (e["x0"] - XMIN) * per_day, (e["x1"] - XMIN) * per_day
                left = M.ARROW_PT / 2.0 if e["cut_left"] else 0.0
                right = M.ARROW_PT / 2.0 if e["cut_right"] else 0.0
                g0, g1 = a - left, b + right
            w = size(str(e["number"]))[0]
            cands = []
            for side in ("right", "left"):
                n0 = g1 + NUMBER_PAD_PT if side == "right" else g0 - NUMBER_PAD_PT - w
                n1 = n0 + w
                if n0 < -1e-9 or n1 > width + 1e-9:
                    continue
                if n0 - M.LINE_CLEAR_PT < x_decl < n1 + M.LINE_CLEAR_PT:
                    continue                      # a number does not stand on the line of the declaration
                cands.append((side, n0, n1, min(g0, n0), max(g1, n1)))
            if not cands:
                raise C.Stop("the number of mark %s finds no place beside its mark" % e["number"])
            lane, hit = 0, None
            while hit is None:
                for side, n0, n1, x0, x1 in cands:
                    if all(p["lane"] != lane or x1 + EXTENT_GAP_PT <= p["x0"] or p["x1"] + EXTENT_GAP_PT <= x0
                           for p in placed):
                        hit = (side, n0, n1, x0, x1)
                        break
                if hit is None:
                    lane += 1
            e["lane"] = {"lane": lane, "side": hit[0], "n0": hit[1], "n1": hit[2], "x0": hit[3], "x1": hit[4],
                         "g0": g0, "g1": g1}
            placed.append({"lane": lane, "x0": hit[3], "x1": hit[4]})
        n_lanes = 1 + max(p["lane"] for p in placed)
        lane_h = body + LANE_GAP_PT
        top = n_lanes * lane_h - LANE_GAP_PT
        plan["lanes"] = {"lanes": n_lanes, "body_pt": body, "lane_height_pt": lane_h, "ink_height_of_a_number_pt": ink}
        plan["numbers"] = [{"text": str(e["number"]), "numbers": [e["number"]]} for e in ev]
        plan["why"] = plan["why"] + "; every mark has its own number and a lane"
        plan["strip_height_pt"] = M._snap_up(top + M.STRIP_TOP_PT)
        plan["stack_pt"] = M._snap_up(M.STRIP_GAP_PT) + plan["strip_height_pt"]
        room, room_total = plan["room"], plan["room_total_pt"]
        if plan["stack_pt"] > room_total + 1e-6:
            raise C.Stop("the strip needs %.1f pt, the figure allows %.1f pt" % (plan["stack_pt"], room_total))
        s = plan["stack_pt"]
        unit = M.HATCH_STEP_PX * 72.0 / C.DPI
        grow = min(math.ceil(s / unit - 1e-9) * unit, room["growth_of_the_figure_pt"])
        red_ab = min(s - grow, room["from_the_gap_between_a_and_b_pt"])
        red_bc = max(0.0, s - grow - red_ab)
        plan["layout"] = {"growth_of_the_figure_pt": grow, "taken_from_the_gap_between_a_and_b_pt": red_ab,
                          "taken_from_the_gap_between_b_and_c_pt": red_bc,
                          "height_of_the_figure_in": M.HEIGHT + grow / 72.0}
        STATE["plan"] = plan
        return plan

    def draw_in_lanes(fig, axs, plan, XMIN, XMAX, DECLARED):
        H = plan["strip_height_pt"]
        L = plan["lanes"]
        axs.set_ylim(0, H)                                # one unit of the strip is one point
        axs.patch.set_visible(False)
        for s in ("left", "right", "top", "bottom"):
            axs.spines[s].set_visible(False)              # no line along the lower edge
        axs.tick_params(axis="both", which="both", bottom=False, left=False, labelbottom=False, labelleft=False)
        decl, = axs.plot([DECLARED, DECLARED], [0, H], color=C.GREY_REF, ls=(0, (4, 3)), lw=0.9, zorder=1.5)
        per_day = plan["per_day_pt"]
        kw = dict(color=M.STRIP_COLOUR, clip_on=False, zorder=4)
        marks, texts = [], []
        for e in plan["events"]:
            la = e["lane"]
            yc = la["lane"] * L["lane_height_pt"] + L["body_pt"] / 2.0
            arts = []
            if e["mark"] == "triangle":
                m, = axs.plot([e["x0"]], [yc], marker="v", ms=M.MARK_PT, mec="none", ls="none", **kw)
                arts.append(("triangle", m))
            elif e["mark"] == "cross":
                m, = axs.plot([e["x0"]], [yc], marker="x", ms=M.MARK_PT, mew=M.CROSS_LW_PT, ls="none", **kw)
                arts.append(("cross", m))
            else:
                if e["mark"] == "duration":
                    ln, = axs.plot([e["x0"], e["x1"]], [yc, yc], lw=BAR_LW_PT, solid_capstyle="butt", **kw)
                    arts.append(("bar", ln))
                else:
                    ln, = axs.plot([e["x0"], e["x1"]], [yc, yc], lw=M.SPAN_LW_PT, ls=(0, M.DOTS),
                                   dash_capstyle="butt", **kw)
                    arts.append(("dotted line", ln))
                for x, cut, head in ((e["x0"], e["cut_left"], "<"), (e["x1"], e["cut_right"], ">")):
                    if cut:
                        m, = axs.plot([x], [yc], marker=head, ms=M.ARROW_PT, mec="none", ls="none", **kw)
                        arts.append(("arrow head", m))
            marks.append({"row": e["row"], "number": e["number"], "mark": e["mark"], "first": e["first"],
                          "last": e["last"], "x0": e["x0"], "x1": e["x1"], "anchor": e["anchor"],
                          "cut_left": e["cut_left"], "cut_right": e["cut_right"], "artists": arts})
            right = la["side"] == "right"
            x_text = XMIN + (la["n0"] if right else la["n1"]) / per_day
            t = axs.text(x_text, yc - L["ink_height_of_a_number_pt"] / 2.0, str(e["number"]),
                         ha="left" if right else "right", va="baseline", fontsize=C.SIZES[2],
                         color=M.STRIP_COLOUR, zorder=5)
            texts.append({"numbers": [e["number"]], "text": str(e["number"]), "artist": t})
        STATE["axs"] = axs
        return {"marks": marks, "texts": texts, "stems": [], "declaration": decl, "words": False}

    M.plan_of_the_strip = plan_in_lanes
    M.draw_strip = draw_in_lanes
    return M


def boxes_in_pixels(fig, axs, plan, dpi):
    """For every mark: the box of its glyph and of its number in pixels of the PNG (left, top, right, bottom)."""
    pos = axs.get_position()
    w_in, h_in = fig.get_size_inches()
    k = dpi / 72.0
    left, bottom, height_px = pos.x0 * w_in * dpi, pos.y0 * h_in * dpi, h_in * dpi
    L = plan["lanes"]
    out = []
    for e in plan["events"]:
        la = e["lane"]
        y0 = la["lane"] * L["lane_height_pt"]
        y1 = y0 + L["body_pt"]
        row = lambda y: height_px - (bottom + y * k)
        out.append({"number": e["number"], "mark": e["mark"], "first": e["first"], "last": e["last"],
                    "lane": la["lane"], "side": la["side"],
                    "glyph": [left + la["g0"] * k, row(y1), left + la["g1"] * k, row(y0)],
                    "number_box": [left + la["n0"] * k, row(y1), left + la["n1"] * k, row(y0)]})
    return {"events": out, "strip": {"left": left, "right": left + pos.width * w_in * dpi,
                                     "top": height_px - (bottom + plan["strip_height_pt"] * k), "bottom": height_px - bottom},
            "height_px": height_px}


def run(code_dir, values_path, outdir, stem):
    M = install(load(code_dir))
    os.makedirs(outdir, exist_ok=True)
    fig, extra = M.draw(values_path)
    plan, axs = STATE["plan"], STATE["axs"]
    boxes = boxes_in_pixels(fig, axs, plan, M.C.DPI)
    saved = M.C.save_figure(fig, stem, outdir)
    M.plt_close(fig)
    return {"saved": {k: v for k, v in saved.items()}, "boxes": boxes,
            "plan": {"rule": plan["rule"], "lanes": plan["lanes"], "strip_height_pt": plan["strip_height_pt"],
                     "stack_pt": plan["stack_pt"], "layout": plan["layout"], "room_total_pt": plan["room_total_pt"],
                     "events": [{"number": e["number"], "mark": e["mark"], "first": e["first"], "last": e["last"],
                                 "lane": e["lane"]["lane"], "side": e["lane"]["side"]} for e in plan["events"]]},
            "height_in": extra.get("height_in")}
