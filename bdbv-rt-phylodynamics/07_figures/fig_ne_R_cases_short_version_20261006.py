"""fig_ne_R_cases_short_version_20261006.py. Figure 1 of the short version of the article.
It is called by draw_figure_1.py.

The figure is drawn by the program fig_ne_R_cases_20261006.py (folder code_main_figures_20261006/) from its
file of values, through fig_ne_R_cases_lanes_20261006.py, which draws the marks of the events in lanes.
What this program adds to that picture or changes in it:
  (a') panel a also draws the effective population size of two analyses with BEAST X, from the rows that the file
      of values fig_program_trajectories_20261006_values.csv marks as drawn, in the colour that the conventions give to BEAST X; solid for the
      first analysis named, dashed for the second (changed after a viewing of the picture: DOTTED for the first, so that
      the median of the primary analysis shows between the dots). The median is a stair line; the bounds of the 95 % HPD interval
      are thin horizontal lines over each interval, without joins. The axis is that of the program. A bound below the
      lowest power of ten of the axis (the room that the program keeps for the bar of the root date) is not drawn.
  (a'') ONE key with four entries stands above the frame of panel a, to the
      right of the title of the panel; the key of the program inside the frame is removed. A word of the program inside
      the frame that touches a line or the band of an analysis is moved to the left to the first free place.
  (w) the wide form: with fig_width the figure is drawn that many
      inches wide, and the left and the right margin of the grid keep their inches. artists_of(fig) lists what the
      axes hold, for the comparison of two drawings.
  (b) the lanes and the sides of the numbers of the marks are found by a search over all assignments: fewest lanes,
      then fewest numbers to the left of their mark; the distances are those of fig_ne_R_cases_lanes_20261006.py. Where
      more lanes than max_lanes are needed the program stops.
  (b') strip_scale: the marks of the strip, the widths of their lines,
      their numbers and the distances between them are drawn at this share of their size. Nothing else is scaled.
Nothing else is changed: the three panels, their sizes and what they hold are drawn by fig_ne_R_cases_20261006.py.
No function prints.
"""
import csv
import importlib.util
import math
import os

MAX_LANES = 2
BEAST_MEDIAN_LW = 1.2          # thinner than the median of the primary analysis (1.7), so that both stay visible
BEAST_BOUND_LW = 0.6
BEAST_DASH = (0, (4.0, 2.0))   # the second analysis: dashed, in widths of the line
BEAST_DOT_PT = (1.2, 1.4)      # the first analysis: dotted; drawn and open, in points
E_NE = "effective population size of an interval"
KEY_MARGIN_PT = 1.0            # the free margin around the key above panel a
WHITE_FROM = 250               # a pixel is free where its three channels are at least this


def style(n, lw):
    """the line style of the n-th analysis with BEAST X for a line of the width lw: the first is dotted, so that the
    median of the primary analysis shows between the dots where the two lie on each other; the second is dashed"""
    return (0, (BEAST_DOT_PT[0] / lw, BEAST_DOT_PT[1] / lw)) if n == 0 else BEAST_DASH


def load_lanes(path):
    spec = importlib.util.spec_from_file_location("fig_ne_R_cases_lanes", path)
    W = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(W)
    return W


def search(cands, n_lanes, gap):
    """cands: for every mark the list of its places (side, n0, n1, x0, x1) in points. Returns for every mark
    (lane, index of the place) so that in a lane two extents are at least `gap` apart, or None. Among all
    assignments: fewest numbers to the left, then the lowest lanes, then right before left, in the order of the marks."""
    n = len(cands)
    best = [None, None]
    cur = []

    def free(lane, x0, x1):
        for (l2, k2), c2 in zip(cur, cands):
            if l2 == lane and not (x1 + gap <= c2[k2][3] or c2[k2][4] + gap <= x0):
                return False
        return True

    def key():
        return (sum(1 for (l, k), c in zip(cur, cands) if c[k][0] == "left"), tuple(l for l, _ in cur),
                tuple(0 if c[k][0] == "right" else 1 for (l, k), c in zip(cur, cands)))

    def rec(i):
        if i == n:
            k = key()
            if best[0] is None or k < best[1]:
                best[0], best[1] = list(cur), k
            return
        for lane in range(n_lanes):
            for k, c in enumerate(cands[i]):
                if free(lane, c[3], c[4]):
                    cur.append((lane, k))
                    rec(i + 1)
                    cur.pop()

    rec(0)
    return best[0]


class scaled_strip:
    """while it is open, the sizes of the strip (marks, widths of their lines, numbers, distances) are multiplied by s"""
    NAMES_M = ("MARK_PT", "CROSS_LW_PT", "SPAN_LW_PT", "ARROW_PT")
    NAMES_W = ("BAR_LW_PT", "NUMBER_PAD_PT", "EXTENT_GAP_PT", "LANE_GAP_PT")

    def __init__(self, M, W, s):
        self.M, self.W, self.s = M, W, float(s)

    def __enter__(self):
        M, W, s = self.M, self.W, self.s
        self.old = [(M, n, getattr(M, n)) for n in self.NAMES_M] + [(W, n, getattr(W, n)) for n in self.NAMES_W]
        self.sizes = M.C.SIZES
        if s != 1.0:
            for obj, n, v in self.old:
                setattr(obj, n, v * s)
            M.C.SIZES = tuple(self.sizes[:2]) + (self.sizes[2] * s,)
        return self

    def __exit__(self, *a):
        for obj, n, v in self.old:
            setattr(obj, n, v)
        self.M.C.SIZES = self.sizes
        return False


def places(M, W, plan, xmin_t, declared_t):
    """for every mark of a plan: the places of its number and the extent of its glyph, by the rule of the program of
    the lanes, with the sizes that the two modules hold when this is called"""
    per_day, width = plan["per_day_pt"], plan["geometry"]["axes_width_pt"]
    XMIN = M.D(xmin_t)
    x_decl = (M.D(declared_t) - XMIN) * per_day
    size = M._measure()
    out, glyphs = [], []
    for e in plan["events"]:
        if e["mark"] in ("triangle", "cross"):
            g0, g1 = e["x_pt"] - M.MARK_PT / 2.0, e["x_pt"] + M.MARK_PT / 2.0
        else:
            a, b = (e["x0"] - XMIN) * per_day, (e["x1"] - XMIN) * per_day
            g0 = a - (M.ARROW_PT / 2.0 if e["cut_left"] else 0.0)
            g1 = b + (M.ARROW_PT / 2.0 if e["cut_right"] else 0.0)
        w = size(str(e["number"]))[0]
        c = []
        for side in ("right", "left"):
            n0 = g1 + W.NUMBER_PAD_PT if side == "right" else g0 - W.NUMBER_PAD_PT - w
            n1 = n0 + w
            if n0 < -1e-9 or n1 > width + 1e-9:
                continue
            if n0 - M.LINE_CLEAR_PT < x_decl < n1 + M.LINE_CLEAR_PT:
                continue
            c.append((side, n0, n1, min(g0, n0), max(g1, n1)))
        if not c:
            raise M.C.Stop("the number of mark %s finds no place beside its mark" % e["number"])
        out.append(c)
        glyphs.append((g0, g1))
    return out, glyphs


def install_search(M, W, max_lanes, state, strip_scale=1.0):
    C = M.C
    in_lanes = M.plan_of_the_strip                    # after W.install(M): the plan in lanes
    drawn_in_lanes = M.draw_strip                     # after W.install(M): the drawing in lanes

    def plan_by_search(mark_rows, number_rows, xmin_t, xmax_t, declared_t):
        plan = in_lanes(mark_rows, number_rows, xmin_t, xmax_t, declared_t)
        ev = plan["events"]
        state["lanes_of_the_plan_of_1403"] = [(e["number"], e["lane"]["lane"], e["lane"]["side"]) for e in ev]
        full, glyphs_full = places(M, W, plan, xmin_t, declared_t)
        for e, c, (g0, g1) in zip(ev, full, glyphs_full):
            la = e["lane"]
            if abs(g0 - la["g0"]) > 1e-9 or abs(g1 - la["g1"]) > 1e-9 or not any(
                    s == la["side"] and abs(n0 - la["n0"]) < 1e-9 and abs(n1 - la["n1"]) < 1e-9 for s, n0, n1, x0, x1 in c):
                raise C.Stop("the places of mark %s at full size are not those of the plan in lanes" % e["number"])
        with scaled_strip(M, W, strip_scale):
            cands, glyphs = places(M, W, plan, xmin_t, declared_t)
            ink = W.ink_height_of_the_digits(C)
            body = max(M.MARK_PT, ink)
            gap, lane_gap = W.EXTENT_GAP_PT, W.LANE_GAP_PT
            state["sizes"] = {"scale": float(strip_scale), "mark_pt": M.MARK_PT, "cross_lw_pt": M.CROSS_LW_PT,
                              "span_lw_pt": M.SPAN_LW_PT, "arrow_pt": M.ARROW_PT, "bar_lw_pt": W.BAR_LW_PT,
                              "number_pad_pt": W.NUMBER_PAD_PT, "extent_gap_pt": gap, "lane_gap_pt": lane_gap,
                              "size_of_the_numbers_pt": C.SIZES[2], "ink_height_of_a_number_pt": ink, "body_pt": body}
        state["places"], state["gap"] = cands, gap
        best, n_l = None, 0
        while best is None and n_l < max_lanes:
            n_l += 1
            best = search(cands, n_l, gap)
        state["searched_lanes"] = n_l
        if best is None:
            raise C.Stop("no assignment of the %d marks to at most %d lanes" % (len(ev), max_lanes))
        for e, (lane, k), c, (g0, g1) in zip(ev, best, cands, glyphs):
            side, n0, n1, x0, x1 = c[k]
            e["lane"] = {"lane": lane, "side": side, "n0": n0, "n1": n1, "x0": x0, "x1": x1, "g0": g0, "g1": g1}
        n_lanes = 1 + max(l for l, _ in best)
        lane_h = body + lane_gap
        top = n_lanes * lane_h - lane_gap
        plan["lanes"] = {"lanes": n_lanes, "body_pt": body, "lane_height_pt": lane_h, "ink_height_of_a_number_pt": ink}
        plan["why"] = plan["why"] + "; the marks stand in at most %d lanes" % max_lanes
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
        W.STATE["plan"] = plan
        return plan

    def draw_scaled(fig, axs, plan, XMIN, XMAX, DECLARED):
        with scaled_strip(M, W, strip_scale):
            return drawn_in_lanes(fig, axs, plan, XMIN, XMAX, DECLARED)

    M.plan_of_the_strip = plan_by_search
    M.draw_strip = draw_scaled
    return M


def beast_rows(path, analyses, stop):
    """analyses: [(analysis, panel of Figure S2), ...]; the rows that the file of values of Figure S2 marks as drawn"""
    rows = list(csv.DictReader(open(path, newline="", encoding="utf-8")))
    out = []
    for key, panel in analyses:
        g = [r for r in rows if r.get("analysis") == key and r["panel"] == panel and r["element"].startswith(E_NE)
             and not r["drawn_as"].lstrip().lower().startswith("not drawn")]
        g.sort(key=lambda r: r["x"])
        if not g:
            raise stop("file of values of Figure S2: no drawn row of %s in panel %s" % (key, panel))
        for a, b in zip(g[:-1], g[1:]):
            if a["x_end"] != b["x"]:
                raise stop("file of values of Figure S2: the intervals of %s are not contiguous" % key)
        out.append({"analysis": key, "panel": panel, "rows": g})
    return out


def add_beast(M, fig, sets):
    """draws into panel a the medians (stair lines) and the bounds (horizontal lines over each interval, without
    joins) of the analyses with BEAST X; a bound below the lowest power of ten of the axis is not drawn"""
    C = M.C
    axa = fig.interventions["axa"]
    y0, y1 = axa.get_ylim()
    x0, x1 = axa.get_xlim()
    floor = C.decade_ceil(y0)
    drawn = []
    for n, s in enumerate(sets):
        g = s["rows"]
        edges = [M.D(r["x"]) for r in g] + [M.D(g[-1]["x_end"])]
        vals = {c: [float(r[c]) for r in g] for c in ("value", "lower", "upper")}
        if not (x0 <= edges[0] and edges[-1] <= x1):
            raise C.Stop("an interval of %s lies outside the time axis" % s["analysis"])
        if not (floor <= min(vals["value"]) and max(vals["value"]) <= y1):
            raise C.Stop("a median of %s lies outside the axis of the full figure or in the room of the root date" % s["analysis"])
        if max(vals["upper"] + vals["lower"]) > y1:
            raise C.Stop("a bound of %s lies above the axis of the full figure" % s["analysis"])
        arts = {"value": axa.stairs(vals["value"], edges, baseline=None, color=C.GREEN, lw=BEAST_MEDIAN_LW,
                                    ls=style(n, BEAST_MEDIAN_LW), zorder=3.6)}
        cut = {}
        for col in ("lower", "upper"):
            xs, ys, out = [], [], []
            for i, v in enumerate(vals[col]):
                if v < floor:
                    out.append(i)
                    continue
                xs += [edges[i], edges[i + 1], float("nan")]
                ys += [v, v, float("nan")]
            arts[col], = axa.plot(xs, ys, color=C.GREEN, lw=BEAST_BOUND_LW, ls=style(n, BEAST_BOUND_LW), solid_capstyle="butt",
                                  dash_capstyle="butt", zorder=3.4, scalex=False, scaley=False)
            cut[col] = out
        drawn.append({"analysis": s["analysis"], "edges": edges, "values": vals, "artists": arts, "not_drawn": cut,
                      "days": [(r["x"], r["x_end"]) for r in g]})
    if (y0, y1) != axa.get_ylim() or (x0, x1) != axa.get_xlim():
        raise C.Stop("the limits of panel a changed while the lines of BEAST X were drawn")
    return axa, floor, drawn


def key_above(M, fig, axa, labels):
    """ONE key above the frame of panel a, to the right of the title of the panel, in two columns: the two entries of
    the existing key with its symbols, and one entry for each analysis with BEAST X. The existing key inside
    the frame is removed. The lower edge of the key stands KEY_MARGIN_PT above the upper edge of the frame, its right
    edge at the right end of the tick labels of the right-hand axis. Where the key with its margin does not fit below
    the upper edge of the picture, the upper edge of the frame is lowered by whole points, and the letter and the
    title of the panel keep their place. Stops if, in the picture drawn without the key, a pixel that is not white
    stands inside the box of the key widened by the margin."""
    import numpy as np
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch
    from matplotlib.text import Annotation
    C = M.C
    old = axa.get_legend()
    labels_of_the_track = [t.get_text() for t in old.get_texts()]
    if len(labels_of_the_track) != 2 or len(labels) != 4:
        raise C.Stop("the key of panel a of the full figure has %d entries (2 expected); %d labels are given (4 expected)"
                     % (len(labels_of_the_track), len(labels)))
    old.remove()
    twins = [a for a in fig.axes if a is not axa and a.get_position().bounds == axa.get_position().bounds]
    heads = [t for t in axa.texts if isinstance(t, Annotation) and tuple(t.xy) == (0, 1)
             and t.xycoords == "axes fraction"]
    if len(heads) != 2:
        raise C.Stop("panel a has %d annotations at its upper left corner, 2 expected (letter and title)" % len(heads))
    handles = [
        (Line2D([], [], color=C.BLUE, lw=1.7), Patch(facecolor=C.blend(C.BLUE, C.BLUE_ALPHA), edgecolor="none")),
        (Line2D([], [], color=C.blend(C.BLUE, 0.55), lw=1.4, ls=(0, (2.5, 1.5))),
         Patch(facecolor="none", hatch="//////", edgecolor=C.blend(C.BLUE, 0.45), linewidth=0)),
        Line2D([], [], color=C.GREEN, lw=BEAST_MEDIAN_LW, ls=style(0, BEAST_MEDIAN_LW)),
        Line2D([], [], color=C.GREEN, lw=BEAST_MEDIAN_LW, ls=style(1, BEAST_MEDIAN_LW)),
    ]
    k = fig.dpi / 72.0
    m = KEY_MARGIN_PT * k

    def place():
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()        # the picture without the key
        hpx, wpx = buf.shape[:2]
        ab = axa.get_window_extent()
        right = max([t.get_window_extent().x1 for a in twins for t in a.get_yticklabels()
                     if t.get_text().strip() and t.get_visible()] + [ab.x1])
        leg = axa.legend(handles, list(labels), ncols=2, loc="lower right",
                         bbox_to_anchor=(right / wpx, (ab.y1 + m) / hpx), bbox_transform=fig.transFigure,
                         fontsize=C.SIZES[2], handlelength=2.2, borderaxespad=0.0, borderpad=0.0, labelspacing=0.3,
                         columnspacing=1.5, frameon=False, handler_map={tuple: M.HandlerOver()})
        fig.canvas.draw()
        return leg, leg.get_window_extent(), buf, hpx, wpx, ab

    leg, bb, buf, hpx, wpx, ab = place()
    head_boxes = [t.get_window_extent() for t in heads]
    lowered = 0
    over = bb.y1 + m - hpx
    if over > 1e-6:
        lowered = int(math.ceil(over / k - 1e-9))
        leg.remove()
        d = lowered / 72.0 / fig.get_figheight()
        for a in [axa] + twins:
            p = a.get_position()
            a.set_position([p.x0, p.y0, p.width, p.height - d])
        for t in heads:
            x, y = t.xyann
            t.xyann = (x, y + lowered)
        leg, bb, buf, hpx, wpx, ab = place()
        moved = max(abs(a.y0 - b.y0) + abs(a.x0 - b.x0) for a, b in zip(head_boxes, [t.get_window_extent() for t in heads]))
        if bb.y1 + m - hpx > 1e-6 or moved > 0.51:
            raise C.Stop("after the frame was lowered by %d pt the key does not fit, or the title has moved (%.2f px)"
                         % (lowered, moved))
    r0, r1 = int(math.floor(hpx - bb.y1 - m)), int(math.ceil(hpx - bb.y0 + m))
    c0, c1 = int(math.floor(bb.x0 - m)), int(math.ceil(bb.x1 + m))
    if r0 < 0 or c0 < 0 or c1 > wpx or bb.y0 < ab.y1:
        raise C.Stop("the key above panel a does not lie inside the picture and above the frame")
    ink = int((buf[r0:r1, c0:c1].min(axis=2) < WHITE_FROM).sum())
    if ink:
        raise C.Stop("in the picture without the key, %d pixels that are not white stand in the box of the key widened by %g pt"
                     % (ink, KEY_MARGIN_PT))
    s = C.DPI / fig.dpi
    title_right = max(b.x1 for b in [t.get_window_extent() for t in heads])
    info = {"labels": list(labels), "labels_of_the_track": labels_of_the_track, "size_pt": C.SIZES[2],
            "margin_pt": KEY_MARGIN_PT, "the_frame_is_lowered_by_pt": lowered, "dpi_of_the_check": float(fig.dpi),
            "box_px": [bb.x0 * s, (hpx - bb.y1) * s, bb.x1 * s, (hpx - bb.y0) * s],              # left, top, right, bottom
            "width_pt": bb.width / k, "height_pt": bb.height / k, "from_the_title_pt": (bb.x0 - title_right) / k,
            "from_the_upper_edge_pt": (hpx - bb.y1) / k, "above_the_frame_pt": (bb.y0 - ab.y1) / k,
            "title_box_px": [min(b.x0 for b in head_boxes) * s, (hpx - max(b.y1 for b in head_boxes)) * s,
                             max(b.x1 for b in head_boxes) * s, (hpx - min(b.y0 for b in head_boxes)) * s]}
    return leg, info


def _shapes(fig):
    """the lines and the band of the three analyses in panel a, in display units: (name, x0, y0, x1, y1); a line has
    no area"""
    S = fig.short_version
    T = S["axa"].transData.transform
    shapes = []

    def flat(name, xa_, xb_, v):
        (xa, ya), (xb, _) = T((xa_, v)), T((xb_, v))
        shapes.append((name, xa, ya, xb, ya))

    def stair(name, edges, vals):
        for i, v in enumerate(vals):
            flat(name, edges[i], edges[i + 1], v)
            if i + 1 < len(vals):
                (xb, ya), (_, yb) = T((edges[i + 1], v)), T((edges[i + 1], vals[i + 1]))
                shapes.append((name, xb, min(ya, yb), xb, max(ya, yb)))

    d = sorted(S["delphy"], key=lambda r: r["x"])
    ed = [r["x"] for r in d] + [d[-1]["x_end"]]
    stair("median of the primary analysis", ed, [r["value"] for r in d])
    for r in d:
        (xa, ya), (xb, yb) = T((r["x"], r["lower"])), T((r["x_end"], r["upper"]))
        shapes.append(("band of the primary analysis", xa, ya, xb, yb))
    for s in S["beast"]:
        stair("median of %s" % s["analysis"], s["edges"], s["values"]["value"])
        for col, word in (("lower", "lower bound"), ("upper", "upper bound")):
            for i, v in enumerate(s["values"][col]):
                if i not in s["not_drawn"][col]:
                    flat("%s of %s" % (word, s["analysis"]), s["edges"][i], s["edges"][i + 1], v)
    return shapes


def _hit(shapes, box, pad, ab, pad_x=None):
    """the names of the shapes that, cut at the frame, come within pad of the box (x0, y0, x1, y1); pad_x: another
    distance to the left and to the right"""
    px = pad if pad_x is None else pad_x
    return sorted(set(name for name, xa, ya, xb, yb in shapes
                      if max(xa, ab.x0) <= box[2] + px and min(xb, ab.x1) >= box[0] - px
                      and max(ya, ab.y0) <= box[3] + pad and min(yb, ab.y1) >= box[1] - pad
                      and max(xa, ab.x0) <= min(xb, ab.x1) and max(ya, ab.y0) <= min(yb, ab.y1)))


def _words(fig, rend, ab):
    """the visible words of panel a whose box reaches into the frame"""
    out = []
    for t in fig.short_version["axa"].texts:
        if t.get_text().strip() and t.get_visible():
            b = t.get_window_extent(rend)
            if b.x1 > ab.x0 and b.x0 < ab.x1 and b.y1 > ab.y0 and b.y0 < ab.y1:
                out.append((t, b))
    return out


def clear_the_words(M, fig, pad_pt=1.0, clear_pt=3.0, clear_x_pt=24.0, step_pt=1.0, most_pt=150.0):
    """A word of the original figure inside panel a that touches a line or the band of an analysis (pad_pt) is moved to the
    LEFT, in steps of step_pt, to the first place at which no line or band comes within clear_pt of its box above or
    below nor within clear_x_pt to its left or right at the height of the box, and at which, in the picture drawn
    without the word, no pixel that is not white stands in its box widened by pad_pt."""
    import numpy as np
    C = M.C
    axa = fig.short_version["axa"]
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    ab = axa.get_window_extent(rend)
    k = fig.dpi / 72.0
    s = C.DPI / fig.dpi
    shapes = _shapes(fig)
    moved = []
    for t, b in _words(fig, rend, ab):
        hit = _hit(shapes, (b.x0, b.y0, b.x1, b.y1), pad_pt * k, ab)
        if not hit:
            continue
        t.set_visible(False)
        fig.canvas.draw()
        buf = np.asarray(fig.canvas.buffer_rgba())[..., :3].copy()        # the picture without the word
        t.set_visible(True)
        hpx = buf.shape[0]
        found, n = None, 0
        while found is None and (n + 1) * step_pt <= most_pt:
            n += 1
            dx = n * step_pt * k
            box = (b.x0 - dx, b.y0, b.x1 - dx, b.y1)
            if box[0] - clear_pt * k < ab.x0:
                break
            if _hit(shapes, box, clear_pt * k, ab, pad_x=clear_x_pt * k):
                continue
            r0, r1 = int(math.floor(hpx - box[3] - pad_pt * k)), int(math.ceil(hpx - box[1] + pad_pt * k))
            c0, c1 = int(math.floor(box[0] - pad_pt * k)), int(math.ceil(box[2] + pad_pt * k))
            if (buf[r0:r1, c0:c1].min(axis=2) < WHITE_FROM).any():
                continue
            found = n
        if found is None:
            raise C.Stop("the word %r touches %s and finds no free place within %g pt to its left"
                         % (t.get_text(), hit, most_pt))
        tr = t.get_transform()
        p = tr.transform(t.get_position())
        q = tr.inverted().transform((p[0] - found * step_pt * k, p[1]))
        t.set_position((float(q[0]), float(q[1])))
        fig.canvas.draw()
        b2 = t.get_window_extent(rend)
        if abs((b.x0 - b2.x0) - found * step_pt * k) > 0.51 or abs(b.y0 - b2.y0) > 0.51:
            raise C.Stop("the word %r does not stand at its planned place" % t.get_text())
        moved.append({"word": t.get_text(), "touched": hit, "moved_to_the_left_pt": found * step_pt,
                      "box_px_before": [b.x0 * s, (hpx - b.y1) * s, b.x1 * s, (hpx - b.y0) * s],
                      "box_px_after": [b2.x0 * s, (hpx - b2.y1) * s, b2.x1 * s, (hpx - b2.y0) * s]})
    return moved


def touches(fig, pad_pt=1.0):
    """the key and the words inside panel a that touch a line of an analysis or the band of the primary analysis:
    [{"what": ..., "touches": [...]}, ...] for all of them (an empty list of touches: free)"""
    S = fig.short_version
    axa = S["axa"]
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    ab = axa.get_window_extent(rend)
    pad = pad_pt * fig.dpi / 72.0
    shapes = _shapes(fig)
    kb = S["legend"].get_window_extent(rend)
    out = [{"what": "the key", "touches": _hit(shapes, (kb.x0, kb.y0, kb.x1, kb.y1), pad, ab)}]
    for t, b in _words(fig, rend, ab):
        out.append({"what": t.get_text(), "touches": _hit(shapes, (b.x0, b.y0, b.x1, b.y1), pad, ab)})
    return out


class wide_figure:
    """while it is open, the figure is drawn width_in inches wide: the width of the conventions, the width with which
    its function makes a figure, the left and the right margin of the grid of the original figure, which keep their inches, and
    (w') the place of the names of the vertical axes, which keep their distance in inches from the frame.
    width_in None: nothing is changed."""

    def __init__(self, M, width_in):
        self.M, self.w = M, (None if width_in is None else float(width_in))
        self.info = None

    def __enter__(self):
        M, C = self.M, self.M.C
        self.old = (C.FIG_WIDTH, C.new_figure, M.GRID, M.YLAB_X)
        w0, make = C.FIG_WIDTH, C.new_figure
        if self.w is not None:
            w = self.w
            axes0 = w0 * (M.GRID["right"] - M.GRID["left"])
            C.FIG_WIDTH = w
            C.new_figure = lambda height, width=None: make(height, w if width is None else width)
            M.GRID = dict(M.GRID, left=M.GRID["left"] * w0 / w, right=1.0 - (1.0 - M.GRID["right"]) * w0 / w)
            M.YLAB_X = M.YLAB_X * axes0 / (w * (M.GRID["right"] - M.GRID["left"]))
        self.info = {"width_in": float(C.FIG_WIDTH), "width_of_the_conventions_in": float(w0),
                     "grid_left": float(M.GRID["left"]), "grid_right": float(M.GRID["right"]),
                     "left_margin_in": float(M.GRID["left"] * C.FIG_WIDTH),
                     "right_margin_in": float((1.0 - M.GRID["right"]) * C.FIG_WIDTH), "ylab_x": float(M.YLAB_X),
                     "names_of_the_vertical_axes_left_of_the_frame_in": float(
                         -M.YLAB_X * C.FIG_WIDTH * (M.GRID["right"] - M.GRID["left"]))}
        return self

    def __exit__(self, *a):
        self.M.C.FIG_WIDTH, self.M.C.new_figure, self.M.GRID, self.M.YLAB_X = self.old
        return False


def artists_of(fig, strip_axes=()):
    """what the axes of a figure hold, for a comparison of two drawings: for every axes its limits, scales, names and
    tick labels, and for every line, patch, collection, word and key its kind, its looks, its numbers and, for the two
    directions, whether the numbers are in the units of the data"""
    import numpy as np
    from matplotlib import colors as mc

    def col(c):
        try:
            a = np.atleast_2d(mc.to_rgba_array(c))
            return [mc.to_hex(x, keep_alpha=True) for x in a[:4]] + ([len(a)] if len(a) > 4 else [])
        except Exception:
            return [str(c)]

    def arr(a):
        a = np.asarray(a, dtype=float).reshape(-1, 2)
        return [[float(v) if np.isfinite(v) else None for v in row] for row in a]

    def unit(tr, ax):
        try:
            x, y = tr.contains_branch_seperately(ax.transData)
            return [bool(x), bool(y)]
        except Exception:
            return [False, False]

    out = []
    for n, ax in enumerate(fig.axes):
        d = {"axes": n, "strip": bool(any(ax is s for s in strip_axes if s is not None)),
             "xlim": [float(v) for v in ax.get_xlim()], "ylim": [float(v) for v in ax.get_ylim()],
             "xscale": ax.get_xscale(), "yscale": ax.get_yscale(), "xlabel": ax.get_xlabel(), "ylabel": ax.get_ylabel(),
             "xticklabels": [t.get_text() for t in ax.get_xticklabels() if t.get_visible()],
             "yticklabels": [t.get_text() for t in ax.get_yticklabels() if t.get_visible()], "held": []}
        for a in ax.lines:
            d["held"].append({"kind": "line", "unit": unit(a.get_transform(), ax), "numbers": arr(a.get_xydata()),
                              "looks": [col(a.get_color()), float(a.get_linewidth()), str(a.get_linestyle()),
                                        str(a.get_marker()), float(a.get_markersize()), col(a.get_markerfacecolor()),
                                        float(a.get_zorder()), bool(a.get_visible())]})
        for a in ax.patches:
            v = a.get_patch_transform().transform(a.get_path().vertices)
            d["held"].append({"kind": type(a).__name__, "unit": unit(a.get_data_transform(), ax), "numbers": arr(v),
                              "looks": [col(a.get_facecolor()), col(a.get_edgecolor()), float(a.get_linewidth()),
                                        str(a.get_linestyle()), str(a.get_hatch()), float(a.get_zorder()),
                                        bool(a.get_visible())]})
        for a in ax.collections:
            d["held"].append({"kind": type(a).__name__, "unit": unit(a.get_transform(), ax),
                              "numbers": [x for p in a.get_paths() for x in arr(p.vertices) + [[None, None]]],
                              "offsets": arr(a.get_offsets()), "offset_unit": unit(a.get_offset_transform(), ax),
                              "looks": [col(a.get_facecolor()), col(a.get_edgecolor()),
                                        [float(x) for x in np.atleast_1d(a.get_linewidth())][:4], str(a.get_hatch()),
                                        float(a.get_zorder()), bool(a.get_visible())]})
        for a in ax.texts:
            note = hasattr(a, "xyann")
            d["held"].append({"kind": "word", "text": a.get_text(),
                              "unit": ([str(a.xycoords) == "data"] * 2) if note else unit(a.get_transform(), ax),
                              "numbers": arr([a.xy]) if note else arr([a.get_position()]),
                              "set_off": arr([a.xyann]) if note else [],
                              "looks": [col(a.get_color()), float(a.get_fontsize()), str(a.get_fontweight()),
                                        str(a.get_ha()), str(a.get_va()), float(a.get_rotation()),
                                        bool(a.get_visible())]})
        leg = ax.get_legend()
        if leg is not None:
            d["held"].append({"kind": "key", "text": [t.get_text() for t in leg.get_texts()], "unit": [False, False],
                              "numbers": [], "looks": [[float(t.get_fontsize()) for t in leg.get_texts()]]})
        out.append(d)
    return out


def run(code_dir, lanes_path, values_path, s2_values_path, outdir, stem, analyses, labels, max_lanes=MAX_LANES,
        strip_scale=1.0, fig_width=None):
    W = load_lanes(lanes_path)
    M = W.install(W.load(code_dir))
    with wide_figure(M, fig_width) as wide:
        return _run(W, M, wide, values_path, s2_values_path, outdir, stem, analyses, labels, max_lanes, strip_scale)


def _run(W, M, wide, values_path, s2_values_path, outdir, stem, analyses, labels, max_lanes, strip_scale):
    state = {}
    M = install_search(M, W, max_lanes, state, strip_scale)
    os.makedirs(outdir, exist_ok=True)
    sets = beast_rows(s2_values_path, analyses, M.C.Stop)
    fig, extra = M.draw(values_path)
    axa, floor, drawn = add_beast(M, fig, sets)
    delphy = [r for r in M.C.read_csv(values_path) if r["element"] == M.E_STAIR]
    fig.short_version = {"axa": axa, "beast": drawn, "legend": None,
                         "delphy": [{"x": M.D(r["x"]), "x_end": M.D(r["x_end"]), "value": float(r["value"]),
                                     "lower": float(r["lower"]), "upper": float(r["upper"])} for r in delphy]}
    moved = clear_the_words(M, fig)
    leg, key = key_above(M, fig, axa, labels)
    fig.short_version["legend"] = leg
    over = touches(fig)
    plan, axs = W.STATE["plan"], W.STATE["axs"]
    boxes = W.boxes_in_pixels(fig, axs, plan, M.C.DPI)
    pos = axa.get_position()
    w_in, h_in = fig.get_size_inches()
    dpi = M.C.DPI
    panel_a = {"left": pos.x0 * w_in * dpi, "right": pos.x1 * w_in * dpi, "top": (1 - pos.y1) * h_in * dpi,
               "bottom": (1 - pos.y0) * h_in * dpi, "xlim": list(axa.get_xlim()), "ylim": list(axa.get_ylim()),
               "yscale": axa.get_yscale(), "height_px": h_in * dpi, "width_px": w_in * dpi,
               "lowest_power_of_ten": floor}
    back = []
    for s in drawn:
        v, e, _ = s["artists"]["value"].get_data()
        d = {"analysis": s["analysis"], "days": s["days"],
             "value": {"values": [float(x) for x in v], "edges": [float(x) for x in e],
                       "colour": s["artists"]["value"].get_edgecolor(),
                       "linewidth": float(s["artists"]["value"].get_linewidth()),
                       "linestyle": str(s["artists"]["value"].get_linestyle())},
             "not_drawn": {c: [s["days"][i] for i in s["not_drawn"][c]] for c in ("lower", "upper")}}
        for col in ("lower", "upper"):
            a = s["artists"][col]
            xs, ys = [float(x) for x in a.get_xdata()], [float(y) for y in a.get_ydata()]
            d[col] = {"segments": [[xs[i], xs[i + 1], ys[i]] for i in range(0, len(xs), 3)], "colour": a.get_color(),
                      "linewidth": float(a.get_linewidth()), "linestyle": str(a.get_linestyle()),
                      "flat": bool(all(ys[i] == ys[i + 1] for i in range(0, len(ys), 3)))}
        back.append(d)
    held = artists_of(fig, [axs])
    saved = M.C.save_figure(fig, stem, outdir)
    M.plt_close(fig)
    return {"saved": {k: v for k, v in saved.items()}, "figure": wide.info, "held": held, "boxes": boxes,
            "panel_a": panel_a, "artists": back,
            "key": key, "words_moved": moved, "touches": over,
            "keys_or_words_over_data_by_the_track": extra.get("keys_or_words_over_data"),
            "search": {"lanes_of_the_plan_of_1403": state["lanes_of_the_plan_of_1403"],
                       "lanes_searched_until_an_assignment": state["searched_lanes"],
                       "places_per_mark": [len(c) for c in state["places"]],
                       "one_lane_is_possible": bool(search(state["places"], 1, state["gap"]) is not None)},
            "strip": state["sizes"],
            "plan": {"rule": plan["rule"], "lanes": plan["lanes"], "strip_height_pt": plan["strip_height_pt"],
                     "stack_pt": plan["stack_pt"], "layout": plan["layout"], "room_total_pt": plan["room_total_pt"],
                     "events": [{"number": e["number"], "mark": e["mark"], "first": e["first"], "last": e["last"],
                                 "lane": e["lane"]["lane"], "side": e["lane"]["side"]} for e in plan["events"]]},
            "height_in": extra.get("height_in"),
            "date_number_of": {d: M.D(d) for s in drawn for pair in s["days"] for d in pair}}
