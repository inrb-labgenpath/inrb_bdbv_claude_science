"""fig_program_trajectories_20261006.py

Figure S2 of the supplement (effective population size over time, Delphy and BEAST X, pair by pair), drawn
from the table of the trajectories of the comparison of the two programs of the round of October 2026
(program_comparison_trajectories_20261002.csv).

    python fig_program_trajectories_20261006.py INPUTS_SUPPLEMENT.json INPUTS_THIRD_TASK.json [OUTDIR]
    python fig_program_trajectories_20261006.py --draw-only VALUES.csv [OUTDIR]

The rules of drawing are those of fig_program_trajectories_20260926.py (the earlier generator, which is NOT
run; its typed constants are read from its source), with the changes of this round:
  panels     four, as the figure stands: the pairs A, B_estimated, C680 and C515 of the table, in its order;
             the pair 'B' of the table is not drawn
  summaries  pair C680: the summary cut at the saved states is drawn as thinner lines, in the lighter tone
             that Figure S1 gives it; pair A: the table holds one row with the state 'not complete' and no
             trajectory of the cut summary: nothing is drawn for it
  dagger     behind 'BEAST X' in the panel of an analysis whose chains were continued from a saved state
  R4 (f)     the cell meets_convergence_criterion is copied; an analysis with a cell that begins with 'no'
             would be marked with a double dagger behind its name in the panel
  axes       R4 (h): limits by a rule from the values of this round
  titles     built from cells of the table of values of the comparison
The file of values is written BEFORE the figure is drawn, and the figure is drawn from it and from nothing
else.
"""
import datetime as dt
import math
import os
import re
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import supp_fig_common_20261006 as S    # noqa: E402

C = S.C
import matplotlib as mpl               # noqa: E402
import matplotlib.dates as mdates      # noqa: E402
from matplotlib.lines import Line2D    # noqa: E402
from matplotlib.patches import Patch   # noqa: E402

NAME = "fig_program_trajectories" + S.SUFFIX
VALUES_FILE = NAME + "_values.csv"
GENERATOR = "fig_program_trajectories_20260926.py"
TABLE = "program_comparison_trajectories_20261002.csv"
FIRST = S.FIRST_SUMMARY

E_STEP = "effective population size of an interval (median and 95 % HPD interval)"
E_DECLARED = "declaration of the outbreak"
E_AXIS_FIRST = "time axis, first day"
E_AXIS_END = "time axis, end"
E_NAME = "name of a program in a panel"
PANELS = "abcd"
TONE_OF_THE_CUT = 0.45              # the lighter tone of Figure S1 of this round: share of the colour of BEAST X
LW = {"Delphy": 1.6, "BEAST X": 1.6, "bound": 0.7, "cut": 0.9, "cut bound": 0.45}   # earlier generator, lines 75, 79, 80; cut: PROPOSAL
TITLE_PATTERN = ("Setting {letter}: {words}", 41)      # the form of the titles of the earlier generator, with its line
KEY_WORDS = {"Delphy": ("Delphy: median and 95 % HPD interval", 93),
             "BEAST X": ("BEAST X: median (thick line) and bounds of the 95 % HPD interval (thin lines)", 94)}
WORD_DECLARED = ("declaration", 92)
Y_LABEL = (r"$N_e\tau$ (years)", 88)
CORNERS = ("upper left", "lower right", "upper right", "lower left")      # places of the names in a panel, in this order (PROPOSAL)
LINE_IN = 9.6 / 72.0 * 10.0 / 8.0   # height of a line of a title of 10 pt


def build_values(I):
    T = I.table(TABLE).need("pair", "program", "analysis", "genomes", "chains_used", "states_used", "n_parameters",
                            "interval", "older_edge", "newer_edge", "Ne_tau_years_median", "Ne_tau_years_hpd_lo",
                            "Ne_tau_years_hpd_hi", "state_of_the_table", "meets_convergence_criterion", "source_table",
                            "source_version", "note_on_the_analysis", "summary_of_the_analysis_of_BEAST_X", "note_on_the_pair",
                            "last_collection_date")
    V = C.Values(NAME, I)
    k = S.constants_of(I, GENERATOR, ["XMIN", "XMAX", "DECLARED", "YMIN0", "YMAX0", "PW", "PH", "LEFT", "GAPX", "TOP", "GAPY", "BOTTOM"])
    lines = I.lines(GENERATOR)
    for words, line in list(KEY_WORDS.values()) + [WORD_DECLARED, (TITLE_PATTERN[0].split("{")[0], TITLE_PATTERN[1])]:
        if words not in lines[line - 1]:
            raise C.Stop(f"{GENERATOR}: line {line} does not hold the words {words!r}")
    if Y_LABEL[0].replace("\\", "\\\\") not in lines[Y_LABEL[1] - 1] and Y_LABEL[0] not in lines[Y_LABEL[1] - 1]:
        raise C.Stop(f"{GENERATOR}: line {Y_LABEL[1]} does not hold the label of the axis")
    cont = S.continued_chains(I, T, TABLE)
    words = S.words_of_the_pairs(I)
    pairs = []
    for r in T.rows:
        if r["pair"] not in pairs:
            pairs.append(r["pair"])
    drawn_pairs = [p for p in pairs if p not in S.PAIRS_NOT_DRAWN]
    if len(drawn_pairs) != len(PANELS):
        raise C.Stop(f"{TABLE}: {len(drawn_pairs)} pairs to draw; the figure has {len(PANELS)} panels")
    panel_of = dict(zip(drawn_pairs, PANELS))
    # ---- time axis and declaration
    first_txt, first_line, _ = C.pattern_in_document(I, C.CONVENTIONS, r"Time axes of time-series panels: (\d+ \w+ \d{4}) to", "first day of the time axis")
    XMIN = C.day_of_words(first_txt, "first day of the time axis")
    decl_txt, decl_line, _ = C.pattern_in_document(I, C.CONVENTIONS, r"Declaration of the outbreak: (\d+ \w+ \d{4})", "declaration of the outbreak")
    DECLARED = C.day_of_words(decl_txt, "declaration")
    # ---- every row of the table: drawn, or not drawn with its reason
    units = {}
    for i, r in enumerate(T.rows):
        key = (r["pair"], r["summary_of_the_analysis_of_BEAST_X"], r["program"])
        units.setdefault(key, r["analysis"])
        if units[key] != r["analysis"]:
            raise C.Stop(f"{TABLE}: {key} holds two analyses")
    states, drawn = {}, []
    for i, r in enumerate(T.rows):
        p, prog, s = r["pair"], r["program"], r["summary_of_the_analysis_of_BEAST_X"]
        reason = ""
        complete = r["state_of_the_table"] == "complete"
        if p in S.PAIRS_NOT_DRAWN:
            reason = S.PAIRS_NOT_DRAWN[p]
        elif not complete:
            reason = f"the row holds no trajectory; its cell state_of_the_table reads '{r['state_of_the_table']}'"
        elif prog == "Delphy" and s != FIRST:
            reason = ("the row of Delphy stands in the table once for each summary of BEAST X of the pair; it is drawn once, "
                      "with the summary of the analysis")
        elif not C.is_day(r["newer_edge"]):
            raise C.Stop(f"{TABLE}: row {i}: newer_edge {r['newer_edge']!r} is no day")
        elif C.to_day(r["newer_edge"]) <= XMIN:
            reason = ("the interval ends on or before the first day of the time axis (earlier generator, line 53: such an "
                      "interval is not drawn)")
        if complete:
            for c in ("Ne_tau_years_median", "Ne_tau_years_hpd_lo", "Ne_tau_years_hpd_hi"):
                if not C.is_number(r[c]) or C.num(r[c]) <= 0:
                    raise C.Stop(f"{TABLE}: row {i}, column {c}: {r[c]!r} cannot be drawn on a logarithmic axis")
        states[i] = "not drawn: " + reason if reason else "drawn"
        panel = panel_of.get(p, "none")
        label = f"interval {int(float(r['interval']))}" if C.is_number(r["interval"]) else ""
        is_open = r["older_edge"] == "open"
        if reason:
            V.add(panel, E_STEP, label=label, value=r["Ne_tau_years_median"], lower=r["Ne_tau_years_hpd_lo"],
                  upper=r["Ne_tau_years_hpd_hi"], x="" if (is_open or not complete) else r["older_edge"],
                  x_end=r["newer_edge"], unit="years (Ne tau)" if complete else "",
                  drawn_as=f"{S.NOT_DRAWN}: {reason}", table=TABLE, row=i,
                  columns=("value=Ne_tau_years_median; lower=Ne_tau_years_hpd_lo; upper=Ne_tau_years_hpd_hi; "
                           + ("" if is_open else "x=older_edge; ") + "x_end=newer_edge") if complete else "",
                  note=(f"{prog}; {s}; pair {p} of the table"
                        + ("; the interval is open towards the past" if is_open else "")
                        + (f"; note of the pair in the table: {r['note_on_the_pair']}" if not complete else "")),
                  analysis=r["analysis"], meets=r["meets_convergence_criterion"])
            continue
        older = XMIN if is_open else C.to_day(r["older_edge"])
        cut = is_open or older < XMIN
        how = ("Delphy: line (median) and band (95 % HPD interval) in the colour of Delphy" if prog == "Delphy" else
               "BEAST X: thick line (median) between two thin lines (bounds of the 95 % HPD interval) in the colour of BEAST X"
               if s == FIRST else
               "BEAST X, second summary: thinner line (median) between two thinner lines (bounds of the 95 % HPD interval) in a "
               "lighter tone of the colour of BEAST X")
        V.add(panel, E_STEP, label=label, value=r["Ne_tau_years_median"], lower=r["Ne_tau_years_hpd_lo"],
              upper=r["Ne_tau_years_hpd_hi"], x=max(older, XMIN), x_end=r["newer_edge"], unit="years (Ne tau)",
              drawn_as=how, table=TABLE, row=i,
              columns="value=Ne_tau_years_median; lower=Ne_tau_years_hpd_lo; upper=Ne_tau_years_hpd_hi; "
                      + ("" if cut else "x=older_edge; ") + "x_end=newer_edge",
              note=f"{prog}; {s}; pair {p} of the table"
                   + (f"; the interval starts before the axis (older_edge {r['older_edge']}); drawn from the first day of the axis" if cut else ""),
              analysis=r["analysis"], meets=r["meets_convergence_criterion"])
        drawn.append(V.rows[-1])
    if not drawn:
        raise C.Stop(f"{TABLE}: no row to draw")
    # ---- axes (R4 h)
    XMAX = max(C.to_day(r["x_end"]) for r in drawn)
    last = [r for r in drawn if C.to_day(r["x_end"]) == XMAX][0]
    V.add("a to d", E_AXIS_FIRST, label="", value=first_txt, x=XMIN, unit="date", drawn_as="left end of the time axis",
          table=C.CONVENTIONS, row=first_line, columns=r"value=pattern Time axes of time-series panels: (\d+ \w+ \d{4}) to",
          note=f"R4 (h), PROPOSED rule: the first day that the conventions name for time axes (as in Figure 1 of this round); "
               f"the earlier generator holds XMIN as a typed day ({k['XMIN'][0]}, line {k['XMIN'][1]})")
    V.add("a to d", E_AXIS_END, label="", value=last["x_end"], x=XMAX, unit="date", drawn_as="right end of the time axis",
          table=TABLE, row=last["source_row_0based"], columns="value=newer_edge",
          note=f"R4 (h), PROPOSED rule: the latest newer edge of the intervals drawn; the earlier generator holds XMAX as a typed "
               f"day ({k['XMAX'][0]}, line {k['XMAX'][1]})")
    if XMIN.year != XMAX.year:
        raise C.Stop("the time axis spans two years; the label of the axis prints one")
    lo = min(C.num(r["lower"]) for r in drawn)
    hi = max(C.num(r["upper"]) for r in drawn)
    r_lo = [r for r in drawn if C.num(r["lower"]) == lo][0]
    r_hi = [r for r in drawn if C.num(r["upper"]) == hi][0]
    YMIN, YMAX = C.decade_floor(lo), C.decade_ceil(hi)
    rule = (f"R4 (h), PROPOSED rule: the smallest range between whole powers of ten that holds every bound drawn, the same in all "
            f"panels; the earlier generator holds a smallest range as typed numbers (YMIN0, YMAX0 = {k['YMIN0'][0]}, {k['YMAX0'][0]}, "
            f"line {k['YMIN0'][1]}) and widens it by whole powers of ten (lines 61 and 62)")
    V.add("a to d", S.E_AXIS_RULE, label="axis of the effective population size", x=repr(YMIN), x_end=repr(YMAX),
          unit="years (Ne tau)", drawn_as="limits of the axis (logarithmic); labels at the whole powers of ten",
          note=rule + f"; the smallest bound drawn stands in row {r_lo['source_row_0based']} of {TABLE} (column Ne_tau_years_hpd_lo), "
                      f"the largest in row {r_hi['source_row_0based']} (column Ne_tau_years_hpd_hi)")
    V.add("a to d", E_DECLARED, label=WORD_DECLARED[0], value=decl_txt, x=DECLARED, unit="date",
          drawn_as="vertical dashed line in every panel; the word in panel a", table=C.CONVENTIONS, row=decl_line,
          columns=r"value=pattern Declaration of the outbreak: (\d+ \w+ \d{4})",
          note=f"the day is read from the conventions (the earlier generator holds it as a typed day, {k['DECLARED'][0]}, line "
               f"{k['DECLARED'][1]}); the word is that of the earlier generator (line {WORD_DECLARED[1]}); the conventions name "
               "the label 'outbreak declared'")
    # ---- titles, names of the programs, key
    VT = I.table(S.VALUES_OF_THE_COMPARISON)
    for p in drawn_pairs:
        if p not in words:
            raise C.Stop(f"{S.VALUES_OF_THE_COMPARISON}: no rows of the pair {p}")
        w = words[p]
        for prog in ("Delphy", "BEAST X"):
            if w[prog + "_analysis"] != units.get((p, FIRST, prog)):
                raise C.Stop(f"pair {p}, {prog}: the two tables of the comparison name different analyses")
        V.add(panel_of[p], S.E_TITLE, label=TITLE_PATTERN[0].format(letter=w["letter"], words=w["words"]), value=w["Delphy_cell"],
              x=p, table=S.VALUES_OF_THE_COMPARISON, row=w["Delphy_row"],
              columns="value=what_differs_from_the_primary_analysis_of_Delphy; x=pair", drawn_as="title of the panel",
              note="PROPOSED, built from cells (message, C 1): the form of the titles of the earlier generator (line "
                   f"{TITLE_PATTERN[1]}: 'Setting', the letter, a colon), with the capital letters at the head of the cell pair "
                   "and the cell what_differs_from_the_primary_analysis_of_Delphy of the row of Delphy up to the first of "
                   "' (', ' that ', ' instead of '")
        for prog in ("Delphy", "BEAST X"):
            a = units[(p, FIRST, prog)]
            i0 = [i for i, r in enumerate(T.rows) if r["pair"] == p and r["program"] == prog and r["summary_of_the_analysis_of_BEAST_X"] == FIRST][0]
            mark = S.DAGGER if a in cont else ""
            no = sorted({r["meets_convergence_criterion"] for r in drawn if r["analysis"] == a and r["panel"] == panel_of[p]
                         and C.marked_by_R4f(r["meets_convergence_criterion"])})
            mark += S.DOUBLE_DAGGER if no else ""
            V.add(panel_of[p], E_NAME, label=prog + mark, value=T.rows[i0]["program"], x=p, table=TABLE, row=i0,
                  columns="value=program; x=pair", drawn_as="name of the program inside the panel, in the colour of the program",
                  note="PROPOSED (message, C 4: the label of an analysis carries the dagger wherever the figure shows the "
                       "analysis; the earlier generator names the programs in the key only): the cell program"
                       + (f"; dagger: {cont[a]['chains_continued']} of {cont[a]['chains']} chains of the analysis were continued "
                          f"from a saved state ({cont[a]['read_from']})" if a in cont else "")
                       + (f"; double dagger (R4 f): cells {no}" if no else ""),
                  analysis=a)
    seconds = []
    for r in drawn:
        if "second summary" in r["drawn_as"] and r["panel"] not in [x[0] for x in seconds]:
            seconds.append((r["panel"], T.rows[int(r["source_row_0based"])]["summary_of_the_analysis_of_BEAST_X"], int(r["source_row_0based"])))
    for what, (kw, line) in KEY_WORDS.items():
        V.add("key", S.E_KEY, label=kw, drawn_as={"Delphy": "band and line in the colour of Delphy",
                                                   "BEAST X": "thick line in the colour of BEAST X"}[what],
              note=f"carried over (earlier generator, line {line})")
    for panel, s, i0 in seconds:
        V.add("key", S.E_KEY, label=f"BEAST X: {s} (panel {panel})", value=s, x=panel, table=TABLE, row=i0,
              columns="value=summary_of_the_analysis_of_BEAST_X", drawn_as="thinner line in a lighter tone of the colour of BEAST X",
              note="PROPOSED, built from the cell (message, C 3); the letter of the panel is that of the pair of the row")
    info = {"pairs_of_the_table": pairs, "pairs_drawn": drawn_pairs, "panel_of_the_pair": panel_of,
            "pairs_not_drawn": {p: S.PAIRS_NOT_DRAWN[p] for p in pairs if p in S.PAIRS_NOT_DRAWN},
            "analyses_with_chains_continued_from_a_saved_state": cont,
            "states_of_the_rows_of_the_table": states,
            "layout_of_the_earlier_generator": {n: {"value": v[0], "line": v[1]} for n, v in k.items()},
            "time_axis": [XMIN.isoformat(), XMAX.isoformat()], "axis_of_the_effective_population_size": [YMIN, YMAX]}
    return V, info


def checks_of_the_table(I, V, info):
    T = I.table(TABLE)
    by_state = {}
    for i, s in info["states_of_the_rows_of_the_table"].items():
        key = "drawn" if s == "drawn" else s.split(":", 1)[1].strip()[:110]
        by_state[key] = by_state.get(key, 0) + 1
    twice = {"compared": 0, "equal": 0, "not_equal": []}
    first = {}
    for i, r in enumerate(T.rows):
        if r["program"] != "Delphy":
            continue
        key = (r["pair"], r["analysis"], r["interval"])
        if key in first:
            twice["compared"] += 1
            same = all(T.rows[first[key]][c] == r[c] for c in ("Ne_tau_years_median", "Ne_tau_years_hpd_lo", "Ne_tau_years_hpd_hi",
                                                                "older_edge", "newer_edge", "meets_convergence_criterion"))
            twice["equal"] += same
            if not same:
                twice["not_equal"].append([first[key], i])
        else:
            first[key] = i
    # R2: the cells of the table against the tables of intervals that are among the inputs. The tables of intervals hold the
    # natural logarithm (columns ln_Ne_tau_...), which the table of the trajectories holds beside the value in years.
    r2 = {"rows_of_the_table_whose_source_table_is_among_the_inputs_in_the_version_named": 0, "cells_compared": 0,
          "cells_equal_as_text": 0, "cells_equal_as_numbers_and_not_as_text": 0, "not_equal": [],
          "rows_whose_source_table_is_not_among_the_inputs": 0, "rows_without_a_source_table": 0,
          "rows_found_by_no_key_or_by_several": [], "columns_compared": "lnNe_median = ln_Ne_tau_median; lnNe_hpd_lo = "
          "ln_Ne_tau_hpd95_lower; lnNe_hpd_hi = ln_Ne_tau_hpd95_upper; older_edge; newer_edge; meets_convergence_criterion",
          "tables": {}}
    for i, r in enumerate(T.rows):
        name = r["source_table"]
        if not name:
            r2["rows_without_a_source_table"] += 1
            continue
        if not (name in I.by_file and I.by_file[name]["version"] == r["source_version"] and I.allowed(name)):
            r2["rows_whose_source_table_is_not_among_the_inputs"] += 1
            continue
        E = I.table(name).need("analysis", "interval", "older_edge", "newer_edge", "ln_Ne_tau_median", "ln_Ne_tau_hpd95_lower",
                               "ln_Ne_tau_hpd95_upper", "meets_convergence_criterion")
        hits = [j for j, e in enumerate(E.rows) if e["analysis"] == r["analysis"] and C.is_number(e["interval"])
                and C.is_number(r["interval"]) and float(e["interval"]) == float(r["interval"])]
        if len(hits) != 1:
            r2["rows_found_by_no_key_or_by_several"].append([i, len(hits)])
            continue
        r2["rows_of_the_table_whose_source_table_is_among_the_inputs_in_the_version_named"] += 1
        r2["tables"][name] = r2["tables"].get(name, 0) + 1
        for a, b in (("lnNe_median", "ln_Ne_tau_median"), ("lnNe_hpd_lo", "ln_Ne_tau_hpd95_lower"), ("lnNe_hpd_hi", "ln_Ne_tau_hpd95_upper"),
                     ("older_edge", "older_edge"), ("newer_edge", "newer_edge"), ("meets_convergence_criterion", "meets_convergence_criterion")):
            r2["cells_compared"] += 1
            x, y = r[a], E.rows[hits[0]][b]
            if x == y:
                r2["cells_equal_as_text"] += 1
            elif C.is_number(x) and C.is_number(y) and float(x) == float(y):
                r2["cells_equal_as_numbers_and_not_as_text"] += 1
            else:
                r2["not_equal"].append([i, a, name, hits[0], b])
    rows = V.rows
    drawn = [r for r in rows if r["element"] == E_STEP and not r["drawn_as"].startswith(S.NOT_DRAWN)]
    return {"rows_of_the_table": len(T.rows), "rows_of_the_table_by_state": dict(sorted(by_state.items())),
            "rows_of_the_table_drawn": len(drawn),
            "intervals_drawn_by_panel_and_program": {
                p: {q: sum(1 for r in drawn if r["panel"] == p and r["drawn_as"].startswith(q)) for q in ("Delphy", "BEAST X: ", "BEAST X, second summary")}
                for p in sorted({r["panel"] for r in drawn})},
            "analyses_drawn_by_panel": {p: sorted({r["analysis"] for r in drawn if r["panel"] == p}) for p in sorted({r["panel"] for r in drawn})},
            "rows_of_Delphy_that_stand_twice_in_the_table": twice,
            "R2_cells_of_the_table_against_the_tables_of_intervals_among_the_inputs": r2,
            "columns_of_the_table_that_are_drawn": ["Ne_tau_years_median", "Ne_tau_years_hpd_lo", "Ne_tau_years_hpd_hi", "older_edge", "newer_edge"]}


# =============================================================================================
def free_corner(ax, fig, texts_of, steps, lines_x):
    """The first corner of CORNERS in which the names of the programs stand on no band and no line of the panel.
    steps: (x0, x1, lower, upper) in data units of everything drawn; lines_x: x of vertical lines."""
    r = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    pad = 4.0 * fig.dpi / 72.0
    box = ax.get_window_extent(r)
    for corner in CORNERS:
        up, left = corner.startswith("upper"), corner.endswith("left")
        x = box.x0 + pad if left else box.x1 - pad
        y = box.y1 - pad if up else box.y0 + pad
        arts, yy = [], y
        order = texts_of if up else list(reversed(texts_of))
        ok = True
        for (words, colour) in order:
            t = ax.text(*inv.transform((x, yy)), words, color=colour, fontsize=C.SIZES[2], ha="left" if left else "right",
                        va="top" if up else "bottom", zorder=6)
            arts.append(t)
            b = t.get_window_extent(r)
            yy = b.y0 - 1.5 * fig.dpi / 72.0 if up else b.y1 + 1.5 * fig.dpi / 72.0
            (x0, y0), (x1, y1) = inv.transform((b.x0 - pad / 2, b.y0 - pad / 2)), inv.transform((b.x1 + pad / 2, b.y1 + pad / 2))
            for (s0, s1, lo, hi) in steps:
                if s0 < x1 and s1 > x0 and lo < y1 and hi > y0:
                    ok = False
            for lx in lines_x:
                if x0 <= lx <= x1:
                    ok = False
        if ok:
            return corner, arts
        for t in arts:
            t.remove()
    return None, []


def draw(values_path):
    rows = C.read_csv(values_path)
    C.style()

    def of(element, panel=None):
        return [r for r in rows if r["element"] == element and (panel is None or r["panel"] == panel)]

    def one(element):
        x = of(element)
        if len(x) != 1:
            raise C.Stop(f"file of values: {len(x)} rows of the element '{element}'")
        return x[0]
    D = lambda t: mdates.date2num(C.to_day(t))        # noqa: E731
    XMIN, XMAX = D(one(E_AXIS_FIRST)["x"]), D(one(E_AXIS_END)["x"])
    year = C.to_day(one(E_AXIS_FIRST)["x"]).year
    DECLARED = D(one(E_DECLARED)["x"])
    YMIN, YMAX = float(one(S.E_AXIS_RULE)["x"]), float(one(S.E_AXIS_RULE)["x_end"])
    PW, PH, LEFT, GAPX, TOP, GAPY, BOTTOM = 2.78, 1.72, 0.72, 0.56, 0.42, 0.62, 0.86
    W = C.FIG_WIDTH
    titles = {r["panel"]: r for r in of(S.E_TITLE)}
    panels = sorted(titles)
    wrapped = {p: S.wrap(titles[p]["label"], PW - 15.0 / 72.0 - 0.02, C.SIZES[0], hang="") for p in panels}
    extra_top = [max(wrapped[p].count("\n") for p in panels[2 * k:2 * k + 2]) * LINE_IN for k in (0, 1)]
    keys = of(S.E_KEY)
    key_w = W - LEFT - 0.35 - 0.05
    key_labels = [S.wrap(e["label"], key_w, C.SIZES[2], hang="") for e in keys]
    key_lines = sum(x.count("\n") + 1 for x in key_labels)
    bottom = BOTTOM + max(0, key_lines - 2) * (9.6 / 72.0 + 0.4 * 8.0 / 72.0)
    height = TOP + extra_top[0] + 2 * PH + GAPY + extra_top[1] + bottom
    fig = C.new_figure(height)
    AX, names, corners, steps_of = {}, {}, {}, {}
    for kk, p in enumerate(panels):
        row, col = divmod(kk, 2)
        top = TOP + extra_top[0] + row * (PH + GAPY + extra_top[1])
        ax = fig.add_axes([(LEFT + col * (PW + GAPX)) / W, 1 - (top + PH) / height, PW / W, PH / height])
        ax.set_gid(p)
        AX[p] = ax
        ax.set_yscale("log")
        ax.set_ylim(YMIN, YMAX)
        ax.set_xlim(XMIN, XMAX)
        occupied = []
        st = [r for r in of(E_STEP, p) if not r["drawn_as"].startswith(S.NOT_DRAWN)]
        groups = []
        for r in st:
            key = (r["analysis"], r["drawn_as"])
            if key not in groups:
                groups.append(key)
        for (a, how) in groups:
            g = sorted([r for r in st if (r["analysis"], r["drawn_as"]) == (a, how)], key=lambda r: C.to_day(r["x"]))
            for u, v in zip(g[:-1], g[1:]):
                if C.to_day(u["x_end"]) != C.to_day(v["x"]):
                    raise C.Stop(f"panel {p}, {a}: the intervals are not contiguous")
            e = [D(r["x"]) for r in g] + [D(g[-1]["x_end"])]
            med = [float(r["value"]) for r in g]
            lo = [float(r["lower"]) for r in g]
            hi = [float(r["upper"]) for r in g]
            occupied += [(e[j], e[j + 1], lo[j], hi[j]) for j in range(len(g))]
            if how.startswith("Delphy"):
                ax.stairs(hi, e, baseline=lo, fill=True, facecolor=C.BLUE, alpha=C.BLUE_ALPHA, edgecolor="none", zorder=2)
                ax.stairs(med, e, baseline=None, color=C.BLUE, lw=LW["Delphy"], zorder=3)
            elif how.startswith("BEAST X: "):
                ax.stairs(med, e, baseline=None, color=C.GREEN, lw=LW["BEAST X"], zorder=4)
                for b in (lo, hi):
                    ax.stairs(b, e, baseline=None, color=C.GREEN, lw=LW["bound"], zorder=3.5)
            else:
                tone = C.blend(C.GREEN, TONE_OF_THE_CUT)
                ax.stairs(med, e, baseline=None, color=tone, lw=LW["cut"], zorder=4.5)
                for b in (lo, hi):
                    ax.stairs(b, e, baseline=None, color=tone, lw=LW["cut bound"], zorder=4.4)
        steps_of[p] = occupied
        ax.axvline(DECLARED, color=C.GREY_REF, ls=(0, (4, 3)), lw=0.9, zorder=1.5)
        e0, e1 = int(round(math.log10(YMIN))), int(round(math.log10(YMAX)))
        yt = [10.0 ** j for j in range(e0, e1 + 1)]
        ax.set_yticks(yt)
        ax.set_yticklabels([(f"{v:g}" if 0.01 <= v <= 1000 else r"$10^{%d}$" % int(round(math.log10(v)))) for v in yt] if col == 0 else [])
        ax.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=[0.1 * j for j in range(2, 10)], numticks=12))
        ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
        xt = [dt.date(year, m, 1) for m in range(1, 13) if XMIN <= mdates.date2num(dt.date(year, m, 1)) <= XMAX]
        ax.set_xticks([mdates.date2num(x) for x in xt])
        ax.set_xticklabels([C.MONTHS[x.month - 1] for x in xt])
        if col == 0:
            ax.set_ylabel(Y_LABEL[0])
        if row == 1:
            ax.set_xlabel(str(year))
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        S.panel_head(ax, p, wrapped[p])
    fig.canvas.draw()
    word = AX[panels[0]].annotate(one(E_DECLARED)["label"], xy=(DECLARED, YMAX), xycoords="data", xytext=(-3, -2),
                                  textcoords="offset points", ha="right", va="top", fontsize=C.SIZES[2], color="0.35")
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    for p in panels:
        ax = AX[p]
        texts = []
        for r in of(E_NAME, p):
            c = C.BLUE if r["value"] == "Delphy" else C.GREEN
            texts.append((r["label"], c))
        occ = list(steps_of[p])
        if p == panels[0]:
            b = word.get_window_extent(rr)
            (x0, y0), (x1, y1) = ax.transData.inverted().transform((b.x0, b.y0)), ax.transData.inverted().transform((b.x1, b.y1))
            occ.append((x0, x1, y0, y1))
        corner, arts = free_corner(ax, fig, texts, occ, [DECLARED])
        if corner is None:
            raise C.Stop(f"panel {p}: no corner of the panel is free for the names of the programs")
        corners[p] = corner
        names[p] = arts
    handles = []
    for e, lab in zip(keys, key_labels):
        d = e["drawn_as"]
        if d.startswith("band"):
            handles.append(Patch(facecolor=C.blend(C.BLUE, C.BLUE_ALPHA), edgecolor=C.BLUE, lw=1.2, label=lab))
        elif d.startswith("thick"):
            handles.append(Line2D([0], [0], color=C.GREEN, lw=LW["BEAST X"], label=lab))
        else:
            handles.append(Line2D([0], [0], color=C.blend(C.GREEN, TONE_OF_THE_CUT), lw=LW["cut"], label=lab))
    leg = fig.legend(handles=handles, loc="lower left", bbox_to_anchor=(LEFT / W, 0.02 / height), frameon=False,
                     fontsize=C.SIZES[2], ncol=1, handlelength=1.8, labelspacing=0.4, borderaxespad=0)
    fig.canvas.draw()
    key_box = leg.get_window_extent(rr)
    lowest = min(t.get_window_extent(rr).y0 for p in panels[2:] for t in [AX[p].xaxis.label])
    extra = {"height_in": round(height, 4), "corner_of_the_names_by_panel": corners,
             "titles_on_more_than_one_line": [wrapped[p].replace("\n", " / ") for p in panels if "\n" in wrapped[p]],
             "key_below_the_labels_of_the_time_axis": bool(key_box.y1 <= lowest),
             "key_inside_the_figure": bool(key_box.x1 <= fig.bbox.x1 and key_box.y0 >= 0),
             "layout_of_the_earlier_generator_carried_over": {"PW": PW, "PH": PH, "LEFT": LEFT, "GAPX": GAPX, "TOP": TOP, "GAPY": GAPY,
                                                              "BOTTOM": BOTTOM},
             "widths_of_lines_pt": LW}
    return fig, extra


def main(argv):
    if len(argv) >= 3 and argv[1] == "--draw-only":
        return S.draw_only(argv, draw, NAME)
    if len(argv) < 3:
        raise C.Stop("usage: python fig_program_trajectories_20261006.py INPUTS_SUPPLEMENT.json INPUTS_THIRD_TASK.json [OUTDIR]")
    outdir = argv[3] if len(argv) > 3 else "."
    os.makedirs(outdir, exist_ok=True)
    I = S.inputs_of(argv, 1)
    V, info = build_values(I)
    k = info["layout_of_the_earlier_generator"]
    typed = {"PW": 2.78, "PH": 1.72, "LEFT": 0.72, "GAPX": 0.56, "TOP": 0.42, "GAPY": 0.62, "BOTTOM": 0.86}
    for n, v in typed.items():
        if abs(float(k[n]["value"]) - v) > 1e-12:
            raise C.Stop(f"{GENERATOR}: {n} = {k[n]['value']}; this program holds {v}")
    values_path = V.save(os.path.join(outdir, VALUES_FILE))       # the file of values is written first
    fig, extra = draw(values_path)                                 # ... and the figure is drawn from it
    out, ok = S.standard(I, fig, values_path, outdir, NAME, os.path.basename(__file__), extra)
    out["info"] = {kk: v for kk, v in info.items() if kk != "states_of_the_rows_of_the_table"}
    out["table"] = t = checks_of_the_table(I, V, info)
    more = (not t["rows_of_Delphy_that_stand_twice_in_the_table"]["not_equal"]
            and not t["R2_cells_of_the_table_against_the_tables_of_intervals_among_the_inputs"]["not_equal"]
            and not t["R2_cells_of_the_table_against_the_tables_of_intervals_among_the_inputs"]["rows_found_by_no_key_or_by_several"]
            and sum(t["rows_of_the_table_by_state"].values()) == t["rows_of_the_table"]
            and extra["key_below_the_labels_of_the_time_axis"] and extra["key_inside_the_figure"]
            and all(c.get("the_two_counts_are_equal", True) for c in info["analyses_with_chains_continued_from_a_saved_state"].values()))
    return S.finish(out, ok, outdir, NAME, more_ok=more,
                    said={"rows_of_the_table_by_state": t["rows_of_the_table_by_state"], "corners": extra["corner_of_the_names_by_panel"]})


if __name__ == "__main__":
    sys.exit(main(sys.argv))
