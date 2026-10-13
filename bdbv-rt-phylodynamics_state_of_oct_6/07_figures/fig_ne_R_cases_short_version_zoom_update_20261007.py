"""fig_ne_R_cases_short_version_zoom_update_20261007.py. Figure 1 of the article with a ZOOM PANEL d below panel c.
The constant R_TOP_D, the end of the axis of R of panel d, is 6.0: the next whole number above the largest upper bound
drawn in the range of panel d (a weekly pair has an upper bound just above 5). Just before the picture is saved, run()
calls the function root_date_measured of the program of the event legend, which measures on the picture the distance
between the dates of the root and the band of the staircase and the distance between the bar of the root and the time
axis, and stops the run before a picture is saved if one of them is below the smallest distance; what it returns goes
into the record. The name of the picture is an argument of run() (stem), not a constant.

What it does. It follows the program fig_ne_R_cases_short_version_event_legend_20261007.py (called P here) step
by step and calls its functions for every step that does not change:
  1. panels a, b and c are drawn as P.draw draws them, up to the names P1 to P3 and the room between the panels
     (P.draw, lines 466 to 481);
  2. NEW. The strip of numbered marks above panel b is taken away, and no legend of the events is drawn. The head of
     panel b is set above its frame as the conventions set the heads of a and c (panel_head). Panel a moves down by
     SHIFT_A_PX pixels of the PNG, and the picture grows by GROWTH_PX pixels below panel c. Panel d is drawn there
     from its file of values, which this program writes from tables of the analyses (values_of_panel_d);
  3. the words of panel a are cleared and the key above panel a is placed (P.draw, lines 482 to 485); the names of
     the periods are set in bold and the keys of b and c are set on white (P.run, lines 579 and 580);
  4. the checks of P that do not concern the legend of the events, and the checks of panel d;
  5. the picture is saved by the function of the conventions (PNG at 300 dpi and PDF), and its two parts
     (panels a to c; panel d) are saved as PNG, each the rows of the whole picture.
No number of panels a to c is touched: they are drawn from the file of values of the figure of panels a to c.

Panel d (reproduction number between 12 April and 15 July 2026, with the response events):
  - primary analysis: every pair of adjacent Skygrid intervals whose two mid-points lie in the range, as a band
    over the time between the two mid-points (thick line: median; height: 95 % HPD interval), in the colour of the
    points of panel b, which are the same estimates; numbers of primary_pairs_sensitivity_20261007.csv;
  - analysis with intervals of 7.2 days (delphy_sens_gridWeekly32): every pair of adjacent intervals that meet in
    the range, as a diamond with its 95 % HPD interval at the time where the two intervals meet; numbers of the
    series 'reproduction number between successive mid-points' of trajectories_delphy_sens_gridWeekly32_20261007.csv;
    open where the posterior standard deviation of the growth rate is 0.90 or more of the prior one (table computed
    with the functions of analysis.py and rt_lib.py);
  - the line R = 1, the declaration of the outbreak, the names of the periods with brackets;
  - under the time axis of panel b a thin bracket marks the range that panel d shows;
  - the nine response events of the checked table at their dates, each with its number and its short name.
Every date stands at the start of its day, as in panels a to c; the estimates stand at their times on the grid of
their analysis (rt_lib.Grid), and these times are held against the dates of the tables.

SHIFT_A_PX is a multiple of 25 pixels: the period of the hatching of panel c, which the backends anchor at the top
of the picture (constant HATCH_STEP_PX of the original figure program), so that the hatching keeps its place.
GROWTH_PX less SHIFT_A_PX is a multiple of 3 pixels of the PNG, that is a whole number of the pixels of the drawing
(200 dpi) on which the programs of panels a to c make their checks and place the key above panel a: otherwise panels a to c
stand on another pixel grid than in the picture without panel d, and the key of panel a is placed elsewhere.
A thin guiding line shows the "zoom in" between the time axis of panel c and the whole range of panel d:
  - a thin grey line joins the range of panel d on the time axis of panel c to the two ends of the time axis of panel d
    (zoom_guides): straight down from the axis of c past the names of the months, slanting to the upper corners of the
    room of panel d, and across the whole width of panel d above its head. The range is that of the row of values of
    the bracket under panel b; the program stops if a piece of the line comes near a text.
  - With the line (GUIDES = True) the two parts of the picture are NOT cut: the line crosses the row of the cut. With
    GUIDES = False the picture is saved in its two parts.
"""
import csv
import datetime as dt
import importlib.util
import math
import os
import re
import sys

SHIFT_A_PX = 175                  # panel a moves down by this (42.0 pt): then a stands 45.3 pt above b; b stands 44.6 pt above c
GROWTH_PX = 1225                  # the picture grows by this below panel c (294.0 pt)
HEIGHT_LIMIT_IN = 12.2            # the largest height for the saving of this figure (the conventions hold 8.5 in)
WINDOW = ("2026-04-12", "2026-07-15")     # the range of panel d
R_TOP_D = 6.0                     # the axis of R of panel d ends here
FRAME_D_PT = 132.0                # height of the frame of panel d
PRIMARY_GRID = (20, 0.8)          # parameters and cut-off (years) of the Skygrid of the primary analysis
WEEKLY_GRID = (32, 0.613699)      # the same of the analysis delphy_sens_gridWeekly32 (command lines of its chains)
SERIES_W = "reproduction number between successive mid-points"
MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")
TICKS_D = ("2026-04-15", "2026-05-01", "2026-05-15", "2026-06-01", "2026-06-15", "2026-07-01", "2026-07-15")
MINOR_TICKS_D = ("2026-04-22", "2026-05-08", "2026-05-22", "2026-06-08", "2026-06-22", "2026-07-08")
HEAD_B = "Reproduction number"
HEAD_D = "Reproduction number and response events, 12 April \u2013 15 July"
KEY_D_TITLE = "Pairs of adjacent intervals: median, 95 % HPD"
KEY_D = ("primary analysis (intervals of 15.4 days)", "weekly grid (intervals of 7.2 days)",
         "the same, posterior SD / prior SD \u2265 0.90")
BAND_ALPHA = 0.40                 # the band of a pair of the primary analysis: the colour of the pairs of panel b at this share
GREY_W = "0.35"                   # the colour of the analysis with intervals of 7.2 days
MS_W = 3.4                        # size of its symbol
# above the frame of panel d; every length in points above the upper edge of the frame
RULER_LINE_PT = 5.0               # the line of a bracket of a period
RULER_TICK_PT = 2.5               # its tick at an end, downwards
RULER_INSET_PT = 1.5              # a bracket ends this far from the day at which its period ends
NAME_BASE_PT = 8.5                # lower edge of the name of a period
EVENT_BASE_PT = 27.0              # the edge on which a triangle stands
CROSS_UP_PT = 6.5                 # a cross stands this much higher, so that marks of neighbouring days do not touch
SPAN_Y_PT = -2.6                  # middle of a bar and of a dotted line: under that edge
ROW0_PT = EVENT_BASE_PT + 18.5    # middle of the lowest row of names: above the crosses
ROW_PITCH_PT = 10.0               # from row to row
LABEL_PAD_PT = 2.5                # between a thin line and its name
STEM_LW_PT = 0.6                  # the thin line from a mark to its name
ARROW_PT = 3.0                    # arrow head at an end of a bracket that is cut at the axis
HEAD_GAP_PT = 7.0                 # from the highest name to the head of the panel
TEXT_HALF_PT = 4.0                # half the height of a name of 8 pt
ZOOM_WORDS = "range of panel d"   # under the time axis of panel b
ZOOM_LINE_PT = 8.0                # the line of the mark of the range, under the lower edge of the frame of panel b
ZOOM_TICK_PT = 2.5                # its tick at an end, upwards
ZOOM_TEXT_GAP_PT = 1.5            # from the line to its words
ZOOM_COLOUR = "0.35"
# ----- addition of 8 October 2026: the guiding line between the time axis of panel c and the two ends of panel d
GUIDES = True                     # False: the picture of 7 October, with its two parts
GUIDE_COLOUR = ZOOM_COLOUR        # the grey of the bracket under panel b
GUIDE_LW_PT = 0.6
GUIDE_BELOW_NAMES_PT = 3.0        # the straight pieces end this far below the names of the months of panel c
GUIDE_ABOVE_HEAD_PT = 4.5         # the piece across panel d stands this far above the head of panel d
GUIDE_LEAST_FALL_PT = 20.0        # least height of the slanting pieces
GUIDE_CLEAR_OF_TEXT_PT = 1.5      # no piece comes nearer than this to the box of a text

E_AX0, E_AX1 = "time axis of panel d, first day", "time axis of panel d, end"
E_TICK = "time axis of panel d, tick with its label"
E_DECL = "declaration of the outbreak"
E_RAX = "limits of an axis by rule"
E_PERIOD = "name of a period, printed above the panel"
E_PAIR_P = "R of a pair of adjacent Skygrid intervals, primary analysis, over the time between the two mid-points"
E_PAIR_W = "R of a pair of adjacent Skygrid intervals, analysis with intervals of 7.2 days"
E_EVENT = "intervention, mark above panel d with its number and short name"
E_KEY = "entry of a key"
E_HEAD = "head of a panel"
E_ZOOM = "range of panel d, marked under the time axis of panel b"
COLUMNS = ["figure", "panel", "element", "label", "value", "lower", "upper", "x", "x_end", "unit", "drawn_as",
           "source_table", "source_version_id", "source_row_0based", "source_columns", "note", "analysis",
           "quantity_id", "delivered", "meets_convergence_criterion"]


class Stop(SystemExit):
    pass


def load(path, name):
    spec = importlib.util.spec_from_file_location(name, path)
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def num(text):
    """a day 'YYYY-MM-DD' (its start) or a time 'YYYY-MM-DDTHH:MM:SS' as a number of the time axis"""
    import matplotlib.dates as mdates
    return float(mdates.date2num(dt.datetime.fromisoformat(text)))


def iso(ts):
    return ts.strftime("%Y-%m-%dT%H:%M:%S")


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


# ---------------------------------------------------------------------------------------------------------------
# the file of values of panel d
# ---------------------------------------------------------------------------------------------------------------
def values_of_panel_d(inp, words, names, stem):
    """the rows of the file of values of panel d. inp: for each of 'primary_pairs', 'weekly_series',
    'weekly_computed', 'events', 'values_abc' a dict with 'path', 'file' and 'version'; 'pipeline_dir'."""
    if inp["pipeline_dir"] not in sys.path:
        sys.path.insert(0, inp["pipeline_dir"])
    import rt_lib as L
    rows = []

    def add(element, **kw):
        r = {c: "" for c in COLUMNS}
        r.update(figure=stem, panel="d", element=element, delivered="yes")
        for k, v in kw.items():
            if k not in r:
                raise Stop("no column %r" % k)
            r[k] = "" if v is None else str(v)
        rows.append(r)
        return r

    def src(key, row=None):
        return dict(source_table=inp[key]["file"], source_version_id=inp[key]["version"],
                    source_row_0based="" if row is None else row)

    w0, w1 = WINDOW
    W0, W1 = dt.datetime.fromisoformat(w0), dt.datetime.fromisoformat(w1)
    abc = read_csv(inp["values_abc"]["path"])

    def one_abc(element):
        got = [(k, r) for k, r in enumerate(abc) if r["element"] == element]
        if len(got) != 1:
            raise Stop("the file of values of panels a to c holds %d rows %r" % (len(got), element))
        return got[0]

    k_first, first = one_abc("time axis, first day")
    add(E_AX0, label=first["label"], x=w0, unit="date",
        drawn_as="left end of the time axis of panel d; the year is printed below the panel",
        source_table="constant WINDOW of this program (the range); the year: %s, row %d"
                     % (inp["values_abc"]["file"], k_first), source_version_id=inp["values_abc"]["version"])
    add(E_AX1, x=w1, unit="date", drawn_as="right end of the time axis of panel d",
        source_table="constant WINDOW of this program (the range)")
    for t in TICKS_D:
        d = dt.date.fromisoformat(t)
        add(E_TICK, label="%d %s" % (d.day, MONTHS[d.month - 1]), x=t, unit="date", drawn_as="tick with its label")
    for t in MINOR_TICKS_D:
        add(E_TICK, label="", x=t, unit="date", drawn_as="small tick without a label")
    add(E_RAX, label="axis of R, panel d", x=0, x_end=R_TOP_D, unit="R",
        drawn_as="limits of the axis of R; ticks at the whole numbers",
        note="the largest upper bound drawn is below the end of the axis (checked when the rows are written)")
    k_decl, decl = one_abc(E_DECL)
    add(E_DECL, label=decl["label"], value=decl["value"], x=decl["x"], unit="date",
        drawn_as="vertical dashed line through the panel up to its words above the panel",
        note="row %d of the file of values of panels a to c, which gives as its source: %s, row %s"
             % (k_decl, decl["source_table"], decl["source_row_0based"]), **src("values_abc", k_decl))
    per = [(k, r) for k, r in enumerate(abc) if r["element"] == "name of a period, printed above the panel"]
    if len(per) != len(names):
        raise Stop("the file of values of panels a to c holds %d names of periods" % len(per))
    for (k, r), name in zip(per, names):
        if r["label"][1:] != name[1:]:
            raise Stop("the periods are not in their order: %s, %s" % (r["label"], name))
        add(E_PERIOD, label=name, value=name, x=r["x"], x_end=r["x_end"], unit="date",
            drawn_as="name above a bracket over the part of the period that lies in the range; an arrow head where "
                     "the period goes on beyond the axis",
            note="row %d of the file of values of panels a to c (there: %s), which gives as its source: %s, row %s"
                 % (k, r["label"], r["source_table"], r["source_row_0based"]), analysis=r["analysis"],
            **src("values_abc", k))

    # ----- the pairs of the primary analysis
    comp = read_csv(inp["weekly_computed"]["path"])
    last = sorted(set(r["most_recent_collection_date"] for r in comp))
    if len(last) != 1:
        raise Stop("the most recent collection date is not one: %s" % last)
    last = last[0]
    P = read_csv(inp["primary_pairs"]["path"])
    g0 = L.Grid.from_cutoff(PRIMARY_GRID[0], PRIMARY_GRID[1], last)
    if len(P) != PRIMARY_GRID[0] - 1:
        raise Stop("the table of the primary pairs holds %d rows" % len(P))
    n_p = 0
    for k, r in enumerate(P):
        m = re.fullmatch(r"(\d+) to (\d+)", r["pair"])
        i = int(m.group(2))
        if int(m.group(1)) != i + 1:
            raise Stop("pair %r" % r["pair"])
        meet, older, newer = L.dec2ts(g0.t_anchor - i * g0.D), L.dec2ts(g0.mids[i]), L.dec2ts(g0.mids[i - 1])
        if (str(meet.date()) != r["date_on_which_the_two_intervals_meet"]
                or str(newer.date()) != r["mid_point_of_the_newer_interval"]
                or (r["mid_point_of_the_older_interval"] and str(older.date()) != r["mid_point_of_the_older_interval"])):
            raise Stop("the grid of the primary analysis does not give the dates of pair %s of its table" % r["pair"])
        if older.to_pydatetime() >= W0 and newer.to_pydatetime() <= W1:
            if r["sd_ratio_of_0_90_or_more"] != "no":
                raise Stop("pair %s of the primary analysis is an open symbol in panel b; panel d has no form for it" % r["pair"])
            if float(r["R_hpd95_upper"]) >= R_TOP_D:
                raise Stop("an upper bound reaches the end of the axis")
            add(E_PAIR_P, label="pair %s, intervals meet on %s" % (r["pair"], r["date_on_which_the_two_intervals_meet"]),
                value=r["R_median"], lower=r["R_hpd95_lower"], upper=r["R_hpd95_upper"], x=iso(older), x_end=iso(newer),
                unit="R", drawn_as="band (95 % HPD) and thick line (median) from the mid-point of the older to the "
                                   "mid-point of the newer interval",
                source_columns="value=R_median; lower=R_hpd95_lower; upper=R_hpd95_upper; x, x_end: the mid-points of "
                               "intervals %d and %d on the grid of the analysis (%d parameters, cut-off %s years, anchored "
                               "at the most recent genome, %s), whose days are the columns mid_point_of_the_older_interval "
                               "and mid_point_of_the_newer_interval" % (i + 1, i, PRIMARY_GRID[0], PRIMARY_GRID[1], last),
                note="the intervals meet at %s; posterior SD / prior SD of the growth rate: %s (column sd_ratio)"
                     % (iso(meet), r["sd_ratio"]), analysis=r["analysis"],
                meets_convergence_criterion=r["meets_convergence_criterion"], **src("primary_pairs", k))
            n_p += 1

    # ----- the pairs of the analysis with intervals of 7.2 days
    T = read_csv(inp["weekly_series"]["path"])
    S = [(k, r) for k, r in enumerate(T) if r["series"] == SERIES_W]
    by_pair = {}
    for k, r in enumerate(comp):
        by_pair[r["pair"]] = (k, r)
    gw = L.Grid.from_cutoff(WEEKLY_GRID[0], WEEKLY_GRID[1], last)
    if len(S) != WEEKLY_GRID[0] - 1 or len(comp) != len(S):
        raise Stop("the series holds %d rows, the computed table %d" % (len(S), len(comp)))
    n_w = n_open = 0
    for k, r in S:
        m = re.fullmatch(r"(\d+) to (\d+)", r["interval"])
        i = int(m.group(2))
        dec = gw.t_anchor - i * gw.D
        if abs(dec - float(r["date_decimal"])) > 1e-9:
            raise Stop("the grid of the analysis with intervals of 7.2 days does not give the time of pair %s" % r["interval"])
        kc, c = by_pair[r["interval"]]
        for a, b in (("R_median", "median"), ("R_hpd95_lower", "hpd95_lower"), ("R_hpd95_upper", "hpd95_upper"),
                     ("ess_pooled", "ess_pooled"), ("rhat_split", "rhat_split")):
            if abs(float(c[a]) - float(r[b])) > 1e-12:
                raise Stop("pair %s: %s of the computed table is not %s of the stored series" % (r["interval"], a, b))
        if c["date_on_which_the_two_intervals_meet"] != r["date"] or c["states_used"] != r["states_used"]:
            raise Stop("pair %s: day or states differ between the computed table and the stored series" % r["interval"])
        meet = L.dec2ts(dec)
        if str(meet.date()) != r["date"]:
            raise Stop("pair %s: the time on the grid is not on the day of the table" % r["interval"])
        if W0 <= meet.to_pydatetime() <= W1:
            if r["note"].strip():
                raise Stop("pair %s in the range carries a note: %s" % (r["interval"], r["note"]))
            if float(r["hpd95_upper"]) >= R_TOP_D:
                raise Stop("an upper bound reaches the end of the axis")
            is_open = c["sd_ratio_of_0_90_or_more"] == "yes"
            if is_open != (float(c["sd_ratio"]) >= 0.90):
                raise Stop("pair %s: the mark of the ratio does not fit the ratio" % r["interval"])
            add(E_PAIR_W, label="pair %s, intervals meet on %s" % (r["interval"], r["date"]), value=r["median"],
                lower=r["hpd95_lower"], upper=r["hpd95_upper"], x=iso(meet), unit="R",
                drawn_as=("open" if is_open else "filled") + " grey diamond (median) with a line (95 % HPD) at the time "
                         "where the two intervals meet" + ("; open: posterior SD / prior SD of 0.90 or more" if is_open else ""),
                source_columns="value=median; lower=hpd95_lower; upper=hpd95_upper; x: date_decimal as a time (the grid of "
                               "the analysis: %d parameters, cut-off %s years, anchored at the most recent genome, %s)"
                               % (WEEKLY_GRID[0], WEEKLY_GRID[1], last),
                note="posterior SD / prior SD of the growth rate: %s (table %s, version %s, row %d, column sd_ratio; the "
                     "stored series does not hold it)" % (c["sd_ratio"], inp["weekly_computed"]["file"],
                                                          inp["weekly_computed"]["version"], kc),
                analysis=r["analysis"], meets_convergence_criterion=r["meets_convergence_criterion"],
                **src("weekly_series", k))
            n_w += 1
            n_open += int(is_open)

    # ----- the response events
    E = read_csv(inp["events"]["path"])
    marks_abc = [r for r in abc if r["element"] == "intervention, mark in the strip above panel b"]
    if len(E) != len(words) or len(marks_abc) != len(E):
        raise Stop("%d events, %d short names, %d marks of the delivered figure" % (len(E), len(words), len(marks_abc)))
    for k, r in enumerate(E):
        n = int(r["order"])
        if n != k + 1:
            raise Stop("the events are not in their order")
        if r["first_day"] == r["last_day"]:
            mark = "cross" if r["kind"] == "attack on the response" else "triangle"
            how = {"cross": "cross on the edge above the panel, at the day",
                   "triangle": "triangle on the edge above the panel, tip at the day"}[mark]
        elif "not known" in r["meaning_of_the_period"]:
            mark, how = "dotted line", "dotted line under the edge above the panel, between the two days: the day is not known"
        else:
            mark, how = "bar", "bar under the edge above the panel, from the first to the last day"
        old = marks_abc[k]
        if not old["drawn_as"].startswith(mark) or (old["x"], old["x_end"]) != (r["first_day"], r["last_day"]):
            raise Stop("event %d: mark or days differ from the delivered figure (%s; %s to %s)"
                       % (n, old["drawn_as"][:30], old["x"], old["x_end"]))
        if not (W0 <= dt.datetime.fromisoformat(r["first_day"]) and dt.datetime.fromisoformat(r["last_day"]) <= W1):
            raise Stop("event %d lies outside the range" % n)
        add(E_EVENT, label="%d %s" % (n, words[k]), value=n, x=r["first_day"], x_end=r["last_day"], unit="date",
            drawn_as=mark + ": " + how + "; a thin line to the number and the short name",
            source_columns="x=first_day; x_end=last_day; the mark by kind (%s)%s; the short name: WORDS of "
                           "fig_ne_R_cases_short_version_event_legend_20261007.py"
                           % (r["kind"], (" and meaning_of_the_period (%s)" % r["meaning_of_the_period"])
                              if r["meaning_of_the_period"] else ""),
            note="%s | sources: %s | row %s of the checked table of events"
                 % (r["words_for_the_legend"], r["sources"], r["row_of_the_checked_table"]), **src("events", k))

    add(E_KEY, label=KEY_D_TITLE, drawn_as="title of the key inside panel d")
    for lab, how in zip(KEY_D, ("band with a thick line", "filled grey diamond on a line", "open grey diamond on a line")):
        add(E_KEY, label=lab, drawn_as="entry of the key inside panel d: " + how)
    add(E_HEAD, label="d", value=HEAD_D, drawn_as="letter and head of panel d")
    add(E_ZOOM, label=ZOOM_WORDS, x=w0, x_end=w1, unit="date",
        drawn_as="thin bracket under the time axis of panel b with its words",
        source_table="constant WINDOW of this program (the range)")
    rows[-1]["panel"] = "b"
    return rows, {"primary_pairs_in_the_range": n_p, "weekly_pairs_in_the_range": n_w, "of_them_open": n_open,
                  "most_recent_collection_date": last, "primary_grid": {"K": g0.K, "D_days": g0.D * L.YEAR_DAYS},
                  "weekly_grid": {"K": gw.K, "D_days": gw.D * L.YEAR_DAYS}}


def write_values(rows, path):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow(r)


# ---------------------------------------------------------------------------------------------------------------
# step 2: the strip goes, the head of panel b, the room, panel d
# ---------------------------------------------------------------------------------------------------------------
def styles_of_the_strip(strip):
    """the looks of the marks of the figure of panels a to c, read from its drawing"""
    import matplotlib as mpl
    seen = {}
    for m in strip["marks"]:
        for kind, a in m["artists"]:
            d = {"lw": float(a.get_linewidth()), "ms": float(a.get_markersize()), "mew": float(a.get_markeredgewidth()),
                 "marker": str(a.get_marker()), "colour": mpl.colors.to_hex(a.get_color())}
            seen.setdefault(kind, []).append(d)
    out = {}
    for kind, ds in seen.items():
        if any(d != ds[0] for d in ds):
            raise Stop("the marks %r of the strip do not all look alike" % kind)
        out[kind] = ds[0]
    if set(out) != {"triangle", "cross", "bar", "dotted line"}:
        raise Stop("the strip holds the kinds %s" % sorted(out))
    if len(set(d["colour"] for d in out.values())) != 1:
        raise Stop("the marks of the strip have more than one colour")
    return out


def strip_away_and_room(fig, C, dpi):
    """the strip of the marks is removed with its head; the head of panel b is set by the conventions; panel a moves
    down by SHIFT_A_PX; the picture grows by GROWTH_PX below panel c. Returns the record and the place of panel d."""
    from matplotlib.patches import Rectangle
    I = fig.interventions
    axb, axc, axs = I["axb"], I["axc"], I["axes"]
    axa = fig.short_version["axa"]
    own = list(fig.texts) + list(fig.legends) + list(fig.lines) + list(fig.patches) + list(fig.images) \
        + [a for a in fig.artists if not (isinstance(a, Rectangle) and a.get_data_transform() == fig.transFigure
                                          and [round(float(v), 9) for v in a.get_bbox().bounds] == [0.0, 0.0, 1.0, 1.0])]
    if own:
        raise Stop("the figure holds artists of its own; they would not move: %s" % [type(a).__name__ for a in own])
    if SHIFT_A_PX % 25 or (GROWTH_PX - SHIFT_A_PX) % 3:
        raise Stop("SHIFT_A_PX must be a multiple of 25, and GROWTH_PX less SHIFT_A_PX a multiple of 3")
    heads = [t.get_text() for t in axs.texts if hasattr(t, "xyann")]
    if heads != ["b", HEAD_B]:
        raise Stop("the strip holds the heads %s" % heads)
    looks = styles_of_the_strip(I["strip"])
    n_marks = len(I["strip"]["marks"])
    axs.remove()
    if any(ax is axs for ax in fig.axes):
        raise Stop("the strip is still an axes of the figure")
    new_heads = C.panel_head(axb, "b", HEAD_B)
    follows = []
    for ax in fig.axes:
        if ax.get_axes_locator() is not None:
            twins = [o for o in ax._twinned_axes.get_siblings(ax) if o is not ax and o.get_axes_locator() is None]
            if len(twins) != 1:
                raise Stop("an axes of the figure is placed by a locator and is not the twin of one axes")
            follows.append(ax)
    w_in, h1 = [float(v) for v in fig.get_size_inches()]
    rows1 = int(h1 * dpi)
    rows2 = rows1 - SHIFT_A_PX + GROWTH_PX
    h2 = rows2 / float(dpi)
    while int(h2 * dpi) < rows2:
        h2 = math.nextafter(h2, math.inf)
    ya, yb, yc = axa.get_position().y0, axb.get_position().y0, axc.get_position().y0
    if not ya > yb > yc:
        raise Stop("panels a, b and c do not stand in this order from above")
    before = [(ax, [float(v) for v in ax.get_position().bounds]) for ax in fig.axes]
    fig.set_size_inches(w_in, h2)
    moved = []
    for ax, (x0, y0, w, h) in before:
        if any(ax is f for f in follows):
            moved.append("with its twin")
            continue
        up_px = GROWTH_PX - (SHIFT_A_PX if y0 >= ya - 1e-9 else 0)
        ax.set_position([x0, (y0 * h1 + up_px / float(dpi)) / h2, w, h * h1 / h2])
        moved.append(up_px)
    xb0, _, wb, _ = [float(v) for v in before[[ax for ax, _ in before].index(axb)][1]]
    bottom_pt = yc * h1 * 72.0                               # panel d stands as far above the lower edge as panel c did
    place = [xb0, bottom_pt / (h2 * 72.0), wb, FRAME_D_PT / (h2 * 72.0)]
    return {"marks_removed": n_marks, "heads_removed": heads, "heads_of_b": [t.get_text() for t in new_heads],
            "height_before_in": h1, "height_in": h2, "rows_before": rows1, "rows": rows2,
            "rows_of_panels_a_to_c": rows1 - SHIFT_A_PX, "axes_moved_up_px": moved,
            "SHIFT_A_PX": SHIFT_A_PX, "GROWTH_PX": GROWTH_PX, "bottom_of_panel_d_pt": bottom_pt,
            "looks_of_the_marks": looks}, place, looks


def range_mark_under_b(fig, rows, C):
    """a thin bracket with its words under the time axis of panel b over the range of panel d"""
    axb = fig.interventions["axb"]
    r = [x for x in rows if x["element"] == E_ZOOM]
    if len(r) != 1:
        raise Stop("%d rows %r" % (len(r), E_ZOOM))
    r = r[0]
    x0, x1 = num(r["x"]), num(r["x_end"])
    lo, hi = axb.get_xlim()
    if not lo <= x0 < x1 <= hi:
        raise Stop("the range of panel d does not lie on the time axis of panel b")
    fig.canvas.draw()
    hb = axb.get_window_extent().height * 72.0 / fig.dpi
    tr = axb.get_xaxis_transform()
    ya, yb_ = -(ZOOM_LINE_PT - ZOOM_TICK_PT) / hb, -ZOOM_LINE_PT / hb
    ln, = axb.plot([x0, x0, x1, x1], [ya, yb_, yb_, ya], color=ZOOM_COLOUR, lw=0.8, transform=tr, clip_on=False,
                   zorder=3, solid_capstyle="butt", solid_joinstyle="miter")
    t = axb.text(0.5 * (x0 + x1), -(ZOOM_LINE_PT + ZOOM_TEXT_GAP_PT) / hb, r["label"], transform=tr, ha="center",
                 va="top", fontsize=C.SIZES[2], color=ZOOM_COLOUR, clip_on=False)
    return {"row_of_values": r, "line": ln, "text": t}


def _crosses(p, q, box):
    """True if the straight piece from p to q meets the box (x0, y0, x1, y1)"""
    (xa, ya), (xb, yb) = p, q
    x0, y0, x1, y1 = box
    t0, t1 = 0.0, 1.0
    dx, dy = xb - xa, yb - ya
    for pp, qq in ((-dx, xa - x0), (dx, x1 - xa), (-dy, ya - y0), (dy, y1 - ya)):
        if pp == 0:
            if qq < 0:
                return False
        else:
            t = qq / float(pp)
            if pp < 0:
                if t > t1:
                    return False
                t0 = max(t0, t)
            else:
                if t < t0:
                    return False
                t1 = min(t1, t)
    return True


def zoom_guides(fig, axc, axd, rows, heads, dpi_png):
    """ONE thin line that joins the range of panel d on the time axis of panel c to the two ends of the time axis
    of panel d: from the axis of c at the first day of the range straight down past the names of the months, slanting
    to the left end of the frame of panel d, across the whole width of panel d above its head, and back up, slanting
    and then straight, to the axis of c at the end of the range. The range is that of the row of values of the bracket
    under panel b. Stops if a piece comes nearer than GUIDE_CLEAR_OF_TEXT_PT to a text of the figure.
    Returns the line and its record (pixels of the PNG; rows counted from above)."""
    from matplotlib.lines import Line2D
    from matplotlib.text import Text
    r = [x for x in rows if x["element"] == E_ZOOM]
    if len(r) != 1:
        raise Stop("%d rows %r" % (len(r), E_ZOOM))
    x0, x1 = num(r[0]["x"]), num(r[0]["x_end"])
    lo, hi = axc.get_xlim()
    if not lo <= x0 < x1 <= hi:
        raise Stop("the range of panel d does not lie on the time axis of panel c")
    d0, d1 = axd.get_xlim()
    if abs(d0 - x0) > 1e-9 or abs(d1 - x1) > 1e-9:
        raise Stop("the time axis of panel d does not run over the range that is marked")
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    k = fig.dpi / 72.0
    months = [t for t in axc.get_xticklabels() if t.get_visible() and t.get_text().strip()]
    if not months or not heads:
        raise Stop("panel c has no names on its time axis, or panel d has no head")
    bc, bd = axc.get_window_extent(rend), axd.get_window_extent(rend)
    xa = [float(axc.transData.transform((x, 0.0))[0]) for x in (x0, x1)]
    left, right = float(bd.x0), float(bd.x1)
    y_axis = float(bc.y0)
    y_down = min(float(t.get_window_extent(rend).y0) for t in months) - GUIDE_BELOW_NAMES_PT * k
    y_across = max(float(t.get_window_extent(rend).y1) for t in heads) + GUIDE_ABOVE_HEAD_PT * k
    if y_down - y_across < GUIDE_LEAST_FALL_PT * k:
        raise Stop("between the names of the months of panel c and the head of panel d the slanting pieces have %.1f pt"
                   % ((y_down - y_across) / k))
    if not left < xa[0] < xa[1] < right:
        raise Stop("the range on the axis of panel c does not lie between the two ends of the frame of panel d")
    pts = [(xa[0], y_axis), (xa[0], y_down), (left, y_across), (right, y_across), (xa[1], y_down), (xa[1], y_axis)]
    xy = fig.transFigure.inverted().transform(pts)
    ln = Line2D([float(v) for v in xy[:, 0]], [float(v) for v in xy[:, 1]], transform=fig.transFigure,
                color=GUIDE_COLOUR, lw=GUIDE_LW_PT, solid_capstyle="butt", solid_joinstyle="round", zorder=1.5)
    fig.add_artist(ln)
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    pad = GUIDE_CLEAR_OF_TEXT_PT * k
    shown = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    crossed = []
    for t in shown:
        b = t.get_window_extent(rend)
        box = (float(b.x0) - pad, float(b.y0) - pad, float(b.x1) + pad, float(b.y1) + pad)
        for i in range(len(pts) - 1):
            if _crosses(pts[i], pts[i + 1], box):
                crossed.append({"text": t.get_text(), "piece": i + 1})
    if crossed:
        raise Stop("the guiding line comes too near to texts: %s" % crossed)
    s = dpi_png / float(fig.dpi)
    H = float(fig.bbox.height)
    fall = y_down - y_across
    return ln, {"dates": [r[0]["x"], r[0]["x_end"]],
                "points_px_from_the_upper_left": [[p[0] * s, (H - p[1]) * s] for p in pts],
                "x_on_the_axis_of_c_px": [v * s for v in xa], "ends_of_the_frame_of_d_px": [left * s, right * s],
                "row_of_the_axis_of_c": (H - y_axis) * s, "row_below_the_names": (H - y_down) * s,
                "row_across": (H - y_across) * s, "fall_of_the_slanting_pieces_pt": fall / k,
                "angles_deg": [math.degrees(math.atan2(fall, xa[0] - left)), math.degrees(math.atan2(fall, right - xa[1]))],
                "texts_checked": len(shown), "texts_crossed": crossed, "colour": GUIDE_COLOUR, "lw_pt": GUIDE_LW_PT}


def draw_panel_d(fig, place, rows, C, M, looks, n_pairs_b):
    """panel d from the rows of its file of values alone. Returns the axes and what was drawn, row by row."""
    import matplotlib as mpl
    from matplotlib.lines import Line2D
    from matplotlib.patches import Patch, Rectangle
    from matplotlib.transforms import offset_copy
    I = fig.interventions
    axa, axb, axc = fig.short_version["axa"], I["axb"], I["axc"]

    def of(element):
        return [r for r in rows if r["element"] == element]

    def one(element):
        got = of(element)
        if len(got) != 1:
            raise Stop("%d rows %r" % (len(got), element))
        return got[0]

    axd = fig.add_axes(place)
    X0, X1 = num(one(E_AX0)["x"]), num(one(E_AX1)["x"])
    rax = [r for r in of(E_RAX) if r["label"] == "axis of R, panel d"][0]
    y_top = float(rax["x_end"])
    axd.set_xlim(X0, X1)
    axd.set_ylim(float(rax["x"]), y_top)
    axd.set_yticks(list(range(0, int(round(y_top)) + 1)))
    for s in ("top", "right"):
        axd.spines[s].set_visible(False)
    axd.spines["left"].set_position(("outward", 5))
    ticks = of(E_TICK)
    major = [r for r in ticks if r["label"]]
    axd.set_xticks([num(r["x"]) for r in major])
    axd.set_xticklabels([r["label"] for r in major])
    axd.set_xticks([num(r["x"]) for r in ticks if not r["label"]], minor=True)
    size_ticks = float(axc.get_xticklabels()[0].get_fontsize())
    axd.tick_params(axis="both", which="major", labelsize=size_ticks)
    axd.set_xlabel(one(E_AX0)["label"], fontsize=axc.xaxis.label.get_fontsize())
    axd.set_ylabel(axb.yaxis.label.get_text(), fontsize=axb.yaxis.label.get_fontsize())
    axd.yaxis.set_label_coords(float(axb.yaxis.label.get_position()[0]), 0.5)
    fig.canvas.draw()
    frame_pt = axd.get_window_extent().height * 72.0 / fig.dpi
    if abs(frame_pt - FRAME_D_PT) > 1e-6:
        raise Stop("the frame of panel d is %.6f pt high" % frame_pt)
    per_day = axd.get_window_extent().width * 72.0 / fig.dpi / (X1 - X0)
    tr = axd.get_xaxis_transform()                      # x: the time axis; y: the frame, 0 below and 1 above

    def yf(pt):
        return 1.0 + pt / FRAME_D_PT

    drawn = {"primary": [], "weekly": [], "events": [], "periods": [], "names": [], "stems": []}
    # ----- the line R = 1, as panel b draws it
    line1 = [l for l in axb.lines if list(l.get_ydata()) == [1.0, 1.0]]
    if len(line1) != 1:
        raise Stop("panel b holds %d lines at R = 1" % len(line1))
    axd.axhline(1.0, color=line1[0].get_color(), lw=line1[0].get_linewidth(), zorder=line1[0].get_zorder())
    # ----- the pairs of the primary analysis: the fill and the thick line of the periods of panel b, without an outline
    med_b = [l for l in axb.lines if abs(l.get_linewidth() - 2.6) < 1e-9 and len(l.get_xdata()) == 2
             and l.get_ydata()[0] == l.get_ydata()[1]]
    fill_b = [p for p in axb.patches if isinstance(p, Rectangle) and p.get_linewidth() == 0.0
              and mpl.colors.to_hex(p.get_facecolor()) == mpl.colors.to_hex(C.BLUE)]
    if len(med_b) != 3 or len(fill_b) != 3:
        raise Stop("panel b holds %d thick lines and %d blue fills of periods" % (len(med_b), len(fill_b)))
    bars_b = [l for l in axb.lines if len(l.get_xdata()) == 2 and l.get_xdata()[0] == l.get_xdata()[1]
              and l.get_linestyle() == "-"]
    pale = sorted(set(mpl.colors.to_hex(l.get_color()) for l in bars_b))
    if len(bars_b) != n_pairs_b or len(pale) != 1:
        raise Stop("panel b holds %d bars of pairs in the colours %s; %d of one colour expected" % (len(bars_b), pale, n_pairs_b))
    fc, lw_med, col_med = mpl.colors.to_rgba(pale[0], BAND_ALPHA), med_b[0].get_linewidth(), pale[0]
    for r in of(E_PAIR_P):
        a, b = num(r["x"]), num(r["x_end"])
        lo, hi, v = float(r["lower"]), float(r["upper"]), float(r["value"])
        band = Rectangle((a, lo), b - a, hi - lo, facecolor=fc, edgecolor="none", linewidth=0.0, zorder=1.6)
        axd.add_patch(band)
        ln, = axd.plot([a, b], [v, v], color=col_med, lw=lw_med, zorder=1.9, solid_capstyle="butt")
        drawn["primary"].append((r, band, ln))
    # ----- the pairs of the analysis with intervals of 7.2 days
    for r in of(E_PAIR_W):
        x, lo, hi, v = num(r["x"]), float(r["lower"]), float(r["upper"]), float(r["value"])
        is_open = r["drawn_as"].startswith("open")
        bar, = axd.plot([x, x], [lo, hi], color=GREY_W, lw=0.8, zorder=2, solid_capstyle="butt")
        mk, = axd.plot([x], [v], marker="D", ms=MS_W, mfc="white" if is_open else GREY_W, mec=GREY_W, mew=0.8,
                       ls="none", zorder=2.2)
        drawn["weekly"].append((r, bar, mk))
    # ----- the periods: brackets with names above the frame
    name_b = [t for t in axb.texts if re.fullmatch(r"P\d", t.get_text())]
    brk_b = [l for l in axb.lines if len(l.get_xdata()) == 4]
    if len(name_b) != 3 or len(brk_b) != 3:
        raise Stop("panel b holds %d names and %d brackets of periods" % (len(name_b), len(brk_b)))
    inset = RULER_INSET_PT / per_day
    for r in of(E_PERIOD):
        a, b = num(r["x"]), num(r["x_end"])
        xa, xb = max(a, X0), min(b, X1)
        if not xb > xa:
            raise Stop("period %s does not reach into the range" % r["label"])
        cut_l, cut_r = a < X0, b > X1
        xs, ys = [], []
        if not cut_l:
            xs += [xa + inset, xa + inset]
            ys += [yf(RULER_LINE_PT - RULER_TICK_PT), yf(RULER_LINE_PT)]
        else:
            xs += [xa]
            ys += [yf(RULER_LINE_PT)]
        if not cut_r:
            xs += [xb - inset, xb - inset]
            ys += [yf(RULER_LINE_PT), yf(RULER_LINE_PT - RULER_TICK_PT)]
        else:
            xs += [xb]
            ys += [yf(RULER_LINE_PT)]
        kw = dict(color=brk_b[0].get_color(), transform=tr, clip_on=False, zorder=3)
        ln, = axd.plot(xs, ys, lw=brk_b[0].get_linewidth(), **kw)
        heads = []
        for x, cut, head in ((xa, cut_l, "<"), (xb, cut_r, ">")):
            if cut:
                h, = axd.plot([x], [yf(RULER_LINE_PT)], marker=head, ms=ARROW_PT, mec="none", ls="none", **kw)
                heads.append(h)
        t = axd.text(0.5 * (xa + xb), yf(NAME_BASE_PT), r["label"], transform=tr, ha="center", va="bottom",
                     fontsize=name_b[0].get_fontsize(), color=name_b[0].get_color(), fontweight="bold", clip_on=False)
        drawn["periods"].append((r, ln, heads, t, (xa, xb, cut_l, cut_r)))
    # ----- the events and the declaration: marks on one edge, names in rows on both sides
    col = looks["triangle"]["colour"]
    items = []
    for r in of(E_EVENT):
        x0, x1 = num(r["x"]), num(r["x_end"])
        items.append({"row_of_values": r, "mark": r["drawn_as"].split(":")[0], "x0": x0, "x1": x1,
                      "anchor": 0.5 * (x0 + x1), "text": r["label"], "colour": col})
    decl = one(E_DECL)
    words_a = [t for t in axa.texts if t.get_text() == decl["label"]]
    line_b = [l for l in axb.lines if l.get_linestyle() == "--" and len(l.get_xdata()) == 2
              and l.get_xdata()[0] == l.get_xdata()[1]]
    if len(words_a) != 1 or len(line_b) != 1:
        raise Stop("panel a holds %d texts %r; panel b holds %d dashed vertical lines" % (len(words_a), decl["label"], len(line_b)))
    if abs(float(line_b[0].get_xdata()[0]) - num(decl["x"])) > 1e-9:
        raise Stop("the dashed line of panel b does not stand at the day of the declaration")
    items.append({"row_of_values": decl, "mark": "declaration", "x0": num(decl["x"]), "x1": num(decl["x"]),
                  "anchor": num(decl["x"]), "text": decl["label"], "colour": words_a[0].get_color()})
    items.sort(key=lambda e: e["anchor"])
    half = (len(items) + 1) // 2
    size = C.SIZES[2]
    for j, e in enumerate(items):
        e["side"], e["row"] = ("left", j + 1) if j < half else ("right", len(items) - j)
        y_row = ROW0_PT + (e["row"] - 1) * ROW_PITCH_PT
        kw = dict(color=e["colour"], transform=tr, clip_on=False, zorder=4)
        arts = []
        if e["mark"] == "triangle":
            m, = axd.plot([e["x0"]], [yf(EVENT_BASE_PT + looks["triangle"]["ms"] / 2)], marker="v",
                          ms=looks["triangle"]["ms"], mec="none", ls="none", **kw)
            arts.append(m)
            foot = EVENT_BASE_PT + looks["triangle"]["ms"]
        elif e["mark"] == "cross":
            m, = axd.plot([e["x0"]], [yf(EVENT_BASE_PT + CROSS_UP_PT + looks["cross"]["ms"] / 2)], marker="x",
                          ms=looks["cross"]["ms"], mew=looks["cross"]["mew"], ls="none", **kw)
            arts.append(m)
            foot = EVENT_BASE_PT + CROSS_UP_PT + looks["cross"]["ms"] + 1.0
        elif e["mark"] == "bar":
            m, = axd.plot([e["x0"], e["x1"]], [yf(EVENT_BASE_PT + SPAN_Y_PT)] * 2, lw=looks["bar"]["lw"],
                          solid_capstyle="butt", **kw)
            arts.append(m)
            foot = EVENT_BASE_PT + 0.5
        elif e["mark"] == "dotted line":
            m, = axd.plot([e["x0"], e["x1"]], [yf(EVENT_BASE_PT + SPAN_Y_PT)] * 2, lw=looks["dotted line"]["lw"],
                          ls=(0, M.DOTS), dash_capstyle="butt", **kw)
            arts.append(m)
            foot = EVENT_BASE_PT + 0.5
        elif e["mark"] == "declaration":
            m, = axd.plot([e["anchor"], e["anchor"]], [0.0, yf(y_row)], color=line_b[0].get_color(), ls=(0, (4, 3)),
                          lw=line_b[0].get_linewidth(), transform=tr, clip_on=False, zorder=line_b[0].get_zorder())
            arts.append(m)
            foot = None
        else:
            raise Stop("no mark %r" % e["mark"])
        if foot is not None:
            st, = axd.plot([e["anchor"], e["anchor"]], [yf(foot), yf(y_row)], color=e["colour"], lw=STEM_LW_PT,
                           transform=tr, clip_on=False, zorder=3, solid_capstyle="butt")
            drawn["stems"].append(st)
        right = e["side"] == "right"
        t = axd.text(e["anchor"], yf(y_row), e["text"], ha="left" if right else "right", va="center", fontsize=size,
                     color=e["colour"], clip_on=False,
                     transform=offset_copy(tr, fig=fig, x=LABEL_PAD_PT if right else -LABEL_PAD_PT, y=0.0, units="points"))
        drawn["names"].append(t)
        drawn["events"].append((e["row_of_values"], e["mark"], arts, t, e["side"], e["row"]))
    # ----- the key and the head
    key_rows = of(E_KEY)
    title, labels = key_rows[0]["label"], [r["label"] for r in key_rows[1:]]
    handles = [(Patch(facecolor=fc, edgecolor="none"), Line2D([], [], color=col_med, lw=lw_med, solid_capstyle="butt")),
               Line2D([], [], color=GREY_W, lw=0.8, marker="D", ms=MS_W, mfc=GREY_W, mec=GREY_W, mew=0.8),
               Line2D([], [], color=GREY_W, lw=0.8, marker="D", ms=MS_W, mfc="white", mec=GREY_W, mew=0.8)]
    leg = axd.legend(handles, labels, loc="upper right", bbox_to_anchor=(1.0, 1.0), frameon=False,
                     fontsize=axb.get_legend().get_texts()[0].get_fontsize(), handlelength=2.2, borderaxespad=0.2,
                     labelspacing=0.35, handler_map={tuple: M.HandlerOver()}, title=title,
                     title_fontsize=axb.get_legend().get_texts()[0].get_fontsize(), alignment="left")
    head = one(E_HEAD)
    top_row = max(e[5] for e in drawn["events"])
    heads = C.panel_head(axd, head["label"], head["value"],
                         dy_pt=ROW0_PT + (top_row - 1) * ROW_PITCH_PT + TEXT_HALF_PT + HEAD_GAP_PT)
    drawn.update(key=leg, heads=heads, per_day_pt=per_day, frame_pt=frame_pt, rows_of_names=top_row)
    return axd, drawn


def read_back(drawn):
    """what stands in the drawing of panel d, held against the rows of the file of values; a list of faults"""
    import matplotlib as mpl
    bad = []

    def near(a, b, what, tol=1e-9):
        if abs(float(a) - float(b)) > tol:
            bad.append("%s: drawn %r, file of values %r" % (what, float(a), float(b)))

    for r, band, ln in drawn["primary"]:
        a, b = num(r["x"]), num(r["x_end"])
        near(band.get_x(), a, r["label"] + ", left edge of the band")
        near(band.get_x() + band.get_width(), b, r["label"] + ", right edge of the band")
        near(band.get_y(), r["lower"], r["label"] + ", lower bound")
        near(band.get_y() + band.get_height(), r["upper"], r["label"] + ", upper bound")
        near(ln.get_xdata()[0], a, r["label"] + ", line from")
        near(ln.get_xdata()[1], b, r["label"] + ", line to")
        near(ln.get_ydata()[0], r["value"], r["label"] + ", median")
        near(ln.get_ydata()[1], r["value"], r["label"] + ", median")
    for r, bar, mk in drawn["weekly"]:
        x = num(r["x"])
        near(bar.get_xdata()[0], x, r["label"] + ", time of the line")
        near(bar.get_xdata()[1], x, r["label"] + ", time of the line")
        near(bar.get_ydata()[0], r["lower"], r["label"] + ", lower bound")
        near(bar.get_ydata()[1], r["upper"], r["label"] + ", upper bound")
        near(mk.get_xdata()[0], x, r["label"] + ", time of the symbol")
        near(mk.get_ydata()[0], r["value"], r["label"] + ", median")
        is_open = mpl.colors.to_hex(mk.get_markerfacecolor()) == "#ffffff"
        if is_open != r["drawn_as"].startswith("open"):
            bad.append("%s: the symbol is %s" % (r["label"], "open" if is_open else "filled"))
    for r, mark, arts, t, side, row in drawn["events"]:
        xs = list(arts[0].get_xdata())
        want = [num(r["x"])] if mark in ("triangle", "cross") else [num(r["x"]), num(r["x"])] if mark == "declaration" \
            else [num(r["x"]), num(r["x_end"])]
        if len(xs) != len(want):
            bad.append("%s: %d points drawn" % (r["label"], len(xs)))
        else:
            for a, b in zip(xs, want):
                near(a, b, r["label"] + ", day of the mark")
        if t.get_text() != r["label"]:
            bad.append("%s: the name drawn is %r" % (r["label"], t.get_text()))
    for r, ln, heads, t, (xa, xb, cut_l, cut_r) in drawn["periods"]:
        if t.get_text() != r["label"]:
            bad.append("period %s: the name drawn is %r" % (r["label"], t.get_text()))
        near(t.get_position()[0], 0.5 * (xa + xb), "period %s, middle of the name" % r["label"])
        if len(heads) != int(cut_l) + int(cut_r):
            bad.append("period %s: %d arrow heads" % (r["label"], len(heads)))
    return bad


# ---------------------------------------------------------------------------------------------------------------
# the run
# ---------------------------------------------------------------------------------------------------------------
def run(event_program, short_program, code_dir, lanes_path, values_path, s2_values_path, outdir, stem, analyses,
        labels, inp):
    from matplotlib.text import Text
    P = load(event_program, "fig1_event_legend_stored")
    S = P.load(short_program, "fig1_short_version_stored")
    # ----- step 1: P.draw, lines 466 to 481
    W = S.load_lanes(lanes_path)
    M = W.install(W.load(code_dir))
    C = M.C
    state = {}
    S.install_search(M, W, P.MAX_LANES, state, P.STRIP_SCALE)
    sets = S.beast_rows(s2_values_path, analyses, C.Stop)
    fig, extra = M.draw(values_path)
    tick_size = P.tick_labels(fig, C.SIZES[0])
    root = P.root_date_down(fig, P.ROOT_DOWN_PT)
    axa, floor, drawn_beast = S.add_beast(M, fig, sets)
    delphy = [r for r in C.read_csv(values_path) if r["element"] == M.E_STAIR]
    fig.short_version = {"axa": axa, "beast": drawn_beast, "legend": None,
                         "delphy": [{"x": M.D(r["x"]), "x_end": M.D(r["x_end"]), "value": float(r["value"]),
                                     "lower": float(r["lower"]), "upper": float(r["upper"])} for r in delphy]}
    room = P.periods_and_room(fig, C.DPI)
    I = fig.interventions
    axb, axc = I["axb"], I["axc"]
    # ----- step 2
    rows, counts = values_of_panel_d(inp, list(P.WORDS), list(P.NAMES), stem)
    os.makedirs(outdir, exist_ok=True)
    values_d = os.path.join(outdir, stem + "_values_panel_d.csv")
    write_values(rows, values_d)
    rows = read_csv(values_d)                              # the drawing is made from the file as it is written
    new_room, place, looks = strip_away_and_room(fig, C, C.DPI)
    n_pairs_b = sum(1 for r in C.read_csv(values_path) if r["element"] == "R of a pair of adjacent Skygrid intervals")
    axd, drawn = draw_panel_d(fig, place, rows, C, M, looks, n_pairs_b)
    zoom = range_mark_under_b(fig, rows, C)
    # ----- step 3: P.draw, lines 482 to 485; P.run, lines 572 to 580
    moved = S.clear_the_words(M, fig)
    leg, key = S.key_above(M, fig, axa, labels)
    fig.short_version["legend"] = leg
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    k = fig.dpi / 72.0
    numb = sorted([t for t in axb.texts if re.match(r"^\d+\.\d+ \(", t.get_text())],
                  key=lambda t: t.get_window_extent(rend).x0)
    nb = [t.get_window_extent(rend) for t in numb]
    gaps = [(nb[i + 1].x0 - nb[i].x1) / k for i in range(len(nb) - 1)]
    if len(numb) != 3 or min(gaps) < P.LEAST_BETWEEN_PERIOD_LABELS_PT:
        raise Stop("the labels of the period numbers (%d) stand %s pt apart" % (len(numb), [round(g, 1) for g in gaps]))
    names = P.names_in_bold(fig)
    keys = P.keys_on_white(S, M, fig)
    # ----- the guiding line between the time axis of panel c and the two ends of panel d (nothing moves after this)
    guide_line, guides = zoom_guides(fig, axc, axd, rows, drawn["heads"], C.DPI) if GUIDES else (None, None)
    # ----- step 4: the checks
    fig.canvas.draw()
    overlap = P.texts_that_overlap(fig)
    white = {"key of panel b": P.not_white_under(fig, [axb.get_legend()], S.WHITE_FROM),
             "key of panel c": P.not_white_under(fig, [axc.get_legend()], S.WHITE_FROM),
             "key above panel a": P.not_white_under(fig, [fig.short_version["legend"]], S.WHITE_FROM),
             "key of panel d": P.not_white_under(fig, [drawn["key"]], S.WHITE_FROM)}
    own_texts = list(drawn["names"]) + [p[3] for p in drawn["periods"]] + list(drawn["heads"]) + [zoom["text"]] \
        + [t for t in axb.texts if hasattr(t, "xyann") and t.get_text() in ("b", HEAD_B)]
    under = []
    for t in own_texts:
        w = P.not_white_under(fig, [t], S.WHITE_FROM, pad_pt=0.5)
        under.append({"text": t.get_text(), "pixels_not_white": w["pixels_not_white"],
                      "inside": w["box_widened_inside_the_picture"], "box_pt": w["box_pt"]})
    back = read_back(drawn)
    zx = [float(v) for v in zoom["line"].get_xdata()]
    if zx != [num(zoom["row_of_values"]["x"])] * 2 + [num(zoom["row_of_values"]["x_end"])] * 2 \
            or zoom["text"].get_text() != zoom["row_of_values"]["label"]:
        back.append("the mark of the range under panel b does not stand as its row gives it")
    from matplotlib import font_manager as fm
    fonts = {}
    for t in fig.findobj(Text):
        if t.get_visible() and t.get_text().strip():
            fonts.setdefault(os.path.basename(fm.findfont(t.get_fontproperties())), set()).add(str(t.get_fontweight()))
    fonts = {k_: sorted(v) for k_, v in fonts.items()}
    touches_at_the_end = S.touches(fig)
    texts = P.boxes_of_the_texts(fig)
    shown = [t for t in fig.findobj(Text) if t.get_visible() and t.get_text().strip()]
    sizes_of_the_texts = sorted(set(float(t.get_fontsize()) for t in shown))
    rend = fig.canvas.get_renderer()
    pos = {n_: [float(v) for v in a.get_position().bounds] for n_, a in (("a", axa), ("b", axb), ("c", axc), ("d", axd))}
    heads_pt = {}
    for ax_, n_ in ((axa, "a"), (axb, "b"), (axc, "c"), (axd, "d")):
        hs = [t for t in ax_.texts if hasattr(t, "xyann") and tuple(t.xy) == (0, 1) and t.xycoords == "axes fraction"]
        heads_pt[n_] = [{"text": t.get_text(), "box_pt": [float(v) / k for v in t.get_window_extent(rend).extents]} for t in hs]
    placed = [{"label": r["label"], "mark": mark, "side": side, "row": row,
               "box_of_the_name_pt": [float(v) / k for v in t.get_window_extent(rend).extents]}
              for r, mark, arts, t, side, row in drawn["events"]]
    # ----- step 5: the saving
    root_measured = P.root_date_measured(fig) if hasattr(P, "root_date_measured") else None
    limit_of_the_conventions = C.MAX_HEIGHT
    C.MAX_HEIGHT = max(C.MAX_HEIGHT, HEIGHT_LIMIT_IN)
    try:
        saved = C.save_figure(fig, stem, outdir)
    finally:
        C.MAX_HEIGHT = limit_of_the_conventions
    saved["limit_of_the_conventions_in"] = float(limit_of_the_conventions)
    saved["limit_for_this_figure_in"] = float(HEIGHT_LIMIT_IN)
    size_in = [float(v) for v in fig.get_size_inches()]
    M.plt_close(fig)
    from PIL import Image
    parts = {}
    with Image.open(saved["png"]) as im:
        wpx, hpx = im.size
        cut = new_room["rows_of_panels_a_to_c"]
        if hpx != new_room["rows"] or not 0 < cut < hpx:
            raise Stop("the saved picture has %d rows, %d expected" % (hpx, new_room["rows"]))
        for name, box in (() if guides is not None else (("panels_a_to_c", (0, 0, wpx, cut)), ("panel_d", (0, cut, wpx, hpx)))):
            path = os.path.join(outdir, "%s_%s.png" % (stem, name))
            im.crop(box).save(path, dpi=(C.DPI, C.DPI))
            parts[name] = {"png": path, "rows": [box[1], box[3]], "width_px": wpx, "height_px": box[3] - box[1]}
    return {"saved": dict(saved), "parts": parts, "guides": guides, "size_in": size_in, "values_of_panel_d": values_d,
            "rows_of_values": len(rows), "counts": counts, "room_of_the_stored_program": room, "room": new_room,
            "place_of_panel_d": place, "key_above_a": key, "words_moved": moved, "names": names, "keys": keys,
            "tick_labels_pt": tick_size, "root_date": root, "root_date_measured": root_measured,
            "period_labels": {"texts": [t.get_text() for t in numb], "gaps_pt": gaps},
            "overlap": overlap, "white": white, "white_under_the_texts_of_panel_d_and_the_head_of_b": under,
            "read_back": back, "touches_at_the_end": touches_at_the_end, "texts": texts,
            "sizes_of_the_texts": sizes_of_the_texts, "font_files_by_weight": fonts, "positions_of_the_axes": pos,
            "mark_of_the_range_under_b": {"words": zoom["text"].get_text(),
                                          "box_pt": [float(v) / k for v in zoom["text"].get_window_extent(rend).extents]}, "heads_pt": heads_pt,
            "events_placed": placed, "per_day_pt": drawn["per_day_pt"], "frame_of_panel_d_pt": drawn["frame_pt"],
            "rows_of_names": drawn["rows_of_names"],
            "drawn": {"primary": len(drawn["primary"]), "weekly": len(drawn["weekly"]), "events": len(drawn["events"]),
                      "periods": len(drawn["periods"]), "stems": len(drawn["stems"])},
            "keys_or_words_over_data_by_the_track": extra.get("keys_or_words_over_data"),
            "strip_of_the_stored_program": state.get("sizes"),
            "constants": {"SHIFT_A_PX": SHIFT_A_PX, "GROWTH_PX": GROWTH_PX, "HEIGHT_LIMIT_IN": HEIGHT_LIMIT_IN,
                          "WINDOW": list(WINDOW), "R_TOP_D": R_TOP_D, "FRAME_D_PT": FRAME_D_PT,
                          "PRIMARY_GRID": list(PRIMARY_GRID), "WEEKLY_GRID": list(WEEKLY_GRID), "TICKS_D": list(TICKS_D),
                          "MINOR_TICKS_D": list(MINOR_TICKS_D), "HEAD_B": HEAD_B, "HEAD_D": HEAD_D,
                          "KEY_D_TITLE": KEY_D_TITLE, "KEY_D": list(KEY_D), "GREY_W": GREY_W, "MS_W": MS_W,
                          "BAND_ALPHA": BAND_ALPHA, "CROSS_UP_PT": CROSS_UP_PT, "ZOOM_WORDS": ZOOM_WORDS,
                          "ZOOM_LINE_PT": ZOOM_LINE_PT, "ZOOM_TICK_PT": ZOOM_TICK_PT,
                          "ZOOM_TEXT_GAP_PT": ZOOM_TEXT_GAP_PT, "ZOOM_COLOUR": ZOOM_COLOUR,
                          "RULER_LINE_PT": RULER_LINE_PT, "RULER_TICK_PT": RULER_TICK_PT, "RULER_INSET_PT": RULER_INSET_PT,
                          "NAME_BASE_PT": NAME_BASE_PT, "EVENT_BASE_PT": EVENT_BASE_PT, "SPAN_Y_PT": SPAN_Y_PT,
                          "ROW0_PT": ROW0_PT, "ROW_PITCH_PT": ROW_PITCH_PT, "LABEL_PAD_PT": LABEL_PAD_PT,
                          "STEM_LW_PT": STEM_LW_PT, "ARROW_PT": ARROW_PT, "HEAD_GAP_PT": HEAD_GAP_PT,
                          "TEXT_HALF_PT": TEXT_HALF_PT,
                          "GUIDES": GUIDES, "GUIDE_COLOUR": GUIDE_COLOUR, "GUIDE_LW_PT": GUIDE_LW_PT,
                          "GUIDE_BELOW_NAMES_PT": GUIDE_BELOW_NAMES_PT, "GUIDE_ABOVE_HEAD_PT": GUIDE_ABOVE_HEAD_PT,
                          "GUIDE_LEAST_FALL_PT": GUIDE_LEAST_FALL_PT, "GUIDE_CLEAR_OF_TEXT_PT": GUIDE_CLEAR_OF_TEXT_PT}}
