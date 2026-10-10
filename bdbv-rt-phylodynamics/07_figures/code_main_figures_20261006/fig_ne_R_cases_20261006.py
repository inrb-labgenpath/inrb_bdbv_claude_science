"""fig_ne_R_cases_20261006.py
Figure 1 of the article (effective population size, reproduction number, confirmed cases), drawn from the
tables of the round of October 2026, with the table of the response events, which are shown in a strip above
panel b. The figure for the article; its picture holds no banner. This is the program of 5 October
(fig_ne_R_cases_20261005b.py) with changed names of the files and a changed sentence for the legend; the
drawing of the events is that of 5 October, and the rest of the drawing is that of 4 October.

    python fig_ne_R_cases_20261006.py INPUTS.json OUTDIR --interventions NAME.csv [--sources KEYS.csv]   the figure of the article
    python fig_ne_R_cases_20261006.py INPUTS.json [OUTDIR]                           the figure without the strip
    python fig_ne_R_cases_20261006.py INPUTS.json OUTDIR --interventions NAME.csv [--sources KEYS.csv] --preview
    python fig_ne_R_cases_20261006.py --draw-only VALUES.csv OUTDIR                  draws from a file of values alone

NAME.csv and KEYS.csv are entries of INPUTS.json ('file name -> path'); no path stands in this program, and
no day, label, kind or source of a table is typed here. Nothing depends on the number of the rows or on their
days. The name of a table ends with _<tag>_<YYYYMMDD>.csv; 'INTERNAL_', 'interventions_' and 'figure1b_'
before the tag are left out of the names of the files that are written. A file of a test with fixtures is
called INTERNAL_TEST_WITH_FIXTURES_<...>.csv; every label of it begins with 'FIXTURE ', and its picture
carries two banners. A file of another name that holds the word FIXTURE stops the program. The version of a
table is that of the table of inputs, or of the entry 'NAME.csv#version' of INPUTS.json.

THE TABLE has one row per intervention, in one of two forms, told by its header. A fault stops the program;
nothing is repaired. What is read, and the column that holds it:
                                   first form                  second form
    first day (YYYY-MM-DD)         first_day                   first_day
    last day                       last_day                    last_day
    what a period means            meaning_of_the_period *     meaning_of_the_period
    kind of the intervention       kind *                      kind
    short label (48 characters)    label                       short_label
    words for the legend           words_for_the_legend        words_for_the_legend
    source: number of a reference  source_reference_number     -
    source: citation               source_citation             -
    source: keys                   source_keys *               sources
    order of the rows              -                           order (counts the rows from 1)
    row of the table it rests on   -                           row_of_the_checked_table
    note (not read by the drawing) note *                      -                        (* optional column)
last_day empty or equal to first_day: one day. A later last_day: a period, and meaning_of_the_period says
what it is; a meaning that holds 'not known' or 'not a duration' is two days between which the day lies (an
interval of uncertainty), every other meaning is a duration. One day has no meaning. Exactly one of the
three sources is filled; keys are separated by '; ' and are printed as the table has them. With --sources
the keys are looked up in the column 'key' of KEYS.csv, and a key that it does not hold fails a check.
<<BR>> in a label fixes the breaks of its lines.

THE MARKS stand at the lower edge of a strip above panel b that shares the time axis, in ONE colour that the
three panels do not use:
    one day                       a triangle with its tip at the day; a cross where the kind begins with
                                  'attack'
    a duration                    a bar UNDER the lower edge from the first to the last day, a tick at each end
    two days between which        a dotted line under the lower edge between the two days, without ticks
    the day lies
    an end beyond the time axis   is cut at the axis and carries an arrow head
The kind is shown by the mark of one day only; the label and the sentence for the legend say what a period
is. Rule 1 (default, 'labels in words'): a thin line rises at the first day to the label, and the label
stands in words beside the line, in at most three lines, to the right or to the left, at one of at most three
heights; a search takes the fewest heights at which no label touches another label, a line of another mark
or the line of the declaration. Rule 2 (fallback, 'numbered marks'), where rule 1 finds no place in the room
that the height of the figure allows: numbers above the marks, by first day and then by the order of the
rows (numbers that would touch are printed as one text, '2, 3' or '4' to '6' with a dash); the words stand in
the sentence for the legend. The checks say which rule was used and why and, where numbers were used, how
many heights the labels in words would need. The line of the declaration runs through the strip as it runs
through the panel; the declaration is no row of the table. Nothing is drawn inside panel b: no mark can hide
a rectangle, a point, the line R = 1 or the line of the declaration. An intervention that lies wholly outside
the time axis is not drawn; it is held in the file of values and listed in the checks.
Panels a, b and c keep their size and their content (a check compares them pixel by pixel with the picture
without interventions). The strip takes its room from the height of the figure (up to the maximum of the
conventions), then from the gap between a and b, then from the gap between b and c; the title of panel b
stands above the strip.

WHAT IS WRITTEN. WITH the table of the events (--interventions; neither a test nor a preview), THE FIGURE OF
THE ARTICLE:
    fig_ne_R_cases_20261006_values.csv                  (BEFORE anything is drawn)
    fig_ne_R_cases_20261006.png/.pdf                    (drawn from the file of values and from nothing else;
                                                         no banner)
    fig_ne_R_cases_20261006_checks.json
    fig_ne_R_cases_20261006_legend_sentence_events.txt  (the sentence for the legend, in the marked form of the
                                                         template; written where an event is drawn)
A table of the article that has no row gives the same files without the strip and without the sentence (a
test or a preview with a table that has no row writes the files of the run without a table).
WITHOUT a table, the figure without the strip, which the checks and the regression compare:
    INTERNAL_fig_ne_R_cases_without_events_20261006_values.csv, INTERNAL_fig_ne_R_cases_without_events_20261006
    .png/.pdf and INTERNAL_fig_ne_R_cases_without_events_20261006_checks.json.
With --preview <stem> = fig_ne_R_cases_<tag>_20261006, every name begins with INTERNAL_PREVIEW_ and the picture
carries two banners; for a test with fixtures <stem> = fig_ne_R_cases_<tag>_20261006 and every name begins
with INTERNAL_TEST_WITH_FIXTURES_ (two banners); both write <...>_legend_sentence_events.txt.
The rows of 4 October stand unchanged at the head of the file of values (but for the column figure, which
holds the name of this figure); the rows of the events follow and name the table, its version, the row and
the column of every cell. With a table that has rows the program returns 1 where a check fails (the run
without a table returns 0).
THE SENTENCE FOR THE LEGEND. Where the rule of the marks is 'numbered marks' the sentence gives for each mark
its number, the words for the legend, the day or the days, and the sources, WITHOUT the short label, which
the figure then does not show; where the rule is 'labels in words' it gives the label, the words, the day or
the days and the sources. The explanation of the marks names the kinds and the forms that the table holds and
no other.

The rules of drawing are those of fig_ne_R_cases_20260924.py (the earlier generator, which is NOT run and of
which no table is read), with the changes of 4 October (R4); every rule is listed in
INTERNAL_main_figures_mapping_of_inputs_20261006.csv, which holds the rules of the strip as well. No value
is computed (R1).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import main_fig_common_20261006 as C          # sets the thread variables before numpy is imported

import datetime as dt
import hashlib
import json
import math
import re
import tempfile

import numpy as np
import matplotlib as mpl
import matplotlib.dates as mdates
from matplotlib.patches import Rectangle, Patch
from matplotlib.lines import Line2D
from matplotlib.legend_handler import HandlerTuple

NAME = "fig_ne_R_cases" + C.SUFFIX
# the figure without the strip (a run without a table): compared by the checks and the regression
NAME_WITHOUT = "fig_ne_R_cases_without_events" + C.SUFFIX
VALUES_FILE = "INTERNAL_" + NAME_WITHOUT + "_values.csv"
FIGURE_STEM = "INTERNAL_" + NAME_WITHOUT
CHECKS_FILE = "INTERNAL_" + NAME_WITHOUT + "_checks.json"
LEGEND_SUFFIX = "_legend_sentence_events.txt"
ARTICLE_NAMES = {"values": NAME + "_values.csv", "figure": NAME, "checks": NAME + "_checks.json",
                 "legend": NAME + LEGEND_SUFFIX}
# words that the name of a file of the article does not hold (compared without regard to case)
NOT_IN_A_NAME_OF_THE_ARTICLE = ("INTERNAL", "PROVISIONAL", "PREVIEW", "FIXTURE", "TEST")

# ---- the table of the events (task of 5 October; the figure of the article since 6 October) -------
PROGRAM_DATE = C.SUFFIX
MARKER = "INTERNAL_TEST_WITH_FIXTURES"
PREVIEW = "INTERNAL_PREVIEW"
FIXTURE_WORD = "FIXTURE"
TABLE_NAME = re.compile(r"^(?:(" + MARKER + r")_)?(?:INTERNAL_)?(?:interventions_|figure1b_)?([A-Za-z0-9_]+?)"
                        r"_(\d{8})\.csv$")
# the two forms of the table: what the program reads -> the column of the file that holds it
FORMS = {
    "first form": {
        "told_by": "label",
        "required": {"first_day": "first_day", "last_day": "last_day", "label": "label",
                     "words": "words_for_the_legend", "reference": "source_reference_number",
                     "citation": "source_citation"},
        "optional": {"keys": "source_keys", "kind": "kind", "meaning": "meaning_of_the_period", "note": "note"}},
    "second form": {
        "told_by": "short_label",
        "required": {"order": "order", "first_day": "first_day", "last_day": "last_day",
                     "meaning": "meaning_of_the_period", "kind": "kind", "label": "short_label",
                     "words": "words_for_the_legend", "keys": "sources", "origin": "row_of_the_checked_table"},
        "optional": {}}}
SOURCES_KEY_COLUMN = "key"                 # column of the table of sources that holds the keys
KIND_WITH_A_CROSS = re.compile(r"^attack\b", re.I)          # a kind that begins with 'attack': a cross
DAY_NOT_KNOWN = re.compile(r"not a duration|not known", re.I)   # meaning of a period: the day is not known
LABEL_MAX, WORDS_MAX = 48, 200             # characters of a label ('about 45'), of the words and of a citation
LINES_MAX, LINE_MAX, LINE_MAX_ONE = 3, 26, 18   # lines of a label; characters of a line; of a label of one line
UNIT_OF_THE_LEGEND = "P022"                # unit of the template that holds the legend of Figure 1
# the strip above panel b; every length in points (1/72 inch)
STRIP_COLOUR = "#762a83"                   # a PROPOSAL: one colour that panels a, b and c do not use
COLOUR_DISTANCE_MIN = 40.0                 # of 255: smallest difference to a colour of the panels, conventions
LINESPACING = 1.1
STRIP_GAP_PT = 3.0                         # from the upper edge of panel b to the lower edge of the strip
STRIP_TOP_PT = 1.5                         # above the highest text of the strip
MARK_PT = 4.5                              # height and width of a triangle and of a cross, upper end of a tick
CROSS_LW_PT = 1.1                          # width of the lines of a cross
SPAN_Y_PT = -1.2                           # middle of a bar and of a dotted line: UNDER the lower edge
SPAN_LW_PT = 1.2                           # width of a bar and of a dotted line
DOTS = (1.0, 1.3)                          # a dot and the gap after it, in widths of the line
ARROW_PT = 3.0                             # arrow head at an end that is cut at the axis
CAP_LW_PT, STEM_LW_PT = 1.0, 0.6           # tick at an end of a bar; thin line and lower edge of the strip
LEVEL0_PT = 6.5                            # lower edge of a label of the lowest height, above the lower edge
LEVELS_MAX = 3                             # 'two or three heights'
LABEL_PAD_PT = 2.5                         # between the thin line and its label
LABEL_GAP_PT = 2.0                         # smallest distance between two labels
LINE_CLEAR_PT = 1.5                        # smallest distance between a label and a line that is not its own
NUMBER_Y_PT = MARK_PT + 1.5                # lower edge of a number (fallback rule)
GAP_GIVES_PT = 11.5                        # what the gap between two panels gives to the strip at most
HATCH_STEP_PX = 25                         # the figure grows in steps of 25 pixels of the PNG (6 pt): the period
#                                            of the hatching '////' of panel c, which the backends anchor at the
#                                            top of the picture; so the hatching of panel c keeps its place
NODES_MAX = 500000                         # placements that the search may try
HEIGHTS_ASKED_ABOUT = 12                   # where numbers are used: up to how many heights the report asks about
TIP_ROW_PT = 0.9                           # row of pixels in which the tips of the triangles are read back
PIXEL_TOLERANCE = 1.6                      # pixels of the PNG (300 dpi): 0.4 pt, less than a quarter of a day
PIXEL_TOLERANCE_CROSS = 2.6                # of a cross: its lines move by one pixel from one row of pixels
#                                            to the next
PIXEL_TOLERANCE_HEAD = 3.2                 # at the point of an arrow head (an end that is cut at the axis):
#                                            the point fades over about two pixels
PANEL_PAD_PT = 5.0                         # below a panel: its ticks
TITLE_BLOCK_PT = 21.0                      # above a panel: its letter and title
BANNER_TOP = "TEST WITH FIXTURES \u2013 invented dates and labels"
BANNER_BOTTOM = "TEST WITH FIXTURES \u2013 not for the article"
BANNERS = {"test": (BANNER_TOP, BANNER_BOTTOM),
           "preview": ("PREVIEW \u2013 estimates provisional; events as the user chose them",
                       "PREVIEW \u2013 not the figure of the article")}
BANNER_COLOUR = "#b00020"
BANNER_PAD_PT = 2.0                        # around a banner: left out where panels are compared pixel by pixel
DRAWN = {"triangle": "triangle on the lower edge of the strip above panel b, tip at the day",
         "cross": "cross on the lower edge of the strip above panel b, at the day",
         "duration": "bar under the lower edge of the strip above panel b, from the first to the last day, a "
                     "tick at each end",
         "day not known": "dotted line under the lower edge of the strip above panel b, between the two days: the day "
                   "is not known"}
NOT_DRAWN_OUTSIDE = "not drawn: outside the time axis"

# elements of the file of values (the drawing finds its rows by these words)
E_AXIS_FIRST = "time axis, first day"
E_AXIS_END = "time axis, end"
E_DECLARED = "declaration of the outbreak"
E_STAIR = "effective population size, staircase"
E_GENOMES = "genomes of the analysed set per Skygrid interval"
E_ROOT = "date of the root (tMRCA), bar and printed dates"
E_FADE = "limit of the hatched part (median date of the root)"
E_PAIR = "R of a pair of adjacent Skygrid intervals"
E_PERIOD = "R of a period, posterior"
E_PERIOD_TEXT = "R of a period, value printed above the panel"
E_PERIOD_NAME = "name of a period, printed above the panel"
E_PRIOR_SIM = "R of a period, smoothing prior alone, quantiles of simulated trajectories"
E_PRIOR_SD = "R of a period, smoothing prior alone, range from the standard deviation"
E_CASES = "confirmed cases per week, by date of report"
E_HATCH = "hatching of a week whose split by province is pooled"
E_PROVISIONAL = "dashed outline of the provisional week"
E_KEY = "entry of a key"
E_AXIS_RULE = "limits of an axis by rule"
E_INT_MARK = "intervention, mark in the strip above panel b"
E_INT_NUMBER = "intervention, number of the mark"
E_INT_WORDS = "intervention, words for the legend"
E_INT_SOURCE = "intervention, source"
E_INT_KIND = "intervention, kind"
E_INT_MEANING = "intervention, meaning of the period"
E_INT_ORDER = "intervention, order in the table"
E_INT_ORIGIN = "intervention, row of the checked table of events"
E_BANNER = "banner of a test with fixtures or of a preview"

# positions that the earlier generator holds as typed numbers; the rule of each is in the mapping
HEIGHT = 7.9
GRID = dict(height_ratios=[2.45, 2.85, 1.75], left=0.105, right=0.925, top=0.962, bottom=0.062, hspace=0.27)
YLAB_X = -0.083
GENOME_AXIS_OVER_SPINE = 300.0 / 125.0     # earlier: axis to 300 where the spine ends at 125
STRIP_SHARE = 0.80 / 4.0                   # earlier: strip of 0.80 above an R axis that ends at 4
NE_ROOM_BELOW = 0.004 / 0.01               # earlier: axis from 0.004 where the lowest decade is 0.01
PROV_TEXT_SHARE = 705.0 / 900.0            # earlier: word 'provisional' at 705 of an axis to 900
PROV_LINE_SHARE = 690.0 / 900.0
PROV_GAP_SHARE = 14.0 / 900.0
PROVINCES = [("ituri", "Ituri", C.ORANGE_ITURI), ("nord_kivu", "Nord-Kivu", C.ORANGE_NK),
             ("other_provinces", "other provinces", C.ORANGE_OTHER)]


# =============================================================================================
# 1. the file of values
# =============================================================================================
def primary_analysis(I):
    """The analysis that the LATEST specification names for the quantities PRIM.* of the tokens (the primary
    analysis)."""
    t, _ = C.tokens_by_the_latest_specification(I)
    t.need("quantity_id", "specification: key_analysis")
    names = sorted({r["specification: key_analysis"] for r in t.rows if r["quantity_id"].startswith("PRIM.")
                    and r["specification: key_analysis"]})
    if len(names) != 1:
        raise C.Stop(f"the latest specification names {len(names)} analyses for the quantities PRIM.* of the "
                     f"tokens: {names}")
    return names[0]


def token_text(I, quantity_id):
    """(text, 0-based row) of a quantity that is unchanged by design, from the table of tokens."""
    t = I.table(C.TABLE_OF_TOKENS).need("quantity_id", "printed_where_the_quantity_is_unchanged_by_design",
                                         "status_of_the_quantity_in_the_specification")
    rows = [i for i, r in enumerate(t.rows) if r["quantity_id"] == quantity_id]
    texts = sorted({t.rows[i]["printed_where_the_quantity_is_unchanged_by_design"] for i in rows})
    if not rows or len(texts) != 1 or not texts[0]:
        raise C.Stop(f"table of tokens: {quantity_id}: no single text 'unchanged by design' ({texts})")
    return texts[0], rows[0]


def quantity_ids(I, analysis, kind, period, end_used, columns, variant=None, row=None):
    """quantity_id(s) of the LATEST specification whose key (definition for the program, op 'cell') names this
    analysis, this kind of quantity, the day(s) of this period, this end and these columns; '' where the
    specification holds none. The latest specification names a period by its first day (and by its last day
    where the key holds one), not by its name: `row` is the row of the table that holds the quantity, and
    its cells period_start and period_end are compared with the key."""
    t = I.table(C.SPECIFICATION).need("quantity_id", "definition_for_the_program")
    if period is not None and row is None:
        raise C.Stop(f"quantity_ids: {analysis}, {kind}, {period}: the row of the table is needed")
    out = []
    for r in t.rows:
        d = r["definition_for_the_program"].strip()
        if not d.startswith("{"):
            continue
        j = json.loads(d)
        w = j.get("where") or {}
        if j.get("op") != "cell" or w.get("analysis") != analysis or "kind" not in w:
            continue
        k = C.kind_of_the_key(w["kind"], r["quantity_id"])
        if k["kind"] != kind or k.get("variant") != variant:
            continue
        if period is not None:
            if w.get("period_start") != row["period_start"]:
                continue
            if "period_end" in w and w["period_end"] != row["period_end"]:
                continue
        elif "period_start" in w or "period_end" in w:
            continue
        if "end_used" in w and w["end_used"] != end_used:
            continue
        if list(j.get("columns", [])) != list(columns):
            continue
        out.append(r["quantity_id"])
    return "; ".join(out)


def build_values(I):
    V = C.Values(NAME, I)
    info = {}
    primary = primary_analysis(I)
    f_est = I.estimates_table_of(primary)
    if f_est is None:
        raise C.Stop(f"the primary analysis {primary} has no delivered table of estimates")
    est = I.table(f_est).need("analysis", "quantity", "period_start", "period_end", "end_used", "median",
                              "hpd95_lower", "hpd95_upper", "meets_convergence_criterion", "genomes",
                              "chains_run", "chains_retained", "states_used")
    info["primary_analysis"] = primary
    level_hpd, row_level_hpd = token_text(I, "LEVEL.hpd95")
    level_prior, row_level_prior = token_text(I, "LEVEL.prior95")
    open_from_text, row_open = token_text(I, "FIG1.open_rule")
    axis_token, row_axis_token = token_text(I, "FIG.axis_0201")
    T = C.TABLE_OF_TOKENS
    COL_T = "printed_where_the_quantity_is_unchanged_by_design"

    # ---- time axis (R4 a) and declaration --------------------------------------------------------
    first_txt, first_line, _ = C.pattern_in_document(
        I, C.CONVENTIONS, r"Time axes of time-series panels: (\d+ \w+ \d{4}) to", "first day of the time axis")
    XMIN = C.day_of_words(first_txt, "first day of the time axis")
    if C.day_of_words(axis_token, "FIG.axis_0201") != XMIN:
        raise C.Stop(f"first day of the axis: conventions {first_txt!r}, table of tokens {axis_token!r}")
    V.add("a; b; c", E_AXIS_FIRST, label=str(XMIN.year), value=first_txt, x=XMIN, unit="date",
          drawn_as="left end of the time axis; the year is printed below panel c", table=C.CONVENTIONS,
          row=first_line, columns=r"value=pattern Time axes of time-series panels: (\d+ \w+ \d{4}) to",
          note=f"the table of tokens gives the same day for FIG.axis_0201 ({axis_token}); line of the document, 0-based",
          quantity_id="FIG.axis_0201", delivered="yes")
    cases = I.table(C.CASES).need("week_ending_sunday", "total_drc", "ituri", "nord_kivu", "other_provinces",
                                  "weeks_pooled_for_province_split", "report_dates_in_week", "complete_week")
    weeks = [C.to_day(r["week_ending_sunday"], "week_ending_sunday") for r in cases.rows]
    if any((b - a).days != 7 for a, b in zip(weeks[:-1], weeks[1:])) or any(w.weekday() != 6 for w in weeks):
        raise C.Stop("case series: the weeks are not consecutive weeks that end on a Sunday")
    last = len(weeks) - 1
    XMAX = weeks[last] + dt.timedelta(days=1)
    V.add("a; b; c", E_AXIS_END, label="", value=cases.rows[last]["week_ending_sunday"], x=XMAX, unit="date",
          drawn_as="right end of the time axis", table=C.CASES, row=last, columns="value=week_ending_sunday",
          note="R4 (a): x is the day after the last day of the last week of the case series", delivered="yes")
    if XMIN.year != XMAX.year:
        raise C.Stop("the time axis spans two years; the label of the axis prints one")
    decl_txt, decl_line, _ = C.pattern_in_document(
        I, C.CONVENTIONS, r"Declaration of the outbreak: (\d+ \w+ \d{4})", "declaration of the outbreak")
    DECLARED = C.day_of_words(decl_txt, "declaration")
    V.add("a; b; c", E_DECLARED, label="outbreak declared", value=decl_txt, x=DECLARED, unit="date",
          drawn_as="vertical dashed line; the words are printed in panel a", table=C.CONVENTIONS, row=decl_line,
          columns=r"value=pattern Declaration of the outbreak: (\d+ \w+ \d{4})",
          note="line of the document, 0-based", delivered="yes")
    info.update(xmin=XMIN, xmax=XMAX, declared=DECLARED)

    # ---- a. effective population size --------------------------------------------------------------
    iv = I.table(C.T_INTERVALS_S).need("analysis", "interval", "older_edge", "newer_edge", "genomes_collected",
                                       "Ne_tau_years_median", "Ne_tau_years_hpd95_lower",
                                       "Ne_tau_years_hpd95_upper", "meets_convergence_criterion", "genomes")
    if {r["analysis"] for r in iv.rows} != {primary}:
        raise C.Stop(f"{C.T_INTERVALS_S}: holds rows of another analysis than {primary}")
    r_root = C.find_row(est, "date of the root", primary, "date of the most recent common ancestor", variant="date")
    root = est.rows[r_root]
    TM, TM_LO, TM_HI = (C.to_day(root[c], "date of the root") for c in ("median", "hpd95_lower", "hpd95_upper"))
    FADE = TM
    closed = [(i, r) for i, r in enumerate(iv.rows) if r["older_edge"] != "open"]
    closed.sort(key=lambda t: C.to_day(t[1]["older_edge"]))
    for (i, a), (j, b) in zip(closed[:-1], closed[1:]):
        if a["newer_edge"] != b["older_edge"]:
            raise C.Stop(f"{C.T_INTERVALS_S}: intervals are not contiguous (rows {i}, {j})")
    n_by_interval = sum(int(r["genomes_collected"]) for r in iv.rows)
    genomes = sorted({r["genomes"] for r in est.rows if r["analysis"] == primary})
    if len(genomes) != 1 or int(genomes[0]) != n_by_interval:
        raise C.Stop(f"genomes of the analysis {genomes} and sum over the intervals {n_by_interval} differ")
    info["genomes"] = genomes[0]
    info["genomes_summed_over_the_intervals"] = n_by_interval
    visible = [(i, r) for i, r in closed if C.to_day(r["newer_edge"]) > XMIN]
    for i, r in visible:
        older, newer = C.to_day(r["older_edge"]), C.to_day(r["newer_edge"])
        cut = older < XMIN
        V.add("a", E_STAIR, label=f"Skygrid interval {int(r['interval'])}", value=r["Ne_tau_years_median"],
              lower=r["Ne_tau_years_hpd95_lower"], upper=r["Ne_tau_years_hpd95_upper"],
              x=max(older, XMIN), x_end=newer, unit="years (Ne tau)",
              drawn_as="median line, 95 % HPD band" + ("; hatched before the median date of the root"
                                                        if older < FADE else ""),
              table=C.T_INTERVALS_S, row=i,
              columns="value=Ne_tau_years_median; lower=Ne_tau_years_hpd95_lower; upper=Ne_tau_years_hpd95_upper; "
                      + ("" if cut else "x=older_edge; ") + "x_end=newer_edge",
              note=("interval starts before the axis (older_edge " + r["older_edge"] + "); drawn from the first day of the axis"
                    if cut else ""),
              analysis=primary, delivered="yes", meets=r["meets_convergence_criterion"])
    for i, r in visible:
        older, newer = C.to_day(r["older_edge"]), C.to_day(r["newer_edge"])
        cut = older < XMIN
        V.add("a", E_GENOMES, label=f"Skygrid interval {int(r['interval'])}", value=r["genomes_collected"],
              x=max(older, XMIN), x_end=newer, unit="genomes", drawn_as="grey bar, right-hand axis",
              table=C.T_INTERVALS_S, row=i,
              columns="value=genomes_collected; " + ("" if cut else "x=older_edge; ") + "x_end=newer_edge",
              note="interval starts before the axis" if cut else "", analysis=primary, delivered="yes",
              meets="")
    root_text = f"tMRCA {C.fmt_date(TM)} ({C.fmt_date(TM_LO)} {C.NDASH} {C.fmt_date(TM_HI)})"
    qid_root = quantity_ids(I, primary, "date of the most recent common ancestor", None, root["end_used"],
                            ["median", "hpd95_lower", "hpd95_upper"], variant="date")
    V.add("a", E_ROOT, label=root_text, value=root["median"], lower=root["hpd95_lower"], upper=root["hpd95_upper"],
          unit="date", drawn_as="horizontal bar with a diamond at the median; dates printed above the bar"
          + ("; arrow head: the lower bound lies before the first day of the axis" if TM_LO < XMIN else ""),
          table=f_est, row=r_root, columns="value=median; lower=hpd95_lower; upper=hpd95_upper",
          note="label of the row as found: " + root["quantity"] + "; dates printed as the earlier generator prints them",
          analysis=primary, quantity_id=qid_root, delivered="yes", meets=root["meets_convergence_criterion"])
    V.add("a", E_FADE, label="limit of the hatched part", value=root["median"], x=FADE, unit="date",
          drawn_as="change from the hatched to the filled band", table=f_est, row=r_root, columns="value=median",
          note="the earlier generator took the first date of its table of trajectories at which half of the posterior "
               "samples have their root (column fraction_of_samples_with_root_before_date); no table of this round "
               "holds that column; the median date of the root is read from the table of estimates, as the task and "
               "the legend of the figure name the limit",
          analysis=primary, quantity_id=qid_root, delivered="yes", meets=root["meets_convergence_criterion"])
    lows = [C.num(r["Ne_tau_years_hpd95_lower"]) for _, r in visible]
    highs = [C.num(r["Ne_tau_years_hpd95_upper"]) for _, r in visible]
    ne_top = C.decade_ceil(max(highs))
    ne_floor = C.decade_floor(min(lows))
    V.add("a", E_AXIS_RULE, label="axis of Ne tau (logarithmic)", x=repr(NE_ROOM_BELOW * ne_floor), x_end=repr(ne_top),
          unit="years (Ne tau)", drawn_as="limits of the left-hand axis; ticks at the powers of ten inside",
          note=f"rule: upper limit = smallest power of ten at or above the largest upper bound drawn; lower limit = "
               f"{NE_ROOM_BELOW:g} x the largest power of ten at or below the smallest lower bound drawn (room for the bar of the root date)")
    gmax = max(int(r["genomes_collected"]) for _, r in visible)
    g_spine = max(25, int(math.ceil(gmax / 25.0 - 1e-9)) * 25)
    V.add("a", E_AXIS_RULE, label="axis of the genomes per interval", x="0", x_end=str(g_spine), unit="genomes",
          drawn_as="right-hand axis: the spine ends at x_end; ticks at the multiples of 50 on the spine",
          note=f"rule: the spine ends at the smallest multiple of 25 at or above the largest count drawn; the axis "
               f"runs to {GENOME_AXIS_OVER_SPINE:g} x the end of the spine, so that the bars keep to the lower part of the panel")
    V.add("a", E_KEY, label=f"median and {level_hpd} HPD", value=level_hpd, drawn_as="entry of the key of panel a",
          table=T, row=row_level_hpd, columns="value=" + COL_T, quantity_id="LEVEL.hpd95", delivered="yes",
          note="unchanged by design")
    V.add("a", E_KEY, label="before the median tMRCA", drawn_as="entry of the key of panel a")
    V.add("a", E_KEY, label="genomes per interval", drawn_as="words at the right-hand axis")

    # ---- b. reproduction number ----------------------------------------------------------------------
    periods = I.table(C.PERIODS).need("period", "start", "end", "stored_name", "status")
    prior = I.table(C.T_PRIOR).need("analysis", "quantity", "period_start", "period_end", "end_used",
                                    "prior_R_lower_from_the_standard_deviation",
                                    "prior_R_upper_from_the_standard_deviation",
                                    "prior_R_2_5_per_cent_quantile_simulated",
                                    "prior_R_97_5_per_cent_quantile_simulated", "range_that_governs")
    fmt_of = {}
    for r in I.table(T).need("format", "token").rows:
        fmt_of.setdefault(r["quantity_id"], (r["format"], r["token"]))
    known_ids = set(fmt_of)
    bounds_b = []
    for P in ("P1", "P2", "P3"):
        rp = periods.one(f"period {P}", period=P)
        name = periods.rows[rp]["stored_name"]
        end_used = "own" if periods.rows[rp]["end"] == "END" else "fixed date"
        i = C.find_row(est, f"R of {P}", primary, "reproduction number", period=P, end_used=end_used)
        r = est.rows[i]
        if r["period_start"] != periods.rows[rp]["start"] or (end_used == "fixed date"
                                                              and r["period_end"] != periods.rows[rp]["end"]):
            raise C.Stop(f"{f_est}: row {i}: the dates of {P} differ from {C.PERIODS}")
        qid = quantity_ids(I, primary, "reproduction number", P, end_used, ["median", "hpd95_lower", "hpd95_upper"],
                           row=r)
        m, lo, hi = (C.num(r[c], f"R of {P}") for c in ("median", "hpd95_lower", "hpd95_upper"))
        bounds_b += [lo, hi]
        V.add("b", E_PERIOD, label=name, value=r["median"], lower=r["hpd95_lower"], upper=r["hpd95_upper"],
              x=r["period_start"], x_end=r["period_end"], unit="R",
              drawn_as="blue rectangle (95 % HPD), thick line (median), bracket above the panel",
              table=f_est, row=i,
              columns="value=median; lower=hpd95_lower; upper=hpd95_upper; x=period_start; x_end=period_end",
              note=f"label of the row as found: {r['quantity']}; end_used = {r['end_used']}", analysis=primary,
              quantity_id=qid, delivered="yes", meets=r["meets_convergence_criterion"])
        text = C.fmt_ci(m, lo, hi)
        first_qid = qid.split("; ")[0] if qid else ""
        if first_qid in fmt_of:
            tok, err = C.print_token(I, fmt_of[first_qid][0], [r["median"], r["hpd95_lower"], r["hpd95_upper"]])
            same = "equal" if tok == text else "NOT equal"
            fnote = (f"printed as the earlier generator prints it (2 decimals); format_token with the format of the "
                     f"token {fmt_of[first_qid][1]} of the same quantity gives {tok or err!r}: {same}")
        else:
            fnote = "printed as the earlier generator prints it (2 decimals); the table of tokens holds no format of this quantity"
        V.add("b", E_PERIOD_TEXT, label=text, value=r["median"], lower=r["hpd95_lower"], upper=r["hpd95_upper"],
              x=r["period_start"], x_end=r["period_end"], unit="R", drawn_as="text in the strip above the panel",
              table=f_est, row=i,
              columns="value=median; lower=hpd95_lower; upper=hpd95_upper; x=period_start; x_end=period_end",
              note=fnote, analysis=primary, quantity_id=qid, delivered="yes",
              meets=r["meets_convergence_criterion"])
        V.add("b", E_PERIOD_NAME, label=name, value=name, x=r["period_start"], x_end=r["period_end"],
              drawn_as="text at the top of the strip above the panel", table=C.PERIODS, row=rp,
              columns="value=stored_name", note=f"period {P}; x and x_end as in the row of the estimate",
              analysis=primary, delivered="yes")
        j = C.find_row(prior, f"prior range of {P}", primary, "reproduction number", period=P, end_used=end_used)
        q = prior.rows[j]
        if (q["period_start"], q["period_end"]) != (r["period_start"], r["period_end"]):
            raise C.Stop(f"{C.T_PRIOR}: row {j}: the dates of {P} differ from the table of estimates")
        bounds_b += [C.num(q["prior_R_2_5_per_cent_quantile_simulated"]),
                     C.num(q["prior_R_97_5_per_cent_quantile_simulated"])]
        V.add("b", E_PRIOR_SIM, label=name, lower=q["prior_R_2_5_per_cent_quantile_simulated"],
              upper=q["prior_R_97_5_per_cent_quantile_simulated"], x=q["period_start"], x_end=q["period_end"],
              unit="R", drawn_as=f"light-grey rectangle with dashed outline ({level_prior} range)",
              table=C.T_PRIOR, row=j,
              columns="lower=prior_R_2_5_per_cent_quantile_simulated; upper=prior_R_97_5_per_cent_quantile_simulated; "
                      "x=period_start; x_end=period_end",
              note=C.NOTE_SIMULATED + "; carried over from "
                   "the earlier generator (R3, R4 d); the table has no column meets_convergence_criterion",
              analysis=primary, delivered="yes")
        V.add("b", E_PRIOR_SD, label=name, lower=q["prior_R_lower_from_the_standard_deviation"],
              upper=q["prior_R_upper_from_the_standard_deviation"], x=q["period_start"], x_end=q["period_end"],
              unit="R", drawn_as="not drawn", table=C.T_PRIOR, row=j,
              columns="lower=prior_R_lower_from_the_standard_deviation; upper=prior_R_upper_from_the_standard_deviation; "
                      "x=period_start; x_end=period_end",
              note=f"R4 (d): the range that Table 1 prints (range_that_governs = {q['range_that_governs']!r}); held "
                   "beside the range drawn; the table has no column meets_convergence_criterion",
              analysis=primary, quantity_id=(f"PRIM.{name}.prior" if f"PRIM.{name}.prior" in known_ids else ""),
              delivered="yes")
    pairs = I.table(C.T_PAIRS).need("analysis", "pair", "date_on_which_the_two_intervals_meet", "R_median",
                                    "R_hpd95_lower", "R_hpd95_upper", "sd_ratio", "sd_ratio_of_0_90_or_more",
                                    "meets_convergence_criterion")
    column_limit = re.fullmatch(r"sd_ratio_of_(\d)_(\d+)_or_more", "sd_ratio_of_0_90_or_more")
    limit_of_the_column = f"{column_limit.group(1)}.{column_limit.group(2)}"
    if limit_of_the_column != open_from_text:
        raise C.Stop(f"limit of the open symbols: column {limit_of_the_column}, table of tokens {open_from_text}")
    OPEN_FROM = float(open_from_text)
    flags = []
    drawn_pairs = 0
    for i, r in enumerate(pairs.rows):
        if r["analysis"] != primary:
            raise C.Stop(f"{C.T_PAIRS}: row {i} is of another analysis")
        flag = r["sd_ratio_of_0_90_or_more"]
        if flag not in ("yes", "no"):
            raise C.Stop(f"{C.T_PAIRS}: row {i}: sd_ratio_of_0_90_or_more is {flag!r}")
        flags.append((flag == "yes") == (C.num(r["sd_ratio"]) >= OPEN_FROM))
        d = C.to_day(r["date_on_which_the_two_intervals_meet"], "pair")
        if not (XMIN <= d <= XMAX):
            continue
        drawn_pairs += 1
        bounds_b += [C.num(r["R_hpd95_lower"]), C.num(r["R_hpd95_upper"])]
        V.add("b", E_PAIR, label=f"pair {r['pair']}, intervals meet on {r['date_on_which_the_two_intervals_meet']}",
              value=r["R_median"], lower=r["R_hpd95_lower"], upper=r["R_hpd95_upper"], x=d, unit="R",
              drawn_as=("open" if flag == "yes" else "filled") + " pale-blue point, 95 % HPD", table=C.T_PAIRS, row=i,
              columns="value=R_median; lower=R_hpd95_lower; upper=R_hpd95_upper; x=date_on_which_the_two_intervals_meet",
              note=f"sd_ratio_of_0_90_or_more = {flag} (decides the symbol); sd_ratio = {r['sd_ratio']}",
              analysis=primary, delivered="yes", meets=r["meets_convergence_criterion"])
    if not all(flags):
        raise C.Stop(f"{C.T_PAIRS}: the column sd_ratio_of_0_90_or_more and sd_ratio >= {open_from_text} differ in "
                     f"{flags.count(False)} rows")
    info["pairs_in_the_table"] = len(pairs)
    info["pairs_drawn"] = drawn_pairs
    info["pairs_flag_equals_sd_ratio_rule"] = f"{sum(flags)} of {len(flags)}"
    r_top = max(1.0, float(math.ceil(max(bounds_b) - 1e-9)))
    V.add("b", E_AXIS_RULE, label="axis of R", x="0", x_end=repr(r_top), unit="R",
          drawn_as="limits of the axis of R; ticks at the whole numbers",
          note=f"rule: from 0 to the smallest whole number at or above the largest bound drawn in the panel "
               f"(largest bound {max(bounds_b)!r}); the strip above the axis takes {STRIP_SHARE:g} x the length of the axis")
    V.add("b", E_KEY, label=f"window: median and {level_hpd} HPD", value=level_hpd,
          drawn_as="entry of the key of panel b", table=T, row=row_level_hpd, columns="value=" + COL_T,
          quantity_id="LEVEL.hpd95", delivered="yes", note="unchanged by design")
    V.add("b", E_KEY, label=f"window: {level_prior} range, smoothing prior alone", value=level_prior,
          drawn_as="entry of the key of panel b", table=T, row=row_level_prior, columns="value=" + COL_T,
          quantity_id="LEVEL.prior95", delivered="yes", note="unchanged by design")
    V.add("b", E_KEY, label=f"pair of adjacent intervals, {level_hpd} HPD", value=level_hpd,
          drawn_as="entry of the key of panel b", table=T, row=row_level_hpd, columns="value=" + COL_T,
          quantity_id="LEVEL.hpd95", delivered="yes", note="unchanged by design")
    V.add("b", E_KEY, label=f"the same, posterior SD / prior SD \u2265 {open_from_text}", value=open_from_text,
          drawn_as="entry of the key of panel b", table=T, row=row_open, columns="value=" + COL_T,
          quantity_id="FIG1.open_rule", delivered="yes",
          note="unchanged by design; the column sd_ratio_of_0_90_or_more of the table of the pairs carries the same limit in its name")

    # ---- c. confirmed cases -----------------------------------------------------------------------------
    for i, r in enumerate(cases.rows):
        s = sum(C.num(r[c]) for c, _, _ in PROVINCES)
        if abs(s - C.num(r["total_drc"])) > 0.11:
            raise C.Stop(f"{C.CASES}: row {i}: the provinces add up to {s}, the national count is {r['total_drc']}")
    counts = [int(r["report_dates_in_week"]) for r in cases.rows]
    fewest = [i for i, n in enumerate(counts) if n == min(counts)]
    provisional = [i for i in fewest if i == last]
    if len(provisional) != 1:
        raise C.Stop(f"provisional week: the rule of the earlier generator gives {len(provisional)} weeks "
                     f"(weeks with the fewest report dates: rows {fewest}; last week: row {last})")
    pv = provisional[0]
    m = re.search(r"The week ending (\d+ \w+) rests on (\d+) report dates[^.]*?it is to be marked as provisional in figures",
                  I.text(C.CASES_README))
    if not m:
        raise C.Stop(f"{C.CASES_README}: the sentence on the provisional week was not found")
    readme_day, readme_n = m.group(1), m.group(2)
    info["provisional_week_by_the_rule"] = cases.rows[pv]["week_ending_sunday"]
    info["provisional_week_report_dates"] = cases.rows[pv]["report_dates_in_week"]
    info["provisional_week_named_by_the_README"] = readme_day
    info["provisional_week_report_dates_by_the_README"] = readme_n
    info["rule_gives_the_week_of_the_README"] = (C.fmt_date(weeks[pv]) == readme_day
                                                 and cases.rows[pv]["report_dates_in_week"] == readme_n)
    info["weeks_with_the_fewest_report_dates"] = [cases.rows[i]["week_ending_sunday"] for i in fewest]
    info["complete_week_values"] = sorted({r["complete_week"] for r in cases.rows})
    for i, r in enumerate(cases.rows):
        pooled = int(r["weeks_pooled_for_province_split"])
        for col, lab, _ in PROVINCES:
            V.add("c", E_CASES, label=f"{lab}, week ending {r['week_ending_sunday']}", value=r[col],
                  x=weeks[i] - dt.timedelta(days=6), x_end=weeks[i], unit="cases",
                  drawn_as="stacked bar" + (f"; hatched (province split pooled over {pooled} weeks)" if pooled > 1 else "")
                  + ("; dashed outline, provisional" if i == pv else ""),
                  table=C.CASES, row=i, columns=f"value={col}; x_end=week_ending_sunday",
                  note=f"x = first day of the week (Monday); national count of the week {r['total_drc']}; "
                       f"report dates in the week {r['report_dates_in_week']}", delivered="yes")
        if pooled > 1:
            V.add("c", E_HATCH, label=f"week ending {r['week_ending_sunday']}", value=r["total_drc"],
                  lower=r["weeks_pooled_for_province_split"], x=weeks[i] - dt.timedelta(days=6), x_end=weeks[i],
                  unit="cases", drawn_as="white hatching over the whole bar (height: national count of the week)",
                  table=C.CASES, row=i,
                  columns="value=total_drc; lower=weeks_pooled_for_province_split; x_end=week_ending_sunday",
                  note="lower holds the number of weeks pooled (more than 1: hatched)", delivered="yes")
    V.add("c", E_PROVISIONAL, label="provisional", value=cases.rows[pv]["total_drc"],
          lower=cases.rows[pv]["report_dates_in_week"], x=weeks[pv] - dt.timedelta(days=6), x_end=weeks[pv],
          unit="cases", drawn_as="dashed black outline of the bar; the word 'provisional' with a leader line",
          table=C.CASES, row=pv, columns="value=total_drc; lower=report_dates_in_week; x_end=week_ending_sunday",
          note=f"R4 (b): rule of the earlier generator (the week with the fewest report dates, if it is the last week; "
               f"exactly one week): row {pv}; the README names the week ending {readme_day} with {readme_n} report dates; "
               f"lower holds the report dates of the week; the column complete_week reads "
               f"{', '.join(info['complete_week_values'])} in all {len(cases)} weeks and does not decide",
          delivered="yes")
    pw = sorted({int(r["weeks_pooled_for_province_split"]) for r in cases.rows
                 if int(r["weeks_pooled_for_province_split"]) > 1})
    for _, lab, _ in PROVINCES:
        V.add("c", E_KEY, label=lab, drawn_as="entry of the key of panel c")
    if pw:
        i_lo = min(i for i, r in enumerate(cases.rows) if int(r["weeks_pooled_for_province_split"]) == pw[0])
        i_hi = min(i for i, r in enumerate(cases.rows) if int(r["weeks_pooled_for_province_split"]) == pw[-1])
        V.add("c", E_KEY, label=f"province split pooled, {pw[0]}{C.NDASH}{pw[-1]} weeks",
              value=cases.rows[i_lo]["weeks_pooled_for_province_split"], drawn_as="entry of the key of panel c",
              table=C.CASES, row=i_lo, columns="value=weeks_pooled_for_province_split",
              note=f"smallest and largest number of weeks pooled (above 1) in the column; the largest stands in row {i_hi}",
              delivered="yes")
    cmax = max(C.num(r["total_drc"]) for r in cases.rows)
    c_top = (math.floor(cmax / 100.0 + 1e-9) + 1) * 100.0
    V.add("c", E_AXIS_RULE, label="axis of the cases per week", x="0", x_end=repr(c_top), unit="cases",
          drawn_as="limits of the axis; ticks at the multiples of 200 below the limit",
          note=f"rule: from 0 to the smallest multiple of 100 above the largest weekly national count; the word "
               f"'provisional' stands at {PROV_TEXT_SHARE:.4f} x the limit, its leader line ends at {PROV_LINE_SHARE:.4f} x the limit")
    # R4 (f): a value whose cell meets_convergence_criterion begins with 'no' (and is none of 'not assessed',
    # 'not evaluated', 'not applicable') is drawn like the others and marked with a double dagger
    for r in V.rows:
        if C.marked_by_R4f(r["meets_convergence_criterion"]) and r["drawn_as"] != "not drawn":
            r["drawn_as"] += "; double dagger"
            if r["element"] in (E_PERIOD_TEXT, E_ROOT):
                r["label"] += " " + C.DOUBLE_DAGGER
    return V, info


# =============================================================================================
# 1b. the optional table of interventions (task of 5 October); no value of it is typed here
# =============================================================================================
def strict_day(text, what):
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", text):
        raise C.Stop(f"{what}: the cell {text!r} is not a day written YYYY-MM-DD")
    try:
        return dt.date.fromisoformat(text)
    except ValueError:
        raise C.Stop(f"{what}: the cell {text!r} is not a day of the calendar")


def table_of_interventions(I, name):
    """Registers the table that the parameter names. Returns {table, test, tag, version, version_from}."""
    m = TABLE_NAME.match(name)
    if not m:
        raise C.Stop("the name of a table of interventions ends with _<YYYYMMDD>.csv (interventions_<tag>_<YYYYMMDD>"
                     f".csv, or another given name), with {MARKER}_ before it for a test with "
                     f"fixtures; the parameter names {name!r}")
    if not I.has(name):
        raise C.Stop(f"the JSON of the inputs holds no entry, or no file, for the table of interventions {name}")
    test = m.group(1) is not None
    if name in I.by_file:
        version, where = I.by_file[name]["version"], "the table of inputs"
    else:
        version = I.paths.get(name + "#version", "")
        where = f"the entry '{name}#version' of the JSON of the inputs" if version else \
            "none: file of a test with fixtures"
        if not version and not test:
            raise C.Stop(f"{name}: no version: neither the table of inputs nor an entry '{name}#version' of the "
                         "JSON of the inputs names one")
        I.by_file[name] = {"file": name, "version": version, "used_for": "Figure 1, panel b",
                           "what_it_is": "table of interventions", "analysis": "", "state_of_the_file": "",
                           "note": ""}
    return {"table": name, "test": test, "tag": m.group(2), "version": version, "version_from": where}


def lines_of_a_label(label):
    """The lines that the user fixed with <<BR>> (mark of the template), or None."""
    return [p.strip() for p in label.split("<<BR>>")] if "<<BR>>" in label else None


def form_of(t):
    """Which of the two forms the header of a table has, and which column holds what."""
    hits = [f for f, d in FORMS.items() if d["told_by"] in t.columns]
    if len(hits) != 1:
        raise C.Stop(f"{t.name}: the header holds neither the column 'label' (first form) nor the column "
                     "'short_label' (second form), or it holds both")
    d = FORMS[hits[0]]
    known = dict(list(d["required"].items()) + list(d["optional"].items()))
    missing = [c for c in d["required"].values() if c not in t.columns]
    unknown = [c for c in t.columns if c not in known.values()]
    if missing or unknown:
        raise C.Stop(f"{t.name} ({hits[0]}): columns missing: {missing}; columns not known: {unknown}; the columns "
                     f"are {list(d['required'].values())} and, optional, {list(d['optional'].values())}")
    return hits[0], {k: c for k, c in known.items() if c in t.columns}


def span_of(first, last, meaning, where):
    """'day', 'duration' or 'day not known' (the day is not known and lies between the two days)."""
    if last is None:
        if meaning:
            raise C.Stop(f"{where}: one day, and the cell meaning_of_the_period is not empty")
        return "day"
    if not meaning:
        raise C.Stop(f"{where}: a period (last_day after first_day) says in the column meaning_of_the_period what "
                     "it is: a duration, or two days between which the day lies; the cell is empty or the "
                     "column is missing")
    return "day not known" if DAY_NOT_KNOWN.search(meaning) else "duration"


def mark_of(span, kind):
    """The mark by which an intervention is drawn."""
    if span != "day":
        return span
    return "cross" if KIND_WITH_A_CROSS.search(kind) else "triangle"


def read_interventions(I, name, test, fresh=False):
    """The rows of the table, checked. Every fault stops the program; nothing is repaired.
    Returns (rows, form, columns); fresh=True reads the file again."""
    t = C.Table(name, I.path(name), I.version(name)) if fresh else I.table(name)
    form, cols = form_of(t)
    out, seen = [], set()
    for k, r in enumerate(t.rows):
        w = f"{name}, row {k}"
        g = {key: r[c] for key, c in cols.items()}
        first = strict_day(g["first_day"], w + ", first_day")
        last = strict_day(g["last_day"], w + ", last_day") if g["last_day"] and g["last_day"] != g["first_day"] \
            else None
        if last is not None and last < first:
            raise C.Stop(f"{w}: last_day {g['last_day']} lies before first_day {g['first_day']}")
        for key, longest in (("label", LABEL_MAX), ("words", WORDS_MAX), ("citation", WORDS_MAX),
                             ("keys", WORDS_MAX), ("kind", WORDS_MAX), ("meaning", WORDS_MAX)):
            if key not in g:
                continue
            s = g[key]
            if key in ("label", "words") and not s:
                raise C.Stop(f"{w}: the cell {cols[key]} is empty")
            plain = s.replace("<<BR>>", " ") if key == "label" else s
            if s != s.strip() or re.search(r"[\n\r\t\[\]]|<<|>>|  ", plain) or \
                    (key in ("label", "words") and plain.endswith(".")):
                raise C.Stop(f"{w}, {cols[key]}: {s!r} holds a line break, a tab, a square bracket, a mark of a "
                             "token, two spaces in a row, a space at an end or a full stop at its end")
            if len(plain) > longest:
                raise C.Stop(f"{w}, {cols[key]}: {len(plain)} characters, at most {longest}")
        fixed = lines_of_a_label(g["label"])
        if fixed is not None and (len(fixed) > LINES_MAX or any(not p for p in fixed)):
            raise C.Stop(f"{w}, label: <<BR>> gives {len(fixed)} lines (at most {LINES_MAX}), or an empty line")
        filled = [key for key in ("reference", "citation", "keys") if g.get(key)]
        if len(filled) != 1:
            raise C.Stop(f"{w}: exactly one source is filled ("
                         + ", ".join(cols[key] for key in ("reference", "citation", "keys") if key in cols) + ")")
        if filled[0] == "reference" and not re.fullmatch(r"[1-9]\d{0,3}", g["reference"]):
            raise C.Stop(f"{w}, {cols['reference']}: {g['reference']!r} is not the number of a reference")
        if filled[0] == "keys" and not all(p and p == p.strip() and ";" not in p for p in g["keys"].split("; ")):
            raise C.Stop(f"{w}, {cols['keys']}: {g['keys']!r}: keys are separated by '; '")
        if "order" in g and g["order"] != str(k + 1):
            raise C.Stop(f"{w}, {cols['order']}: {g['order']!r}; the column counts the rows from 1")
        if "origin" in g and not re.fullmatch(r"\d{1,6}", g["origin"]):
            raise C.Stop(f"{w}, {cols['origin']}: {g['origin']!r} is not the number of a row")
        span = span_of(first, last, g.get("meaning", ""), w)
        if test and not (g["label"].startswith(FIXTURE_WORD + " ") and FIXTURE_WORD in g["words"]
                         and (filled[0] == "reference" or FIXTURE_WORD in g[filled[0]])):
            raise C.Stop(f"{w}: the file is a test with fixtures; its label begins with '{FIXTURE_WORD} ' and its "
                         f"words, its citation and its keys hold '{FIXTURE_WORD}'")
        if not test and any(FIXTURE_WORD.lower() in c.lower() for c in r.values()):
            raise C.Stop(f"{w}: a cell holds '{FIXTURE_WORD}', and the name of the file does not begin with {MARKER}")
        key3 = (g["first_day"], g["last_day"] if last is not None else "", g["label"])
        if key3 in seen:
            raise C.Stop(f"{w}: the same days and label stand in an earlier row")
        seen.add(key3)
        out.append({"row": k, "first": first, "last": last, "span": span, "kind": g.get("kind", ""),
                    "mark": mark_of(span, g.get("kind", "")), "meaning": g.get("meaning", ""),
                    "label": g["label"], "words": g["words"], "source": g[filled[0]], "source_key": filled[0],
                    "cells": g})
    return out, form, cols


def on_the_axis(rows, xmin, xmax):
    """Marks which rows are drawn, where a period is cut, and numbers the rows that are drawn by their first
    day on the axis, then by the order of the rows (the thin line and the number stand at the first day)."""
    drawn = []
    for r in rows:
        if r["last"] is None:
            r["drawn"], r["cut_left"], r["cut_right"] = xmin <= r["first"] <= xmax, False, False
        else:
            r["drawn"] = r["last"] >= xmin and r["first"] <= xmax
            r["cut_left"], r["cut_right"] = r["drawn"] and r["first"] < xmin, r["drawn"] and r["last"] > xmax
        if r["drawn"]:
            a = max(r["first"], xmin)
            b = a if r["last"] is None else min(r["last"], xmax)
            r["x0"], r["x1"] = a, b
            r["anchor_days"] = float((a - xmin).days)
            drawn.append(r)
    for n, r in enumerate(sorted(drawn, key=lambda r: (r["anchor_days"], r["row"])), 1):
        r["number"] = n
    return drawn


def cut_words(r):
    return "; ".join(w for w, c in (("begins before the first day of the axis: arrow head there", r["cut_left"]),
                                    ("ends after the end of the axis: arrow head there", r["cut_right"])) if c)


def add_interventions(V, I, table, info, preview=False):
    """Adds the rows of the interventions to the file of values (after the rows of 4 October, which stay as
    they are). Returns the summary for the checks."""
    name, test = table["table"], table["test"]
    rows, form, cols = read_interventions(I, name, test)
    xmin, xmax = info["xmin"], info["xmax"]
    drawn = on_the_axis(rows, xmin, xmax)
    patterns = C.word_patterns(I)
    words = []
    said = "; test with fixtures: invented, never for the article" if test else ""
    for r in rows:
        g, k = r["cells"], r["row"]
        how = DRAWN[r["mark"]] + ("; " + cut_words(r) if r["drawn"] and cut_words(r) else "") if r["drawn"] \
            else NOT_DRAWN_OUTSIDE
        V.add("b", E_INT_MARK, label=g["label"], x=g["first_day"], x_end=g["last_day"], unit="date", drawn_as=how,
              table=name, row=k,
              columns=f"label={cols['label']}; x={cols['first_day']}; x_end={cols['last_day']}",
              note=("the time axis runs from %s to %s" % (xmin.isoformat(), xmax.isoformat())) + said,
              delivered="yes")
        if r["drawn"]:
            V.add("b", E_INT_NUMBER, label=str(r["number"]), x=g["first_day"], x_end=g["last_day"], unit="date",
                  drawn_as="number of the mark; printed only where the labels in words cannot be placed", row=k,
                  note="by rule, no cell: the interventions on the time axis are numbered by their first day "
                       "on the axis, then by the order of the rows; source_row_0based names the row of the table")
        V.add("b", E_INT_WORDS, label=g["words"], x=g["first_day"], unit="date",
              drawn_as="not drawn: printed in the sentence for the legend", table=name, row=k,
              columns=f"label={cols['words']}; x={cols['first_day']}", delivered="yes")
        V.add("b", E_INT_SOURCE, label=r["source"], x=g["first_day"], unit="date",
              drawn_as="not drawn: printed in the sentence for the legend", table=name, row=k,
              columns=f"label={cols[r['source_key']]}; x={cols['first_day']}",
              note={"reference": "number of a reference of the article", "citation": "citation",
                    "keys": "keys of sources, printed as the table has them"}[r["source_key"]], delivered="yes")
        for key, element, how2 in (("kind", E_INT_KIND, "not drawn as a text: decides the mark of one day (a kind "
                                                       "that begins with 'attack': cross; any other: triangle)"),
                                   ("meaning", E_INT_MEANING, "not drawn as a text: decides how a period is drawn "
                                                             "(bar: duration; dotted line: the day is not known)"),
                                   ("order", E_INT_ORDER, "not drawn"), ("origin", E_INT_ORIGIN, "not drawn")):
            if key in cols:
                V.add("b", element, label=g[key], x=g["first_day"], unit="date", drawn_as=how2, table=name, row=k,
                      columns=f"label={cols[key]}; x={cols['first_day']}", delivered="yes")
        for key in ("label", "words", "citation"):
            for h in C.words_found(g.get(key, ""), patterns):
                words.append((k, cols[key], h["found"], h["kind"]))
    mode = "test" if test else ("preview" if preview else "")
    if mode and rows:
        for text, where in zip(BANNERS[mode], ("top", "bottom")):
            V.add("picture", E_BANNER, label=text, drawn_as=f"banner at the {where} right of the picture",
                  note=f"the table of interventions is a file of a test ({MARKER})" if test
                  else "preview: the estimates are provisional, the events are those of the table")
    return {"table": name, "version": table["version"], "version_from": table["version_from"], "test": test,
            "preview": bool(preview), "tag": table["tag"], "rows": len(rows), "on_the_time_axis": len(drawn),
            "form_of_the_table": form, "columns_of_the_table": cols,
            "not_drawn_outside_the_time_axis": [
                {"row": r["row"], "label": r["label"], "first_day": r["cells"]["first_day"],
                 "last_day": r["cells"]["last_day"]} for r in rows if not r["drawn"]],
            "periods_cut_at_the_axis": [{"row": r["row"], "label": r["label"], "cut": cut_words(r)}
                                        for r in rows if r["drawn"] and cut_words(r)],
            "words_of_R4j_in_the_cells": words}


# =============================================================================================
# 1c. the strip of the interventions: where the labels stand (decided before the figure is made)
# =============================================================================================
def wraps_of(label):
    """Candidate lines of a label: the lines fixed with <<BR>>, or the balanced splits at spaces into 1, 2 and
    3 lines whose longest line holds at most LINE_MAX characters (fewest lines first)."""
    fixed = lines_of_a_label(label)
    if fixed is not None:
        return ["\n".join(fixed)]
    words = label.split(" ")
    out = []
    for k in range(1, LINES_MAX + 1):
        if k > len(words):
            break
        best = None
        for cuts in _cuts(len(words), k):
            lines = [" ".join(words[a:b]) for a, b in zip((0,) + cuts, cuts + (len(words),))]
            score = (max(len(s) for s in lines), sum(len(s) ** 2 for s in lines))
            if best is None or score < best[0]:
                best = (score, lines)
        if best[0][0] <= (LINE_MAX_ONE if k == 1 else LINE_MAX):
            out.append("\n".join(best[1]))
    if not out:                                 # a word longer than a line: the split with the most lines
        out.append("\n".join(best[1]))
    return out


def _cuts(n, k):
    import itertools
    return list(itertools.combinations(range(1, n), k - 1))


def _measure():
    """Width and height of a text of the strip in points, measured with the renderer of a figure of its own."""
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    f = Figure(figsize=(2, 2), dpi=200)
    FigureCanvasAgg(f)
    r = f.canvas.get_renderer()
    cache = {}

    def size(text):
        if text not in cache:
            t = f.text(0, 0, text, fontsize=C.SIZES[2], linespacing=LINESPACING)
            bb = t.get_window_extent(r)
            t.remove()
            cache[text] = (bb.width * 72.0 / f.dpi, bb.height * 72.0 / f.dpi)
        return cache[text]
    return size


def geometry_of_4_october():
    """Rectangles of the three panels of 4 October in inches from the bottom of the figure."""
    from matplotlib.figure import Figure
    from matplotlib.gridspec import GridSpec
    tmp = Figure(figsize=(C.FIG_WIDTH, HEIGHT))
    gs = GridSpec(3, 1, figure=tmp, **GRID)
    pos = [gs[i].get_position(tmp) for i in range(3)]
    return {"x0": pos[0].x0, "w": pos[0].width,
            "a": (pos[0].y0 * HEIGHT, pos[0].y1 * HEIGHT), "b": (pos[1].y0 * HEIGHT, pos[1].y1 * HEIGHT),
            "c": (pos[2].y0 * HEIGHT, pos[2].y1 * HEIGHT), "axes_width_pt": pos[0].width * C.FIG_WIDTH * 72.0}


def _snap_up(pt):
    """Points, rounded up to a whole number of pixels of the PNG, so that a panel that moves keeps its pixels."""
    px = math.ceil(pt * C.DPI / 72.0 - 1e-9)
    return px * 72.0 / C.DPI


def room_for_the_strip():
    steps = math.floor((C.MAX_HEIGHT - HEIGHT) * C.DPI / HATCH_STEP_PX + 1e-9)
    return {"growth_of_the_figure_pt": steps * HATCH_STEP_PX * 72.0 / C.DPI,
            "from_the_gap_between_a_and_b_pt": _snap_up(GAP_GIVES_PT),
            "from_the_gap_between_b_and_c_pt": _snap_up(GAP_GIVES_PT)}


def mark_of_a_row(drawn_as):
    for mark, text in DRAWN.items():
        if drawn_as == text or drawn_as.startswith(text + ";"):
            return mark
    raise C.Stop(f"file of values: 'drawn_as' of an intervention is not known: {drawn_as!r}")


def events_of_the_values(mark_rows, number_rows, XMIN, XMAX, per_day):
    """The interventions that are drawn, from the rows of the file of values."""
    number = {r["source_row_0based"]: int(r["label"]) for r in number_rows}
    ev = []
    for r in mark_rows:
        mark = mark_of_a_row(r["drawn_as"])
        first = D(r["x"])
        last = D(r["x_end"]) if r["x_end"] and r["x_end"] != r["x"] else None
        if (last is None) != (mark in ("triangle", "cross")):
            raise C.Stop(f"file of values: row {r['source_row_0based']}: 'drawn_as' and the days do not agree")
        if last is None and not (XMIN <= first <= XMAX) or last is not None and (last < XMIN or first > XMAX):
            raise C.Stop(f"file of values: the intervention of row {r['source_row_0based']} is drawn and lies "
                         "outside the time axis")
        a = max(first, XMIN)
        b = a if last is None else min(last, XMAX)
        anchor = a
        ev.append({"row": r["source_row_0based"], "number": number[r["source_row_0based"]], "label": r["label"],
                   "mark": mark, "first": r["x"], "last": r["x_end"] if last is not None else "",
                   "x0": a, "x1": b, "anchor": anchor, "foot": 0.0 if mark == "day not known" else MARK_PT,
                   "cut_left": first < XMIN, "cut_right": last is not None and last > XMAX,
                   "x_pt": (anchor - XMIN) * per_day})
    ev.sort(key=lambda e: e["number"])
    if [e["number"] for e in ev] != list(range(1, len(ev) + 1)) or \
            any(a["x_pt"] > b["x_pt"] + 1e-9 for a, b in zip(ev, ev[1:])):
        raise C.Stop("file of values: the numbers of the interventions do not run from left to right")
    return ev


def plan_of_the_strip(mark_rows, number_rows, xmin_t, xmax_t, declared_t):
    """Decides, from the rows of the file of values, which rule is used and where every text stands.
    Rule 1 (default): labels in words beside a thin line, at as few heights as possible (at most LEVELS_MAX).
    Rule 2 (fallback): numbers above the marks, where rule 1 finds no place without overlap in the room that
    the height of the figure allows. All positions in points; x from the first day of the time axis."""
    geo = geometry_of_4_october()
    size = _measure()
    XMIN, XMAX = D(xmin_t), D(xmax_t)
    per_day = geo["axes_width_pt"] / (XMAX - XMIN)
    ev = events_of_the_values(mark_rows, number_rows, XMIN, XMAX, per_day)
    room = room_for_the_strip()
    room_total = sum(room.values())
    # room_total: what the strip as a whole may take above panel b (the gap under it, the labels, the edge
    # above them). h_max: what of it the LABELS may take, measured from the lower edge of the strip (one pixel
    # is kept for the rounding to whole pixels)
    h_max = room_total - _snap_up(STRIP_GAP_PT) - STRIP_TOP_PT - 72.0 / C.DPI
    x_decl = (D(declared_t) - XMIN) * per_day
    width = geo["axes_width_pt"]

    two, three = size("Xg\nXg")[1] + LABEL_GAP_PT, size("Xg\nXg\nXg")[1] + LABEL_GAP_PT
    eps = 1e-6
    # the heights that are tried, fewest first: the lower edges of the labels above the lower edge of the strip;
    # 'two' leaves room for a label of two lines below, 'three' for a label of three lines
    grids = [g for g in ([LEVEL0_PT], [LEVEL0_PT, LEVEL0_PT + two], [LEVEL0_PT, LEVEL0_PT + three],
                         [LEVEL0_PT, LEVEL0_PT + two, LEVEL0_PT + 2 * two]) if len(g) <= LEVELS_MAX]

    def candidates(e, grid):
        out = []
        for text in wraps_of(e["label"]):
            w, h = size(text)
            for level, y0 in enumerate(grid):
                if y0 + h > h_max + 1e-9:
                    continue
                for side in ("right", "left"):
                    x0 = e["x_pt"] + LABEL_PAD_PT if side == "right" else e["x_pt"] - LABEL_PAD_PT - w
                    if x0 < -1e-9 or x0 + w > width + 1e-9:
                        continue
                    if x0 - LINE_CLEAR_PT < x_decl < x0 + w + LINE_CLEAR_PT:
                        continue                    # a label does not stand on the line of the declaration
                    out.append({"text": text, "lines": text.count("\n") + 1, "level": level, "side": side,
                                "x0": x0, "x1": x0 + w, "y0": y0, "y1": y0 + h, "w": w, "h": h})
        out.sort(key=lambda c: (c["level"], c["lines"], c["side"] != "right"))
        return out

    def conflict(ei, ci, ej, cj):
        if ci["x0"] < cj["x1"] + LABEL_GAP_PT - eps and cj["x0"] < ci["x1"] + LABEL_GAP_PT - eps and \
                ci["y0"] < cj["y1"] + LABEL_GAP_PT - eps and cj["y0"] < ci["y1"] + LABEL_GAP_PT - eps:
            return True
        for (ea, ca, eb, cb) in ((ei, ci, ej, cj), (ej, cj, ei, ci)):       # the line of b through the label of a
            if ca["x0"] - LINE_CLEAR_PT < eb["x_pt"] < ca["x1"] + LINE_CLEAR_PT and \
                    cb["y1"] > ca["y0"] - LABEL_GAP_PT + eps:
                return True
        return False

    def search_in(cands, chosen, i=0):
        nonlocal nodes
        if i == len(ev):
            return True
        for c in cands[i]:
            nodes += 1
            if nodes > NODES_MAX:
                return False
            if all(not conflict(ev[i], c, ev[j], chosen[j]) for j in range(i)):
                chosen.append(c)
                if search_in(cands, chosen, i + 1):
                    return True
                chosen.pop()
        return False

    found, nodes, tried = None, 0, []
    for grid in grids:
        cands = [candidates(e, grid) for e in ev]
        tried.append({"heights": len(grid), "lower_edges_of_the_labels_pt": [round(y, 2) for y in grid],
                      "interventions_without_a_place": sum(1 for c in cands if not c)})
        if any(not c for c in cands):
            continue
        chosen = []
        ok = search_in(cands, chosen)
        tried[-1]["placed"] = bool(ok)
        if ok:
            found = chosen
            break
        if nodes > NODES_MAX:
            tried[-1]["search_stopped_after_nodes"] = nodes
            break
    plan = {"events": ev, "per_day_pt": per_day, "room": room, "room_total_pt": room_total,
            "room_of_the_labels_pt": h_max, "tried": tried, "nodes_searched": nodes, "geometry": geo}
    if found is None:
        # what the labels in words would need, were the height of the figure no limit (for the report only)
        need, spent = None, nodes
        h_limit, h_max = h_max, float("inf")
        for n in range(LEVELS_MAX + 1, HEIGHTS_ASKED_ABOUT + 1):
            grid = [LEVEL0_PT + k * two for k in range(n)]
            cands = [candidates(e, grid) for e in ev]
            if any(not c for c in cands):
                continue
            chosen, nodes = [], 0
            if search_in(cands, chosen):
                top = max(c["y1"] for c in chosen)
                stack = _snap_up(STRIP_GAP_PT) + _snap_up(top + STRIP_TOP_PT)
                need = {"heights": n,
                        "for_the_strip_above_panel_b_pt": round(stack, 2),
                        "for_the_strip_above_panel_b_the_figure_allows_pt": round(room_total, 2),
                        "of_it_for_the_labels_pt": round(top, 2),
                        "of_it_for_the_labels_the_figure_allows_pt": round(h_limit, 2),
                        "height_of_the_figure_with_panels_of_this_size_in":
                            round(HEIGHT + (stack - room["from_the_gap_between_a_and_b_pt"]
                                            - room["from_the_gap_between_b_and_c_pt"]) / 72.0, 2),
                        "largest_height_of_the_conventions_in": C.MAX_HEIGHT}
                break
        h_max, nodes = h_limit, spent
        plan["what_the_labels_in_words_would_need"] = need if need is not None else {
            "heights": f"more than {HEIGHTS_ASKED_ABOUT}, or the search was stopped"}
    if found is not None:
        for e, c in zip(ev, found):
            e["place"] = c
        plan["rule"] = "labels in words"
        plan["heights_used"] = len({c["y0"] for c in found})
        plan["why"] = f"every label has a place without overlap at {plan['heights_used']} height(s)"
        top = max(c["y1"] for c in found)
    else:
        plan["rule"] = "numbered marks"
        plan["heights_used"] = 0
        plan["why"] = (f"the labels in words find no place without overlap at up to {LEVELS_MAX} heights: above "
                       f"panel b the figure allows {room_total:.2f} pt for the strip, of which {h_max:.2f} pt for "
                       "the labels (the rest is the gap under the strip and the edge above the labels); numbers "
                       "are printed and the words stand in the sentence for the legend")
        groups = [[e] for e in ev]
        while True:
            boxes = []
            for g in groups:
                nums = [e["number"] for e in g]
                text = str(nums[0]) if len(nums) == 1 else (f"{nums[0]}, {nums[1]}" if len(nums) == 2
                                                            else f"{nums[0]}{C.NDASH}{nums[-1]}")
                w, h = size(text)
                c = 0.5 * (g[0]["x_pt"] + g[-1]["x_pt"])
                if c - w / 2 - LINE_CLEAR_PT < x_decl < c + w / 2 + LINE_CLEAR_PT:
                    # a number does not stand on the line of the declaration: it moves to the side of its mark
                    c = x_decl - LINE_CLEAR_PT - w / 2 if c <= x_decl else x_decl + LINE_CLEAR_PT + w / 2
                c = min(max(c, w / 2), width - w / 2)
                boxes.append({"text": text, "x0": c - w / 2, "x1": c + w / 2, "centre": c, "w": w, "h": h})
            hit = next((i for i in range(len(groups) - 1)
                        if boxes[i]["x1"] + LABEL_GAP_PT > boxes[i + 1]["x0"]), None)
            if hit is None:
                break
            groups[hit:hit + 2] = [groups[hit] + groups[hit + 1]]
        plan["numbers"] = [{"text": b["text"], "centre_pt": b["centre"], "x0": b["x0"], "x1": b["x1"], "h": b["h"],
                            "numbers": [e["number"] for e in g]} for g, b in zip(groups, boxes)]
        top = NUMBER_Y_PT + max(b["h"] for b in boxes)
    plan["strip_height_pt"] = _snap_up(top + STRIP_TOP_PT)
    plan["stack_pt"] = _snap_up(STRIP_GAP_PT) + plan["strip_height_pt"]
    if plan["stack_pt"] > room_total + 1e-6:
        raise C.Stop(f"the strip needs {plan['stack_pt']:.1f} pt, the figure allows {room_total:.1f} pt")
    s = plan["stack_pt"]
    unit = HATCH_STEP_PX * 72.0 / C.DPI
    grow = min(math.ceil(s / unit - 1e-9) * unit, room["growth_of_the_figure_pt"])
    red_ab = min(s - grow, room["from_the_gap_between_a_and_b_pt"])       # below 0: the gap grows
    red_bc = max(0.0, s - grow - red_ab)
    plan["layout"] = {"growth_of_the_figure_pt": grow, "taken_from_the_gap_between_a_and_b_pt": red_ab,
                      "taken_from_the_gap_between_b_and_c_pt": red_bc,
                      "height_of_the_figure_in": HEIGHT + grow / 72.0}
    return plan


def figure_with_strip(plan):
    """The figure with the strip above panel b. The three panels keep the size of 4 October; panel c stays where
    it is, panel b moves down only if the gap to c gives room, panel a moves up with the figure."""
    geo, lay = plan["geometry"], plan["layout"]
    H1 = lay["height_of_the_figure_in"]
    red_bc = lay["taken_from_the_gap_between_b_and_c_pt"] / 72.0
    fig = C.new_figure(H1)

    def rect(y0, y1):
        return [geo["x0"], y0 / H1, geo["w"], (y1 - y0) / H1]
    a0, a1 = geo["a"]
    b0, b1 = geo["b"][0] - red_bc, geo["b"][1] - red_bc
    c0, c1 = geo["c"]
    grow = lay["growth_of_the_figure_pt"] / 72.0
    axa = fig.add_axes(rect(a0 + grow, a1 + grow))
    axb = fig.add_axes(rect(b0, b1), sharex=axa)
    axc = fig.add_axes(rect(c0, c1), sharex=axa)
    s0 = b1 + _snap_up(STRIP_GAP_PT) / 72.0
    axs = fig.add_axes(rect(s0, s0 + plan["strip_height_pt"] / 72.0), sharex=axa)
    return fig, axa, axb, axc, axs


def draw_strip(fig, axs, plan, XMIN, XMAX, DECLARED):
    """Marks on the lower edge of the strip; rule 1: thin line and label in words; rule 2: numbers."""
    H = plan["strip_height_pt"]
    axs.set_ylim(0, H)                               # one unit of the strip is one point
    axs.patch.set_visible(False)
    for s in ("left", "right", "top"):
        axs.spines[s].set_visible(False)
    axs.spines["bottom"].set_color(STRIP_COLOUR)
    axs.spines["bottom"].set_linewidth(STEM_LW_PT)
    axs.tick_params(axis="both", which="both", bottom=False, left=False, labelbottom=False, labelleft=False)
    axs.spines["bottom"].set_zorder(2.0)
    decl, = axs.plot([DECLARED, DECLARED], [0, H], color=C.GREY_REF, ls=(0, (4, 3)), lw=0.9, zorder=1.5)
    per_day = plan["per_day_pt"]
    words = plan["rule"] == "labels in words"
    marks, texts, stems = [], [], []
    kw = dict(color=STRIP_COLOUR, clip_on=False, zorder=4)
    for e in plan["events"]:
        arts = []
        if e["mark"] == "triangle":
            m, = axs.plot([e["x0"]], [MARK_PT / 2], marker="v", ms=MARK_PT, mec="none", ls="none", **kw)
            arts.append(("triangle", m))
        elif e["mark"] == "cross":
            m, = axs.plot([e["x0"]], [MARK_PT / 2], marker="x", ms=MARK_PT, mew=CROSS_LW_PT, ls="none", **kw)
            arts.append(("cross", m))
        else:
            if e["mark"] == "duration":
                ln, = axs.plot([e["x0"], e["x1"]], [SPAN_Y_PT, SPAN_Y_PT], lw=SPAN_LW_PT, solid_capstyle="butt", **kw)
                arts.append(("bar", ln))
            else:
                ln, = axs.plot([e["x0"], e["x1"]], [SPAN_Y_PT, SPAN_Y_PT], lw=SPAN_LW_PT, ls=(0, DOTS),
                               dash_capstyle="butt", **kw)
                arts.append(("dotted line", ln))
            for x, cut, head in ((e["x0"], e["cut_left"], "<"), (e["x1"], e["cut_right"], ">")):
                if cut:
                    m, = axs.plot([x], [SPAN_Y_PT], marker=head, ms=ARROW_PT, mec="none", ls="none", **kw)
                    arts.append(("arrow head", m))
                elif e["mark"] == "duration":
                    m, = axs.plot([x, x], [SPAN_Y_PT - SPAN_LW_PT / 2, MARK_PT], lw=CAP_LW_PT,
                                  solid_capstyle="butt", **kw)
                    arts.append(("tick", m))
        marks.append({"row": e["row"], "number": e["number"], "mark": e["mark"], "first": e["first"],
                      "last": e["last"], "x0": e["x0"], "x1": e["x1"], "anchor": e["anchor"],
                      "cut_left": e["cut_left"], "cut_right": e["cut_right"], "artists": arts})
        if words:
            c = e["place"]
            st, = axs.plot([e["anchor"], e["anchor"]], [e["foot"], c["y1"]], color=STRIP_COLOUR, lw=STEM_LW_PT,
                           zorder=3, solid_capstyle="butt", clip_on=False)
            stems.append((e["row"], st))
            x_text = e["anchor"] + (LABEL_PAD_PT if c["side"] == "right" else -LABEL_PAD_PT) / per_day
            t = axs.text(x_text, c["y1"], c["text"], ha="left" if c["side"] == "right" else "right", va="top",
                         fontsize=C.SIZES[2], linespacing=LINESPACING, color=STRIP_COLOUR, zorder=5)
            texts.append({"row": e["row"], "number": e["number"], "text": c["text"], "side": c["side"],
                          "height": c["level"] + 1, "lines": c["lines"], "artist": t})
    if not words:
        for g in plan["numbers"]:
            x = XMIN + g["centre_pt"] / per_day
            t = axs.text(x, NUMBER_Y_PT, g["text"], ha="center", va="bottom", fontsize=C.SIZES[2],
                         linespacing=LINESPACING, color=STRIP_COLOUR, zorder=5)
            texts.append({"numbers": g["numbers"], "text": g["text"], "artist": t})
    return {"marks": marks, "texts": texts, "stems": stems, "declaration": decl, "words": words}


# =============================================================================================
# 1d. the sentence for the legend and the checks of the strip
# =============================================================================================
def interventions_of_the_table(I, table, xmin, xmax):
    """The table read AGAIN from its file (not from the file of values): what must stand in the figure."""
    rows, form, cols = read_interventions(I, table["table"], table["test"], fresh=True)
    on_the_axis(rows, xmin, xmax)
    out = []
    for r in rows:
        e = {"row": str(r["row"]), "mark": r["mark"], "span": r["span"], "kind": r["kind"],
             "meaning": r["meaning"], "first": r["cells"]["first_day"],
             "last": r["cells"]["last_day"] if r["last"] is not None else "",
             "last_cell": r["cells"]["last_day"], "on_the_axis": r["drawn"], "label": r["label"],
             "words": r["words"], "source": r["source"], "cut_left": r["cut_left"], "cut_right": r["cut_right"]}
        if r["drawn"]:
            e.update(x0=r["x0"].isoformat(), x1=r["x1"].isoformat(), number=r["number"],
                     anchor=D(xmin.isoformat()) + r["anchor_days"])
        out.append(e)
    return out


def sentence_for_the_legend(I, items, rule):
    """One sentence in the form of the template. items: the interventions that are drawn, in the order of their
    numbers: [{number, mark, kind, label, words, first, last, source, cut_left, cut_right}]."""
    tok = I.table(C.TABLE_OF_TOKENS).need("unit_of_text", "kind_of_format", "format", "token")

    def form(kind):
        hits = [r for r in tok.rows if r["unit_of_text"] == UNIT_OF_THE_LEGEND and r["kind_of_format"] == kind]
        if not hits:
            raise C.Stop(f"table of tokens: the unit {UNIT_OF_THE_LEGEND} holds no token of the kind {kind!r}")
        return hits[0]
    f_day, f_range = form("date"), form("range of dates")

    def printed(fmt, values, e):
        text, err = C.print_token(I, fmt["format"], values)
        if err:
            raise C.Stop(f"date of the intervention of row {e['row']}: {err}")
        return text
    parts = []
    for e in items:
        if e["mark"] in ("triangle", "cross"):
            when = printed(f_day, [e["first"]], e)
        elif e["mark"] == "duration":
            when = printed(f_range, [[e["first"], e["last"]]], e)
        else:
            when = f"between {printed(f_day, [e['first']], e)} and {printed(f_day, [e['last']], e)}"
        if rule == "numbered marks":
            # number, words for the legend, day or days, sources; the short label, which the figure does not
            # show, is left out
            parts.append(f"{e['number']}, {e['words']}, {when} [{e['source']}]")
        else:
            parts.append(f"{e['label'].replace('<<BR>>', ' ')}: {e['words']}, {when} [{e['source']}]")
    key = []
    for mark, shape in (("triangle", "triangle"), ("cross", "cross")):
        kinds = sorted({e["kind"] for e in items if e["mark"] == mark and e["kind"]})
        if any(e["mark"] == mark for e in items):
            key.append(f"{shape}: " + (", ".join(kinds) + ", " if kinds else "") + "one day")
    if any(e["mark"] == "duration" for e in items):
        key.append("bar under the line, a tick at each end: from the first to the last day")
    if any(e["mark"] == "day not known" for e in items):
        key.append("dotted line under the line: the day is not known and lies between the two days")
    if any(e["cut_left"] or e["cut_right"] for e in items):
        key.append("arrow head: the period continues beyond the time axis")
    head = ("Numbered marks above the panel" if rule == "numbered marks" else "Strip above the panel") \
        + f": events at their dates ({'; '.join(key)}): "
    return head + "; ".join(parts) + ".", {"date": f_day["token"], "range of dates": f_range["token"]}


def explanation_against_the_table(sentence, expected):
    """The explanation of the marks (the words in brackets at the head of the sentence) against the table read
    again: it names every kind and every form that the table holds on the time axis, and no other."""
    m = re.match(r"^[^:]*: events at their dates \((.*?)\): ", sentence)
    if not m:
        return {"ok": False, "note": "the head of the sentence was not found"}
    said = m.group(1)
    forms = {"triangle": "triangle: ", "cross": "cross: ",
             "duration": "bar under the line, a tick at each end: ",
             "day not known": "dotted line under the line: "}
    held = sorted({e["mark"] for e in expected})
    named = sorted(k for k, w in forms.items() if w in said)
    kinds_held = sorted({e["kind"] for e in expected if e["mark"] in ("triangle", "cross") and e["kind"]})
    kinds_named = []
    for mark in ("triangle", "cross"):
        k = re.search(r"(?:^|; )" + mark + r": (.*?)one day", said)
        if k:
            kinds_named += [x.strip() for x in k.group(1).split(",") if x.strip()]
    arrow_held = any(e["cut_left"] or e["cut_right"] for e in expected)
    arrow_named = "arrow head: " in said
    return {"explanation": said, "forms_that_the_table_holds": held, "forms_named": named,
            "kinds_of_the_marks_of_one_day_that_the_table_holds": kinds_held,
            "kinds_named": sorted(set(kinds_named)),
            "a_period_is_cut_at_the_axis": arrow_held, "the_arrow_head_is_named": arrow_named,
            "ok": held == named and kinds_held == sorted(set(kinds_named)) and arrow_held == arrow_named}


def _runs(mask_row):
    """[(first column, last column)] of the runs of True in a row of a mask."""
    idx = np.flatnonzero(mask_row)
    if idx.size == 0:
        return []
    cut = np.flatnonzero(np.diff(idx) > 1)
    starts = np.concatenate(([idx[0]], idx[cut + 1]))
    ends = np.concatenate((idx[cut], [idx[-1]]))
    return [(int(a), int(b)) for a, b in zip(starts, ends)]


def asked_of_the_pixels(expected, col, s, words):
    """What the table asks of three rows of pixels of the strip: [(row of pixels, first column, last column,
    tolerance at the first, tolerance at the last, what)]. Lengths in pixels; the rows are TIP_ROW_PT and half a
    mark above the lower edge of the strip, and the middle of the bars under it."""
    out = []
    t, tc = PIXEL_TOLERANCE, PIXEL_TOLERANCE_CROSS
    for e in expected:
        r = f"row {e['row']}"
        if e["mark"] == "triangle":
            for row, wide in (("tip", TIP_ROW_PT), ("mid", MARK_PT / 2)):
                out.append((row, col(e["x0"]) - wide * s / 2, col(e["x0"]) + wide * s / 2, t, t,
                            f"{r}: triangle at {e['x0']}"))
        elif e["mark"] == "cross":
            half = CROSS_LW_PT * s / math.sqrt(2.0)
            off = (MARK_PT / 2 - TIP_ROW_PT) * s
            out.append(("tip", col(e["x0"]) - off - half, col(e["x0"]) - off + half, tc, tc,
                        f"{r}: cross at {e['x0']}, left arm"))
            out.append(("tip", col(e["x0"]) + off - half, col(e["x0"]) + off + half, tc, tc,
                        f"{r}: cross at {e['x0']}, right arm"))
            out.append(("mid", col(e["x0"]) - half, col(e["x0"]) + half, tc, tc, f"{r}: cross at {e['x0']}"))
        else:
            lo, hi, t_lo, t_hi = col(e["x0"]), col(e["x1"]), t, t
            if e["cut_left"]:
                lo, t_lo = lo - ARROW_PT * s / 2, PIXEL_TOLERANCE_HEAD
            elif e["mark"] == "duration":
                lo -= CAP_LW_PT * s / 2
            if e["cut_right"]:
                hi, t_hi = hi + ARROW_PT * s / 2, PIXEL_TOLERANCE_HEAD
            elif e["mark"] == "duration":
                hi += CAP_LW_PT * s / 2
            if e["mark"] == "duration":
                out.append(("span", lo, hi, t_lo, t_hi, f"{r}: bar from {e['x0']} to {e['x1']}"))
                for x, cut in ((e["x0"], e["cut_left"]), (e["x1"], e["cut_right"])):
                    if not cut:
                        for row in ("tip", "mid"):
                            out.append((row, col(x) - CAP_LW_PT * s / 2, col(x) + CAP_LW_PT * s / 2, t, t,
                                        f"{r}: tick at {x}"))
            else:
                on, period = DOTS[0] * SPAN_LW_PT * s, (DOTS[0] + DOTS[1]) * SPAN_LW_PT * s
                a, b, k = col(e["x0"]), col(e["x1"]), 0
                if e["cut_left"]:
                    out.append(("span", lo, a + ARROW_PT * s / 2, t_lo, t, f"{r}: arrow head at {e['x0']}"))
                while a + k * period < b - 1e-9:
                    out.append(("span", a + k * period, min(a + k * period + on, b), t, t,
                                f"{r}: dot {k + 1} of the dotted line from {e['x0']} to {e['x1']}"))
                    k += 1
                if e["cut_right"]:
                    out.append(("span", b - ARROW_PT * s / 2, hi, t, t_hi, f"{r}: arrow head at {e['x1']}"))
                for row in ("tip", "mid") if words else ():
                    out.append((row, a - STEM_LW_PT * s / 2, a + STEM_LW_PT * s / 2, t, t,
                                f"{r}: thin line at {e['x0']}"))
    return out


def marks_read_from_the_png(png_path, strip_box, xmin, xmax, expected, words=True):
    """Reads the marks back from the PIXELS of the saved picture. strip_box: (x0, y0, w, h) of the strip in
    fractions of the figure. In three rows of pixels the runs of pixels in the colour of the strip are set
    against the stretches that the table asks for: each stretch (marks that touch: their stretches as one)
    must be one run whose two ends lie within the tolerance of the ends that are asked for, and no other run
    may be there. The tolerance is wider for a cross and at the point of an arrow head; every tolerance and
    the largest difference found under it are written to the checks."""
    from PIL import Image
    with Image.open(png_path) as im:
        a = np.asarray(im.convert("RGB")).astype(int)
    Hpx, Wpx = a.shape[:2]
    s = C.DPI / 72.0
    mask = ((a[:, :, 0] - a[:, :, 1]) > 25) & ((a[:, :, 2] - a[:, :, 1]) > 25)
    x0_px, w_px = strip_box[0] * Wpx, strip_box[2] * Wpx
    base = strip_box[1] * Hpx                      # pixels from the bottom

    def col(day_text):
        return x0_px + (D(day_text) - xmin) / (xmax - xmin) * w_px

    def row(pt):
        return int(math.floor(Hpx - (base + pt * s)))
    asked = asked_of_the_pixels(expected, col, s, words)
    names = {PIXEL_TOLERANCE: "usual", PIXEL_TOLERANCE_CROSS: "of a cross",
             PIXEL_TOLERANCE_HEAD: "at the point of an arrow head"}
    out = {"rows_of_pixels": {}, "not_found": [], "not_asked_for": [], "outside_the_tolerance": [], "read": [],
           "tolerances_px": {"usual": PIXEL_TOLERANCE, "of a cross": PIXEL_TOLERANCE_CROSS,
                             "at the point of an arrow head": PIXEL_TOLERANCE_HEAD},
           "ends_compared_under_each_tolerance": {v: 0 for v in names.values()},
           "largest_difference_px_under_each_tolerance": {v: 0.0 for v in names.values()},
           "one_day_px": round(w_px / (xmax - xmin), 3)}
    for name, pt in (("tip", TIP_ROW_PT), ("mid", MARK_PT / 2.0), ("span", SPAN_Y_PT)):
        r = row(pt)
        runs = _runs(mask[r, :])
        mine = sorted((x for x in asked if x[0] == name), key=lambda x: x[1])
        groups = []
        for x in mine:                               # marks that touch are read as one stretch
            if groups and x[1] <= groups[-1]["hi"] + 1.0:
                g = groups[-1]
                if x[2] > g["hi"]:
                    g["hi"], g["t_hi"] = x[2], x[4]
                g["what"].append(x[5])
            else:
                groups.append({"lo": x[1], "hi": x[2], "t_lo": x[3], "t_hi": x[4], "what": [x[5]]})
        out["rows_of_pixels"][name] = {"row_from_the_top": r, "points_from_the_lower_edge_of_the_strip": pt,
                                       "runs_in_the_picture": len(runs), "stretches_asked_for": len(mine),
                                       "stretches_after_joining_marks_that_touch": len(groups)}
        used = set()
        for g in groups:
            what = " + ".join(g["what"])
            hit = [k for k, (p, q) in enumerate(runs)
                   if p <= g["hi"] + g["t_hi"] and q + 1 >= g["lo"] - g["t_lo"]]
            if not hit:
                if g["hi"] - g["lo"] < 1.5 and len(g["what"]) == 1:
                    continue                         # a stretch of less than 1.5 pixels need not be seen
                out["not_found"].append(f"{name}: {what}: asked for at columns {g['lo']:.1f} to {g['hi']:.1f}")
                continue
            used.update(hit)
            p, q = runs[hit[0]][0], runs[hit[-1]][1]
            holes = sum(runs[b][0] - runs[a][1] - 1 for a, b in zip(hit, hit[1:]))
            d_lo, d_hi = p - g["lo"], (q + 1) - g["hi"]
            centre = 0.5 * (p + q + 1)
            entry = {"row_of_pixels": name, "what": what, "asked_for_columns": [round(g["lo"], 2),
                                                                                round(g["hi"], 2)],
                     "run_of_the_picture": [int(p), int(q) + 1], "difference_px": [round(d_lo, 2), round(d_hi, 2)],
                     "difference_days": [round(d_lo / w_px * (xmax - xmin), 3),
                                         round(d_hi / w_px * (xmax - xmin), 3)],
                     "tolerance_px": [g["t_lo"], g["t_hi"]], "marks_joined": len(g["what"]),
                     "pixels_without_the_colour_inside": int(holes),
                     "centre_read_as_day": mdates.num2date(
                         xmin + (centre - x0_px) / w_px * (xmax - xmin)).strftime("%Y-%m-%d %H:%M")}
            out["read"].append(entry)
            for d, tol in ((d_lo, g["t_lo"]), (d_hi, g["t_hi"])):
                out["ends_compared_under_each_tolerance"][names[tol]] += 1
                out["largest_difference_px_under_each_tolerance"][names[tol]] = round(max(
                    out["largest_difference_px_under_each_tolerance"][names[tol]], abs(d)), 2)
            if abs(d_lo) > g["t_lo"] or abs(d_hi) > g["t_hi"] or holes > (2 if len(g["what"]) > 1 else 0):
                out["outside_the_tolerance"].append(
                    f"{name}: {what}: asked for columns {g['lo']:.1f} to {g['hi']:.1f}, found {p} to {q + 1}"
                    + (f", {holes} pixels inside without the colour" if holes else ""))
        for k, (p, q) in enumerate(runs):
            if k not in used:
                out["not_asked_for"].append(f"{name}: columns {p} to {q}")
    out["largest_difference_px"] = max(out["largest_difference_px_under_each_tolerance"].values())
    out["ok"] = not out["not_found"] and not out["not_asked_for"] and not out["outside_the_tolerance"]
    return out


def marks_read_from_the_artists(strip, expected, xmin, xmax):
    """Reads the marks back from the artists of the figure (place, shape, kind of line) and sets them against
    the table."""
    out = {"marks_in_the_figure": len(strip["marks"]), "marks_asked_for": len(expected), "not_equal": []}
    by_row = {m["row"]: m for m in strip["marks"]}
    stems = dict(strip["stems"])

    def days(art):
        return [mdates.num2date(x).strftime("%Y-%m-%d") for x in np.atleast_1d(art.get_xdata())]

    def shape(art):
        dotted = art.get_linestyle() not in ("-", "solid", "None", "none", "", " ")
        return {"marker": art.get_marker() if art.get_marker() not in ("None", "", " ", None) else "",
                "line": "" if art.get_linestyle() in ("None", "none", "", " ") else ("dotted" if dotted
                                                                                    else "solid")}
    forms = {"triangle": {"marker": "v", "line": ""}, "cross": {"marker": "x", "line": ""},
             "bar": {"marker": "", "line": "solid"}, "dotted line": {"marker": "", "line": "dotted"},
             "tick": {"marker": "", "line": "solid"}}
    for e in expected:
        m = by_row.pop(e["row"], None)
        if m is None:
            out["not_equal"].append(f"row {e['row']}: no mark in the figure")
            continue
        got = [(what, days(art)) for what, art in m["artists"]]
        for what, art in m["artists"]:
            ys = [float(y) for y in np.atleast_1d(art.get_ydata())]
            want_y = {"triangle": [MARK_PT / 2], "cross": [MARK_PT / 2], "bar": [SPAN_Y_PT] * 2,
                      "dotted line": [SPAN_Y_PT] * 2, "arrow head": [SPAN_Y_PT],
                      "tick": [SPAN_Y_PT - SPAN_LW_PT / 2, MARK_PT]}[what]
            if len(ys) != len(want_y) or any(abs(a - b) > 1e-9 for a, b in zip(ys, want_y)):
                out["not_equal"].append(f"row {e['row']}: the {what} does not stand at its height")
            if what in forms and shape(art) != forms[what]:
                out["not_equal"].append(f"row {e['row']}: the {what} is drawn as {shape(art)}")
            if what == "arrow head" and art.get_marker() not in ("<", ">"):
                out["not_equal"].append(f"row {e['row']}: the arrow head is drawn as {art.get_marker()!r}")
        if e["mark"] in ("triangle", "cross"):
            want = [(e["mark"], [e["x0"]])]
        else:
            want = [("bar" if e["mark"] == "duration" else "dotted line", [e["x0"], e["x1"]])]
            for x, cut in ((e["x0"], e["cut_left"]), (e["x1"], e["cut_right"])):
                if cut:
                    want.append(("arrow head", [x]))
                elif e["mark"] == "duration":
                    want.append(("tick", [x, x]))
        if got != want:
            out["not_equal"].append(f"row {e['row']}: figure {got}, table {want}")
        st = stems.get(e["row"])
        if strip["words"]:
            if st is None:
                out["not_equal"].append(f"row {e['row']}: no thin line")
            elif any(abs(x - e["anchor"]) > 1e-6 for x in st.get_xdata()):
                out["not_equal"].append(f"row {e['row']}: the thin line does not stand at its place")
    for row in by_row:
        out["not_equal"].append(f"row {row}: a mark in the figure that the table does not ask for")
    out["ok"] = not out["not_equal"]
    return out


def texts_of_the_strip_checked(fig, strip, axs):
    """No label on another label, on a line of another mark, on the line of the declaration or on a mark; every
    label inside the strip; every line beside its own label."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    hits = []
    sb = axs.get_window_extent(r)
    boxes = [(t, t["artist"].get_window_extent(r)) for t in strip["texts"]]
    px = fig.dpi / 72.0
    for i, (t, bb) in enumerate(boxes):
        name = t["text"].replace("\n", " ")
        if bb.x0 < sb.x0 - 0.5 or bb.x1 > sb.x1 + 0.5 or bb.y0 < sb.y0 - 0.5 or bb.y1 > sb.y1 + 0.5:
            hits.append(f"'{name}' leaves the strip")
        for u, bu in boxes[i + 1:]:
            if bb.x0 < bu.x1 and bu.x0 < bb.x1 and bb.y0 < bu.y1 and bu.y0 < bb.y1:
                hits.append(f"'{name}' on '{u['text']}'".replace("\n", " "))
        lines = [(f"line of row {row}", st) for row, st in strip["stems"] if row != t.get("row")]
        lines.append(("line of the declaration", strip["declaration"]))
        for m in strip["marks"]:
            lines += [(f"{what} of row {m['row']}", art) for what, art in m["artists"]]
        for what, art in lines:
            ba = art.get_window_extent(r)
            pad = 0.5 * art.get_linewidth() * px if what.startswith("line") else 0.0
            if bb.x0 < ba.x1 + pad and ba.x0 - pad < bb.x1 and bb.y0 < ba.y1 and ba.y0 < bb.y1:
                hits.append(f"'{name}' on the {what}")
    return {"texts": len(boxes), "overlaps": hits, "ok": not hits}


def panels_against_the_figure_without_interventions(png_with, png_without, plan, banners_px):
    """Panels a, b and c of the picture with the strip, pixel by pixel against the picture without interventions.
    The rows of pixels are those of the three blocks; a block that moved is compared where it stands now."""
    from PIL import Image
    with Image.open(png_with) as im:
        A = np.asarray(im.convert("RGBA")).copy()
    with Image.open(png_without) as im:
        B = np.asarray(im.convert("RGBA")).copy()
    s = C.DPI / 72.0
    geo, lay = plan["geometry"], plan["layout"]
    grow = int(round(lay["growth_of_the_figure_pt"] * s))
    red_bc = int(round(lay["taken_from_the_gap_between_b_and_c_pt"] * s))
    if A.shape[0] - B.shape[0] != grow or A.shape[1] != B.shape[1]:
        return {"ok": False, "why": f"sizes {A.shape[:2]} and {B.shape[:2]}, growth asked for {grow} px"}
    for (x0, y0, x1, y1) in banners_px:              # a banner of a test is no part of a panel
        A[max(0, y0):y1, max(0, x0):x1] = 255
        yb0, yb1 = (y0, y1) if y1 <= A.shape[0] // 2 else (y0 - grow, y1 - grow)
        B[max(0, yb0):max(0, yb1), max(0, x0):x1] = 255
    HB = B.shape[0]

    def rows_from_bottom(y_in):                     # row of pixels (from the top of B) of a height in inches
        return int(round(HB - y_in * C.DPI))
    pad = int(round(PANEL_PAD_PT * s))
    a_bot = rows_from_bottom(geo["a"][0]) + pad
    b_top, b_bot = rows_from_bottom(geo["b"][1]), rows_from_bottom(geo["b"][0]) + pad
    c_top = rows_from_bottom(geo["c"][1]) - int(round(TITLE_BLOCK_PT * s))
    blocks = {"panel a (with its title, to below its time axis)": ((0, a_bot), 0),
              "panel b (from its upper edge to below its time axis)": ((b_top, b_bot), grow + red_bc),
              "panel c (with its title, to the lower edge of the picture)": ((c_top, HB), grow)}
    out = {"blocks": {}, "ok": True}
    for name, ((r0, r1), shift) in blocks.items():
        a_part, b_part = A[r0 + shift:r1 + shift], B[r0:r1]
        same = a_part.shape == b_part.shape and bool((a_part == b_part).all())
        n = int((a_part != b_part).any(axis=2).sum()) if a_part.shape == b_part.shape else -1
        out["blocks"][name] = {"rows_of_pixels": r1 - r0, "moved_down_by_px": shift, "pixels_that_differ": n,
                               "identical": same}
        out["ok"] = out["ok"] and same
    return out


# =============================================================================================
# 2. the drawing (from the file of values and from nothing else)
# =============================================================================================
def D(text):
    return mdates.date2num(dt.datetime.combine(dt.date.fromisoformat(text), dt.time()))


class HandlerOver(HandlerTuple):
    """Patch with the line drawn over it (earlier generator, class HandlerOver)."""

    def create_artists(self, legend, orig_handle, xdescent, ydescent, width, height, fontsize, trans):
        out = []
        for h in orig_handle[::-1]:
            handler = legend.get_legend_handler(legend.get_legend_handler_map(), h)
            out += handler.create_artists(legend, h, xdescent, ydescent, width, height, fontsize, trans)
        return out


def draw(values_path):
    rows = C.read_csv(values_path)
    font = C.style()

    def of(element, panel=None):
        return [r for r in rows if r["element"] == element and (panel is None or r["panel"] == panel)]

    def one(element, label=None):
        hits = [r for r in rows if r["element"] == element and (label is None or r["label"] == label)]
        if len(hits) != 1:
            raise C.Stop(f"file of values: {element} / {label}: expected 1 row, found {len(hits)}")
        return hits[0]

    def key(panel):
        return [r["label"] for r in of(E_KEY, panel)]

    XMIN_T, XMAX_T = one(E_AXIS_FIRST)["x"], one(E_AXIS_END)["x"]
    XMIN, XMAX = D(XMIN_T), D(XMAX_T)
    DECLARED = D(one(E_DECLARED)["x"])
    data_artists = {"a": [], "b": [], "c": []}

    drawn_int = [r for r in of(E_INT_MARK) if not r["drawn_as"].startswith("not drawn")]
    plan = plan_of_the_strip(drawn_int, of(E_INT_NUMBER), XMIN_T, XMAX_T, one(E_DECLARED)["x"]) \
        if drawn_int else None
    if plan is None:                           # no intervention to draw: the figure of 4 October
        fig = C.new_figure(HEIGHT)
        gs = fig.add_gridspec(3, 1, **GRID)
        axa = fig.add_subplot(gs[0])
        axb = fig.add_subplot(gs[1], sharex=axa)
        axc = fig.add_subplot(gs[2], sharex=axa)
        axs = None
    else:                                      # the same three panels, and the strip above panel b
        fig, axa, axb, axc, axs = figure_with_strip(plan)
    d0, d1 = dt.date.fromisoformat(XMIN_T), dt.date.fromisoformat(XMAX_T)
    month_starts = []
    y, m = d0.year, d0.month
    while dt.date(y, m, 1) <= d1:
        if dt.date(y, m, 1) >= d0:
            month_starts.append(dt.date(y, m, 1))
        y, m = (y + 1, 1) if m == 12 else (y, m + 1)
    axis_R = one(E_AXIS_RULE, "axis of R")
    R_AXIS_MAX = float(axis_R["x_end"])
    STRIP = STRIP_SHARE * R_AXIS_MAX
    for ax in (axa, axb, axc):
        ax.set_xlim(XMIN, XMAX)
        ax.set_xticks([D(d.isoformat()) for d in month_starts])
        ax.set_xticklabels([C.MONTHS[d.month - 1] for d in month_starts])
        ax.spines["left"].set_position(("outward", 5))
        if ax is not axb:
            ax.axvline(DECLARED, color=C.GREY_REF, ls=(0, (4, 3)), lw=0.9, zorder=1.5)
    axb.plot([DECLARED, DECLARED], [0, R_AXIS_MAX], color=C.GREY_REF, ls=(0, (4, 3)), lw=0.9, zorder=1.5)
    for ax in (axa, axb):
        ax.tick_params(labelbottom=False)
    axc.set_xlabel(one(E_AXIS_FIRST)["label"])

    # ----- a. effective population size --------------------------------------------------------------
    st = sorted(of(E_STAIR), key=lambda r: r["x"])
    edges = [D(r["x"]) for r in st] + [D(st[-1]["x_end"])]
    for a, b in zip(st[:-1], st[1:]):
        if a["x_end"] != b["x"]:
            raise C.Stop("file of values: the intervals of the staircase are not contiguous")
    med = [float(r["value"]) for r in st]
    lo = [float(r["lower"]) for r in st]
    hi = [float(r["upper"]) for r in st]
    FADE = D(one(E_FADE)["x"])
    k = int(np.searchsorted(np.array(edges), FADE))          # edges[k-1] < FADE <= edges[k]
    if k >= len(edges):
        raise C.Stop("the median date of the root lies after the last interval")
    if k > 0:
        e_old = edges[:k] + [FADE]
        h1 = axa.stairs(hi[:k], e_old, baseline=lo[:k], fill=True, facecolor="none", hatch="//////",
                        edgecolor=C.blend(C.BLUE, 0.45), linewidth=0, zorder=2)
        h2 = axa.stairs(med[:k], e_old, baseline=None, color=C.blend(C.BLUE, 0.55), lw=1.4,
                        ls=(0, (2.5, 1.5)), zorder=3)
        e_new = [FADE] + edges[k:]
        s1 = axa.stairs(hi[k - 1:], e_new, baseline=lo[k - 1:], fill=True, facecolor=C.BLUE,
                        alpha=C.BLUE_ALPHA, edgecolor="none", zorder=2)
        s2 = axa.stairs(med[k - 1:], e_new, baseline=None, color=C.BLUE, lw=1.7, zorder=3)
    else:
        s1 = axa.stairs(hi, edges, baseline=lo, fill=True, facecolor=C.BLUE, alpha=C.BLUE_ALPHA,
                        edgecolor="none", zorder=2)
        s2 = axa.stairs(med, edges, baseline=None, color=C.BLUE, lw=1.7, zorder=3)
    axis_ne = one(E_AXIS_RULE, "axis of Ne tau (logarithmic)")
    ne_lo, ne_hi = float(axis_ne["x"]), float(axis_ne["x_end"])
    axa.set_yscale("log")
    axa.set_ylim(ne_lo, ne_hi)
    decades = []
    t = C.decade_ceil(ne_lo)
    while t <= ne_hi * (1 + 1e-9):
        decades.append(t)
        t *= 10.0
    axa.set_yticks(decades)
    axa.set_yticklabels([f"{d:g}" for d in decades])
    axa.yaxis.set_minor_locator(mpl.ticker.LogLocator(base=10, subs=np.arange(2, 10) * 0.1, numticks=12))
    axa.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    axa.set_ylabel(r"$N_e\tau$ (years)")
    axa.yaxis.set_label_coords(YLAB_X, 0.5)
    for r in st:
        if "double dagger" in r["drawn_as"]:
            axa.text(0.5 * (D(r["x"]) + D(r["x_end"])), float(r["upper"]), C.DOUBLE_DAGGER, ha="center",
                     va="bottom", fontsize=C.SIZES[2], color=C.BLUE)
    # staircase as boxes of the canvas, for the test 'a key or a text over data'
    stair_boxes = [(D(r["x"]), D(r["x_end"]), float(r["lower"]), float(r["upper"])) for r in st]

    axg = axa.twinx()
    axg.set_zorder(axa.get_zorder() - 1)
    axa.patch.set_visible(False)
    gen = sorted(of(E_GENOMES), key=lambda r: r["x"])
    axg.bar([D(r["x"]) for r in gen], [int(r["value"]) for r in gen],
            width=[D(r["x_end"]) - D(r["x"]) for r in gen], align="edge", color=C.GREY_BAR, edgecolor="white",
            linewidth=0.6, zorder=1)
    axis_g = one(E_AXIS_RULE, "axis of the genomes per interval")
    G_SPINE = float(axis_g["x_end"])
    G_TOP = GENOME_AXIS_OVER_SPINE * G_SPINE
    axg.set_ylim(0, G_TOP)
    axg.set_yticks([v for v in range(0, int(G_SPINE) + 1, 50)])
    axg.spines["right"].set_visible(True)
    axg.spines["right"].set_bounds(0, G_SPINE)
    axg.spines["right"].set_position(("outward", 5))
    axg.spines["left"].set_visible(False)
    axg.tick_params(axis="y", colors="0.25")
    axg.spines["right"].set_color("0.25")
    t_gen = axa.annotate(key("a")[2], xy=(1, G_SPINE / G_TOP), xycoords="axes fraction", xytext=(5, 5),
                         textcoords="offset points", ha="right", va="bottom", fontsize=C.SIZES[1], color="0.25",
                         annotation_clip=False)

    root = one(E_ROOT)
    TM, TM_LO, TM_HI = D(root["value"]), D(root["lower"]), D(root["upper"])
    ytm = 0.075
    tr = mpl.transforms.blended_transform_factory(axa.transData, axa.transAxes)
    x0 = max(TM_LO, XMIN)
    axa.plot([x0, TM_HI], [ytm, ytm], color=C.BLUE, lw=1.3, transform=tr, zorder=4, solid_capstyle="butt",
             clip_on=False)
    axa.plot([TM_HI, TM_HI], [ytm - 0.022, ytm + 0.022], color=C.BLUE, lw=1.3, transform=tr, zorder=4)
    if TM_LO < XMIN:
        axa.plot([XMIN], [ytm], marker="<", ms=5.5, color=C.BLUE, transform=tr, zorder=5, clip_on=False,
                 mec="none")
    else:
        axa.plot([TM_LO, TM_LO], [ytm - 0.022, ytm + 0.022], color=C.BLUE, lw=1.3, transform=tr)
    axa.plot([TM], [ytm], marker="D", ms=5, color=C.BLUE, mec="white", mew=0.6, transform=tr, zorder=6)
    t_root = axa.text(XMIN + 5, ytm + 0.035, root["label"], transform=tr, ha="left", va="bottom",
                      fontsize=C.SIZES[1], color=C.BLUE)
    t_decl = axa.text(DECLARED - 3, 0.965, one(E_DECLARED)["label"], transform=tr, ha="right", va="top",
                      fontsize=C.SIZES[1], color="0.35")
    handles_a = [
        (Line2D([], [], color=C.BLUE, lw=1.7), Patch(facecolor=C.blend(C.BLUE, C.BLUE_ALPHA), edgecolor="none")),
        (Line2D([], [], color=C.blend(C.BLUE, 0.55), lw=1.4, ls=(0, (2.5, 1.5))),
         Patch(facecolor="none", hatch="//////", edgecolor=C.blend(C.BLUE, 0.45), linewidth=0)),
    ]
    leg_a = axa.legend(handles_a, key("a")[:2], loc="upper left", bbox_to_anchor=(0.0, 0.905),
                       fontsize=C.SIZES[1], handlelength=2.2, borderaxespad=0.2, labelspacing=0.35,
                       handler_map={tuple: HandlerOver()})
    C.panel_head(axa, "a", "Effective population size")

    # ----- b. reproduction number ----------------------------------------------------------------------
    axb.set_ylim(0, R_AXIS_MAX + STRIP)
    axb.set_yticks(list(range(0, int(round(R_AXIS_MAX)) + 1)))
    axb.spines["left"].set_bounds(0, R_AXIS_MAX)
    axb.set_ylabel("reproduction number $R$")
    axb.yaxis.set_label_coords(YLAB_X, 0.5 * R_AXIS_MAX / (R_AXIS_MAX + STRIP))
    axb.axhline(1.0, color="0.15", lw=0.7, zorder=1.4)
    PALE = C.blend(C.BLUE, 0.50)
    u = R_AXIS_MAX / 4.0                       # the earlier positions were typed for an axis that ends at 4
    for r in of(E_PAIR):
        is_open = r["drawn_as"].startswith("open")
        x, m_, lo_, hi_ = D(r["x"]), float(r["value"]), float(r["lower"]), float(r["upper"])
        if hi_ > R_AXIS_MAX:
            raise C.Stop("an interval of a pair extends beyond the axis of R")
        ln, = axb.plot([x, x], [lo_, hi_], color=PALE, lw=0.8, zorder=2, solid_capstyle="butt", clip_on=False)
        mk, = axb.plot([x], [m_], marker="o", ms=3.6, mfc="white" if is_open else PALE, mec=PALE, mew=0.8,
                       zorder=2.2, clip_on=False)
        data_artists["b"] += [("interval of a pair", ln), ("symbol of a pair", mk)]
        if "double dagger" in r["drawn_as"]:
            axb.text(x, hi_, C.DOUBLE_DAGGER, ha="center", va="bottom", fontsize=C.SIZES[2], color=PALE)
    prior = {r["label"]: r for r in of(E_PRIOR_SIM)}
    texts = {r["x"]: r for r in of(E_PERIOD_TEXT)}
    names = {r["x"]: r for r in of(E_PERIOD_NAME)}
    for r in of(E_PERIOD):
        start, end = D(r["x"]), D(r["x_end"])
        wdt = end - start
        m_, lo_, hi_ = float(r["value"]), float(r["lower"]), float(r["upper"])
        q = prior[r["label"]]
        qlo, qhi = float(q["lower"]), float(q["upper"])
        axb.add_patch(Rectangle((start, qlo), wdt, qhi - qlo, facecolor=C.blend(C.GREY_PRIOR, 0.45),
                                edgecolor="none", linewidth=0, zorder=1.0))
        p2 = axb.add_patch(Rectangle((start, qlo), wdt, qhi - qlo, facecolor="none", edgecolor="0.35",
                                     linewidth=0.9, linestyle=(0, (3, 2)), zorder=1.3))
        axb.add_patch(Rectangle((start, lo_), wdt, hi_ - lo_, facecolor=mpl.colors.to_rgba(C.BLUE, 0.22),
                                edgecolor="none", linewidth=0, zorder=1.6))
        p4 = axb.add_patch(Rectangle((start, lo_), wdt, hi_ - lo_, facecolor="none", edgecolor=C.BLUE,
                                     linewidth=1.0, zorder=3))
        axb.plot([start, end], [m_, m_], color=C.BLUE, lw=2.6, zorder=1.9, solid_capstyle="butt")
        data_artists["b"] += [("range under the prior alone", p2), ("rectangle of a period", p4)]
        yb = R_AXIS_MAX + 0.13 * u
        axb.plot([start + 0.7, start + 0.7, end - 0.7, end - 0.7], [yb - 0.07 * u, yb, yb, yb - 0.07 * u],
                 color=C.BLUE, lw=0.8, zorder=3, clip_on=False)
        axb.text(start + wdt / 2, yb + 0.06 * u, texts[r["x"]]["label"], ha="center", va="bottom",
                 fontsize=C.SIZES[1], color=C.BLUE)
        axb.text(start + wdt / 2, R_AXIS_MAX + STRIP - 0.02 * u, names[r["x"]]["label"], ha="center", va="top",
                 fontsize=C.SIZES[1], color="0.1")
    handles_b = [
        (Patch(facecolor=mpl.colors.to_rgba(C.BLUE, 0.22), edgecolor=C.BLUE, linewidth=1.0),
         Line2D([], [], color=C.BLUE, lw=2.6)),
        Patch(facecolor=C.blend(C.GREY_PRIOR, 0.45), edgecolor="0.35", linewidth=0.9, linestyle=(0, (3, 2))),
        Line2D([], [], color=PALE, lw=0.8, marker="o", ms=3.6, mfc=PALE, mec=PALE),
        Line2D([], [], color=PALE, lw=0.8, marker="o", ms=3.6, mfc="white", mec=PALE, mew=0.8),
    ]
    leg_b = axb.legend(handles_b, key("b"), loc="upper right", fontsize=C.SIZES[1], handlelength=2.2,
                       bbox_to_anchor=(1.0, (R_AXIS_MAX - 0.05 * u) / (R_AXIS_MAX + STRIP)), borderaxespad=0.0,
                       labelspacing=0.35, handler_map={tuple: HandlerTuple(ndivide=1, pad=0)})
    C.panel_head(axb if axs is None else axs, "b", "Reproduction number")
    strip = None if axs is None else draw_strip(fig, axs, plan, XMIN, XMAX, DECLARED)

    # ----- c. confirmed cases -----------------------------------------------------------------------------
    BARW = 6.2
    cs = of(E_CASES)
    week_ends = sorted({r["x_end"] for r in cs})
    centre = {w: D(w) - 2.5 for w in week_ends}
    bottom = {w: 0.0 for w in week_ends}
    for _, lab, colour in PROVINCES:
        part = {r["x_end"]: float(r["value"]) for r in cs if r["label"].startswith(lab + ", week ending ")}
        if sorted(part) != week_ends:
            raise C.Stop(f"file of values: the weeks of {lab} are not those of the series")
        bars = axc.bar([centre[w] for w in week_ends], [part[w] for w in week_ends], width=BARW,
                       bottom=[bottom[w] for w in week_ends], color=colour, edgecolor="none", zorder=2, label=lab)
        data_artists["c"] += [("bar of the cases", b) for b in bars]
        for w in week_ends:
            bottom[w] += part[w]
    hatch = of(E_HATCH)
    if hatch:
        axc.bar([centre[r["x_end"]] for r in hatch], [float(r["value"]) for r in hatch], width=BARW, color="none",
                edgecolor="white", hatch="////", linewidth=0, zorder=2.5)
    pv = one(E_PROVISIONAL)
    axis_c = one(E_AXIS_RULE, "axis of the cases per week")
    C_TOP = float(axis_c["x_end"])
    pvx, pvy = centre[pv["x_end"]], float(pv["value"])
    axc.bar([pvx], [pvy], width=BARW, facecolor="none", edgecolor="black", linewidth=1.0,
            linestyle=(0, (2.5, 1.8)), zorder=3)
    t_prov = axc.text(XMAX, PROV_TEXT_SHARE * C_TOP, pv["label"], ha="right", va="bottom", fontsize=C.SIZES[1],
                      color="0.1")
    if pvy + PROV_GAP_SHARE * C_TOP >= PROV_LINE_SHARE * C_TOP:
        raise C.Stop("the bar of the provisional week reaches the word 'provisional'")
    axc.plot([pvx, pvx], [pvy + PROV_GAP_SHARE * C_TOP, PROV_LINE_SHARE * C_TOP], color="0.3", lw=0.6, zorder=3)
    axc.set_ylim(0, C_TOP)
    axc.set_yticks([v for v in range(0, int(C_TOP), 200)])
    axc.set_ylabel("cases per week")
    axc.yaxis.set_label_coords(YLAB_X, 0.5)
    kc = key("c")
    handles_c = [Patch(facecolor=c, edgecolor="none") for _, _, c in PROVINCES]
    if len(kc) > len(PROVINCES):
        handles_c.append(Patch(facecolor=C.ORANGE_ITURI, edgecolor="white", hatch="////", linewidth=0))
    leg_c = axc.legend(handles_c, kc, loc="upper left", bbox_to_anchor=(0.0, 1.0), fontsize=C.SIZES[1],
                       handlelength=1.6, borderaxespad=0.2, labelspacing=0.35)
    C.panel_head(axc, "c", "Confirmed cases by date of report")

    # ----- keys and words over data (conventions: 'Legends never cover data') ------------------------------
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    over = []
    # panel a: boxes of the staircase band on the canvas
    stair_px = []
    for (xa, xb, ya, yb_) in stair_boxes:
        p0 = axa.transData.transform((xa, ya))
        p1 = axa.transData.transform((xb, yb_))
        stair_px.append(mpl.transforms.Bbox.from_extents(p0[0], p0[1], p1[0], p1[1]))
    for what, art in (("key of panel a", leg_a), ("dates of the root", t_root), ("words 'outbreak declared'", t_decl),
                      ("words 'genomes per interval'", t_gen)):
        bb = art.get_window_extent(rend)
        for sb in stair_px:
            if bb.overlaps(sb):
                over.append((what, "band of the staircase"))
                break
    tb = t_root.get_window_extent(rend)
    under = [sb for sb in stair_px if sb.x1 > tb.x0 and sb.x0 < tb.x1]
    gap_pt = min((sb.y0 - tb.y1) for sb in under) * 72.0 / fig.dpi if under else None
    bb = leg_b.get_window_extent(rend)
    over += C.boxes_overlap_artists(fig, [("key of panel b", bb)], data_artists["b"])
    over += C.boxes_overlap_artists(fig, [("key of panel c", leg_c.get_window_extent(rend)),
                                          ("word 'provisional'", t_prov.get_window_extent(rend))],
                                    data_artists["c"])
    extra = {"keys_or_words_over_data": sorted(set(over)), "font": font, "height_in": HEIGHT,
             "distance_between_the_dates_of_the_root_and_the_band_pt": gap_pt,
             "hatched_intervals": k, "months_on_the_axis": [C.MONTHS[d.month - 1] for d in month_starts]}
    banners = []
    for r in of(E_BANNER):
        top = "top" in r["drawn_as"]
        banners.append(fig.text(0.995, 0.997 if top else 0.003, r["label"], ha="right",
                                va="top" if top else "bottom", fontsize=C.SIZES[1], fontweight="bold",
                                color=BANNER_COLOUR))
    fixtures = [r for r in of(E_INT_MARK) + of(E_INT_WORDS) if FIXTURE_WORD in r["label"]]
    if bool(fixtures) != (mode_of_the_banners(rows) == "test"):
        raise C.Stop("file of values: a picture with fixtures, and no other, carries the two banners of a test")
    if plan is not None:
        extra["height_in"] = plan["layout"]["height_of_the_figure_in"]
        extra["strip"] = {"rule": plan["rule"], "why": plan["why"], "heights_used": plan["heights_used"]}
    if plan is not None or banners:
        fig.interventions = {"axes": axs, "axa": axa, "axb": axb, "axc": axc, "strip": strip, "plan": plan,
                             "banners": banners}
    return fig, extra


# =============================================================================================
# 3. run
# =============================================================================================
def main_of_4_october(argv, names=None):
    """The run of the program of 4 October, line by line (no table of the events, or a table with no row).
    names: the names of the files; without them the names of the figure without the strip."""
    names = names or {"values": VALUES_FILE, "figure": FIGURE_STEM, "checks": CHECKS_FILE}
    if len(argv) >= 3 and argv[1] == "--draw-only":
        outdir = argv[3] if len(argv) > 3 else "."
        os.makedirs(outdir, exist_ok=True)
        fig, extra = draw(argv[2])
        info = C.save_figure(fig, stem_of_a_file_of_values(argv[2]), outdir)
        print(json.dumps({"drawn_from": os.path.basename(argv[2]), "png": os.path.basename(info["png"]),
                          "keys_or_words_over_data": extra["keys_or_words_over_data"]}))
        return 0
    if len(argv) < 2:
        raise C.Stop("usage: python fig_ne_R_cases_20261006.py INPUTS.json [OUTDIR] [--interventions NAME.csv]")
    outdir = argv[2] if len(argv) > 2 else "."
    os.makedirs(outdir, exist_ok=True)
    I = C.Inputs(argv[1])
    V, info = build_values(I)
    values_path = V.save(os.path.join(outdir, names["values"]))      # the file of values is written first
    fig, extra = draw(values_path)                                    # ... and the figure is drawn from it
    patterns = C.word_patterns(I)
    checks = C.check_figure(fig, patterns)
    saved = C.save_figure(fig, names["figure"], outdir)
    rb = C.read_back(values_path, I)
    rows = C.read_csv(values_path)
    by = {}
    for r in rows:
        by[f"{r['panel']} | {r['element']}"] = by.get(f"{r['panel']} | {r['element']}", 0) + 1
    marked = [f"{r['panel']} | {r['element']} | {r['label']}" for r in rows
              if C.marked_by_R4f(r["meets_convergence_criterion"])]
    saved = {k: (os.path.basename(v) if isinstance(v, str) else v) for k, v in saved.items()}
    out = {"program": os.path.basename(__file__), "figure": saved, "values_file": os.path.basename(values_path),
           "figure_sha256": {"png": C.sha256(os.path.join(outdir, saved["png"])),
                             "pdf": C.sha256(os.path.join(outdir, saved["pdf"]))},
           "values_sha256": C.sha256(values_path), "rows": len(rows), "rows_by_panel_and_element": by,
           "read_back": rb, "rows_not_yet_delivered": sum(r["delivered"] == "no" for r in rows),
           "rows_marked_by_R4f": marked, "rows_not_drawn": sum(r["drawn_as"] == "not drawn" for r in rows),
           "figure_checks": checks, "drawing": extra, "info": info,
           "colours_and_sizes_against_the_conventions": C.colours_against_the_conventions(I),
           "table_of_inputs_version": I.version_of_the_table_of_inputs}
    out["values_marked_by_R4f_with_the_cell_as_it_reads"] = marked_with_their_cells(rows)
    out["not_yet_delivered"] = not_yet_delivered(rows, checks)
    out["banner"] = no_banner(checks, names)
    C.json_dump(out, os.path.join(outdir, names["checks"]))
    ok = (checks["ok"] and not rb["not_equal"] and not extra["keys_or_words_over_data"]
          and out["not_yet_delivered"]["ok"] and out["banner"]["ok"])
    print(json.dumps({"png": saved["png"], "rows": len(rows), "checks_ok": bool(ok),
                      "min_fontsize": checks["min_fontsize"], "read_back_not_equal": len(rb["not_equal"]),
                      "over_data": extra["keys_or_words_over_data"], "overlaps": checks["text_overlaps"],
                      "outside": checks["outside_figure"], "on_spine": checks["text_on_spine"],
                      "words": checks["words_of_R4j"]}))
    return 0


def names_of_the_outputs(table, preview=False):
    """Names of the files of a run with a table of the events: those of the figure of the article, or those
    of a test with fixtures or of a preview."""
    if table["test"]:
        stem, pre = "fig_ne_R_cases_" + table["tag"] + PROGRAM_DATE, MARKER + "_"
    elif preview:
        stem, pre = "fig_ne_R_cases_" + table["tag"] + PROGRAM_DATE, PREVIEW + "_"
    else:
        return dict(ARTICLE_NAMES)
    return {"values": pre + stem + "_values.csv", "figure": pre + stem, "checks": pre + stem + "_checks.json",
            "legend": pre + stem + LEGEND_SUFFIX}


def stem_of_a_file_of_values(path):
    """--draw-only: the name of the picture follows the name of the file of values."""
    name = os.path.basename(path)
    m = re.match(r"^((?:" + MARKER + "_|" + PREVIEW + r"_|INTERNAL_)?fig_ne_R_cases_[A-Za-z0-9_]+)_values\.csv$",
                 name)
    if not m:
        return FIGURE_STEM
    return m.group(1)


def marked_with_their_cells(rows):
    """The values marked by R4 (f), each with its cell meets_convergence_criterion as it reads."""
    return [{"panel": r["panel"], "element": r["element"], "label": r["label"], "analysis": r["analysis"],
             "value": r["value"], "source_table": r["source_table"], "source_row_0based": r["source_row_0based"],
             "meets_convergence_criterion": r["meets_convergence_criterion"]}
            for r in rows if C.marked_by_R4f(r["meets_convergence_criterion"])]


def not_yet_delivered(rows, checks):
    """No row 'not yet delivered' is expected."""
    n = sum(r["delivered"] == "no" for r in rows)
    texts = [t for t, _ in checks["texts"] if C.NOT_YET in t]
    return {"expected": 0, "rows_of_the_file_of_values_with_delivered_no": n,
            "texts_of_the_picture_that_hold_the_words": texts, "ok": n == 0 and not texts}


def no_banner(checks, names, mode=""):
    """The picture of the article holds no banner, and the words INTERNAL, PROVISIONAL and PREVIEW stand
    nowhere in it (compared as written, in capitals). The same words in another case are listed with the
    text that holds them."""
    every = [x for pair in BANNERS.values() for x in pair]
    have = [t for t, _ in checks["texts"] if t in every]
    article = all(n in ARTICLE_NAMES.values() for n in names.values())
    in_names = sorted({w for n in names.values() for w in NOT_IN_A_NAME_OF_THE_ARTICLE if w.lower() in n.lower()})
    out = {"figure_of_the_article": article, "banners_in_the_picture": have,
           "words_of_a_banner_in_a_text_of_the_picture": checks["words_of_a_banner"],
           "the_same_words_in_another_case": checks["words_of_a_banner_in_another_case"],
           "words_in_the_names_of_the_files": in_names if article else "not compared: no file of the article"}
    out["ok"] = (not have and not checks["words_of_a_banner"] and (not article or not in_names)) if not mode \
        else True
    return out


def mode_of_the_banners(rows):
    """'test', 'preview' or '' from the rows of the banners of a file of values."""
    texts = tuple(r["label"] for r in rows if r["element"] == E_BANNER)
    if not texts:
        return ""
    for mode, pair in BANNERS.items():
        if texts == tuple(pair):
            return mode
    raise C.Stop("file of values: the rows of the banners are not those of a test or of a preview")


def box_in_the_png(fig, artist):
    """(x0, y0, x1, y1) of an artist in pixels of the saved PNG, rows counted from the top."""
    dpi = fig.dpi
    fig.set_dpi(C.DPI)                          # the box is measured at the resolution of the PNG
    fig.canvas.draw()
    bb = artist.get_window_extent(fig.canvas.get_renderer())
    fig.set_dpi(dpi)
    fig.canvas.draw()
    f = 1.0
    H = int(round(fig.get_figheight() * C.DPI))
    pad = int(round(BANNER_PAD_PT * C.DPI / 72.0))
    return (int(math.floor(bb.x0 * f)) - pad, int(math.floor(H - bb.y1 * f)) - pad,
            int(math.ceil(bb.x1 * f)) + pad, int(math.ceil(H - bb.y0 * f)) + pad)


def items_of(rows):
    """The interventions that are drawn, in the order of their numbers, from the rows of the file of values."""
    def by(element):
        return {r["source_row_0based"]: r["label"] for r in rows if r["element"] == element}
    numbers, words, sources, kinds = by(E_INT_NUMBER), by(E_INT_WORDS), by(E_INT_SOURCE), by(E_INT_KIND)
    out = []
    for r in rows:
        k = r["source_row_0based"]
        if r["element"] != E_INT_MARK or k not in numbers:
            continue
        mark = mark_of_a_row(r["drawn_as"])
        out.append({"row": k, "number": int(numbers[k]), "mark": mark, "kind": kinds.get(k, ""),
                    "label": r["label"], "words": words[k], "source": sources[k], "first": r["x"],
                    "last": r["x_end"] if mark in ("duration", "day not known") else "",
                    "cut_left": "begins before" in r["drawn_as"], "cut_right": "ends after" in r["drawn_as"]})
    return sorted(out, key=lambda e: e["number"])


def register_sources(I, name):
    """Registers the table that resolves the keys of the sources; a fault stops the program before a file is
    written. Returns the keys of the table."""
    if not I.has(name):
        raise C.Stop(f"the JSON of the inputs holds no entry, or no file, for the table of sources {name}")
    if name not in I.by_file:
        version = I.paths.get(name + "#version", "")
        if not version and not name.startswith(MARKER + "_"):
            raise C.Stop(f"{name}: no version: neither the table of inputs nor an entry '{name}#version' of the "
                         "JSON of the inputs names one")
        I.by_file[name] = {"file": name, "version": version, "used_for": "Figure 1, panel b",
                           "what_it_is": "table of the sources of the interventions", "analysis": "",
                           "state_of_the_file": "", "note": ""}
    t = I.table(name).need(SOURCES_KEY_COLUMN)
    if not t.rows:
        raise C.Stop(f"{name}: the table of sources holds no row")
    return [r[SOURCES_KEY_COLUMN] for r in t.rows]


def table_of_sources(I, name, used):
    """The keys of the sources of the interventions against the table that resolves them."""
    keys = register_sources(I, name)
    return {"table": name, "version": I.version(name), "keys_of_the_table": len(keys),
            "keys_named_twice": sorted({k for k in keys if keys.count(k) > 1}),
            "keys_of_the_interventions": sorted(used),
            "keys_that_the_table_does_not_hold": sorted(set(used) - set(keys)),
            "keys_of_the_table_that_no_intervention_names": sorted(set(keys) - set(used)),
            "ok": not (set(used) - set(keys)) and len(set(keys)) == len(keys)}


def strip_against_the_panels(fig, strip, axs, axa, axb):
    """Nothing of the strip stands on panel b or on panel a (so that no mark hides a rectangle, a point, the
    line R = 1 or the line of the declaration inside the panel)."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    arts = [("baseline of the strip", axs.spines["bottom"]), ("line of the declaration in the strip",
                                                               strip["declaration"])]
    arts += [(f"line of row {row}", st) for row, st in strip["stems"]]
    arts += [(f"text '{t['text']}'".replace("\n", " "), t["artist"]) for t in strip["texts"]]
    for m in strip["marks"]:
        arts += [(f"{what} of row {m['row']}", a) for what, a in m["artists"]]
    hits = []
    pb = axb.get_window_extent(r)
    pa = axa.get_tightbbox(r)
    lo = None
    for what, a in arts:
        bb = a.get_window_extent(r)
        lo = bb.y0 if lo is None else min(lo, bb.y0)
        if bb.y0 < pb.y1 - 0.01 and bb.y1 > pb.y0 and bb.x1 > pb.x0 and bb.x0 < pb.x1:
            hits.append(f"{what} on panel b")
        if bb.y1 > pa.y0 + 0.01 and bb.y0 < pa.y1 and bb.x1 > pa.x0 and bb.x0 < pa.x1:
            hits.append(f"{what} on panel a")
    return {"artists_of_the_strip": len(arts), "on_a_panel": hits,
            "lowest_point_of_the_strip_above_the_upper_edge_of_panel_b_pt":
                None if lo is None else round((lo - pb.y1) * 72.0 / fig.dpi, 2),
            "ok": not hits}


def checks_of_the_interventions(I, table, fig, values_path, png_path, n_head, summary, patterns, names,
                                sources=None):
    """The checks that the task of 5 October adds. Returns (dict of the checks, text of the legend file)."""
    rows = C.read_csv(values_path)
    xmin = dt.date.fromisoformat([r for r in rows if r["element"] == E_AXIS_FIRST][0]["x"])
    xmax = dt.date.fromisoformat([r for r in rows if r["element"] == E_AXIS_END][0]["x"])
    out = {"table": summary["table"], "version": summary["version"], "version_from": summary["version_from"],
           "test_with_fixtures": summary["test"], "preview": summary["preview"],
           "rows_of_the_table": summary["rows"], "on_the_time_axis": summary["on_the_time_axis"],
           "form_of_the_table": {"form": summary["form_of_the_table"],
                                 "what_the_program_reads_from_which_column": summary["columns_of_the_table"],
                                 "forms_that_the_program_knows": {f: dict(list(d["required"].items())
                                                                          + list(d["optional"].items()))
                                                                  for f, d in FORMS.items()}},
           "not_drawn_outside_the_time_axis": summary["not_drawn_outside_the_time_axis"],
           "periods_cut_at_the_axis": summary["periods_cut_at_the_axis"],
           "words_of_R4j_in_the_cells": summary["words_of_R4j_in_the_cells"]}
    info = getattr(fig, "interventions", None) or {}
    strip, plan = info.get("strip"), info.get("plan")
    table_again = interventions_of_the_table(I, table, xmin, xmax)
    expected = [e for e in table_again if e["on_the_axis"]]
    out["kinds_and_marks"] = {
        "rule_of_the_kind": f"a kind that matches {KIND_WITH_A_CROSS.pattern!r} is drawn as a cross on its day, "
                            "every other kind as a triangle; a period shows no kind",
        "rule_of_the_meaning": f"a period whose meaning matches {DAY_NOT_KNOWN.pattern!r} is drawn as a dotted "
                               "line (the day is not known), every other period as a bar with a tick at each end",
        "rows": [{"row": e["row"], "kind": e["kind"], "meaning_of_the_period": e["meaning"],
                  "days": e["span"], "mark": e["mark"] if e["on_the_axis"] else "not drawn"} for e in table_again],
        "marks_by_kind": {k: sorted({e["mark"] for e in expected if e["kind"] == k})
                          for k in sorted({e["kind"] for e in expected})}}
    # 1. the rows of 4 October stand unchanged at the head of the file of values
    with tempfile.TemporaryDirectory() as tmp:
        V0, _ = build_values(I)
        p0 = V0.save(os.path.join(tmp, VALUES_FILE))
        with open(p0, "rb") as fh:
            head = fh.read()
        with open(values_path, "rb") as fh:
            whole = fh.read()
        out["head_of_the_file_of_values"] = {
            "rows_of_4_october": n_head, "rows_built_again_without_a_table": len(V0.rows),
            "sha256_of_the_head": hashlib.sha256(whole[:len(head)]).hexdigest(),
            "sha256_without_a_table": hashlib.sha256(head).hexdigest(),
            "ok": whole[:len(head)] == head and len(V0.rows) == n_head}
        # 2. the rows of the file of values against the table read again
        marks = [r for r in rows if r["element"] == E_INT_MARK]
        faults = []
        if [r["source_row_0based"] for r in marks] != [e["row"] for e in table_again]:
            faults.append("the rows of the marks are not the rows of the table, in its order")
        for r, e in zip(marks, table_again):
            if (r["x"], r["x_end"], r["label"]) != (e["first"], e["last_cell"], e["label"]):
                faults.append(f"row {e['row']}: days or label differ from the table")
            if (not r["drawn_as"].startswith("not drawn")) != e["on_the_axis"]:
                faults.append(f"row {e['row']}: 'drawn_as' says {r['drawn_as']!r}, the table and the axis say "
                              f"{'on' if e['on_the_axis'] else 'outside'} the time axis")
            elif e["on_the_axis"] and mark_of_a_row(r["drawn_as"]) != e["mark"]:
                faults.append(f"row {e['row']}: 'drawn_as' says {mark_of_a_row(r['drawn_as'])}, the kind and the "
                              f"meaning of the table give {e['mark']}")
        for element, key in ((E_INT_KIND, "kind"), (E_INT_MEANING, "meaning"), (E_INT_WORDS, "words"),
                             (E_INT_SOURCE, "source")):
            got = {r["source_row_0based"]: r["label"] for r in rows if r["element"] == element}
            for e in table_again:
                if e["row"] in got and got[e["row"]] != e[key]:
                    faults.append(f"row {e['row']}: the {key} of the file of values differs from the table")
        numbers = {r["source_row_0based"]: int(r["label"]) for r in rows if r["element"] == E_INT_NUMBER}
        if numbers != {e["row"]: e["number"] for e in expected}:
            faults.append("the numbers of the marks are not those of the rule (first day, then order of the rows)")
        order = {r["source_row_0based"]: int(r["label"]) for r in rows if r["element"] == E_INT_ORDER}
        if order:
            by_number = [order[k] for k, n in sorted(numbers.items(), key=lambda x: x[1])]
            out["order_of_the_table"] = {
                "column": summary["columns_of_the_table"].get("order"),
                "order_of_the_rows_that_are_drawn_by_their_number": by_number,
                "the_numbers_follow_the_order_of_the_table": by_number == sorted(by_number),
                "the_numbers_are_the_order_of_the_table": by_number == list(range(1, len(by_number) + 1)),
                "ok": by_number == sorted(by_number)}
        out["rows_against_the_table"] = {"rows_of_marks": len(marks), "faults": faults, "ok": not faults}
        listed = [str(x["row"]) for x in summary["not_drawn_outside_the_time_axis"]]
        drawn_rows = [m["row"] for m in strip["marks"]] if strip else []
        out["outside_the_time_axis"] = {
            "listed": listed, "by_the_table_read_again": [e["row"] for e in table_again if not e["on_the_axis"]],
            "drawn_although_outside": [k for k in listed if k in drawn_rows]}
        out["outside_the_time_axis"]["ok"] = (
            listed == out["outside_the_time_axis"]["by_the_table_read_again"]
            and not out["outside_the_time_axis"]["drawn_although_outside"])
        # 3. to 6. the strip
        if strip is None:
            ok_none = not expected
            out["rule_of_the_marks"] = {"rule": "no strip: no intervention lies on the time axis", "ok": ok_none}
            for k in ("marks_read_from_the_artists", "marks_read_from_the_pixels", "texts_of_the_strip",
                      "strip_against_the_panels"):
                out[k] = {"ok": ok_none, "note": "no strip is drawn"}
        else:
            axs = info["axes"]
            out["rule_of_the_marks"] = {
                "rule": plan["rule"], "why": plan["why"], "heights_used": plan["heights_used"],
                "tried": plan["tried"], "placements_searched": plan["nodes_searched"],
                "room_pt": plan["room"],
                "room_for_the_strip_above_panel_b_the_figure_allows_pt": round(plan["room_total_pt"], 2),
                "of_it_for_the_labels_the_figure_allows_pt": round(plan["room_of_the_labels_pt"], 2),
                "strip_height_pt": round(plan["strip_height_pt"], 3),
                "what_the_labels_in_words_would_need": plan.get("what_the_labels_in_words_would_need"),
                "layout": {k: round(v, 4) for k, v in plan["layout"].items()},
                "labels": [{"row": t.get("row"), "numbers": t.get("numbers"), "text": t["text"],
                            "side": t.get("side"), "height": t.get("height"), "lines": t.get("lines")}
                           for t in strip["texts"]],
                "ok": plan["rule"] in ("labels in words", "numbered marks")
                and (plan["rule"] == "numbered marks") == (plan["heights_used"] == 0)}
            out["marks_read_from_the_artists"] = marks_read_from_the_artists(strip, expected, D(xmin.isoformat()),
                                                                              D(xmax.isoformat()))
            pos = axs.get_position()
            out["marks_read_from_the_pixels"] = marks_read_from_the_png(
                png_path, (pos.x0, pos.y0, pos.width, pos.height), D(xmin.isoformat()), D(xmax.isoformat()),
                expected, strip["words"])
            out["marks_read_from_the_pixels"]["strip_in_fractions_of_the_figure"] = [pos.x0, pos.y0, pos.width,
                                                                                    pos.height]
            out["texts_of_the_strip"] = texts_of_the_strip_checked(fig, strip, axs)
            out["strip_against_the_panels"] = strip_against_the_panels(fig, strip, axs, info["axa"], info["axb"])
            out["colour_of_the_strip"] = colour_of_the_strip_checked(I, fig, info)
        # 7. panels a, b and c against the picture without interventions, pixel by pixel
        fig0, _ = draw(p0)
        saved0 = C.save_figure(fig0, "without", tmp)
        plt_close(fig0)
        banners_px = [box_in_the_png(fig, b) for b in info.get("banners", [])]
        plan0 = plan if plan is not None else {
            "geometry": geometry_of_4_october(),
            "layout": {"growth_of_the_figure_pt": 0.0, "taken_from_the_gap_between_b_and_c_pt": 0.0}}
        out["panels_against_the_picture_without_interventions"] = panels_against_the_figure_without_interventions(
            png_path, saved0["png"], plan0, banners_px)
        out["panels_against_the_picture_without_interventions"]["sha256_of_the_picture_without_interventions"] = \
            C.sha256(saved0["png"])
    # 8. the banners
    mode = "test" if summary["test"] else ("preview" if summary["preview"] else "")
    texts = [t.get_text() for t in C.visible_texts(fig)]
    every = [x for pair in BANNERS.values() for x in pair]
    have = [t for t in texts if t in every]
    prefix = {"test": MARKER + "_", "preview": PREVIEW + "_", "": ""}[mode]
    in_names = sorted({w for n in names.values() for w in NOT_IN_A_NAME_OF_THE_ARTICLE if w.lower() in n.lower()})
    words_in_texts = [(w, t[:60]) for t in texts for w in C.WORDS_OF_A_BANNER
                      if re.search(r"(?<![A-Za-z])" + w + r"(?![A-Za-z])", t)]
    other_case = [(w, t[:60]) for t in texts for w in C.WORDS_OF_A_BANNER
                  if re.search(r"(?<![A-Za-z])" + w + r"(?![A-Za-z])", t, flags=re.I)
                  and not re.search(r"(?<![A-Za-z])" + w + r"(?![A-Za-z])", t)]
    out["banner"] = {"test_with_fixtures": summary["test"], "preview": summary["preview"],
                     "figure_of_the_article": not mode,
                     "banners_in_the_picture": have,
                     "words_of_a_banner_in_a_text_of_the_picture": words_in_texts,
                     "the_same_words_in_another_case": other_case,
                     "every_name_begins_with": prefix,
                     "names_begin_so": all(n.startswith(prefix) for n in names.values()),
                     "words_in_the_names_of_the_files": in_names,
                     "ok": (sorted(have) == sorted(BANNERS[mode]) and all(n.startswith(prefix)
                                                                         for n in names.values())) if mode
                     else (not have and not words_in_texts and not in_names
                           and sorted(names.values()) == sorted(ARTICLE_NAMES.values()))}
    # 9. the keys of the sources
    if sources:
        used = sorted({k for e in table_again if e["source"] and "keys" in summary["columns_of_the_table"]
                       for k in e["source"].split("; ")})
        out["sources_against_their_table"] = table_of_sources(I, sources, used)
    # 10. the sentence for the legend
    legend_text = None
    if expected:
        rule = plan["rule"]
        sentence, tokens = sentence_for_the_legend(I, items_of(rows), rule)
        again, _ = sentence_for_the_legend(I, sorted(expected, key=lambda e: e["number"]), rule)
        head_lines = [
            f"# sentence for the legend of Figure 1, panel b (unit {UNIT_OF_THE_LEGEND} of the template of the "
            "article); written by " + os.path.basename(__file__),
            f"# table of interventions: {summary['table']}"
            + (f" (version {summary['version']})" if summary["version"] else " (no version: file of a test)"),
            f"# rule of the marks: {rule}",
            "# form of the template: unit, tab, text; a source in square brackets (the number of a reference, a "
            "citation, or keys as the table has them); ^{...} superscript, _{...} subscript",
            f"# a day is printed with the format of the token {tokens['date']}, a period with the format of "
            f"the token {tokens['range of dates']}, both of the unit {UNIT_OF_THE_LEGEND}",
            "# where the rule of the marks is 'numbered marks' the sentence holds, for each mark, its number, the "
            "words for the legend, the day or the days and the sources, and not the short label",
            "# to be added to the legend by the track 'Article', after the sentence on panel b; the keys of the "
            "sources are replaced there by the numbers of the references"]
        if summary["test"]:
            head_lines.insert(0, f"# {MARKER}: the interventions are INVENTED; this sentence never enters "
                                 "the article")
        if summary["preview"]:
            head_lines.insert(0, f"# {PREVIEW}: the events are those of the table as given; the keys of the sources "
                                 "are replaced by the numbers of the references when the article is written")
        legend_text = "\n".join(head_lines) + f"\n[{UNIT_OF_THE_LEGEND}]\t{sentence}\n"
        hits = C.words_found(sentence, patterns)
        out["sentence_for_the_legend"] = {
            "file": names["legend"], "sentence": sentence, "characters": len(sentence),
            "equal_to_the_sentence_built_from_the_table_read_again": sentence == again,
            "rule_of_the_marks": rule,
            "every_label_words_and_source_of_the_table_in_it": all(
                (rule == "numbered marks" or e["label"].replace("<<BR>>", " ") in sentence)
                and e["words"] in sentence and f"[{e['source']}]" in sentence for e in expected),
            "parts_as_the_rule_asks": all(
                (f"{e['number']}, {e['words']}, " if rule == "numbered marks"
                 else f"{e['label'].replace('<<BR>>', ' ')}: {e['words']}, ") in sentence
                for e in sorted(expected, key=lambda e: e["number"])),
            "short_labels_that_stand_in_the_sentence": [
                e["label"].replace("<<BR>>", " ") for e in expected
                if e["label"].replace("<<BR>>", " ") + ": " in sentence],
            "kinds_and_forms_that_the_explanation_names": explanation_against_the_table(sentence, expected),
            "interventions_in_it": sentence.count("]"), "interventions_drawn": len(expected),
            "words_of_R4j": [(h["found"], h["kind"]) for h in hits]}
        s = out["sentence_for_the_legend"]
        s["ok"] = (s["equal_to_the_sentence_built_from_the_table_read_again"]
                   and s["every_label_words_and_source_of_the_table_in_it"]
                   and s["parts_as_the_rule_asks"]
                   and (rule != "numbered marks" or not s["short_labels_that_stand_in_the_sentence"])
                   and s["kinds_and_forms_that_the_explanation_names"]["ok"]
                   and s["interventions_in_it"] == s["interventions_drawn"] and not hits)
    else:
        out["sentence_for_the_legend"] = {"file": None, "note": "no intervention is drawn: no sentence is written",
                                          "ok": True}
    out["ok"] = all(v["ok"] for v in out.values() if isinstance(v, dict) and "ok" in v) \
        and not out["words_of_R4j_in_the_cells"]
    return out, legend_text


def colour_of_the_strip_checked(I, fig, info):
    """The colour of the strip is no colour of an artist of panels a, b and c and no colour that the
    conventions name; every artist of the strip has this colour, but the baseline and the declaration."""
    strip = mpl.colors.to_hex(STRIP_COLOUR)
    panels = set()
    for ax in fig.axes:
        if ax is info["axes"]:
            continue
        for a in ax.findobj():
            if not a.get_visible():
                continue
            gets = []
            if isinstance(a, Line2D):                # a line counts if it is drawn, a marker if there is one
                if a.get_linestyle() not in ("None", "", " ") and a.get_linewidth() > 0 \
                        and len(a.get_xdata()) > 1:
                    gets.append("get_color")
                if a.get_marker() not in ("None", "", " ", None) and len(a.get_xdata()):
                    gets += ["get_markerfacecolor", "get_markeredgecolor"]
            elif isinstance(a, mpl.text.Text):
                if a.get_text().strip():
                    gets.append("get_color")
            elif isinstance(a, (mpl.patches.Patch, mpl.collections.Collection)):
                gets += ["get_facecolor", "get_edgecolor"]
            for get in gets:
                try:
                    c = getattr(a, get)()
                    for one in (c if isinstance(c, (list, np.ndarray)) and np.ndim(c) == 2 else [c]):
                        if not isinstance(one, str) or one not in ("none", "None", "auto"):
                            if mpl.colors.to_rgba(one)[3] > 0:
                                panels.add(mpl.colors.to_hex(one))
                except (ValueError, TypeError):
                    pass
    named = sorted({c.lower() for c in re.findall(r"#[0-9a-fA-F]{6}", I.text(C.CONVENTIONS))})
    rgb = np.array(mpl.colors.to_rgb(strip)) * 255
    near = sorted((float(np.abs(np.array(mpl.colors.to_rgb(c)) * 255 - rgb).max()), c)
                  for c in sorted(panels | set(named)))
    own = set()
    for m in info["strip"]["marks"]:
        own |= {mpl.colors.to_hex(a.get_color()) for _, a in m["artists"]}
    own |= {mpl.colors.to_hex(st.get_color()) for _, st in info["strip"]["stems"]}
    own |= {mpl.colors.to_hex(t["artist"].get_color()) for t in info["strip"]["texts"]}
    own.add(mpl.colors.to_hex(info["axes"].spines["bottom"].get_edgecolor()))
    return {"colour": strip, "a_proposal": "yes: the user may name another colour",
            "colours_of_the_artists_of_the_panels": sorted(panels),
            "colours_named_in_the_conventions": named,
            "nearest_of_these": {"colour": near[0][1], "largest_difference_of_a_channel_of_255": near[0][0]},
            "colours_of_the_marks_lines_and_texts_of_the_strip": sorted(own),
            "ok": strip not in panels and strip not in named and own == {strip}
            and near[0][0] >= COLOUR_DISTANCE_MIN}


def plt_close(fig):
    import matplotlib.pyplot as plt
    plt.close(fig)


def main_with_interventions(I, table, outdir, preview=False, sources=None):
    V, info = build_values(I)
    n_head = len(V.rows)
    summary = add_interventions(V, I, table, info, preview)
    names = names_of_the_outputs(table, preview)
    values_path = V.save(os.path.join(outdir, names["values"]))      # the file of values is written first
    fig, extra = draw(values_path)                                    # ... and the figure is drawn from it
    patterns = C.word_patterns(I)
    mode = "test" if table["test"] else ("preview" if preview else "")
    checks = C.check_figure(fig, patterns, banners=BANNERS.get(mode, ()))
    saved = C.save_figure(fig, names["figure"], outdir)
    new, legend_text = checks_of_the_interventions(I, table, fig, values_path, saved["png"], n_head, summary,
                                                   patterns, names, sources)
    rb = C.read_back(values_path, I)
    rows = C.read_csv(values_path)
    legend_path = os.path.join(outdir, names["legend"])
    if legend_text is not None:
        with open(legend_path, "w", encoding="utf-8", newline="\n") as fh:
            fh.write(legend_text)
        new["sentence_for_the_legend"]["sha256"] = C.sha256(legend_path)
    elif os.path.exists(legend_path):
        os.remove(legend_path)
    by = {}
    for r in rows:
        by[f"{r['panel']} | {r['element']}"] = by.get(f"{r['panel']} | {r['element']}", 0) + 1
    marked = [f"{r['panel']} | {r['element']} | {r['label']}" for r in rows
              if C.marked_by_R4f(r["meets_convergence_criterion"])]
    saved = {k: (os.path.basename(v) if isinstance(v, str) else v) for k, v in saved.items()}
    out = {"program": os.path.basename(__file__), "figure": saved, "values_file": os.path.basename(values_path),
           "figure_sha256": {"png": C.sha256(os.path.join(outdir, saved["png"])),
                             "pdf": C.sha256(os.path.join(outdir, saved["pdf"]))},
           "values_sha256": C.sha256(values_path), "rows": len(rows), "rows_of_4_october": n_head,
           "rows_by_panel_and_element": by,
           "read_back": rb, "rows_not_yet_delivered": sum(r["delivered"] == "no" for r in rows),
           "rows_marked_by_R4f": marked, "rows_not_drawn": sum(r["drawn_as"] == "not drawn" for r in rows),
           "figure_checks": checks, "drawing": extra, "info": info, "interventions": new,
           "colours_and_sizes_against_the_conventions": C.colours_against_the_conventions(I),
           "table_of_inputs_version": I.version_of_the_table_of_inputs}
    out["values_marked_by_R4f_with_the_cell_as_it_reads"] = marked_with_their_cells(rows)
    out["not_yet_delivered"] = not_yet_delivered(rows, checks)
    ok = bool(checks["ok"] and not rb["not_equal"] and not extra["keys_or_words_over_data"] and new["ok"]
              and out["not_yet_delivered"]["ok"])
    out["all_checks_hold"] = ok
    C.json_dump(out, os.path.join(outdir, names["checks"]))
    print(json.dumps({"png": saved["png"], "rows": len(rows), "checks_ok": ok,
                      "rule": new["rule_of_the_marks"]["rule"],
                      "min_fontsize": checks["min_fontsize"], "read_back_not_equal": len(rb["not_equal"]),
                      "over_data": extra["keys_or_words_over_data"], "overlaps": checks["text_overlaps"],
                      "outside": checks["outside_figure"], "on_spine": checks["text_on_spine"],
                      "words": checks["words_of_R4j"],
                      "checks_of_the_interventions_that_fail": sorted(
                          k for k, v in new.items() if isinstance(v, dict) and v.get("ok") is False)}))
    plt_close(fig)
    return 0 if ok else 1


def main(argv):
    args = list(argv)
    given = {}
    for flag in ("--interventions", "--sources"):
        if flag in args:
            k = args.index(flag)
            if k + 1 >= len(args) or args[k + 1].startswith("--"):
                raise C.Stop(f"{flag} is followed by the name of the table (an entry of the JSON of the inputs)")
            given[flag] = args[k + 1]
            del args[k:k + 2]
            if os.path.basename(given[flag]) != given[flag]:
                raise C.Stop(f"{flag} takes the NAME of the table, not a path; the path stands in the JSON of "
                             "the inputs")
    preview = "--preview" in args
    if preview:
        args.remove("--preview")
    name, sources = given.get("--interventions"), given.get("--sources")
    if len(args) >= 3 and args[1] == "--draw-only":
        if name is not None or sources is not None or preview:
            raise C.Stop("--draw-only draws from a file of values alone and takes no table")
        rows = C.read_csv(args[2])
        if not any(r["element"] in (E_INT_MARK, E_BANNER) for r in rows):
            return main_of_4_october(args)
        outdir = args[3] if len(args) > 3 else "."
        os.makedirs(outdir, exist_ok=True)
        stem = stem_of_a_file_of_values(args[2])
        mode = mode_of_the_banners(rows)
        if (mode == "test") != stem.startswith(MARKER + "_") or (mode == "preview") != stem.startswith(PREVIEW + "_"):
            raise C.Stop(f"--draw-only: a file of values of a test with fixtures, and no other, has a name that "
                         f"begins with {MARKER}_; a file of values of a preview, and no other, one that begins "
                         f"with {PREVIEW}_")
        fig, extra = draw(args[2])
        info = C.save_figure(fig, stem, outdir)
        print(json.dumps({"drawn_from": os.path.basename(args[2]), "png": os.path.basename(info["png"]),
                          "rule": (extra.get("strip") or {}).get("rule"),
                          "keys_or_words_over_data": extra["keys_or_words_over_data"]}))
        plt_close(fig)
        return 0
    if name is None:
        if sources is not None or preview:
            raise C.Stop("--sources and --preview belong to a run with --interventions")
        return main_of_4_october(args)
    if len(args) < 2:
        raise C.Stop("usage: python fig_ne_R_cases_20261006.py INPUTS.json [OUTDIR] [--interventions NAME.csv "
                     "[--sources NAME.csv] [--preview]]")
    outdir = args[2] if len(args) > 2 else "."
    os.makedirs(outdir, exist_ok=True)
    I = C.Inputs(args[1])
    table = table_of_interventions(I, name)
    if preview and table["test"]:
        raise C.Stop("--preview: a test with fixtures is no preview")
    if sources is not None:
        if sources == name:
            raise C.Stop("--sources names the table of interventions")
        register_sources(I, sources)
    if not read_interventions(I, name, table["test"])[0]:
        # a table with no row: the run without a table. The table of the article gives its files the names of
        # the figure of the article (no strip, no sentence); a test or a preview writes the files of the
        # figure without the strip, as a run without a table does.
        print(json.dumps({"interventions": f"the table {name} holds no row: the run is the run without a table"}))
        if table["test"] or preview:
            return main_of_4_october(args)
        names = dict(ARTICLE_NAMES)
        legend = os.path.join(outdir, names.pop("legend"))
        if os.path.exists(legend):
            os.remove(legend)
        return main_of_4_october(args, names)
    return main_with_interventions(I, table, outdir, preview, sources)


if __name__ == "__main__":
    sys.exit(main(sys.argv))
