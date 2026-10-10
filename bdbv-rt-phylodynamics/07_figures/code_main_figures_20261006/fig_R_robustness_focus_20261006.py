"""fig_R_robustness_focus_20261006.py
Figure 2 of the article (reproduction number of the three periods W1, W2, W3: one row per analysis, three
columns), drawn from the tables of the round of October 2026. The picture holds no banner.

    python fig_R_robustness_focus_20261006.py INPUTS.json [OUTDIR]          builds the file of values, then draws
    python fig_R_robustness_focus_20261006.py --draw-only VALUES.csv OUTDIR draws from a file of values alone

Writes, in this order:
    fig_R_robustness_focus_20261006_values.csv   (BEFORE anything is drawn)
    fig_R_robustness_focus_20261006.png/.pdf     (drawn from the file of values and from nothing else)
    fig_R_robustness_focus_20261006_checks.json

ROWS: the primary analysis; the analyses on the other genome sets (new group, directly below the primary
analysis); the 13 sensitivity analyses in the order of the schedule; the two analyses with other numbers of
cells that the specification names as 'shown in Figure 2'; the other generation times. The rows are read from
the schedules, the specification and the tables; an analysis whose table of estimates is not among the inputs
keeps its row: no such row is expected, and the checks report every such row as an error.
The rules of drawing are those of fig_R_robustness_focus_20260925.py (the earlier generator, which is NOT run;
its lists of labels are READ from its source), with the changes of this round; every rule is listed in
INTERNAL_main_figures_mapping_of_inputs_20261006.csv. No value is computed (R1).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import main_fig_common_20261006 as C          # sets the thread variables before numpy is imported

import ast
import json
import math
import re

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch, Rectangle
from matplotlib.legend_handler import HandlerTuple

from fig_ne_R_cases_20261006 import primary_analysis, token_text, quantity_ids

NAME = "fig_R_robustness_focus" + C.SUFFIX
VALUES_FILE = NAME + "_values.csv"
FIGURE_STEM = NAME
CHECKS_FILE = NAME + "_checks.json"

E_HEAD_NAME = "name of the period at the head of a column"
E_HEAD_DATES = "dates of the period at the head of a column"
E_HEADING = "heading of a group of rows"
E_LABEL = "label of a row"
E_R = "R of the period, median and 95 % HPD"
E_R_COMMON = "R of the last period with the common end"
E_END = "end of the last period, printed at the right of the row"
E_NOT_YET = "words printed in place of the estimates of a row"
E_KEY = "entry of the key"
E_AXIS_RULE = "limits of an axis by rule"
E_TOP_LEFT = "words above the labels of the rows"

PERIODS_DRAWN = ("P1", "P2", "P3")
# layout that the earlier generator holds as typed numbers (inches); listed in the mapping
TOP_IN, BOTTOM_IN = 0.62, 0.92
LEFT_IN, LABEL_IN, GAP_IN, RIGHT_IN = 0.10, 2.76, 0.24, 0.14
COL_IN = (C.FIG_WIDTH - LEFT_IN - LABEL_IN - RIGHT_IN - 2 * GAP_IN) / 3
DATE_PAD_IN = 0.04
END_OFFSET_PT = 3.0
STEP = 0.5                 # R4 (h): limits in steps of 0.5
MARGIN = 0.1               # beyond the first and the last tick
HEADING_NEW_GROUP = "Other thresholds of the called fraction"      # PROPOSAL (group 4 of the task)


# =============================================================================================
# 1. the rows and their labels
# =============================================================================================
def earlier_lists(I):
    """Lists of labels and headings of the earlier generator, READ from its source (not executed)."""
    src = I.text(C.EARLIER_FIG2)
    tree = ast.parse(src)
    out, headings = {}, []
    for node in tree.body:
        if isinstance(node, ast.Assign) and len(node.targets) == 1 and isinstance(node.targets[0], ast.Name):
            n = node.targets[0].id
            if n in ("FOCUS_LABELS", "GRID_SETS", "LAST_DATE_SETS", "SMOOTH_SETS", "CELL_ROWS"):
                out[n] = (ast.literal_eval(node.value), node.lineno, node.end_lineno)
        if isinstance(node, ast.Expr) and isinstance(node.value, ast.Call):
            f = node.value.func
            if isinstance(f, ast.Attribute) and f.attr == "heading" and getattr(f.value, "id", "") == "L" \
                    and node.value.args and isinstance(node.value.args[0], ast.Constant):
                headings.append((node.value.args[0].value, node.lineno))
    for n in ("FOCUS_LABELS", "GRID_SETS", "LAST_DATE_SETS", "SMOOTH_SETS", "CELL_ROWS"):
        if n not in out:
            raise C.Stop(f"{C.EARLIER_FIG2}: the list {n} was not found")
    if len(headings) != 5:
        raise C.Stop(f"{C.EARLIER_FIG2}: {len(headings)} headings found, 5 expected")
    out["HEADINGS"] = headings
    return out


def plan_thresholds(I):
    """Genome set -> words of the plan for its called fraction (table of section 2 of the plan)."""
    out = {}
    for i, line in enumerate(I.lines(C.PLAN)):
        m = re.match(r"^\| (C\d\d) \| ([^|]+?) \|", line)
        if m:
            out[m.group(1)] = (m.group(2).strip(), i)
    if not out:
        raise C.Stop(f"{C.PLAN}: the table of the genome sets was not found")
    return out


def rows_of_the_figure(I):
    """The rows of the figure, in order: list of dicts (group, heading, analysis, kind, label rule, ...)."""
    E = earlier_lists(I)
    H = [h for h, _ in E["HEADINGS"]]
    primary = primary_analysis(I)
    s3 = I.table(C.SCHEDULE_PHASE3).need("analysis", "program", "setting", "genome_set", "track")
    s5 = I.table(C.SCHEDULE_PHASE5).need("group", "analysis", "what_is_varied", "description",
                                         "counterpart_in_the_round_of_24_September", "chains_of_the_analysis")
    rows = [dict(group=1, heading="", analysis=primary, kind="primary")]
    # group 4: the analyses on the other genome sets (same program and setting)
    p = [r for r in s3.rows if r["analysis"] == primary]
    if not p:
        raise C.Stop(f"{C.SCHEDULE_PHASE3}: the primary analysis {primary} is not in the schedule")
    thr = plan_thresholds(I)
    seen = []
    for i, r in enumerate(s3.rows):
        if (r["program"], r["setting"], r["track"]) == (p[0]["program"], p[0]["setting"], p[0]["track"]) \
                and r["analysis"] != primary and r["analysis"] not in seen and r["genome_set"] in thr:
            seen.append(r["analysis"])
            rows.append(dict(group=4, heading=HEADING_NEW_GROUP, analysis=r["analysis"], kind="other set",
                             genome_set=r["genome_set"], schedule=C.SCHEDULE_PHASE3, schedule_row=i))
    # group 2: the 13 sensitivity analyses, in the order of the schedule; checked against section 2 of the specification
    group2 = [g for g in dict.fromkeys(r["group"] for r in s5.rows) if g.startswith("sensitivity")]
    if len(group2) != 1:
        raise C.Stop(f"{C.SCHEDULE_PHASE5}: {len(group2)} groups 'sensitivity ...'")
    order = list(dict.fromkeys(r["analysis"] for r in s5.rows if r["group"] == group2[0]))
    spec_rows = re.findall(r"^\| (\d+) \| (\w+) \| ([^|]+) \| (\w+) \|", I.text(C.PHASE5_SPEC), flags=re.M)
    spec_order = ["delphy_sens_" + a for _, a, _, _ in spec_rows]
    if order != spec_order:
        raise C.Stop("the order of the sensitivity analyses in the schedule differs from section 2 of the specification")
    for a in order:
        i = min(k for k, r in enumerate(s5.rows) if r["analysis"] == a)
        r = s5.rows[i]
        cp = r["counterpart_in_the_round_of_24_September"]
        if cp != dict((("delphy_sens_" + x[1]), x[3]) for x in spec_rows)[a]:
            raise C.Stop(f"counterpart of {a}: schedule {cp!r} differs from the specification")
        if cp in E["FOCUS_LABELS"][0]:
            heading, kind = H[0], "genome set"
        elif cp in E["LAST_DATE_SETS"][0]:
            heading, kind = H[1], "last date"
        elif cp in [g for g, _ in E["GRID_SETS"][0]]:
            heading, kind = H[2], "grid"
        elif cp in E["SMOOTH_SETS"][0]:
            heading, kind = H[2], "smoothing"
        else:
            raise C.Stop(f"{a}: the counterpart {cp!r} has no label in the earlier generator")
        rows.append(dict(group=2, heading=heading, analysis=a, kind=kind, counterpart=cp,
                         description=r["description"], schedule=C.SCHEDULE_PHASE5, schedule_row=i))
    # group 3: the numbers of cells that the specification names as 'shown in Figure 2'
    m = re.search(r"^\| ([^|]*cells[^|]*) \| [^|]* \| shown in Figure 2[^|]*\|", I.text(C.PHASE5_SPEC), flags=re.M)
    if not m:
        raise C.Stop(f"{C.PHASE5_SPEC}: the row 'shown in Figure 2' of section 3 was not found")
    cells = [int(x.replace(",", "")) for x in re.findall(r"([\d,]+) cells", m.group(1) + " cells")]
    cells = list(dict.fromkeys(cells))
    if cells != list(E["CELL_ROWS"][0]):
        raise C.Stop(f"numbers of cells: specification {cells}, earlier generator {E['CELL_ROWS'][0]}")
    for n in cells:
        hit = [k for k, r in enumerate(s5.rows)
               if r["counterpart_in_the_round_of_24_September"] == f"primary | {n} cells"]
        names = list(dict.fromkeys(s5.rows[k]["analysis"] for k in hit))
        if len(names) != 1:
            raise C.Stop(f"{C.SCHEDULE_PHASE5}: {len(names)} analyses with the counterpart 'primary | {n} cells'")
        rows.append(dict(group=3, heading=H[3], analysis=names[0], kind="cells", cells=n,
                         counterpart=f"primary | {n} cells", description=s5.rows[hit[0]]["description"],
                         schedule=C.SCHEDULE_PHASE5, schedule_row=hit[0]))
    # group 5: the other generation times (table without chains, posterior of the primary analysis)
    gt = I.table(C.T_GT).need("analysis", "quantity", "period_start", "period_end", "end_used", "median",
                              "hpd95_lower", "hpd95_upper", "generation_time", "generation_time_mean_days",
                              "generation_time_sd_days", "generation_time_sd_days_as_used",
                              "meets_convergence_criterion")
    for g in dict.fromkeys(r["generation_time"] for r in gt.rows):
        if "reference" in g:
            continue
        rows.append(dict(group=5, heading=H[4], analysis=primary, kind="generation time", generation_time=g))
    return rows, E, thr, primary


def label_of(I, row, est_name, E, thr):
    """(text of the label, value cells for the file of values, how the label was made)."""
    k = row["kind"]
    delivered = est_name is not None
    cell = dict(table="", row="", columns="", value="", upper="")
    n = ""
    if delivered and k != "generation time":
        t = I.table(est_name).need("analysis", "genomes", "most_recent_collection_date", "setting")
        i = C.find_row(t, "R of P1", row["analysis"], "reproduction number", period="P1", end_used="fixed date")
        vals = {r["genomes"] for r in t.rows if r["analysis"] == row["analysis"]}
        if len(vals) != 1:
            raise C.Stop(f"{est_name}: the rows of {row['analysis']} hold {len(vals)} numbers of genomes")
        n = t.rows[i]["genomes"]
        cell = dict(table=est_name, row=i, columns="value=genomes", value=n, upper="")
        if "description" in row and row["description"] not in t.rows[i]["setting"]:
            raise C.Stop(f"{est_name}: the setting of {row['analysis']} does not hold the description of the schedule")
    br = f" ({n})" if n else ""
    if not delivered and k not in ("generation time", "primary", "other set"):
        cell = dict(table=row["schedule"], row=row["schedule_row"], columns="value=description",
                    value=row["description"], upper="")
    if k == "primary":
        return f"Primary analysis{br}", cell, "carried over (earlier generator, line 144)"
    if k == "other set":
        words, line = thr[row["genome_set"]]
        return f"{words}{br}", cell, (f"PROPOSED: words of the plan for the called fraction of the set "
                                      f"({C.PLAN}, line {line}, 0-based); no counterpart in the figure as it stands")
    if k == "genome set":
        return (f"{E['FOCUS_LABELS'][0][row['counterpart']]}{br}", cell,
                f"carried over: words of the label of the counterpart {row['counterpart']} (earlier generator, "
                f"lines {E['FOCUS_LABELS'][1]}-{E['FOCUS_LABELS'][2]})")
    if k == "last date":
        if delivered:
            t = I.table(est_name)
            i = cell["row"]
            day = C.to_day(t.rows[i]["most_recent_collection_date"])
            cell = dict(cell, columns="value=genomes; upper=most_recent_collection_date",
                        upper=t.rows[i]["most_recent_collection_date"])
            how = "carried over (earlier generator, line 153): most recent collection date and genomes of the table of estimates"
        else:
            m = re.search(r"up to (\d+ \w+ \d{4})", row["description"])
            if not m:
                raise C.Stop(f"{row['analysis']}: no date in the description of the schedule")
            day = C.day_of_words(m.group(1))
            how = ("carried over in its words (earlier generator, line 153); the date is that of the description of "
                   "the schedule of phase 5, because the analysis is not delivered (the earlier generator printed the "
                   "most recent collection date of the set); no number of genomes")
        return f"genomes up to {C.fmt_date(day)}{br}", cell, how
    if k == "grid":
        m = re.match(r"(\d+) parameters", row["description"])
        if not m:
            raise C.Stop(f"{row['analysis']}: no number of parameters in the description of the schedule")
        lab = dict(E["GRID_SETS"][0])[row["counterpart"]]
        if delivered:
            cell = dict(table=row["schedule"], row=row["schedule_row"], columns="value=description",
                        value=row["description"], upper="")
        return (f"{lab} ({int(m.group(1))} parameters)", cell,
                f"carried over: words of the counterpart {row['counterpart']} (earlier generator, lines "
                f"{E['GRID_SETS'][1]}-{E['GRID_SETS'][2]} and 158); the number of parameters is read from the "
                "description of the schedule of phase 5 (the earlier generator read the column grid of its table)")
    if k == "smoothing":
        m = re.search(r"double-half time of (\d+(?:\.\d+)?) days", row["description"])
        if not m:
            raise C.Stop(f"{row['analysis']}: no time of the smoothing in the description of the schedule")
        if delivered:
            cell = dict(table=row["schedule"], row=row["schedule_row"], columns="value=description",
                        value=row["description"], upper="")
        return (f"smoothing {float(m.group(1)):g} d", cell,
                "carried over (earlier generator, line 161); the days are read from the description of the schedule "
                "of phase 5 (the earlier generator read the column double_half_time_d of its table)")
    if k == "cells":
        if delivered:
            cell = dict(table=row["schedule"], row=row["schedule_row"], columns="value=description",
                        value=row["description"], upper="")
        return (f"{row['cells']:,} cells", cell,
                "carried over (earlier generator, line 166); the number is that of the row 'shown in Figure 2' of "
                "section 3 of the specification of phase 5 and of the description of the schedule")
    if k == "generation time":
        t = I.table(C.T_GT)
        i = min(j for j, r in enumerate(t.rows) if r["generation_time"] == row["generation_time"])
        r = t.rows[i]
        mean, sd = C.num(r["generation_time_mean_days"]), C.num(r["generation_time_sd_days"])
        if C.num_label(C.num(r["generation_time_sd_days_as_used"])) != C.num_label(sd):
            raise C.Stop(f"{C.T_GT}: row {i}: the two columns of the SD print differently")
        cell = dict(table=C.T_GT, row=i, columns="value=generation_time_mean_days; upper=generation_time_sd_days",
                    value=r["generation_time_mean_days"], upper=r["generation_time_sd_days"])
        return (f"mean {C.num_label(mean)} d, SD {C.num_label(sd)} d", cell,
                "carried over (earlier generator, line 188); mean and SD of the table of the generation times")
    raise C.Stop(f"no rule of the label for the kind {k}")


def text_width_in(text, size, weight="normal"):
    """Width of a text in inches, measured with the renderer on a figure of its own."""
    fig = plt.figure(figsize=(4, 1))
    t = fig.text(0, 0, text, fontsize=size, fontweight=weight)
    fig.canvas.draw()
    w = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
    plt.close(fig)
    return w


# =============================================================================================
# 2. the file of values
# =============================================================================================
def build_values(I):
    C.style()
    V = C.Values(NAME, I)
    info = {}
    rows, E, thr, primary = rows_of_the_figure(I)
    periods = I.table(C.PERIODS).need("period", "start", "end", "stored_name")
    level_hpd, row_level = token_text(I, "LEVEL.hpd95")
    T, COL_T = C.TABLE_OF_TOKENS, "printed_where_the_quantity_is_unchanged_by_design"
    name_of = {}
    for P in PERIODS_DRAWN:
        name_of[P] = periods.rows[periods.one(f"period {P}", period=P)]["stored_name"]
    f_prim = I.estimates_table_of(primary)
    if f_prim is None:
        raise C.Stop(f"the primary analysis {primary} has no delivered table of estimates")
    est_p = I.table(f_prim)
    prim_row = {}
    for P in PERIODS_DRAWN:
        rp = periods.one(f"period {P}", period=P)
        end_used = "own" if periods.rows[rp]["end"] == "END" else "fixed date"
        i = C.find_row(est_p, f"R of {P}", primary, "reproduction number", period=P, end_used=end_used)
        prim_row[P] = (i, est_p.rows[i], end_used)
        r = est_p.rows[i]
        V.add(name_of[P], E_HEAD_NAME, label=name_of[P], value=name_of[P], drawn_as="bold words above the column",
              table=C.PERIODS, row=rp, columns="value=stored_name", note=f"period {P}", delivered="yes")
        d0, d1 = C.to_day(r["period_start"]), C.to_day(r["period_end"])
        V.add(name_of[P], E_HEAD_DATES, label=f"{C.fmt_date(d0)} {C.NDASH} {C.fmt_date(d1)}", x=r["period_start"],
              x_end=r["period_end"], unit="date", drawn_as="words above the column", table=f_prim, row=i,
              columns="x=period_start; x_end=period_end",
              note=f"dates of the row of the primary analysis (end_used = {r['end_used']}); printed as the earlier "
                   "generator prints them", analysis=primary,
              quantity_id=quantity_ids(I, primary, "reproduction number", P, end_used, ["period_start", "period_end"],
                                       row=r),
              delivered="yes")
    V.add("rows", E_TOP_LEFT, label="analysis (number of genomes)", drawn_as="small grey words above the labels")
    end_of_primary = prim_row["P3"][1]["period_end"]

    drawn = {P: [] for P in PERIODS_DRAWN}       # bounds drawn in a column
    ends_printed = []                            # (label of the row, upper bound of the row, text)
    labels_seen, not_delivered, open_rows, marked_rows = set(), [], [], []
    last_heading = None
    for row in rows:
        k = row["kind"]
        if k == "generation time":
            est_name, delivered = C.T_GT, True
        else:
            est_name = I.estimates_table_of(row["analysis"])
            delivered = est_name is not None
        label, cell, how = label_of(I, row, est_name if k != "generation time" else None, E, thr)
        if row["heading"] and row["heading"] != last_heading:
            proposed = row["group"] == 4
            V.add("rows", E_HEADING, label=row["heading"], drawn_as="bold heading of a group",
                  note=("PROPOSED heading of the new group (position and heading are a proposal; the authors decide)"
                        if proposed else "carried over (heading of the earlier generator)"))
            last_heading = row["heading"]
        # the values of the three columns
        vals = {}
        if delivered:
            t = I.table(est_name).need("analysis", "quantity", "period_start", "period_end", "end_used", "median",
                                       "hpd95_lower", "hpd95_upper", "meets_convergence_criterion")
            for P in PERIODS_DRAWN:
                end_used = prim_row[P][2]
                if k == "generation time":
                    hits = [j for j, r in enumerate(t.rows) if r["generation_time"] == row["generation_time"]
                            and r["analysis"] == primary and r["end_used"] == end_used
                            and C.parse_quantity(r["quantity"]) == {"kind": "reproduction number", "period": P,
                                                                    "end_in_label": ""}]
                    if len(hits) != 1:
                        raise C.Stop(f"{C.T_GT}: {row['generation_time']}, {P}: expected 1 row, found {len(hits)}")
                    j = hits[0]
                else:
                    j = C.find_row(t, f"R of {P}", row["analysis"], "reproduction number", period=P, end_used=end_used)
                r = t.rows[j]
                for c in ("median", "hpd95_lower", "hpd95_upper"):
                    C.num(r[c], f"{est_name}, row {j}, {c}")
                if P != "P3" and (r["period_start"], r["period_end"]) != (prim_row[P][1]["period_start"],
                                                                         prim_row[P][1]["period_end"]):
                    raise C.Stop(f"{est_name}: row {j}: the dates of {P} differ from those of the primary analysis")
                vals[P] = (j, r)
        any_marked = any(C.marked_by_R4f(r["meets_convergence_criterion"]) for _, r in vals.values())
        shown = label + (" " + C.DOUBLE_DAGGER if any_marked else "")
        if shown in labels_seen:
            raise C.Stop(f"two rows of the figure have the label {shown!r}")
        labels_seen.add(shown)
        V.add("rows", E_LABEL, label=shown, value=cell["value"], upper=cell["upper"],
              drawn_as=("bold blue label, flush left (primary analysis)" if k == "primary" else "label of the row")
              + ("" if delivered else "; the estimates of the row are not delivered"),
              table=cell["table"], row=cell["row"], columns=cell["columns"], note=how,
              analysis=row["analysis"], delivered="yes" if delivered else "no")
        if not delivered:
            not_delivered.append(row["analysis"])
            V.add(name_of["P2"], E_NOT_YET, label=shown, value="", drawn_as=f"grey words '{C.NOT_YET}' in the row",
                  note="R4 (e): no value stands in for the analysis; the words are printed once, in the middle column",
                  analysis=row["analysis"], delivered="no")
            for P in PERIODS_DRAWN:
                V.add(name_of[P], E_R, label=shown, x=prim_row[P][1]["period_start"], unit="R",
                      drawn_as=f"not drawn: {C.NOT_YET}", note="the table of estimates of the analysis is not among the inputs",
                      analysis=row["analysis"], delivered="no")
            continue
        for P in PERIODS_DRAWN:
            j, r = vals[P]
            other_end = (P == "P3" and r["period_end"] != end_of_primary)
            qid = quantity_ids(I, row["analysis"], "reproduction number", P, r["end_used"],
                               ["median", "hpd95_lower", "hpd95_upper"], row=r) if k != "generation time" else ""
            mark = C.marked_by_R4f(r["meets_convergence_criterion"])
            V.add(name_of[P], E_R, label=shown, value=r["median"], lower=r["hpd95_lower"], upper=r["hpd95_upper"],
                  x=r["period_start"], x_end=r["period_end"], unit="R",
                  drawn_as=("open" if other_end else "filled") + " symbol, horizontal bar"
                  + ("; bold (primary analysis), band and dotted line over all rows" if k == "primary" else "")
                  + ("; double dagger behind the label of the row" if mark else ""),
                  table=est_name, row=j,
                  columns="value=median; lower=hpd95_lower; upper=hpd95_upper; x=period_start; x_end=period_end",
                  note=f"label of the row as found: {r['quantity']}; end_used = {r['end_used']}"
                       + (f"; generation time: {row['generation_time']}" if k == "generation time" else "")
                       + ("; the period ends on another day than in the primary analysis: open symbol" if other_end else ""),
                  analysis=row["analysis"], quantity_id=qid, delivered="yes",
                  meets=r["meets_convergence_criterion"])
            drawn[P] += [C.num(r["hpd95_lower"]), C.num(r["hpd95_upper"])]
            if mark:
                marked_rows.append(f"{shown} | {name_of[P]} | {r['meets_convergence_criterion']}")
            if other_end:
                text = f"to {C.fmt_date(C.to_day(r['period_end']))}"
                V.add(name_of[P], E_END, label=text, value=r["period_end"], x_end=r["period_end"], unit="date",
                      drawn_as="date at the right of the row", table=est_name, row=j,
                      columns="value=period_end; x_end=period_end",
                      note=f"row {shown!r}; the primary analysis ends on {end_of_primary}; printed as the earlier "
                           "generator prints it", analysis=row["analysis"], delivered="yes")
                ends_printed.append((shown, C.num(r["hpd95_upper"]), text))
                open_rows.append(f"{shown}: {text} ({r['period_end']})")
        # the same period with the common end: held, not drawn (R4 c)
        if k != "generation time":
            t = I.table(est_name)
            c = C.rows_found(t, row["analysis"], "reproduction number", period="P3", end_used="common")
            if len(c) == 1:
                r = t.rows[c[0]]
                V.add(name_of["P3"], E_R_COMMON, label=shown, value=r["median"], lower=r["hpd95_lower"],
                      upper=r["hpd95_upper"], x=r["period_start"], x_end=r["period_end"], unit="R",
                      drawn_as="not drawn", table=est_name, row=c[0],
                      columns="value=median; lower=hpd95_lower; upper=hpd95_upper; x=period_start; x_end=period_end",
                      note=f"R4 (c): row with end_used = common, held so that the other form can be drawn; label of "
                           f"the row as found: {r['quantity']}", analysis=row["analysis"],
                      quantity_id=quantity_ids(I, row["analysis"], "reproduction number", "P3", "common",
                                               ["median", "hpd95_lower", "hpd95_upper"], row=r),
                      delivered="yes", meets=r["meets_convergence_criterion"])
            elif len(c) == 0:
                V.add(name_of["P3"], E_R_COMMON, label=shown, unit="R", drawn_as="not drawn",
                      note="no source: the table of estimates of the analysis holds no row with end_used = common",
                      analysis=row["analysis"], delivered="yes")
            else:
                raise C.Stop(f"{est_name}: {len(c)} rows of P3 with the common end")
        else:
            V.add(name_of["P3"], E_R_COMMON, label=shown, unit="R", drawn_as="not drawn",
                  note="no source: the table of the generation times holds no row with end_used = common",
                  analysis=row["analysis"], delivered="yes")
    # ---- limits of the three axes (R4 h) ----------------------------------------------------------
    pm = {P: C.num(prim_row[P][1]["median"]) for P in PERIODS_DRAWN}
    ph = {P: C.num(prim_row[P][1]["hpd95_upper"]) for P in PERIODS_DRAWN}
    for P in PERIODS_DRAWN:
        lo, hi = C.steps_range(drawn[P] + [1.0], STEP)
        added = 0
        if P == "P3" and ends_printed:
            widest = max(text_width_in(t, C.SIZES[2]) for _, _, t in ends_printed)
            while True:
                a0 = 1.0 + (END_OFFSET_PT / 72.0 - widest - DATE_PAD_IN) / COL_IN
                span = (hi + MARGIN) - (lo - MARGIN)
                worst = max(max(u, 1.0, ph[P], pm[P]) for _, u, _ in ends_printed)
                if (worst - (lo - MARGIN)) / span < a0:
                    break
                hi = round(hi + STEP, 10)
                added += 1
                if added > 20:
                    raise C.Stop("no room for the printed ends of the last period")
        V.add(name_of[P], E_AXIS_RULE, label=f"axis of R, column {name_of[P]}", x=repr(round(lo - MARGIN, 10)),
              x_end=repr(round(hi + MARGIN, 10)), unit="R",
              drawn_as=f"limits of the axis; ticks in steps of {STEP:g} from {lo:g} to {hi:g}",
              columns="", note=f"limits by rule, no cell of a table. Rule: first and last tick = smallest range in "
              f"steps of {STEP:g} that holds every bound drawn in the column and the line R = 1; limits {MARGIN:g} "
              f"beyond the first and the last tick"
              + (f"; the last tick was moved up by {added} step(s) so that every printed end has room at the right of its row"
                 if P == "P3" else "") + f"; smallest and largest bound drawn: {min(drawn[P])!r}, {max(drawn[P])!r}")
    keys = ["primary analysis", f"other analysis, median and {level_hpd} HPD",
            f"{name_of['P3']} ends at the date printed", f"primary analysis: {level_hpd} HPD and median"]
    for i, kx in enumerate(keys):
        with_level = level_hpd in kx
        V.add("key", E_KEY, label=kx, value=level_hpd if with_level else "", drawn_as="entry of the key",
              table=T if with_level else "", row=row_level if with_level else "",
              columns=("value=" + COL_T) if with_level else "",
              quantity_id="LEVEL.hpd95" if with_level else "", delivered="yes" if with_level else "",
              note="unchanged by design" if with_level else "")
    info.update(primary_analysis=primary, rows_of_the_figure=len(rows), analyses_not_delivered=not_delivered,
                rows_with_an_open_symbol=open_rows, values_marked_by_R4f=marked_rows,
                end_of_the_primary_analysis=end_of_primary,
                rows_by_group={str(g): sum(1 for r in rows if r["group"] == g) for g in (1, 4, 2, 3, 5)})
    return V, info


# =============================================================================================
# 3. the drawing (from the file of values and from nothing else)
# =============================================================================================
def draw(values_path):
    rows = C.read_csv(values_path)
    font = C.style()
    heads = [r for r in rows if r["element"] == E_HEAD_NAME]
    cols = [r["panel"] for r in heads]
    dates = {r["panel"]: r["label"] for r in rows if r["element"] == E_HEAD_DATES}
    axis = {r["panel"]: r for r in rows if r["element"] == E_AXIS_RULE}
    L = C.RowLayout()
    DATA, flags = {}, {}
    for r in rows:
        if r["element"] == E_HEADING:
            L.heading(r["label"])
        elif r["element"] == E_LABEL:
            prim = r["drawn_as"].startswith("bold blue label")
            if prim:
                L.row(r["label"], r["label"], bold=True, flush=True, size=C.SIZES[1], label_colour=C.BLUE,
                      primary=True)
            else:
                L.row(r["label"], r["label"])
            flags[r["label"]] = dict(primary=prim, delivered=r["delivered"] == "yes")
            DATA[r["label"]] = {}
        elif r["element"] == E_R:
            DATA[r["label"]][r["panel"]] = r
    ends = {}
    for r in rows:
        if r["element"] == E_END:
            m = [k for k in DATA if DATA[k].get(r["panel"], {}).get("source_row_0based") == r["source_row_0based"]
                 and DATA[k][r["panel"]]["source_table"] == r["source_table"]]
            if len(m) != 1:
                raise C.Stop("file of values: a printed end has no single row")
            ends[m[0]] = r
    primary_key = [k for k, f in flags.items() if f["primary"]]
    if len(primary_key) != 1:
        raise C.Stop("file of values: no single primary analysis")
    PRIMARY = DATA[primary_key[0]]

    PITCH = C.pitch_for(L.units, TOP_IN, BOTTOM_IN, preferred=0.17)
    HEIGHT = TOP_IN + L.units * PITCH + BOTTOM_IN
    fig = C.new_figure(HEIGHT)
    axes = {}
    for k, w in enumerate(cols):
        x0 = LEFT_IN + LABEL_IN + k * (COL_IN + GAP_IN)
        ax = fig.add_axes([x0 / C.FIG_WIDTH, BOTTOM_IN / HEIGHT, COL_IN / C.FIG_WIDTH, L.units * PITCH / HEIGHT])
        axes[w] = ax
        ax.set_ylim(L.units - 0.45, -0.55)
        ax.set_yticks([])
        ax.spines["left"].set_visible(False)
        ax.set_xlim(float(axis[w]["x"]), float(axis[w]["x_end"]))
        ax.set_xticks(C.ticks_between(round(float(axis[w]["x"]) + MARGIN, 10),
                                      round(float(axis[w]["x_end"]) - MARGIN, 10), STEP))
        ax.set_xlabel("$R$", labelpad=2)
        ax.axvline(1.0, color="0.15", lw=0.7, zorder=1)
        p = PRIMARY[w]
        ax.axvspan(float(p["lower"]), float(p["upper"]), color=C.BLUE, alpha=0.10, lw=0, zorder=0.5)
        ax.axvline(float(p["value"]), color=C.BLUE, lw=0.7, ls=(0, (1, 1.5)), zorder=0.8)
        ax.text(0.5, 1.0, dates[w], transform=ax.transAxes, ha="center", va="bottom", fontsize=C.SIZES[1])
        ax.annotate(w, xy=(0.5, 1.0), xycoords="axes fraction", xytext=(0, 13), textcoords="offset points",
                    ha="center", va="bottom", fontsize=C.SIZES[0], fontweight="bold")
    TEXT = {}             # (row, column) -> text in whose place the thin line of the row is interrupted
    words = [r for r in rows if r["element"] == E_NOT_YET]
    for it in L.rows():
        d = DATA[it["key"]]
        big = bool(it.get("primary"))
        for w in cols:
            e = d[w]
            if e["drawn_as"].startswith("not drawn"):
                continue
            C.interval(axes[w], it["y"], float(e["value"]), float(e["lower"]), float(e["upper"]), C.BLUE,
                       ms=5.6 if big else 4.2, lw=1.8 if big else 1.1,
                       filled=not e["drawn_as"].startswith("open"))
            lo_, hi_ = axes[w].get_xlim()
            if float(e["lower"]) < lo_ or float(e["upper"]) > hi_:
                raise C.Stop(f"an interval of the row {it['key']!r} extends beyond the axis")
        if it["key"] in ends:
            e = ends[it["key"]]
            TEXT[(it["key"], e["panel"])] = axes[e["panel"]].annotate(
                e["label"], xy=(1.0, it["y"]), xycoords=("axes fraction", "data"), xytext=(END_OFFSET_PT, 0),
                textcoords="offset points", ha="right", va="center", fontsize=C.SIZES[2], color="0.1",
                annotation_clip=False)
    for r in words:
        it = [i for i in L.rows() if i["key"] == r["label"]][0]
        TEXT[(it["key"], r["panel"])] = axes[r["panel"]].annotate(
            C.NOT_YET, xy=(0.5, it["y"]), xycoords=("axes fraction", "data"), ha="center", va="center",
            fontsize=C.SIZES[2], color=C.GREY_NOT_YET, annotation_clip=False, zorder=6)
    # lines of the rows; under a text the line is interrupted, so that no line of a row crosses a text
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    no_room = []
    for w in cols:
        ax = axes[w]
        to_fraction = ax.transAxes.inverted()
        for it in L.rows():
            t = TEXT.get((it["key"], w))
            if t is None:
                ax.axhline(it["y"], color="0.88", lw=0.4, zorder=0.2)
                continue
            bb = t.get_window_extent(rend)
            pad = DATE_PAD_IN * fig.dpi
            a0 = float(to_fraction.transform((bb.x0 - pad, bb.y0))[0])
            a1 = float(to_fraction.transform((bb.x1 + pad, bb.y0))[0])
            if (it["key"]) in ends and ends[it["key"]]["panel"] == w:
                e = DATA[it["key"]][w]
                stops = [float(e["upper"]), 1.0, float(PRIMARY[w]["upper"]), float(PRIMARY[w]["value"])]
                f = [float(to_fraction.transform(ax.transData.transform((v, it["y"])))[0]) for v in stops]
                if not (0.0 < a0 and max(f) < a0):
                    no_room.append(it["key"])
            if (it["key"], w) in TEXT and it["key"] not in ends:
                # the words stand on a white field of the height of the row, so that no line crosses them
                ax.add_patch(Rectangle((a0, it["y"] - 0.5), a1 - a0, 1.0,
                                       transform=mpl.transforms.blended_transform_factory(ax.transAxes, ax.transData),
                                       facecolor="white", edgecolor="none", linewidth=0, zorder=5))
            if a0 > 0.0:
                ax.axhline(it["y"], xmin=0.0, xmax=min(a0, 1.0), color="0.88", lw=0.4, zorder=0.2)
            if a1 < 1.0:
                ax.axhline(it["y"], xmin=max(a1, 0.0), xmax=1.0, color="0.88", lw=0.4, zorder=0.2)
    if no_room:
        raise C.Stop(f"no room for the printed end in the rows {no_room}")
    C.draw_row_labels(fig, axes[cols[0]], L, LEFT_IN)
    top_left = [r["label"] for r in rows if r["element"] == E_TOP_LEFT][0]
    axes[cols[0]].annotate(top_left, xy=(0, 1.0), xycoords=("figure fraction", "axes fraction"),
                           xytext=(LEFT_IN * 72, 0), textcoords="offset points", ha="left", va="bottom",
                           fontsize=C.SIZES[2], color="0.3", annotation_clip=False)
    hk = [Line2D([], [], color=C.BLUE, lw=1.8, marker="o", ms=5.6),
          Line2D([], [], color=C.BLUE, lw=1.1, marker="o", ms=4.2),
          Line2D([], [], color=C.BLUE, lw=1.1, marker="o", ms=4.2, mfc="white", mew=1.0),
          (Patch(facecolor=mpl.colors.to_rgba(C.BLUE, 0.10), edgecolor="none"),
           Line2D([], [], color=C.BLUE, lw=0.7, ls=(0, (1, 1.5))))]
    lk = [r["label"] for r in rows if r["element"] == E_KEY]
    fig.legend(hk, lk, loc="lower left", bbox_to_anchor=(LEFT_IN / C.FIG_WIDTH, 0.06 / HEIGHT),
               fontsize=C.SIZES[2], handlelength=2.0, labelspacing=0.3, borderaxespad=0.0,
               handler_map={tuple: HandlerTuple(ndivide=1, pad=0)})
    # labels of the rows must end before the first column
    fig.canvas.draw()
    rend = fig.canvas.get_renderer()
    left_edge = axes[cols[0]].get_window_extent(rend).x0
    too_long = [t.get_text() for t in fig.findobj(mpl.text.Text)
                if t.get_text() in set(DATA) | {i["text"] for i in L.items if i["kind"] == "heading"}
                and t.get_window_extent(rend).x1 > left_edge - 2]
    extra = {"font": font, "height_in": HEIGHT, "pitch_in": PITCH, "lines": L.units,
             "rows": len(L.rows()), "headings": sum(1 for i in L.items if i["kind"] == "heading"),
             "labels_that_reach_the_first_column": too_long,
             "rows_with_the_words_not_yet_delivered": len(words), "rows_with_a_printed_end": len(ends)}
    return fig, extra


# =============================================================================================
# 4. run
# =============================================================================================
def main(argv):
    if len(argv) >= 3 and argv[1] == "--draw-only":
        outdir = argv[3] if len(argv) > 3 else "."
        os.makedirs(outdir, exist_ok=True)
        fig, extra = draw(argv[2])
        info = C.save_figure(fig, FIGURE_STEM, outdir)
        print(json.dumps({"drawn_from": argv[2], "png": info["png"]}))
        return 0
    if len(argv) < 2:
        raise C.Stop("usage: python fig_R_robustness_focus_20261006.py INPUTS.json [OUTDIR]")
    outdir = argv[2] if len(argv) > 2 else "."
    os.makedirs(outdir, exist_ok=True)
    I = C.Inputs(argv[1])
    V, info = build_values(I)
    values_path = V.save(os.path.join(outdir, VALUES_FILE))          # the file of values is written first
    fig, extra = draw(values_path)                                    # ... and the figure is drawn from it
    patterns = C.word_patterns(I)
    checks = C.check_figure(fig, patterns)
    saved = C.save_figure(fig, FIGURE_STEM, outdir)
    rb = C.read_back(values_path, I)
    rows = C.read_csv(values_path)
    by = {}
    for r in rows:
        by[f"{r['panel']} | {r['element']}"] = by.get(f"{r['panel']} | {r['element']}", 0) + 1
    saved = {k: (os.path.basename(v) if isinstance(v, str) else v) for k, v in saved.items()}
    out = {"program": os.path.basename(__file__), "figure": saved, "values_file": os.path.basename(values_path),
           "figure_sha256": {"png": C.sha256(os.path.join(outdir, saved["png"])),
                             "pdf": C.sha256(os.path.join(outdir, saved["pdf"]))},
           "values_sha256": C.sha256(values_path), "rows": len(rows), "rows_by_panel_and_element": by,
           "read_back": rb,
           "rows_not_yet_delivered": sum(r["delivered"] == "no" for r in rows),
           "not_yet_delivered": {"expected": 0, "analyses_not_delivered": info["analyses_not_delivered"],
                                 "rows_of_the_file_of_values_with_delivered_no": sum(r["delivered"] == "no" for r in rows),
                                 "texts_of_the_picture_that_hold_the_words": [
                                     t for t, _ in checks["texts"] if C.NOT_YET in t],
                                 "ok": not info["analyses_not_delivered"]
                                 and not any(r["delivered"] == "no" for r in rows)
                                 and not any(C.NOT_YET in t for t, _ in checks["texts"])},
           "rows_with_an_open_symbol_and_their_end": [
               {"label": r["label"], "analysis": r["analysis"], "end_of_the_last_period": r["x_end"],
                "end_of_the_primary_analysis": info["end_of_the_primary_analysis"],
                "source_table": r["source_table"], "source_row_0based": r["source_row_0based"],
                "printed_at_the_right_of_the_row": [e["label"] for e in rows if e["element"] == E_END
                                                    and e["source_table"] == r["source_table"]
                                                    and e["source_row_0based"] == r["source_row_0based"]]}
               for r in rows if r["element"] == E_R and r["drawn_as"].startswith("open")],
           "values_marked_by_R4f_with_the_cell_as_it_reads": [
               {"label": r["label"], "panel": r["panel"], "analysis": r["analysis"], "value": r["value"],
                "source_table": r["source_table"], "source_row_0based": r["source_row_0based"],
                "meets_convergence_criterion": r["meets_convergence_criterion"]}
               for r in rows if r["element"] == E_R and C.marked_by_R4f(r["meets_convergence_criterion"])],
           "rows_not_drawn": sum(r["drawn_as"].startswith("not drawn") for r in rows),
           "rows_marked_by_R4f": info["values_marked_by_R4f"], "figure_checks": checks, "drawing": extra,
           "info": info, "table_of_inputs_version": I.version_of_the_table_of_inputs}
    C.json_dump(out, os.path.join(outdir, CHECKS_FILE))
    ok = (checks["ok"] and not rb["not_equal"] and not extra["labels_that_reach_the_first_column"]
          and out["not_yet_delivered"]["ok"])
    print(json.dumps({"png": saved["png"], "rows": len(rows), "checks_ok": bool(ok),
                      "min_fontsize": checks["min_fontsize"], "read_back_not_equal": rb["not_equal"][:5],
                      "overlaps": checks["text_overlaps"], "outside": checks["outside_figure"],
                      "on_spine": checks["text_on_spine"], "words": checks["words_of_R4j"],
                      "too_long": extra["labels_that_reach_the_first_column"], "height_in": extra["height_in"]}))
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
