"""add_case_counts_to_figure_1_values.py - written for this repository; it is not one of the programs that ran in
the analysis.

The repository holds no table of case counts. The file of values of Figure 1 is therefore held without the rows of
panel c (results/fig_ne_R_cases_20261006_values_without_case_counts.csv). This program builds those rows from the
weekly table of confirmed cases, which 06_case_counts_and_response_events/regenerate_case_tables.py writes from the
public source, and writes the complete file of values (results/fig_ne_R_cases_20261006_values.csv), from which
draw_figure_1.py draws.

The rows are built by the rules of the program that wrote the file of values in the analysis
(code_main_figures_20261006/fig_ne_R_cases_20261006.py, function build_values, part 'c. confirmed cases'):
for every week one row for each of Ituri, Nord-Kivu and the other provinces; one row for the hatching of a week whose
split by province is pooled over several weeks; one row for the provisional week (the week with the fewest report
dates, if it is the last week); the entries of the key; and the limits of the axis. That program reads further
working files of the article and cannot be run from the repository alone, which is why its rules are repeated here.

The program compares the SHA-256 of the file it wrote with that of the file of values of the analysis, as it is
held in this repository (three cells of the column note were shortened when the repository was made; no value
changed), and prints whether they are the same. They are the same when the weekly table is that of the analysis
(source repository at commit 8eb57154).

usage, from the top folder of the repository:

    python 07_figures/add_case_counts_to_figure_1_values.py [WEEKLY.csv]

WEEKLY.csv: default results/weekly_confirmed_cases_by_province_20261002.csv
exit status: 0 if the file written is that of the analysis, 3 if it is another, 2 if the program stopped.
"""
import csv
import datetime as dt
import hashlib
import math
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
TOP = os.path.dirname(HERE)
FIGURE = "fig_ne_R_cases_20261006"
CASES = "weekly_confirmed_cases_by_province_20261002.csv"
WITHOUT = os.path.join(TOP, "results", FIGURE + "_values_without_case_counts.csv")
OUT = os.path.join(TOP, "results", FIGURE + "_values.csv")
CASES_README = os.path.join(TOP, "06_case_counts_and_response_events", "README_case_series_20261002.md")
SHA256_OF_THE_ANALYSIS = "a5b33e7a926757bd58fa899bf29f38a0a5737acc8923dd18b3ea239780df4160"

E_AXIS_END = "time axis, end"
E_CASES = "confirmed cases per week, by date of report"
E_HATCH = "hatching of a week whose split by province is pooled"
E_PROVISIONAL = "dashed outline of the provisional week"
E_KEY = "entry of a key"
E_AXIS_RULE = "limits of an axis by rule"
PROVINCES = [("ituri", "Ituri"), ("nord_kivu", "Nord-Kivu"), ("other_provinces", "other provinces")]
PROV_TEXT_SHARE = 705.0 / 900.0            # the word 'provisional' at 705 of an axis to 900
PROV_LINE_SHARE = 690.0 / 900.0
NDASH = "\u2013"
NEEDED = ("week_ending_sunday", "total_drc", "ituri", "nord_kivu", "other_provinces",
          "weeks_pooled_for_province_split", "report_dates_in_week", "complete_week")


class Stop(SystemExit):
    pass


def read(path):
    with open(path, encoding="utf-8", newline="") as f:
        rows = list(csv.reader(f))
    return rows[0], rows[1:]


def case_rows(columns, cases):
    """the rows of panel c, as dictionaries by the columns of the file of values"""
    out = []

    def add(element, **cells):
        row = dict.fromkeys(columns, "")
        row.update(figure=FIGURE, panel="c", element=element, delivered="yes")
        for k, v in cells.items():
            row[k] = v.isoformat() if isinstance(v, dt.date) else str(v)
        out.append(row)

    weeks = [dt.date.fromisoformat(r["week_ending_sunday"]) for r in cases]
    if any((b - a).days != 7 for a, b in zip(weeks[:-1], weeks[1:])) or any(w.weekday() != 6 for w in weeks):
        raise Stop("case series: the weeks are not consecutive weeks that end on a Sunday")
    last = len(weeks) - 1
    for i, r in enumerate(cases):
        s = sum(float(r[c]) for c, _ in PROVINCES)
        if abs(s - float(r["total_drc"])) > 0.11:
            raise Stop(f"{CASES}: row {i}: the provinces add up to {s}, the national count is {r['total_drc']}")
    counts = [int(r["report_dates_in_week"]) for r in cases]
    fewest = [i for i, n in enumerate(counts) if n == min(counts)]
    provisional = [i for i in fewest if i == last]
    if len(provisional) != 1:
        raise Stop(f"provisional week: the rule gives {len(provisional)} weeks "
                   f"(weeks with the fewest report dates: rows {fewest}; last week: row {last})")
    pv = provisional[0]
    m = re.search(r"The week ending (\d+ \w+) rests on (\d+) report dates[^.]*?it is to be marked as provisional in figures",
                  open(CASES_README, encoding="utf-8").read())
    if not m:
        raise Stop(f"{CASES_README}: the sentence on the provisional week was not found")
    readme_day, readme_n = m.group(1), m.group(2)
    complete = sorted({r["complete_week"] for r in cases})
    week = dt.timedelta(days=6)
    for i, r in enumerate(cases):
        pooled = int(r["weeks_pooled_for_province_split"])
        for col, lab in PROVINCES:
            add(E_CASES, label=f"{lab}, week ending {r['week_ending_sunday']}", value=r[col],
                x=weeks[i] - week, x_end=weeks[i], unit="cases",
                drawn_as="stacked bar" + (f"; hatched (province split pooled over {pooled} weeks)" if pooled > 1 else "")
                + ("; dashed outline, provisional" if i == pv else ""),
                source_table=CASES, source_row_0based=i, source_columns=f"value={col}; x_end=week_ending_sunday",
                note=f"x = first day of the week (Monday); national count of the week {r['total_drc']}; "
                     f"report dates in the week {r['report_dates_in_week']}")
        if pooled > 1:
            add(E_HATCH, label=f"week ending {r['week_ending_sunday']}", value=r["total_drc"],
                lower=r["weeks_pooled_for_province_split"], x=weeks[i] - week, x_end=weeks[i], unit="cases",
                drawn_as="white hatching over the whole bar (height: national count of the week)",
                source_table=CASES, source_row_0based=i,
                source_columns="value=total_drc; lower=weeks_pooled_for_province_split; x_end=week_ending_sunday",
                note="lower holds the number of weeks pooled (more than 1: hatched)")
    add(E_PROVISIONAL, label="provisional", value=cases[pv]["total_drc"], lower=cases[pv]["report_dates_in_week"],
        x=weeks[pv] - week, x_end=weeks[pv], unit="cases",
        drawn_as="dashed black outline of the bar; the word 'provisional' with a leader line",
        source_table=CASES, source_row_0based=pv,
        source_columns="value=total_drc; lower=report_dates_in_week; x_end=week_ending_sunday",
        note=f"R4 (b): rule of the earlier generator (the week with the fewest report dates, if it is the last week; "
             f"exactly one week): row {pv}; the README names the week ending {readme_day} with {readme_n} report dates; "
             f"lower holds the report dates of the week; the column complete_week reads "
             f"{', '.join(complete)} in all {len(cases)} weeks and does not decide")
    pw = sorted({int(r["weeks_pooled_for_province_split"]) for r in cases if int(r["weeks_pooled_for_province_split"]) > 1})
    for _, lab in PROVINCES:
        add(E_KEY, label=lab, drawn_as="entry of the key of panel c")
    if pw:
        i_lo = min(i for i, r in enumerate(cases) if int(r["weeks_pooled_for_province_split"]) == pw[0])
        i_hi = min(i for i, r in enumerate(cases) if int(r["weeks_pooled_for_province_split"]) == pw[-1])
        add(E_KEY, label=f"province split pooled, {pw[0]}{NDASH}{pw[-1]} weeks",
            value=cases[i_lo]["weeks_pooled_for_province_split"], drawn_as="entry of the key of panel c",
            source_table=CASES, source_row_0based=i_lo, source_columns="value=weeks_pooled_for_province_split",
            note=f"smallest and largest number of weeks pooled (above 1) in the column; the largest stands in row {i_hi}")
    cmax = max(float(r["total_drc"]) for r in cases)
    c_top = (math.floor(cmax / 100.0 + 1e-9) + 1) * 100.0
    add(E_AXIS_RULE, label="axis of the cases per week", x="0", x_end=repr(c_top), unit="cases",
        drawn_as="limits of the axis; ticks at the multiples of 200 below the limit",
        note=f"rule: from 0 to the smallest multiple of 100 above the largest weekly national count; the word "
             f"'provisional' stands at {PROV_TEXT_SHARE:.4f} x the limit, its leader line ends at {PROV_LINE_SHARE:.4f} x the limit")
    return out, weeks


def main(argv):
    if len(argv) > 2:
        raise Stop("usage: python 07_figures/add_case_counts_to_figure_1_values.py [WEEKLY.csv]")
    weekly = argv[1] if len(argv) == 2 else os.path.join(TOP, "results", CASES)
    if not os.path.exists(weekly):
        raise Stop(f"{weekly} does not exist: run 06_case_counts_and_response_events/regenerate_case_tables.py first")
    head, rows = read(weekly)
    missing = [c for c in NEEDED if c not in head]
    if missing:
        raise Stop(f"{weekly}: columns missing: {missing}")
    cases = [dict(zip(head, r)) for r in rows]
    columns, values = read(WITHOUT)
    el = columns.index("element")
    if any(r[columns.index("panel")] == "c" for r in values):
        raise Stop(f"{WITHOUT} holds rows of panel c")
    new, weeks = case_rows(columns, cases)
    ends = [dict(zip(columns, r)) for r in values if r[el] == E_AXIS_END]
    if len(ends) != 1 or ends[0]["value"] != cases[-1]["week_ending_sunday"] \
            or ends[0]["x"] != (weeks[-1] + dt.timedelta(days=1)).isoformat():
        raise Stop("the end of the time axis in the file of values is not the day after the last week of the weekly table")
    at = [i for i, r in enumerate(values) if r[el].startswith("intervention, ")]
    at = at[0] if at else len(values)
    block = [[r[c] for c in columns] for r in new]
    with open(OUT, "w", encoding="utf-8", newline="") as f:
        csv.writer(f, lineterminator="\n").writerows([columns] + values[:at] + block + values[at:])
    got = hashlib.sha256(open(OUT, "rb").read()).hexdigest()
    same = got == SHA256_OF_THE_ANALYSIS
    print("written:", os.path.relpath(OUT, TOP), "| rows of panel c added:", len(block))
    print("SHA-256:", got)
    print("the file of values of the analysis:", "the same" if same else "NOT the same (" + SHA256_OF_THE_ANALYSIS + ")")
    return 0 if same else 3


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv))
    except Stop as e:
        print("stopped:", e, file=sys.stderr)
        sys.exit(2)
