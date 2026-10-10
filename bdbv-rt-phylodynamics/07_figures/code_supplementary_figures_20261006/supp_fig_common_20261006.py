"""supp_fig_common_20261006.py

Common module of the programs that draw Figures S1 to S5 of the supplement from the tables of the round of
October 2026. It rests on main_fig_common_20261006.py, the common module of the main figures, which is used
UNCHANGED and is imported here as C.

What this module adds:
  Inputs        the table of inputs of the supplement (INTERNAL_supplementary_figures_inputs_20261006.csv)
                and, beside it, the table of inputs of the main figures. A VALUE that a figure draws is read
                from a file of the table of inputs of the supplement and from no other (check 'R1'); the
                files of the table of the main figures are read for rules, names, formats and the template.
  constants_of  typed constants of an earlier generator, read from its SOURCE (the generator is not run)
  Rows          rows and headings with labels of one or more lines (a label that is built from cells of a
                table can be longer than the column of the labels and is then set on two lines)
  read_back     every cell that a row of a file of values names is read again and compared; the rows of the
                one exception to R1 (Figure S3 c, a density that is computed) are counted apart
  standard      the checks that every figure has (sizes, overlaps, words of R4 (j), banner, R1, read-back,
                rows by panel and element), and the file of checks
"""
import ast
import csv
import datetime as dt
import gzip
import json
import math
import os
import re
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "code_main_figures_20261006"))
import main_fig_common_20261006 as C     # noqa: E402  (sets the backend and one thread, as on 4 October)

import matplotlib as mpl                 # noqa: E402
import matplotlib.pyplot as plt          # noqa: E402

SUFFIX = C.SUFFIX
TABLE_OF_INPUTS = "INTERNAL_supplementary_figures_inputs_20261006.csv"
LIST_OF_FIGURES = "article_figures_and_tables_20261002.csv"
TEMPLATE = "INTERNAL_article_template_20261003.txt"
COMPUTED = "computed"                     # first word of drawn_as of a row of the exception to R1 (S3 c)
NOT_DRAWN = "not drawn"                   # first words of drawn_as of a row that the figure does not draw
DAGGER = "\u2020"                         # mark of an analysis with chains continued from a saved state
DOUBLE_DAGGER = "\u2021"                  # mark of R4 (f)
LINE_PT = 9.6                             # height of a line of a label of 8 pt (line spacing 1.2)
IN_CELL = " ~ "                           # in source_columns: 'part=column ~ pattern': the text that the pattern finds in the cell
E_AXIS_RULE = "limits of an axis by rule"
E_KEY = "entry of a key"
E_TITLE = "title of a panel"
E_HEADING = "heading of a group of rows"
E_LABEL = "label of a row"
E_WORDS = "words printed in a panel"


class Inputs(C.Inputs):
    """The two JSON 'file name -> path' and the two tables of inputs. A file that both tables name must
    have the same version in both."""

    def __init__(self, json_of_the_supplement, json_of_the_third_task):
        C.Inputs.__init__(self, json_of_the_third_task)
        self.third = {"by_file": dict(self.by_file), "version": self.version_of_the_table_of_inputs,
                      "table": self.inputs, "paths": dict(self.paths)}
        if not os.path.exists(json_of_the_supplement):
            raise C.Stop(f"the JSON of the inputs of the supplement is missing: {json_of_the_supplement}")
        with open(json_of_the_supplement, encoding="utf-8") as fh:
            own = json.load(fh)
        if TABLE_OF_INPUTS not in own:
            raise C.Stop(f"the JSON of the supplement holds no entry {TABLE_OF_INPUTS}")
        self.paths = dict(self.third["paths"])
        self.paths.update(own)
        self.version_of_the_table_of_inputs = own.get(TABLE_OF_INPUTS + "#version", "")
        t = C.Table(TABLE_OF_INPUTS, self._path(TABLE_OF_INPUTS), self.version_of_the_table_of_inputs)
        t.need("file", "version", "used_for", "what_it_is", "analysis", "state_of_the_file", "note")
        self.supplement = t
        own_rows = {r["file"]: r for r in t.rows}
        if len(own_rows) != len(t.rows):
            raise C.Stop("the table of inputs of the supplement names a file twice")
        for name, r in own_rows.items():
            if name in self.third["by_file"] and self.third["by_file"][name]["version"] != r["version"]:
                raise C.Stop(f"{name}: the two tables of inputs name different versions")
            if name not in own:
                raise C.Stop(f"the JSON of the supplement holds no entry for {name}")
        self.by_file = dict(self.third["by_file"])
        self.by_file.update(own_rows)
        self._tables[TABLE_OF_INPUTS] = t
        # the rows of both tables, those of the supplement first (for the names of analyses)
        both = C.Table.__new__(C.Table)
        both.name, both.path, both.version = TABLE_OF_INPUTS, t.path, t.version
        both.columns = list(t.columns)
        both.rows = list(t.rows) + [r for r in self.third["table"].rows if r["file"] not in own_rows]
        self.inputs = both

    def version(self, name):
        if name == TABLE_OF_INPUTS:
            return self.version_of_the_table_of_inputs
        if name == C.TABLE_OF_INPUTS:
            return self.third["version"]
        return C.Inputs.version(self, name)

    def of_the_supplement(self, name):
        """Whether the table of inputs of the SUPPLEMENT names the file."""
        return any(r["file"] == name for r in self.supplement.rows)

    def table_of(self, name):
        """Name of the table of inputs that names the file (that of the supplement where both do)."""
        if self.of_the_supplement(name):
            return TABLE_OF_INPUTS
        if name in self.third["by_file"]:
            return C.TABLE_OF_INPUTS
        return ""

    def gz_columns(self, name, columns):
        """The cells of the named columns of a compressed table (a posterior file), as TEXT, row by row."""
        self.__dict__.setdefault("_gz", set()).add(name)
        with gzip.open(self._path(name), "rt", encoding="utf-8", newline="") as fh:
            rd = csv.reader(fh)
            head = next(rd)
            missing = [c for c in columns if c not in head]
            if missing:
                raise C.Stop(f"{name}: column(s) missing: {', '.join(missing)}")
            idx = [head.index(c) for c in columns]
            rows = [[rec[i] for i in idx] for rec in rd]
        return head, rows


def inputs_of(argv, first):
    """The two JSON of a command line: INPUTS_OF_THE_SUPPLEMENT.json INPUTS_OF_THE_THIRD_TASK.json."""
    if len(argv) < first + 2:
        raise C.Stop("two JSON of inputs are needed: that of the supplement and that of the third task")
    return Inputs(argv[first], argv[first + 1])


def quantity_ids(inputs, analysis, kind, end_used, columns, variant=None, row=None, kind_as_text=None):
    """quantity_id(s) of the LATEST specification whose key (definition for the program, op 'cell') names this
    analysis, this kind of quantity, the day(s) of the period of `row` (a row of a table with the cells
    period_start and period_end; None for a quantity without a period), this end and these columns; '' where
    the specification holds none (as in fig_ne_R_cases_20261006.py, function quantity_ids). kind_as_text: the
    kind as the key writes it, for a kind that the common module does not know."""
    t = inputs.table(C.SPECIFICATION).need("quantity_id", "definition_for_the_program")
    out = []
    for r in t.rows:
        d = r["definition_for_the_program"].strip()
        if not d.startswith("{"):
            continue
        j = json.loads(d)
        w = j.get("where") or {}
        if j.get("op") != "cell" or w.get("analysis") != analysis or "kind" not in w:
            continue
        if kind_as_text is not None:
            if w["kind"] != kind_as_text:
                continue
        else:
            try:
                k = C.kind_of_the_key(w["kind"], r["quantity_id"])
            except SystemExit:
                continue
            if k["kind"] != kind or k.get("variant") != variant:
                continue
        if row is not None:
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


# ---------------------------------------------------------------------------------------------
# typed constants of an earlier generator, read from its source
# ---------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------
# words of the pairs of the comparison of the programs, built from cells
# ---------------------------------------------------------------------------------------------
VALUES_OF_THE_COMPARISON = "program_comparison_values_20261002.csv"
CHAINS = "chains_beast_sets_20261002.csv"        # table of inputs of the third task: read for the dagger only
FIRST_SUMMARY = "the summary of the analysis"    # cell of the summary that is the estimate of BEAST X
CUTS = (" (", " that ", " instead of ")           # a cell is used up to the first of these
PAIRS_NOT_DRAWN = {"B": "pair 'B' of the table (setting B beside the primary analysis): not drawn as a pair; in the "
                        "article of this round the pair of setting B is the pair 'B_estimated'"}


def cut_cell(cell):
    k = min([cell.index(c) for c in CUTS if c in cell] or [len(cell)])
    return cell[:k].strip()


def items_of(cell):
    return [cut_cell(x.strip()) for x in cell.split(";") if x.strip()]


def continued_chains(inputs, table, name):
    """{analysis: {chains_continued, chains, read_from, ...}} for the analyses of `table` (file `name`) whose
    chains were continued from a saved state: rows of the table of chains whose cell started_again begins
    with 'continued', and cells note_on_the_analysis that say so of every chain. Where the note of a summary
    cut at the saved states names the count of the analysis that it cuts, the two counts are compared."""
    out = {}
    ch = inputs.table(CHAINS).need("analysis", "started_again")
    for a in sorted({r["analysis"] for r in table.rows}):
        rows = ch.where(analysis=a)
        n = sum(1 for i in rows if ch.rows[i]["started_again"].startswith("continued"))
        if n:
            out[a] = {"chains_continued": n, "chains": len(rows), "read_from": f"{CHAINS}, column started_again, "
                      f"rows of {a} whose cell begins with 'continued'", "table": CHAINS}
    for i, r in enumerate(table.rows):
        m = re.search(r"every chain of this analysis was continued from a saved state.*\((\d+) of (\d+) chains\)",
                      r["note_on_the_analysis"])
        if m and r["analysis"] not in out:
            out[r["analysis"]] = {"chains_continued": int(m.group(1)), "chains": int(m.group(2)),
                                  "read_from": f"{name}, column note_on_the_analysis, row {i}", "table": name, "row": i}
    for i, r in enumerate(table.rows):
        m = (re.search(r"not the summary of the analysis (\S+): each of the (\d+) chains that were continued", r["note_on_the_analysis"])
             or re.search(r"not the summary of the analysis (\S+): .*\((\d+) of \d+ chains cut\)", r["note_on_the_analysis"]))
        if m and m.group(1) in out:
            out[m.group(1)]["count_in_the_note_of_the_cut_summary"] = int(m.group(2))
            out[m.group(1)]["the_two_counts_are_equal"] = int(m.group(2)) == out[m.group(1)]["chains_continued"]
    return out


def words_of_the_pairs(inputs):
    """Per pair of the table of values of the comparison, in the order of the table: the capital letters at the
    head of the cell pair, the words of the cell what_differs_from_the_primary_analysis_of_Delphy of the row of
    Delphy up to the first of CUTS, the row of that cell, and the parts of the cell of the row of BEAST X that
    the words of the pair do not say (without 'the program')."""
    T = inputs.table(VALUES_OF_THE_COMPARISON).need("pair", "program", "analysis", "summary_of_the_analysis_of_BEAST_X",
                                                    "what_differs_from_the_primary_analysis_of_Delphy")
    out = {}
    for i, r in enumerate(T.rows):
        if r["summary_of_the_analysis_of_BEAST_X"] != FIRST_SUMMARY:
            continue
        w = out.setdefault(r["pair"], {"letter": re.match(r"[A-Z]+", r["pair"]).group(0)})
        key = "Delphy" if r["program"] == "Delphy" else "BEAST X"
        if key + "_row" not in w:
            w[key + "_row"] = i
            w[key + "_cell"] = r["what_differs_from_the_primary_analysis_of_Delphy"]
            w[key + "_analysis"] = r["analysis"]
    for p, w in out.items():
        if "Delphy_row" not in w or "BEAST X_row" not in w:
            raise C.Stop(f"{VALUES_OF_THE_COMPARISON}: pair {p} lacks a program")
        w["words"] = cut_cell(w["Delphy_cell"])
        w["extra_of_BEAST_X"] = [x for x in items_of(w["BEAST X_cell"]) if x != "the program" and x not in items_of(w["Delphy_cell"])]
    return out


def constants_of(inputs, generator, names):
    """{name: (value, first line, last line)} of the assignments of the module level of the generator (also
    where several names are assigned in one line, separated by ';', or as a tuple 'A, B = 1, 2'). A value that
    is no literal is given as the text of its source."""
    src = inputs.text(generator)
    tree = ast.parse(src)
    out = {}

    def take(target, value):
        if isinstance(target, ast.Name) and target.id in names:
            try:
                v = ast.literal_eval(value)
            except (ValueError, SyntaxError):
                v = ast.get_source_segment(src, value)
                if isinstance(value, ast.Call) and isinstance(value.func, ast.Name) and value.func.id == "dict" and not value.args:
                    try:                                   # dict(left=1.42, top=0.36, ...)
                        v = {kw.arg: ast.literal_eval(kw.value) for kw in value.keywords}
                    except (ValueError, SyntaxError):
                        pass
            out.setdefault(target.id, (v, value.lineno, value.end_lineno))
    for node in tree.body:                 # statements of one line that ';' separates are statements of their own
        if not (isinstance(node, ast.Assign) and len(node.targets) == 1):
            continue
        t = node.targets[0]
        if isinstance(t, ast.Tuple) and isinstance(node.value, ast.Tuple) and len(t.elts) == len(node.value.elts):
            for a, b in zip(t.elts, node.value.elts):
                take(a, b)
        else:
            take(t, node.value)
    missing = [n for n in names if n not in out]
    if missing:
        raise C.Stop(f"{generator}: no assignment found for {', '.join(missing)}")
    return out


def line_of(inputs, generator, words, after=0):
    """1-based number of the first line of the generator (after line `after`) that holds the words."""
    for k, line in enumerate(inputs.lines(generator), 1):
        if k > after and words in line:
            return k
    raise C.Stop(f"{generator}: no line holds {words!r}")


# ---------------------------------------------------------------------------------------------
# labels of one or more lines
# ---------------------------------------------------------------------------------------------
_WIDTHS = {}


def text_width_in(text, size, weight="normal"):
    """Width of a text in inches, measured with the renderer on a figure of its own."""
    key = (text, size, weight)
    if key not in _WIDTHS:
        fig = plt.figure(figsize=(4, 1))
        t = fig.text(0, 0, text, fontsize=size, fontweight=weight)
        fig.canvas.draw()
        _WIDTHS[key] = t.get_window_extent(fig.canvas.get_renderer()).width / fig.dpi
        plt.close(fig)
    return _WIDTHS[key]


def wrap(text, width_in, size, weight="normal", hang="   ", seps=("; ", ", ", " ")):
    """The text on as few lines as its width asks for. A line is broken after '; ' where a line that ends there
    keeps inside width_in, else after ', ', else at a blank - each time at the last such place; a following
    line begins with `hang`."""
    text = " ".join(text.split("\n"))
    if text_width_in(text, size, weight) <= width_in:
        return text
    lines, rest = [], text
    while rest and text_width_in((hang if lines else "") + rest, size, weight) > width_in:
        cut = None
        for sep in seps:
            for m in re.finditer(re.escape(sep), rest):
                if rest[m.end():m.end() + 1] == "%":          # a number and its per cent sign stay on one line
                    continue
                head = rest[:m.end()].rstrip()
                if text_width_in((hang if lines else "") + head, size, weight) <= width_in:
                    cut = m.end()
            if cut is not None:
                break
        if cut is None:
            break
        lines.append((hang if lines else "") + rest[:cut].rstrip())
        rest = rest[cut:]
    lines.append((hang if lines else "") + rest)
    return "\n".join(lines)


STEP_R = 0.25                              # step of the limits of an axis of R (earlier generator of Figure S1, line 89)
TICK_R = (0.5, 1.0, 2.5)                   # labels every 0.5, every 1.0 on a range wider than 2.5 (the same line)


def axis_of_R(values, step=STEP_R, tick=TICK_R, also=(1.0,)):
    """R4 (h), rule for an axis of R: the smallest range in steps of `step` that holds every value given and
    the line R = 1; labels every tick[0], every tick[1] on a range wider than tick[2].
    Returns (lower limit, upper limit, step of the labels, labels)."""
    lo, hi = C.steps_range(list(values) + list(also), step)
    st = tick[0] if hi - lo <= tick[2] + 1e-9 else tick[1]
    t0 = math.ceil(lo / st - 1e-9) * st
    ticks = [round(t0 + j * st, 6) for j in range(int(math.floor((hi - t0) / st + 1e-9)) + 1)]
    return lo, hi, st, ticks


def free_place(ax, fig, boxes_px, steps, lines_x=(), pad_pt=2.0):
    """True if none of the boxes (display units) stands on a band or line of the panel. steps: (x0, x1, lower,
    upper) in data units of everything drawn; lines_x: x of vertical lines, in data units."""
    inv = ax.transData.inverted()
    pad = pad_pt * fig.dpi / 72.0
    for b in boxes_px:
        (x0, y0), (x1, y1) = inv.transform((b.x0 - pad, b.y0 - pad)), inv.transform((b.x1 + pad, b.y1 + pad))
        for (s0, s1, lo, hi) in steps:
            if s0 < x1 and s1 > x0 and lo < y1 and hi > y0:
                return False
        for lx in lines_x:
            if x0 <= lx <= x1:
                return False
    return True


def panel_head(ax, letter, title=None):
    """The letter and the title of a panel as the common module of the main figures sets them; where the title
    has more than one line, the letter stands on the level of its first line."""
    out = C.panel_head(ax, letter, title)
    n = (title or "").count("\n")
    if n and letter:
        x, y = out[0].xyann
        out[0].xyann = (x, y + n * 1.2 * C.SIZES[0])
    return out


class Rows:
    """Headings and rows from the top down; a label of n lines takes 1 + (n - 1) * share units, so that the
    lines of two labels never touch. y of an item: the middle of its lines."""

    def __init__(self, heading_gap=0.45, share=0.85):
        self.items, self.y, self.heading_gap, self.share = [], 0.0, heading_gap, share

    def _add(self, kind, shown, **kw):
        n = shown.count("\n") + 1
        h = 1.0 + (n - 1) * self.share
        item = dict(kind=kind, y=self.y + (h - 1.0) / 2.0, lines=n, **kw)
        self.items.append(item)
        self.y += h
        return item

    def heading(self, text, **kw):
        if self.items:
            self.y += self.heading_gap
        return self._add("heading", text, text=text, **kw)

    def row(self, key, label, **kw):
        return self._add("row", label, key=key, label=label, **kw)

    @property
    def units(self):
        return self.y

    def rows(self):
        return [i for i in self.items if i["kind"] == "row"]


def draw_labels(fig, ax, layout, x_in, indent_in=0.12):
    tr = mpl.transforms.blended_transform_factory(fig.dpi_scale_trans, ax.transData)
    out = []
    for it in layout.items:
        head = it["kind"] == "heading"
        out.append(ax.text(x_in + (0 if head else indent_in), it["y"], it["text"] if head else it["label"],
                           transform=tr, ha="left", va="center", linespacing=1.15,
                           fontsize=C.SIZES[1] if head else C.SIZES[2],
                           fontweight="bold" if head else "normal", color="black", clip_on=False))
    return out


# ---------------------------------------------------------------------------------------------
# read-back and R1
# ---------------------------------------------------------------------------------------------
def read_back(values_path, inputs):
    """Reads every cell that a row of the file of values names AGAIN from its table and compares it with the
    columns of the row (value, lower, upper, x, x_end, label - whichever the row names in source_columns).
    A row of the exception to R1 (drawn_as begins with 'computed') names the column from which it is
    computed and is counted apart: its check is the recomputation in the program of the figure."""
    rows = C.read_csv(values_path)
    out = {"rows": len(rows), "rows_with_a_source": 0, "rows_whose_cells_are_equal": 0, "cells_compared": 0,
           "rows_without_a_source": 0, "rows_computed_by_the_exception_to_R1": 0, "not_equal": [],
           "rows_drawn_with_a_source": 0, "rows_not_drawn_with_a_source": 0,
           "sources_that_the_table_of_inputs_of_the_supplement_does_not_name": []}
    for k, r in enumerate(rows):
        if not r["source_table"]:
            out["rows_without_a_source"] += 1
            if r["value"] or r["lower"] or r["upper"]:
                out["not_equal"].append((k, "a row without a source holds a value"))
            continue
        if r["drawn_as"].startswith(COMPUTED):
            out["rows_computed_by_the_exception_to_R1"] += 1
            if not inputs.of_the_supplement(r["source_table"]):
                out["sources_that_the_table_of_inputs_of_the_supplement_does_not_name"].append((k, r["source_table"]))
            if r["source_version_id"] != inputs.version(r["source_table"]):
                out["not_equal"].append((k, "version id"))
            continue
        out["rows_with_a_source"] += 1
        out["rows_not_drawn_with_a_source" if r["drawn_as"].startswith(NOT_DRAWN) else "rows_drawn_with_a_source"] += 1
        ok = True
        holds_a_value = bool(r["value"] or r["lower"] or r["upper"])
        if holds_a_value and not inputs.of_the_supplement(r["source_table"]):
            out["sources_that_the_table_of_inputs_of_the_supplement_does_not_name"].append((k, r["source_table"]))
        if r["source_version_id"] != inputs.version(r["source_table"]):
            ok = False
            out["not_equal"].append((k, "version id"))
        named = dict(p.split("=", 1) for p in r["source_columns"].split("; ") if "=" in p)
        for part in ("value", "lower", "upper"):
            if r[part] and part not in named:
                ok = False
                out["not_equal"].append((k, f"{part} is filled, no cell is named for it"))
        if r["source_table"].endswith(".csv"):
            t = inputs.table(r["source_table"])
            for part, col in named.items():
                if IN_CELL in col:              # a part of a cell: 'column ~ pattern with one group'
                    col, pat = col.split(IN_CELL, 1)
                    m = re.search(pat, t.cell(int(r["source_row_0based"]), col))
                    cell = m.group(1) if m else None
                    out["parts_of_cells_compared"] = out.get("parts_of_cells_compared", 0) + 1
                else:
                    cell = t.cell(int(r["source_row_0based"]), col)
                out["cells_compared"] += 1
                if cell != r[part]:
                    ok = False
                    out["not_equal"].append((k, f"{part}: file of values {r[part]!r}, cell {cell!r}"))
        else:                                   # a document or a program: the text that the pattern finds in the line
            line = inputs.lines(r["source_table"])[int(r["source_row_0based"])]
            for part, pat in named.items():
                m = re.search(pat[len("pattern "):], line) if pat.startswith("pattern ") else None
                out["cells_compared"] += 1
                if not m or m.group(1) != r[part]:
                    ok = False
                    out["not_equal"].append((k, f"{part}: the pattern does not find {r[part]!r} in the line"))
        out["rows_whose_cells_are_equal"] += ok
    return out


def by_panel_and_element(rows):
    by = {}
    for r in rows:
        k = f"{r['panel']} | {r['element']}"
        by[k] = by.get(k, 0) + 1
    return dict(sorted(by.items()))


def files_named_by_the_specification(inputs, qid):
    """The files that the specification names for a quantity, by where they stand: among the inputs of the
    supplement, among the inputs of the third task only, among neither. Counted by program."""
    spec = inputs.table(C.SPECIFICATION).need("quantity_id", "file")
    rows = [r for r in spec.rows if r["quantity_id"] == qid]
    if len(rows) != 1:
        raise C.Stop(f"specification: {len(rows)} rows {qid}")
    files = [x.strip() for x in rows[0]["file"].split(";") if x.strip()]
    here = [f for f in files if inputs.of_the_supplement(f)]
    third = [f for f in files if f not in here and f in inputs.third["by_file"]]
    none = [f for f in files if f not in here and f not in third]
    return {"files": files, "among_the_inputs_of_the_supplement": here, "among_the_inputs_of_the_third_task_only": third,
            "not_among_the_inputs": none}


def not_among_the_inputs(inputs, names, what):
    """Stops where a file that the program says is not among the inputs is among them. Returns the record."""
    found = [n for n in names if inputs.of_the_supplement(n) or n in inputs.third["by_file"] or n in inputs.paths]
    if found:
        raise C.Stop(f"{what}: the program says that these files are not among the inputs, and they are: {', '.join(found)}")
    return {n: "not among the inputs (neither table of inputs names it)" for n in names}


def files_opened(inputs):
    """The files of the inputs that the program opened, by the table of inputs that names them."""
    names = set(inputs._tables) | set(inputs._texts) | set(inputs.__dict__.get("_gz", ()))
    out = {"named_by_the_table_of_inputs_of_the_supplement": [], "named_by_the_table_of_inputs_of_the_third_task_only": [],
           "the_two_tables_of_inputs": [], "other": []}
    for n in sorted(names):
        if n in (TABLE_OF_INPUTS, C.TABLE_OF_INPUTS):
            out["the_two_tables_of_inputs"].append(n)
        elif inputs.of_the_supplement(n):
            out["named_by_the_table_of_inputs_of_the_supplement"].append(n)
        elif n in inputs.third["by_file"]:
            out["named_by_the_table_of_inputs_of_the_third_task_only"].append(n)
        else:
            out["other"].append(n)
    return out


def standard(inputs, fig, values_path, outdir, stem, program, extra, allow_overlap=(), earlier_allowed=()):
    """The checks that every figure has; writes <stem>_checks.json and returns (checks, ok)."""
    opened = files_opened(inputs)                  # what the program opened to build and to draw the figure
    patterns = C.word_patterns(inputs)
    more = files_opened(inputs)                    # ... and what the check of the words of R4 (j) opens beside
    opened = {"to_build_and_to_draw_the_figure": opened,
              "beside_for_the_names_that_must_not_appear_(R4_j)": {k: [n for n in v if n not in opened[k]] for k, v in more.items()}}
    checks = C.check_figure(fig, patterns, allow_overlap=allow_overlap)
    saved = C.save_figure(fig, stem, outdir)
    rb = read_back(values_path, inputs)
    rows = C.read_csv(values_path)
    saved = {k: (os.path.basename(v) if isinstance(v, str) else v) for k, v in saved.items()}
    marked = [{"panel": r["panel"], "element": r["element"], "label": r["label"], "analysis": r["analysis"],
               "source_table": r["source_table"], "source_row_0based": r["source_row_0based"],
               "drawn": not r["drawn_as"].startswith(NOT_DRAWN),
               "meets_convergence_criterion": r["meets_convergence_criterion"]}
              for r in rows if C.marked_by_R4f(r["meets_convergence_criterion"])]
    other = rb["sources_that_the_table_of_inputs_of_the_supplement_does_not_name"]
    earlier = [[k, name, inputs.table_of(name), rows[k]["element"]] for k, name in other if inputs.table_of(name)]
    r1 = {"rows_with_a_value_whose_source_is_a_file_of_the_table_of_inputs_of_the_supplement":
          sum(1 for r in rows if r["source_table"] and (r["value"] or r["lower"] or r["upper"])
              and not r["drawn_as"].startswith(COMPUTED) and inputs.of_the_supplement(r["source_table"])),
          "rows_with_a_value_whose_source_is_a_file_of_the_table_of_inputs_of_an_earlier_task": earlier,
          "elements_for_which_the_program_reads_a_table_of_an_earlier_task": sorted(earlier_allowed),
          "rows_with_a_value_whose_source_no_table_of_inputs_names": [[k, name] for k, name in other if not inputs.table_of(name)],
          "rows_computed_by_the_exception": rb["rows_computed_by_the_exception_to_R1"],
          "every_value_is_one_cell_of_a_file_of_the_table_of_inputs_of_the_supplement": not other}
    r1["ok"] = (not r1["rows_with_a_value_whose_source_no_table_of_inputs_names"]
                and all(e[3] in earlier_allowed for e in earlier))
    out = {"program": program, "figure": saved, "values_file": os.path.basename(values_path),
           "figure_sha256": {"png": C.sha256(os.path.join(outdir, saved["png"])),
                             "pdf": C.sha256(os.path.join(outdir, saved["pdf"]))},
           "values_sha256": C.sha256(values_path), "rows": len(rows),
           "rows_by_panel_and_element": by_panel_and_element(rows),
           "rows_drawn": sum(not r["drawn_as"].startswith(NOT_DRAWN) for r in rows),
           "rows_not_drawn": sum(r["drawn_as"].startswith(NOT_DRAWN) for r in rows),
           "read_back": rb,
           "R1": r1,
           "values_marked_by_R4f_with_the_cell_as_it_reads": marked,
           "banner": {"words_of_a_banner_in_a_text_of_the_picture": checks.get("words_of_a_banner", []),
                      "the_same_words_in_another_case": checks.get("words_of_a_banner_in_another_case", []),
                      "words_in_the_names_of_the_files": sorted(
                          {w for n in (saved["png"], saved["pdf"], os.path.basename(values_path))
                           for w in ("INTERNAL", "PROVISIONAL", "PREVIEW", "FIXTURE", "TEST") if w.lower() in n.lower()})},
           "files_opened": opened,
           "figure_checks": checks, "drawing": extra,
           "table_of_inputs_version": inputs.version_of_the_table_of_inputs,
           "table_of_inputs_of_the_third_task_version": inputs.third["version"]}
    out["banner"]["ok"] = not out["banner"]["words_of_a_banner_in_a_text_of_the_picture"] \
        and not out["banner"]["words_in_the_names_of_the_files"]
    ok = bool(checks["ok"] and not rb["not_equal"] and out["R1"]["ok"] and out["banner"]["ok"])
    return out, ok


def finish(out, ok, outdir, stem, more_ok=True, said=None):
    out["ok"] = bool(ok and more_ok)
    C.json_dump(out, os.path.join(outdir, stem + "_checks.json"))
    ch = out["figure_checks"]
    line = {"png": out["figure"]["png"], "rows": out["rows"], "checks_ok": out["ok"],
            "min_fontsize": ch["min_fontsize"], "read_back_not_equal": out["read_back"]["not_equal"][:5],
            "overlaps": ch["text_overlaps"][:6], "outside": ch["outside_figure"][:6], "on_spine": ch["text_on_spine"][:6],
            "words": ch["words_of_R4j"][:6], "height_in": out["figure"]["height_in"]}
    line.update(said or {})
    print(json.dumps(line))
    return 0 if out["ok"] else 1


def draw_only(argv, draw, stem):
    """--draw-only VALUES.csv [OUTDIR]: the figure from its file of values and from nothing else."""
    outdir = argv[3] if len(argv) > 3 else "."
    os.makedirs(outdir, exist_ok=True)
    fig, _ = draw(argv[2])
    info = C.save_figure(fig, stem, outdir)
    print(json.dumps({"drawn_from": os.path.basename(argv[2]), "png": os.path.basename(info["png"])}))
    return 0
