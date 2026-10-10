"""main_fig_common_20261006.py
Common module of the programs that draw Figures 1 and 2 and build Tables 1 and 2 of the article
from the tables of the round of October 2026.

Used by
    fig_ne_R_cases_20261006.py
    fig_R_robustness_focus_20261006.py
    main_tables_20261006.py
    build_documents_20261006.py

Input of every program: a JSON 'file name -> path', built from the table of inputs
(INTERNAL_main_figures_inputs_final_20261006.csv) with make_inputs_json_20261004.py. The JSON also holds
the table of inputs itself and, under the key '<name of the table of inputs>#version', its version.
File and key of a quantity are read from the latest specification of the canonical numbers
(INTERNAL_canonical_numbers_specification_20261004.csv), by quantity_id.

Rules that this module carries:
  R1  every value is the content of ONE cell of a table of the inputs; cells are read as TEXT and written
      into the files of values as they stand in the table;
  R2  a row is found by the KIND of the quantity, its period and end_used, never by the label alone;
  R4g a token is printed with format_token of the program of the formats (imported from its path);
  R4j the list FORBIDDEN is READ from the source of the earlier common module (not executed), and the
      internal names of the round are read from the schedules and the tables.
One thread: the variables of the thread pools are set before numpy is imported.
"""
import copy
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS",
           "VECLIB_MAXIMUM_THREADS"):
    os.environ[_v] = "1"
os.environ.setdefault("MPLBACKEND", "Agg")

import ast
import csv
import datetime as dt
import hashlib
import importlib.util
import json
import math
import re
import sys
from decimal import Decimal, getcontext

import numpy as np
import pandas as pd
import matplotlib as mpl

mpl.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.text as mtext
from matplotlib.patches import Rectangle

getcontext().prec = 60
SUFFIX = "_20261006"
TABLE_OF_INPUTS = "INTERNAL_main_figures_inputs_final_20261006.csv"
TABLE_OF_TOKENS = "main_figures_tokens_20261004.csv"
SPECIFICATION = "INTERNAL_canonical_numbers_specification_20261004.csv"   # the latest: it decides file and key
SPECIFICATION_OF_3_OCTOBER = "canonical_numbers_specification_20261003.csv"   # the copy in the table of tokens
TEMPLATE = "INTERNAL_article_template_20261003.txt"
FORMAT_PROGRAM = "format_article_numbers_20261003.py"
CONVENTIONS = "post_figure_conventions_20260924.md"
EARLIER_COMMON = "post_fig_common_20260924.py"
EARLIER_FIG1 = "fig_ne_R_cases_20260924.py"
EARLIER_FIG2 = "fig_R_robustness_focus_20260925.py"
PLAN = "analysis_plan_round_20261002.md"
PHASE5_SPEC = "phase5_specification_20261003.md"
SCHEDULE_PHASE3 = "chain_schedule_round_20261002.csv"
SCHEDULE_PHASE5 = "phase5_schedule_20261003.csv"
PERIODS = "periods_round_20261002.csv"
ENDS_BY_SET = "periods_ends_by_set_20261002.csv"
CASES = "weekly_confirmed_cases_by_province_20261002.csv"
CASES_README = "README_case_series_20261002.md"
T_GT = "generation_times_sensitivity_20261002.csv"
T_PRIOR = "prior_alone_sensitivity_20261002.csv"
T_PERIODS1 = "primary_periods_of_table_1_sensitivity_20261002.csv"
T_PAIRS = "primary_pairs_sensitivity_20261002.csv"
T_INTERVALS_S = "primary_intervals_sensitivity_20261002.csv"
T_LAST = "primary_last_periods_of_other_analyses_sensitivity_20261002.csv"

VALUE_COLUMNS = ["figure", "panel", "element", "label", "value", "lower", "upper", "x", "x_end", "unit",
                 "drawn_as", "source_table", "source_version_id", "source_row_0based", "source_columns",
                 "note", "analysis", "quantity_id", "delivered", "meets_convergence_criterion"]

NOTE_SIMULATED = "simulated; compared within the error of the simulation, not reproduced"
NOT_YET = "not yet delivered"
# these words stand nowhere in a picture (compared as written, in capitals)
WORDS_OF_A_BANNER = ("INTERNAL", "PROVISIONAL", "PREVIEW")
DOUBLE_DAGGER = "\u2021"
NDASH = "\u2013"
MONTHS = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
MONTHS_LONG = ["January", "February", "March", "April", "May", "June", "July", "August", "September",
               "October", "November", "December"]


class Stop(SystemExit):
    """The program stops with a message (a file or a column that it needs is missing, a key finds
    no row or more than one, a check of a rule of drawing fails)."""

    def __init__(self, message):
        super().__init__("STOP: " + message)
        self.message = message


# ---------------------------------------------------------------------------------------------
# inputs
# ---------------------------------------------------------------------------------------------
class Table:
    """A table of the inputs, every cell as the TEXT of the file. Row i is the 0-based data row."""

    def __init__(self, name, path, version):
        self.name, self.path, self.version = name, path, version
        with open(path, newline="", encoding="utf-8") as fh:
            rd = csv.reader(fh)
            try:
                self.columns = next(rd)
            except StopIteration:
                raise Stop(f"{name}: the file is empty")
            self.rows = []
            for i, rec in enumerate(rd):
                if len(rec) != len(self.columns):
                    raise Stop(f"{name}: data row {i} holds {len(rec)} cells, the header {len(self.columns)}")
                self.rows.append(dict(zip(self.columns, rec)))
        if len(set(self.columns)) != len(self.columns):
            raise Stop(f"{name}: a column name stands twice in the header")

    def need(self, *columns):
        missing = [c for c in columns if c not in self.columns]
        if missing:
            raise Stop(f"{self.name}: column(s) missing: {', '.join(missing)}")
        return self

    def cell(self, row, column):
        self.need(column)
        return self.rows[row][column]

    def where(self, **eq):
        self.need(*eq.keys())
        return [i for i, r in enumerate(self.rows) if all(r[k] == v for k, v in eq.items())]

    def one(self, what, **eq):
        hits = self.where(**eq)
        if len(hits) != 1:
            raise Stop(f"{self.name}: {what}: expected 1 row, found {len(hits)} ({eq})")
        return hits[0]

    def __len__(self):
        return len(self.rows)


class Inputs:
    """The JSON 'file name -> path' and the table of inputs (version and analysis of every file)."""

    def __init__(self, json_path):
        if not os.path.exists(json_path):
            raise Stop(f"the JSON of the inputs is missing: {json_path}")
        with open(json_path, encoding="utf-8") as fh:
            self.paths = json.load(fh)
        self.json_path = json_path
        self._tables, self._texts = {}, {}
        if TABLE_OF_INPUTS not in self.paths:
            raise Stop(f"the JSON holds no entry {TABLE_OF_INPUTS}")
        self.version_of_the_table_of_inputs = self.paths.get(TABLE_OF_INPUTS + "#version", "")
        t = Table(TABLE_OF_INPUTS, self._path(TABLE_OF_INPUTS), self.version_of_the_table_of_inputs)
        t.need("file", "version", "used_for", "what_it_is", "analysis", "state_of_the_file", "note")
        self.inputs = t
        self.by_file = {r["file"]: r for r in t.rows}
        if len(self.by_file) != len(t.rows):
            raise Stop("the table of inputs names a file twice")
        self._tables[TABLE_OF_INPUTS] = t

    def _path(self, name):
        if name not in self.paths:
            raise Stop(f"the JSON of the inputs holds no entry for {name}")
        p = self.paths[name]
        if not os.path.exists(p):
            raise Stop(f"file of the inputs not found: {name} -> {p}")
        return p

    def has(self, name):
        return name in self.paths and os.path.exists(self.paths[name])

    def path(self, name):
        return self._path(name)

    def added(self, name):
        """A table that the table of inputs does NOT name and that is added to the inputs for ONE quantity: the
        JSON of the inputs holds it with the entries '<name>#version', '<name>#only_for_quantity_id' and
        '<name>#added_by'. Returns the three, or None."""
        if name in self.by_file or name == TABLE_OF_INPUTS:
            return None
        v, q = self.paths.get(name + "#version"), self.paths.get(name + "#only_for_quantity_id")
        if not v or not q or name not in self.paths:
            return None
        return {"version": v, "only_for_quantity_id": q, "added_by": self.paths.get(name + "#added_by", "")}

    def allowed(self, name, quantity_id=None):
        """Whether the file may be opened: every file of the table of inputs; an added table for its one
        quantity and for no other."""
        if name in self.by_file or name == TABLE_OF_INPUTS:
            return True
        a = self.added(name)
        return bool(a and quantity_id is not None and a["only_for_quantity_id"] == quantity_id)

    def version(self, name):
        if name == TABLE_OF_INPUTS:
            return self.version_of_the_table_of_inputs
        if name not in self.by_file:
            a = self.added(name)
            if a:
                return a["version"]
            raise Stop(f"the table of inputs does not name {name}")
        return self.by_file[name]["version"]

    def audit(self, name):
        return self.by_file[name]["state_of_the_file"] if name in self.by_file else ""

    def table(self, name, quantity_id=None):
        a = self.added(name)
        if a and a["only_for_quantity_id"] != quantity_id:
            raise Stop(f"{name} is among the inputs for the quantity {a['only_for_quantity_id']} only "
                       f"({a['added_by']}); it is not opened for {quantity_id or 'another purpose'}")
        if name not in self._tables:
            self._tables[name] = Table(name, self._path(name), self.version(name))
        return self._tables[name]

    def text(self, name):
        if name not in self._texts:
            with open(self._path(name), encoding="utf-8") as fh:
                self._texts[name] = fh.read()
        return self._texts[name]

    def lines(self, name):
        return self.text(name).split("\n")

    # -- which analyses have their table of estimates among the inputs ----------------------------
    @staticmethod
    def analyses_of_a_row(r):
        """The analyses that a row of the table of inputs names (one, or several separated by '; ')."""
        return [a.strip() for a in r["analysis"].split(";") if a.strip()]

    def estimates_table_of(self, analysis):
        """Name of the table of estimates that the table of inputs gives the analysis, or None."""
        hits = [r["file"] for r in self.inputs.rows
                if analysis in self.analyses_of_a_row(r) and r["what_it_is"].startswith("table of estimates")]
        if len(hits) > 1:
            raise Stop(f"the table of inputs names {len(hits)} tables of estimates for {analysis}")
        if not hits or not self.has(hits[0]):
            return None
        t = self.table(hits[0]).need("analysis")
        return hits[0] if t.where(analysis=analysis) else None

    def chains_table_of(self, analysis):
        hits = [r["file"] for r in self.inputs.rows
                if analysis in self.analyses_of_a_row(r) and r["what_it_is"].startswith("table of chains")]
        return hits[0] if len(hits) == 1 and self.has(hits[0]) else None

    def delivered(self, analysis):
        return self.estimates_table_of(analysis) is not None


def pattern_in_document(inputs, name, pattern, what):
    """The one text that `pattern` (one group) finds in the document: (text, 0-based line)."""
    hits = []
    for i, line in enumerate(inputs.lines(name)):
        for m in re.finditer(pattern, line):
            hits.append((m.group(1), i))
    texts = sorted({h[0] for h in hits})
    if len(texts) != 1:
        raise Stop(f"{name}: {what}: the pattern {pattern!r} finds {len(texts)} distinct texts")
    return hits[0][0], hits[0][1], len(hits)


# ---------------------------------------------------------------------------------------------
# the program of the formats (imported from its path; format_token is called, nothing else)
# ---------------------------------------------------------------------------------------------
_FORMATS = {}


def format_program(inputs):
    if "m" not in _FORMATS:
        spec = importlib.util.spec_from_file_location("format_article_numbers_20261003",
                                                      inputs.path(FORMAT_PROGRAM))
        m = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(m)
        if not hasattr(m, "format_token"):
            raise Stop(f"{FORMAT_PROGRAM} holds no function format_token")
        _FORMATS["m"] = m
    return _FORMATS["m"]


def slots_of(fmt):
    f = json.loads(fmt) if isinstance(fmt, str) else fmt
    return f.get("slots", [])


def print_token(inputs, fmt, values):
    """(text, error) of a token: format_token(format, values) of the program of the formats."""
    m = format_program(inputs)
    try:
        return m.format_token(fmt, values), ""
    except m.FormatError as e:
        return "", "FormatError: " + str(e)
    except Exception as e:                      # reported, never hidden
        return "", f"{type(e).__name__}: {e}"


# ---------------------------------------------------------------------------------------------
# kinds of quantities (R2: a row is found by kind, period and end_used)
# ---------------------------------------------------------------------------------------------
_Q_PERIOD = re.compile(r"^(growth rate|reproduction number|doubling time),?\s+(P\d+)(?:\s+\((own|common) end\))?$")
_Q_TMRCA = re.compile(r"^date of the most recent common ancestor \((date|decimal year)\)$")
_Q_PREC = re.compile(r"^precision of the Skygrid(?: \(estimated\))?$")
_Q_PROB = re.compile(r"^P\((.+?)\)(?:\s+\((own|common) end\))?$")


def parse_quantity(label):
    s = label.strip()
    m = _Q_PERIOD.match(s)
    if m:
        return {"kind": m.group(1), "period": m.group(2), "end_in_label": m.group(3) or ""}
    m = _Q_TMRCA.match(s)
    if m:
        return {"kind": "date of the most recent common ancestor", "variant": m.group(1)}
    if _Q_PREC.match(s):
        return {"kind": "precision of the Skygrid"}
    m = _Q_PROB.match(s)
    if m:
        return {"kind": "probability", "expression": re.sub(r"\s+", " ", m.group(1)).strip(),
                "end_in_label": m.group(2) or ""}
    if s in ("evolutionary rate", "kappa"):
        return {"kind": s}
    return {"kind": "not recognised", "label": s}


# The latest specification (column key_kind and the entry 'kind' of a definition for the program) names the
# KIND of a quantity in words of its own; the period of the quantity stands in the days of the key
# (period_start, period_end), not in the name of the kind. Each kind is read as the kind that parse_quantity
# gives the label of a row. A kind that is not listed here stops the program.
KINDS_OF_THE_LATEST_SPECIFICATION = {
    "reproduction number": {"kind": "reproduction number"},
    "growth rate": {"kind": "growth rate"},
    "doubling time": {"kind": "doubling time"},
    "evolutionary rate": {"kind": "evolutionary rate"},
    "kappa": {"kind": "kappa"},
    "precision of the Skygrid": {"kind": "precision of the Skygrid"},
    "date of the most recent common ancestor, as a date":
        {"kind": "date of the most recent common ancestor", "variant": "date"},
    "date of the most recent common ancestor, as a decimal year":
        {"kind": "date of the most recent common ancestor", "variant": "decimal year"},
    "probability that R is above 1": {"kind": "probability", "expression_pattern": r"R of P\d+ > 1"},
    "probability that R of P2 is below R of P1": {"kind": "probability", "expression": "R of P2 < R of P1"},
}


def kind_of_the_key(kind_in_words, what=""):
    """The kind of the latest specification as the arguments of find_row (kind, variant, expression ...)."""
    if kind_in_words not in KINDS_OF_THE_LATEST_SPECIFICATION:
        raise Stop(f"{SPECIFICATION}: {what}: the kind {kind_in_words!r} is not known to this program")
    return dict(KINDS_OF_THE_LATEST_SPECIFICATION[kind_in_words])


def find_row(table, what, analysis, kind, period=None, end_used=None, period_start=None,
             period_end=None, variant=None, expression=None, expression_pattern=None):
    """0-based row of `table` with this analysis, KIND of quantity, period and end_used. Exactly one
    row, else the program stops. The label of the row is NOT compared as a string."""
    table.need("analysis", "quantity", "period_start", "period_end", "end_used")
    hits = []
    for i, r in enumerate(table.rows):
        if r["analysis"] != analysis:
            continue
        q = parse_quantity(r["quantity"])
        if q["kind"] != kind:
            continue
        if period is not None and q.get("period") != period:
            continue
        if variant is not None and q.get("variant") != variant:
            continue
        if expression is not None and q.get("expression") != expression:
            continue
        if expression_pattern is not None and not re.fullmatch(expression_pattern, q.get("expression", "")):
            continue
        if end_used is not None and r["end_used"] != end_used:
            continue
        if period_start is not None and r["period_start"] != period_start:
            continue
        if period_end is not None and r["period_end"] != period_end:
            continue
        if q.get("end_in_label") and q["end_in_label"] != r["end_used"]:
            raise Stop(f"{table.name}: row {i}: the label {r['quantity']!r} names the end "
                       f"{q['end_in_label']!r}, the column end_used holds {r['end_used']!r}")
        hits.append(i)
    if len(hits) != 1:
        raise Stop(f"{table.name}: {what}: expected 1 row, found {len(hits)} (analysis {analysis}, kind {kind}, "
                   f"period {period}, end_used {end_used}, start {period_start}, end {period_end})")
    return hits[0]


def rows_found(table, analysis, kind, period=None, end_used=None):
    """As find_row, but returns the list of rows (may be empty)."""
    out = []
    for i, r in enumerate(table.rows):
        if r["analysis"] != analysis:
            continue
        q = parse_quantity(r["quantity"])
        if q["kind"] == kind and (period is None or q.get("period") == period) \
                and (end_used is None or r["end_used"] == end_used):
            out.append(i)
    return out


# ---------------------------------------------------------------------------------------------
# the table of tokens with the file and the key of the LATEST specification
# ---------------------------------------------------------------------------------------------
SPEC = "specification: "
# (column of the latest specification, column 'specification: ...' of the table of tokens that holds the copy
# of the specification of 3 October)
COLUMNS_OF_THE_COPY = [
    ("quantity_in_words", "quantity_in_words"), ("track", "track"), ("file", "file"),
    ("file_is_delivered", "file_is_delivered"), ("key_analysis", "key_analysis"), ("key_kind", "key_quantity"),
    ("key_period_start", "key_period_start"), ("key_period_end", "key_period_end"),
    ("key_end_used", "key_end_used"), ("key_column", "key_column"),
    ("key_further_conditions", "key_further_conditions"), ("key_is_complete", "key_is_complete"),
    ("definition_for_the_program", "definition_for_the_program"),
    ("depends_on_step_5", "depends_on_step_5"), ("depends_on_step_6", "depends_on_step_6"),
    ("source_if_step_5_finds_that_the_programs_agree", "source_if_step_5_finds_that_the_programs_agree"),
    ("source_if_step_5_finds_that_they_do_not_agree", "source_if_step_5_finds_that_they_do_not_agree"),
    ("source_if_step_6_takes_the_smoothing_estimated", "source_if_step_6_takes_the_smoothing_estimated"),
    ("source_if_step_6_takes_the_smoothing_fixed", "source_if_step_6_takes_the_smoothing_fixed")]
FILE_OF_THE_KEY = ["file"]
COLUMNS_OF_THE_KEY = ["key_analysis", "key_period_start", "key_period_end", "key_end_used", "key_column",
                      "key_further_conditions"]
COMPARISON_COLUMNS = (["quantity_id", "tokens", "units_of_text", "row_0based_in_the_latest_specification",
                       "file_and_key_equal_to_the_copy", "file_equal", "key_equal", "columns_that_differ"]
                      + [f"{c}: {w}" for c in FILE_OF_THE_KEY + COLUMNS_OF_THE_KEY
                         for w in ("copy of 3 October", "latest specification")]
                      + ["key_quantity: copy of 3 October", "key_kind: latest specification",
                         "text_of_key_kind_equal_to_key_quantity",
                         "definition_for_the_program_equal_as_text",
                         "definition_for_the_program: copy of 3 October",
                         "definition_for_the_program: latest specification",
                         "key_is_complete: copy of 3 October", "key_is_complete: latest specification",
                         "status: copy of 3 October", "status: latest specification", "status_equal",
                         "rule_of_4_October_for_the_program"])
UNCHANGED_BY_DESIGN = "unchanged by design"


def tokens_by_the_latest_specification(I):
    """(table, comparison): the table of tokens of 4 October in which the cells 'specification: ...' of every
    token are those of the LATEST specification, found by quantity_id, and one row of comparison for every
    quantity (latest specification against the copy that the table of tokens holds). The column
    'specification: key_quantity' of the copy is not carried: 'specification: key_kind' takes its place.
    The file of the table of tokens is not changed; a source in the table of tokens keeps its row."""
    if getattr(I, "_tokens_latest", None) is not None:
        return I._tokens_latest
    tok = I.table(TABLE_OF_TOKENS).need("token", "unit_of_text", "quantity_id",
                                        "status_of_the_quantity_in_the_specification",
                                        *[SPEC + c for _, c in COLUMNS_OF_THE_COPY])
    sp = I.table(SPECIFICATION).need("quantity_id", "status", "rule_of_4_October_for_the_program",
                                     *[c for c, _ in COLUMNS_OF_THE_COPY])
    by_id = {}
    for i, r in enumerate(sp.rows):
        if r["quantity_id"] in by_id:
            raise Stop(f"{SPECIFICATION}: the quantity {r['quantity_id']} stands twice")
        by_id[r["quantity_id"]] = i
    merged = copy.copy(tok)
    merged.columns = [c for c in tok.columns if c != SPEC + "key_quantity"] + [
        SPEC + "key_kind", SPEC + "rule_of_4_October_for_the_program", SPEC + "row_0based"]
    merged.rows = []
    comparison = {}
    for r in tok.rows:
        q = r["quantity_id"]
        if q not in by_id:
            raise Stop(f"{SPECIFICATION}: the quantity {q} of the token {r['token']} is missing")
        s = sp.rows[by_id[q]]
        m = {k: v for k, v in r.items() if k != SPEC + "key_quantity"}
        for c_new, c_copy in COLUMNS_OF_THE_COPY:
            m[SPEC + (c_new if c_new == "key_kind" else c_copy)] = s[c_new]
        m[SPEC + "rule_of_4_October_for_the_program"] = s["rule_of_4_October_for_the_program"]
        m[SPEC + "row_0based"] = str(by_id[q])
        status_copy = r["status_of_the_quantity_in_the_specification"]
        if (status_copy == UNCHANGED_BY_DESIGN) != (s["status"] == UNCHANGED_BY_DESIGN):
            raise Stop(f"{q}: the status is {status_copy!r} in the table of tokens and {s['status']!r} in the "
                       "latest specification: the text of a quantity that is unchanged by design stands in the "
                       "table of tokens alone")
        m["status_of_the_quantity_in_the_specification"] = s["status"]
        merged.rows.append(m)
        differ = [c for c in FILE_OF_THE_KEY + COLUMNS_OF_THE_KEY if r[SPEC + c] != s[c]]
        rec = {"quantity_id": q, "tokens": r["token"], "units_of_text": r["unit_of_text"],
               "row_0based_in_the_latest_specification": str(by_id[q]),
               "file_and_key_equal_to_the_copy": "yes" if not differ else "no",
               "file_equal": "yes" if r[SPEC + "file"] == s["file"] else "no",
               "key_equal": "yes" if all(r[SPEC + c] == s[c] for c in COLUMNS_OF_THE_KEY) else "no",
               "columns_that_differ": "; ".join(differ),
               "key_quantity: copy of 3 October": r[SPEC + "key_quantity"],
               "key_kind: latest specification": s["key_kind"],
               "text_of_key_kind_equal_to_key_quantity": "yes" if r[SPEC + "key_quantity"] == s["key_kind"] else "no",
               "definition_for_the_program_equal_as_text":
                   "yes" if r[SPEC + "definition_for_the_program"] == s["definition_for_the_program"] else "no",
               "definition_for_the_program: copy of 3 October": r[SPEC + "definition_for_the_program"],
               "definition_for_the_program: latest specification": s["definition_for_the_program"],
               "key_is_complete: copy of 3 October": r[SPEC + "key_is_complete"],
               "key_is_complete: latest specification": s["key_is_complete"],
               "status: copy of 3 October": status_copy, "status: latest specification": s["status"],
               "status_equal": "yes" if status_copy == s["status"] else "no",
               "rule_of_4_October_for_the_program": s["rule_of_4_October_for_the_program"]}
        for c in FILE_OF_THE_KEY + COLUMNS_OF_THE_KEY:
            rec[f"{c}: copy of 3 October"] = r[SPEC + c]
            rec[f"{c}: latest specification"] = s[c]
        if q in comparison:
            same = {k: v for k, v in comparison[q].items() if k not in ("tokens", "units_of_text")}
            if same != {k: v for k, v in rec.items() if k not in ("tokens", "units_of_text")}:
                raise Stop(f"table of tokens: the copies of the specification differ between the tokens of {q}")
            comparison[q]["tokens"] += "; " + r["token"]
            if r["unit_of_text"] not in comparison[q]["units_of_text"].split("; "):
                comparison[q]["units_of_text"] += "; " + r["unit_of_text"]
        else:
            comparison[q] = rec
    I._tokens_latest = (merged, [comparison[q] for q in sorted(comparison)])
    return I._tokens_latest


def marked_by_R4f(cell):
    """R4 (f): the cell begins with 'no' and is none of 'not assessed', 'not evaluated', 'not applicable'."""
    c = cell.strip().lower()
    return c.startswith("no") and not c.startswith(("not assessed", "not evaluated", "not applicable"))


def is_number(text):
    try:
        return math.isfinite(float(text))
    except (TypeError, ValueError):
        return False


def num(text, what=""):
    if not is_number(text):
        raise Stop(f"{what}: the cell {text!r} is not a number")
    return float(text)


# ---------------------------------------------------------------------------------------------
# dates
# ---------------------------------------------------------------------------------------------
def to_day(text, what=""):
    try:
        return dt.date.fromisoformat(str(text).strip()[:10])
    except ValueError:
        raise Stop(f"{what}: the cell {text!r} is not a day written YYYY-MM-DD")


def is_day(text):
    try:
        dt.date.fromisoformat(str(text).strip())
        return True
    except ValueError:
        return False


def day_of_words(text, what=""):
    """'15 May 2026' or '31 August 2026' -> date (documents of the inputs write days in words)."""
    m = re.fullmatch(r"\s*(\d{1,2}) (\w+) (\d{4})\s*", text)
    if m:
        for names in (MONTHS, MONTHS_LONG):
            if m.group(2) in names:
                return dt.date(int(m.group(3)), names.index(m.group(2)) + 1, int(m.group(1)))
    raise Stop(f"{what}: {text!r} is not a day written 'D Month YYYY'")


def fmt_date(d, year=False):
    """A day as the earlier generators print it in a figure: '13 Feb' (post_fig_common_20260924.py, fmt_date)."""
    s = f"{d.day} {MONTHS[d.month - 1]}"
    return f"{s} {d.year}" if year else s


def fmt_ci(m, lo, hi, nd=2):
    """'1.23 (0.99-1.56)' as the earlier generators print it (post_fig_common_20260924.py, fmt_ci)."""
    return f"{m:.{nd}f} ({lo:.{nd}f}{NDASH}{hi:.{nd}f})"


def num_label(x, nd=1):
    """A number of a label as the earlier generator of Figure 2 prints it (function num)."""
    return f"{x:.{nd}f}".rstrip("0").rstrip(".")


def dec_year(d):
    """Decimal year as the earlier round defines it (dec_year)."""
    a = dt.date(d.year, 1, 1)
    b = dt.date(d.year + 1, 1, 1)
    return d.year + (d - a).days / (b - a).days


# ---------------------------------------------------------------------------------------------
# colours and style (conventions file, binding)
# ---------------------------------------------------------------------------------------------
BLUE = "#27496d"          # Delphy
BLUE_ALPHA = 0.18         # ribbons
GREEN = "#2f8f5b"         # BEAST X (no element of the two figures of this round uses it)
ORANGE_ITURI = "#dd8452"
ORANGE_NK = "#8b3a0e"
ORANGE_OTHER = "#f2c29b"
GREY_PRIOR = "#c8c8c8"    # smoothing prior alone
GREY_REF = "0.5"          # reference lines (declaration)
GREY_BAR = "#dedede"      # genomes per interval (earlier generator of Figure 1, line 197)
GREY_NOT_YET = "0.45"     # the words 'not yet delivered' (R4 e: 'in grey')

FIG_WIDTH = 7.0
MAX_HEIGHT = 8.5
SIZES = (10, 9, 8)
FONT = "Liberation Sans"
LETTER_DX_PT = 0.0
LETTER_DY_PT = 7.0
LETTER_SIZE = 11
DPI = 300


def colours_against_the_conventions(inputs):
    """Every colour of this module that the conventions file names is compared with the file."""
    text = inputs.text(CONVENTIONS)
    want = [("Delphy", BLUE, r"Delphy: blue `(#[0-9a-fA-F]{6})`"),
            ("ribbon alpha", str(BLUE_ALPHA), r"ribbons: same colour, alpha ([0-9.]+)"),
            ("BEAST X", GREEN, r"BEAST X: green `(#[0-9a-fA-F]{6})`"),
            ("Ituri", ORANGE_ITURI, r"Ituri `(#[0-9a-fA-F]{6})`"),
            ("Nord-Kivu", ORANGE_NK, r"Nord-Kivu `(#[0-9a-fA-F]{6})`"),
            ("other provinces", ORANGE_OTHER, r"other provinces `(#[0-9a-fA-F]{6})`"),
            ("smoothing prior alone", GREY_PRIOR, r"light grey `(#[0-9a-fA-F]{6})`"),
            ("declaration", GREY_REF, r"vertical dashed line, colour ([0-9.]+)"),
            ("width of the figure (in)", f"{FIG_WIDTH:.1f}", r"Draw every figure ([0-9.]+) in wide"),
            ("largest height (in)", f"{MAX_HEIGHT:.1f}", r"at most ([0-9.]+) in"),
            ("sizes of the type", str(SIZES), r"apply_figure_style\(sizes=(\([0-9, ]+\))\)"),
            ("smallest size (pt)", str(SIZES[2]), r"below (\d+) pt"),
            ("dpi", str(DPI), r"Save PNG at (\d+) dpi")]
    out = []
    for what, mine, pat in want:
        m = re.search(pat, text)
        found = m.group(1) if m else ""
        out.append({"what": what, "module": mine, "conventions": found,
                    "equal": bool(m) and found.lower() == mine.lower()})
    return out


def style():
    """The settings of apply_figure_style(sizes=(10, 9, 8)), as the earlier common module sets them for a
    stand-alone run (post_fig_common_20260924.py, _fallback_style and style)."""
    import matplotlib.font_manager as fm
    have = {f.name for f in fm.fontManager.ttflist}
    font = FONT if FONT in have else None
    base, secondary, tick = SIZES
    rc = {
        "font.family": "sans-serif", "font.size": base, "axes.labelsize": base,
        "axes.titlesize": base, "legend.fontsize": secondary,
        "xtick.labelsize": tick, "ytick.labelsize": tick, "axes.linewidth": 0.6,
        "xtick.direction": "out", "ytick.direction": "out",
        "xtick.major.size": 3, "ytick.major.size": 3,
        "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "axes.spines.top": False, "axes.spines.right": False,
        "axes.spines.left": True, "axes.spines.bottom": True,
        "axes.grid": False, "legend.frameon": False, "figure.dpi": 200,
        "savefig.dpi": DPI, "savefig.bbox": "tight", "axes.titleweight": "normal",
        "axes.titlelocation": "left", "axes.labelweight": "normal", "lines.linewidth": 1.2,
        "patch.linewidth": 0.6, "pdf.fonttype": 42, "ps.fonttype": 42,
        "figure.facecolor": "white", "savefig.facecolor": "white", "axes.facecolor": "white",
        "hatch.linewidth": 0.5, "mathtext.fontset": "custom" if font else "dejavusans",
        "svg.hashsalt": "main_figures_20261004",
    }
    if font:
        rc["font.sans-serif"] = [font, "DejaVu Sans"]
        rc.update({"mathtext.rm": font, "mathtext.it": f"{font}:italic", "mathtext.bf": f"{font}:bold",
                   "mathtext.cal": font, "mathtext.sf": font, "mathtext.tt": font})
    mpl.rcParams.update(rc)
    return font


def blend(colour, alpha, background="#ffffff"):
    c = np.array(mpl.colors.to_rgb(colour))
    b = np.array(mpl.colors.to_rgb(background))
    return mpl.colors.to_hex(alpha * c + (1 - alpha) * b)


def new_figure(height, width=FIG_WIDTH):
    if height > MAX_HEIGHT + 1e-9:
        raise Stop(f"height {height:.2f} in exceeds the maximum of {MAX_HEIGHT} in")
    fig = plt.figure(figsize=(width, height))
    frame = Rectangle((0, 0), 1, 1, transform=fig.transFigure, facecolor="none", edgecolor="none",
                      linewidth=0, zorder=-100)
    frame.set_gid("full-figure-frame")
    fig.add_artist(frame)
    return fig


def panel_head(ax, letter, title=None, dx_pt=LETTER_DX_PT, dy_pt=LETTER_DY_PT):
    out = []
    if letter:
        out.append(ax.annotate(letter.lower(), xy=(0, 1), xycoords="axes fraction",
                               xytext=(-dx_pt, dy_pt), textcoords="offset points", ha="left",
                               va="bottom", fontweight="bold", fontsize=LETTER_SIZE,
                               annotation_clip=False))
    if title:
        out.append(ax.annotate(title, xy=(0, 1), xycoords="axes fraction",
                               xytext=((-dx_pt + 15.0) if letter else -dx_pt, dy_pt),
                               textcoords="offset points", ha="left", va="bottom",
                               fontsize=SIZES[0], annotation_clip=False))
    return out


# ---------------------------------------------------------------------------------------------
# file of values
# ---------------------------------------------------------------------------------------------
class Values:
    """One row per element drawn or printed. value, lower and upper hold the TEXT of cells; the column
    source_columns says which cell stands in which column ('value=median; lower=hpd95_lower; ...').
    The column delivered holds 'yes' or 'no' in EVERY row: 'no' where the table of estimates of the analysis
    of the row is not among the inputs, 'yes' in every other row, also in a row that names no analysis (a limit
    by rule, a heading, an entry of a key, words). In earlier versions the cell of a row without an analysis
    was empty."""

    def __init__(self, figure, inputs):
        self.figure, self.inputs, self.rows = figure, inputs, []

    def add(self, panel, element, label="", value="", lower="", upper="", x="", x_end="", unit="",
            drawn_as="", table="", row="", columns="", note="", analysis="", quantity_id="",
            delivered="", meets=""):
        def s(v):
            if v is None:
                return ""
            if isinstance(v, (dt.date, dt.datetime)):
                return v.strftime("%Y-%m-%d")
            if isinstance(v, float):
                return repr(v)
            return str(v)
        version = self.inputs.version(table) if table else ""
        self.rows.append({
            "figure": self.figure, "panel": panel, "element": element, "label": label,
            "value": s(value), "lower": s(lower), "upper": s(upper), "x": s(x), "x_end": s(x_end),
            "unit": unit, "drawn_as": drawn_as, "source_table": table, "source_version_id": version,
            "source_row_0based": s(row), "source_columns": columns, "note": note,
            "analysis": analysis, "quantity_id": quantity_id, "delivered": delivered or "yes",
            "meets_convergence_criterion": meets})
        return self.rows[-1]

    def save(self, path):
        write_csv(path, VALUE_COLUMNS, self.rows)
        return path


def write_csv(path, columns, rows):
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=columns, lineterminator="\n", quoting=csv.QUOTE_MINIMAL)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in columns})


def read_csv(path):
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def read_back(values_path, inputs):
    """Reads every cell that a row of the file of values names AGAIN from the table and compares it with
    value, lower, upper (and x, x_end where the row names a cell for them). Returns the counts."""
    rows = read_csv(values_path)
    out = {"rows": len(rows), "rows_with_a_source": 0, "rows_whose_cells_are_equal": 0, "cells_compared": 0,
           "rows_without_a_source": 0, "not_equal": []}
    for k, r in enumerate(rows):
        if not r["source_table"]:
            out["rows_without_a_source"] += 1
            if r["value"] or r["lower"] or r["upper"]:
                out["not_equal"].append((k, "a row without a source holds a value"))
            continue
        out["rows_with_a_source"] += 1
        ok = True
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
                cell = t.cell(int(r["source_row_0based"]), col)
                out["cells_compared"] += 1
                if cell != r[part]:
                    ok = False
                    out["not_equal"].append((k, f"{part}: file of values {r[part]!r}, cell {cell!r}"))
        else:                                   # a document: the text that the pattern finds in the line
            line = inputs.lines(r["source_table"])[int(r["source_row_0based"])]
            for part, pat in named.items():
                m = re.search(pat[len("pattern "):], line) if pat.startswith("pattern ") else None
                out["cells_compared"] += 1
                if not m or m.group(1) != r[part]:
                    ok = False
                    out["not_equal"].append((k, f"{part}: the pattern does not find {r[part]!r} in the line"))
        out["rows_whose_cells_are_equal"] += ok
    return out


# ---------------------------------------------------------------------------------------------
# words that must not appear (R4 j)
# ---------------------------------------------------------------------------------------------
def forbidden_of_the_earlier_module(inputs):
    """The list FORBIDDEN, read from the SOURCE of the earlier common module (the module is not run)."""
    tree = ast.parse(inputs.text(EARLIER_COMMON))
    for node in tree.body:
        if isinstance(node, ast.Assign) and any(getattr(t, "id", "") == "FORBIDDEN" for t in node.targets):
            return [str(x) for x in ast.literal_eval(node.value)], node.lineno, node.end_lineno
    raise Stop(f"{EARLIER_COMMON}: no list FORBIDDEN found")


def internal_names(inputs):
    """Names of analyses and of genome sets of this round, read from the schedules, the table of inputs, the
    specification and the tables. Returns {kind: sorted list}."""
    analyses, sets, tracks = set(), set(), set()
    s3 = inputs.table(SCHEDULE_PHASE3).need("analysis", "genome_set", "track")
    for r in s3.rows:
        analyses.add(r["analysis"])
        sets.add(r["genome_set"].split(" ")[0])
        tracks.add(r["track"])
    s5 = inputs.table(SCHEDULE_PHASE5).need("analysis", "counterpart_in_the_round_of_24_September")
    for r in s5.rows:
        analyses.add(r["analysis"])
    for r in inputs.inputs.rows:
        if r["analysis"]:
            analyses.add(r["analysis"])
        if r["what_it_is"].startswith("table of estimates") and inputs.has(r["file"]):
            t = inputs.table(r["file"])
            for q in t.rows:
                analyses.add(q["analysis"])
                sets.add(q["genome_set"])
        for m in re.finditer(r"track '([^']+)'", r["state_of_the_file"] + " " + r["what_it_is"]):
            tracks.add(m.group(1))
    sp = inputs.table(SPECIFICATION).need("track", "key_analysis")
    for r in sp.rows:
        for x in r["track"].split(";"):
            if x.strip():
                tracks.add(x.strip())
        for x in r["key_analysis"].split(";"):
            if x.strip():
                analyses.add(x.strip())
    tracks.add("Main figures")
    return {"analyses": sorted(a for a in analyses if a), "sets": sorted(s for s in sets if s),
            "tracks": sorted(t for t in tracks if t)}


def word_patterns(inputs):
    """[(pattern, kind)] of R4 (j): the list FORBIDDEN and the internal names of this round."""
    fb, _, _ = forbidden_of_the_earlier_module(inputs)
    pats = [(p, "FORBIDDEN of the earlier common module") for p in fb]
    names = internal_names(inputs)
    for a in names["analyses"]:
        pats.append((r"(?<![A-Za-z0-9_])" + re.escape(a) + r"(?![A-Za-z0-9_])", "name of an analysis"))
    for s in names["sets"]:
        pats.append((r"(?<![A-Za-z0-9_])" + re.escape(s) + r"(?![A-Za-z0-9_])", "name of a genome set"))
    for t in names["tracks"]:
        pats.append((r"(?<![A-Za-z0-9_])" + re.escape(t) + r"(?![A-Za-z0-9_])", "name of a track"))
    return pats


def words_found(text, patterns):
    hits = []
    for pat, kind in patterns:
        for m in re.finditer(pat, text):
            hits.append({"pattern": pat, "kind": kind, "found": m.group(0),
                         "context": text[max(0, m.start() - 30): m.end() + 30].replace("\n", " ")})
    return hits


# ---------------------------------------------------------------------------------------------
# checks of a figure (conventions, 'Verification before saving')
# ---------------------------------------------------------------------------------------------
def visible_texts(fig):
    out = []
    for t in fig.findobj(mtext.Text):
        s = t.get_text()
        if not s or not s.strip() or not t.get_visible():
            continue
        ax = t.axes
        if ax is not None and not ax.get_visible():
            continue
        out.append(t)
    return out


def check_figure(fig, patterns, min_pt=8.0, allow_overlap=(), banners=()):
    """Sizes, texts outside the figure, overlaps of texts, texts on a spine, words of R4 (j), words of a
    banner. banners: the texts of the banners of a test with fixtures or of a preview, which such a picture
    must carry; the figure of the article has none."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    boxes = []
    for t in visible_texts(fig):
        try:
            bb = t.get_window_extent(r)
        except Exception:
            continue
        if bb.width <= 0 or bb.height <= 0:
            continue
        boxes.append((t, bb))
    small = sorted({(round(t.get_fontsize(), 2), t.get_text()[:40]) for t, _ in boxes
                    if t.get_fontsize() < min_pt - 1e-6})
    fb = fig.bbox
    outside = [t.get_text()[:50] for t, bb in boxes
               if bb.x0 < fb.x0 - 0.5 or bb.x1 > fb.x1 + 0.5 or bb.y0 < fb.y0 - 0.5 or bb.y1 > fb.y1 + 0.5]
    overlaps = []
    shrink = 0.6
    for i, (a, ba) in enumerate(boxes):
        for b, bb in boxes[i + 1:]:
            A = mpl.transforms.Bbox.from_extents(ba.x0 + shrink, ba.y0 + shrink, ba.x1 - shrink, ba.y1 - shrink)
            if A.overlaps(bb):
                pair = (a.get_text()[:40], b.get_text()[:40])
                if any(p in pair[0] or p in pair[1] for p in allow_overlap):
                    continue
                overlaps.append(pair)
    spine_hits = []
    for ax in fig.axes:
        if not ax.get_visible():
            continue
        own = set(ax.get_xticklabels(which="both") + ax.get_yticklabels(which="both"))
        for s in ax.spines.values():
            if not s.get_visible():
                continue
            bs = s.get_window_extent(r)
            for t, bt in boxes:
                if t in own or t.axes is not ax:
                    continue
                if bt.overlaps(bs) and bt.width > 0:
                    inner = mpl.transforms.Bbox.from_extents(bt.x0 + 1, bt.y0 + 1, bt.x1 - 1, bt.y1 - 1)
                    if inner.overlaps(bs):
                        spine_hits.append((t.get_text()[:40], s.spine_type))
    words = []
    for t, _ in boxes:
        for h in words_found(t.get_text(), patterns):
            words.append((h["pattern"], h["kind"], t.get_text()[:60]))
    sizes = sorted({round(t.get_fontsize(), 2) for t, _ in boxes})
    banner, in_other_case = [], []
    for t, _ in boxes:
        if t.get_text() in banners:
            continue
        for w in WORDS_OF_A_BANNER:
            if re.search(r"(?<![A-Za-z])" + w + r"(?![A-Za-z])", t.get_text()):
                banner.append((w, t.get_text()[:60]))
            elif re.search(r"(?<![A-Za-z])" + w + r"(?![A-Za-z])", t.get_text(), flags=re.I):
                in_other_case.append((w, t.get_text()[:60]))
    res = {"n_texts": len(boxes),
           "min_fontsize": min((t.get_fontsize() for t, _ in boxes), default=None),
           "font_sizes": sizes, "below_min_size": small, "outside_figure": outside,
           "text_overlaps": overlaps, "text_on_spine": spine_hits, "words_of_R4j": words,
           "words_of_a_banner": banner,
           "words_of_a_banner_in_another_case": in_other_case,
           "texts": sorted({(t.get_text(), round(t.get_fontsize(), 2)) for t, _ in boxes})}
    res["ok"] = not (small or outside or overlaps or spine_hits or words or banner)
    return res


def boxes_overlap_artists(fig, boxes, artists, pad=0.0):
    """Pairs (name of the box, artist) whose extents overlap on the canvas. boxes: [(name, Bbox)]."""
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    hits = []
    for name, bb in boxes:
        B = mpl.transforms.Bbox.from_extents(bb.x0 - pad, bb.y0 - pad, bb.x1 + pad, bb.y1 + pad)
        for what, a in artists:
            try:
                ba = a.get_window_extent(r)
            except Exception:
                continue
            if ba.width <= 0 and ba.height <= 0:
                continue
            A = mpl.transforms.Bbox.from_extents(ba.x0, ba.y0, max(ba.x1, ba.x0 + 0.01), max(ba.y1, ba.y0 + 0.01))
            if B.overlaps(A):
                hits.append((name, what))
    return hits


def save_figure(fig, stem, outdir="."):
    """PNG (300 dpi) and PDF, exactly FIG_WIDTH inches wide (the 'tight' box is the full figure)."""
    png = os.path.join(outdir, stem + ".png")
    pdf = os.path.join(outdir, stem + ".pdf")
    fig.savefig(png, dpi=DPI, bbox_inches="tight", pad_inches=0, facecolor="white",
                metadata={"Software": None})
    fig.savefig(pdf, bbox_inches="tight", pad_inches=0, facecolor="white",
                metadata={"CreationDate": None, "ModDate": None, "Producer": None, "Creator": None})
    from PIL import Image
    with Image.open(png) as im:
        w, h = im.size
    info = {"png": png, "pdf": pdf, "width_px": w, "height_px": h, "width_in": w / DPI,
            "height_in": h / DPI, "dpi": DPI}
    if abs(w - FIG_WIDTH * DPI) > 2:
        raise Stop(f"saved width {w} px is not {FIG_WIDTH} in at {DPI} dpi: an element extends beyond the figure")
    if h / DPI > MAX_HEIGHT + 0.01:
        raise Stop(f"saved height {h / DPI:.2f} in exceeds {MAX_HEIGHT} in")
    return info


# ---------------------------------------------------------------------------------------------
# rows of a forest plot (as the earlier common module lays them out)
# ---------------------------------------------------------------------------------------------
class RowLayout:
    def __init__(self, heading_gap=0.45):
        self.items, self.y, self.heading_gap = [], 0.0, heading_gap

    def heading(self, text, **kw):
        if self.items:
            self.y += self.heading_gap
        self.items.append(dict(kind="heading", y=self.y, text=text, **kw))
        self.y += 1.0
        return self.items[-1]

    def row(self, key, label, **kw):
        self.items.append(dict(kind="row", y=self.y, key=key, label=label, **kw))
        self.y += 1.0
        return self.items[-1]

    @property
    def units(self):
        return self.y

    def rows(self):
        return [i for i in self.items if i["kind"] == "row"]


def pitch_for(units, top_in, bottom_in, preferred=0.17, floor=0.145, max_height=MAX_HEIGHT):
    p = min(preferred, (max_height - top_in - bottom_in) / units)
    if p < floor - 1e-9:
        raise Stop(f"{units:.1f} lines do not fit into {max_height} in at {floor} in per line")
    return p


def draw_row_labels(fig, ax, layout, x_in, indent_in=0.12):
    tr = mpl.transforms.blended_transform_factory(fig.dpi_scale_trans, ax.transData)
    out = []
    for it in layout.items:
        if it["kind"] == "heading":
            out.append(ax.text(x_in, it["y"], it["text"], transform=tr, ha="left", va="center",
                               fontsize=SIZES[1], fontweight="bold", color="black", clip_on=False))
        else:
            out.append(ax.text(x_in + indent_in * (0 if it.get("flush") else 1), it["y"], it["label"],
                               transform=tr, ha="left", va="center", fontsize=it.get("size", SIZES[2]),
                               fontweight="bold" if it.get("bold") else "normal",
                               color=it.get("label_colour", "black"), clip_on=False))
    return out


def interval(ax, y, m, lo, hi, colour, ms=4.2, lw=1.1, filled=True, zorder=3, mew=1.0):
    """Horizontal interval with a symbol at the median (earlier common module, function interval). The axis
    holds every bound by its rule, so that nothing is cut."""
    ax.plot([lo, hi], [y, y], color=colour, lw=lw, solid_capstyle="butt", zorder=zorder)
    ax.plot([m], [y], marker="o", ms=ms, mfc=colour if filled else "white", mec=colour, mew=mew,
            zorder=zorder + 0.2, ls="none")


# ---------------------------------------------------------------------------------------------
# rules for the limits of axes (R4 h); every rule is listed in the mapping
# ---------------------------------------------------------------------------------------------
def steps_range(values, step):
    """Smallest range with limits at multiples of `step` that holds every value."""
    lo = math.floor(min(values) / step + 1e-9) * step
    hi = math.ceil(max(values) / step - 1e-9) * step
    return round(lo, 10), round(hi, 10)


def decade_floor(x):
    return 10.0 ** math.floor(math.log10(x) + 1e-12)


def decade_ceil(x):
    return 10.0 ** math.ceil(math.log10(x) - 1e-12)


def ticks_between(lo, hi, step):
    n = int(round((hi - lo) / step))
    return [round(lo + i * step, 10) for i in range(n + 1)]


def json_dump(obj, path):
    def default(o):
        if isinstance(o, (dt.date, dt.datetime)):
            return o.isoformat()
        if isinstance(o, (np.integer,)):
            return int(o)
        if isinstance(o, (np.floating,)):
            return float(o)
        if isinstance(o, (set, tuple)):
            return list(o)
        return str(o)
    with open(path, "w", encoding="utf-8") as fh:
        json.dump(obj, fh, indent=1, ensure_ascii=False, sort_keys=True, default=default)
        fh.write("\n")
