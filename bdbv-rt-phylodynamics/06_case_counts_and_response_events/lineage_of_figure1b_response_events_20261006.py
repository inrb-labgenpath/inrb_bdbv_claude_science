import pandas as pd
import re

f0 = "INTERNAL_figure1b_response_events_20261005.csv"
f1 = "figure1b_response_events_20261006.csv"

raw = open("inputs/INTERNAL_figure1b_response_events_20261005.csv", encoding="utf-8").read()
E = pd.read_csv("inputs/INTERNAL_figure1b_response_events_20261005.csv", dtype=str, keep_default_na=False)
assert len(E) == 9

i = E.index[E.order == "8"]
assert len(i) == 1
old = "MSF's treatment centre in Mongbwalu opens; the day is not known and lies between these dates"
new = "MSF's treatment centre in Mongbwalu opens"
assert E.loc[i[0], "words_for_the_legend"] == old and raw.count(old) == 1

E2 = E.copy()
E2.loc[i[0], "words_for_the_legend"] = new
assert int((E2 != E).sum().sum()) == 1
raw2 = raw.replace(old, new)
open(f1, "w", encoding="utf-8", newline="").write(raw2)
B = pd.read_csv(f1, dtype=str, keep_default_na=False)
assert B.equals(E2) and list(B.columns) == list(E.columns)