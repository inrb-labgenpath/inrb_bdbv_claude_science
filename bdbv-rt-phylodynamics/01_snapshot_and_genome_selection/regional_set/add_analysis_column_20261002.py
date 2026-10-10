#!/usr/bin/env python3
"""add_analysis_column_20261002.py - writes a table of the earlier form with the name of the analysis of the schedule in a
first column 'analysis_of_the_schedule' (every table of the earlier form of this round has this column in front).
usage: add_analysis_column_20261002.py <analysis> <table in> <table out>"""
import sys
import pandas as pd
analysis, src, dst = sys.argv[1:4]
d = pd.read_csv(src, dtype=str, keep_default_na=False, low_memory=False)
assert 'analysis_of_the_schedule' not in d.columns
d.insert(0, 'analysis_of_the_schedule', analysis)
d.to_csv(dst, index=False)
print(dst, d.shape)
