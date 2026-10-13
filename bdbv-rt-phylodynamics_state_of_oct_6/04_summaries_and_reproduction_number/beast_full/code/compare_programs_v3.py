#!/usr/bin/env python3
"""compare_programs_v3.py - compare_programs_v2.py with FIVE changes, made on 5 October 2026 before any table of the comparison
was saved:
  1. an input table can hold several analyses: the rows of a side are those whose column analysis equals the key
     'analysis_in_the_table' of the side (default: its key 'analysis'); compare_programs_v2.py read every row of a table;
  2. the input tables may name a quantity differently ('reproduction number, P3' with end_used 'own' in one table,
     'reproduction number P3 (own end)' in another): a side can give, in 'quantity_labels', the label of its table for a quantity;
     the label that was read stands in the column source_row, so that the row can be found in the table as it is;
  3. a pair can be computed for two summaries of the analysis of BEAST X (continued chains): the key 'summary' of a pair is
     written to the new column summary_of_the_analysis_of_BEAST_X;
  4. the key 'pair_in_the_form' of a pair is written to the new column pair_in_the_form_of_26_September;
  5. the key 'note_on_the_pair' of a pair is written to the new column note_on_the_pair.
The three new columns stand at the end of both tables.  No value is computed otherwise than in compare_programs_v2.py
(validate_compare_v3.sh compares the two scripts on input that needs none of the changes).
What compare_programs_v2.py was: compare_programs.py with ONE change, made on 4 October 2026 for the continued chains: both
tables have one further column at their end, note_on_the_analysis, copied from the row of the table of
estimates (table of values) or of the table of intervals (table of trajectories) that the row was read from; it is empty where
that table has no such column or its cell is empty.  It marks the rows of an analysis whose chains were continued from a saved
state.  No other column and no value is changed (validate_compare_v2.sh compares the two scripts on the same input).
What compare_programs.py was: tables of the comparison of the two programs, in the form of the tables
program_comparison_values_20260926.csv and program_comparison_trajectories_20260926.csv (their columns and rows are
the form; none of their values is read here).

Input: a file of pairs (JSON).  Every pair names two analyses, one of Delphy and one of BEAST X, each with its table of
estimates and, for the Skygrid, its table of intervals, both in the common form.  The script reads the tables and computes,
for every pair, program and quantity:
  whether the median of this program lies inside the 95 % interval of the other program, and the reverse;
  the difference of the medians (BEAST X minus Delphy) and that difference as a share of the width of the interval of Delphy.
It COMPUTES ONLY.  The reading rule of the comparison is not applied here: the columns of the form that hold a reading
(rule_3 and the two columns on when the rule was fixed) say so in every row.
A quantity that one of the two tables does not hold, or holds as 'not evaluated', gets a row that says so.

Usage: compare_programs_v3.py --pairs pairs.json --periods periods.csv --out-values V.csv --out-trajectories T.csv
"""
import argparse
import json

import pandas as pd

NOT_READ = "not applied here: the value is computed, not read"
VALUE_COLUMNS = ["pair", "program", "analysis", "what_differs_from_the_primary_analysis_of_Delphy", "differs_within_the_pair_in", "genomes",
                 "last_collection_date", "chains_run", "chains_used", "states_used", "quantity", "kind", "period_start", "period_end", "median",
                 "hpd_lo", "hpd_hi", "median_date", "hpd_early_date", "hpd_late_date", "unit", "ess_pooled", "rhat_split",
                 "meets_criteria_of_convergence", "source_table", "source_version", "source_row", "median_inside_interval_of_the_other_program",
                 "median_of_the_other_program_inside_this_interval", "rule_3", "rule_fixed_before_values_seen",
                 "interim_values_displayed_to_the_track_before_the_rule", "difference_of_medians_BEAST_X_minus_Delphy",
                 "difference_BEAST_X_minus_Delphy_as_share_of_the_width_of_the_Delphy_interval", "state_of_the_table",
                 "period", "status_of_the_period", "end_used", "period_ends_on_the_same_date_in_both_programs",
                 "ess_pooled_of_the_reproduction_number_series", "rhat_split_of_the_reproduction_number_series",
                 "meets_convergence_criterion_by_the_reproduction_number_series", "states_with_R_set_to_0", "note_on_the_analysis",
                 "summary_of_the_analysis_of_BEAST_X", "pair_in_the_form_of_26_September", "note_on_the_pair"]
TRAJ_COLUMNS = ["pair", "program", "analysis", "genomes", "last_collection_date", "chains_used", "states_used", "n_parameters", "interval",
                "older_edge_exact", "newer_edge_exact", "lnNe_median", "lnNe_hpd_lo", "lnNe_hpd_hi", "Ne_tau_years_median", "Ne_tau_years_hpd_lo",
                "Ne_tau_years_hpd_hi", "state_of_the_table", "older_edge", "newer_edge", "mid_point", "ess_pooled", "rhat_split",
                "meets_convergence_criterion", "source_table", "source_version", "note_on_the_analysis",
                "summary_of_the_analysis_of_BEAST_X", "pair_in_the_form_of_26_September", "note_on_the_pair"]


def num(x):
    try:
        v = float(x)
        return v if v == v else None
    except Exception:
        return None


def row_of(E, quantity, end_used=None):
    d = E[E["quantity"] == quantity]
    if end_used is not None and len(d) > 1:
        d = d[d["end_used"] == end_used]
    return d.iloc[0] if len(d) == 1 else None


def rows_of_the_analysis(s, key):
    """the rows of the table that belong to the analysis of the side; None where the side names no table"""
    if not s.get(key):
        return None
    t = pd.read_csv(s[key], dtype=str, keep_default_na=False)
    name = s.get("analysis_in_the_table", s["analysis"])
    if "analysis" in t.columns:
        t = t[t["analysis"] == name].reset_index(drop=True)
        assert len(t), (s[key], name)
    return t


def main():
    import math
    ap = argparse.ArgumentParser()
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--periods", required=True)
    ap.add_argument("--out-values", required=True)
    ap.add_argument("--out-trajectories", required=True)
    a = ap.parse_args()
    pairs = json.load(open(a.pairs))["pairs"]
    per = pd.read_csv(a.periods, dtype=str, keep_default_na=False).set_index("period")
    vals, traj = [], []
    for p in pairs:
        side = {}
        for prog in ("Delphy", "BEAST X"):
            s = p[prog]
            side[prog] = dict(s, E=rows_of_the_analysis(s, "estimates"), I=rows_of_the_analysis(s, "intervals"))
        state = "complete" if all(side[x]["E"] is not None for x in side) else "not complete: the table of estimates of " + " and ".join(
            f"{x} ({side[x]['analysis']})" for x in side if side[x]["E"] is None) + " was not available"
        if p["model"] == "skygrid":
            quantities = [(f"R, {k} ({per.loc[k, 'start']} - {per.loc[k, 'end'] if per.loc[k, 'end'] != 'END' else 'end of the last period'})", "R",
                           f"reproduction number, {k}", k) for k in ("P1", "P2", "P3", "P4", "P5", "P6")]
            quantities += [("evolutionary rate", "rate", "evolutionary rate", ""), ("tMRCA", "date", "date of the most recent common ancestor (decimal year)", "")]
        else:
            quantities = [("evolutionary rate", "rate", "evolutionary rate", ""), ("tMRCA", "date", "date of the most recent common ancestor (decimal year)", ""),
                          ("growth rate, constant exponential growth", "growth", "growth rate, constant exponential growth", ""),
                          ("doubling time, constant exponential growth", "doubling", "doubling time, constant exponential growth", ""),
                          ("R, constant exponential growth", "R_constant", "reproduction number, constant exponential growth", "")]
        for label, kind, q, k in quantities:
            got = {}
            for prog in side:
                E = side[prog]["E"]
                lab = side[prog].get("quantity_labels", {})
                r = None if E is None else row_of(E, lab.get(q, q), ("own" if k in ("P3", "P6") else None))
                rd = None if (E is None or kind != "date") else row_of(E, lab.get("date of the most recent common ancestor (date)", "date of the most recent common ancestor (date)"))
                got[prog] = (r, rd)
            same_end = ""
            if k:
                ends = [got[x][0]["period_end"] for x in got if got[x][0] is not None]
                same_end = ("yes" if len(ends) == 2 and ends[0] == ends[1] else ("no" if len(ends) == 2 else ""))
            for prog, other in (("Delphy", "BEAST X"), ("BEAST X", "Delphy")):
                r, rd = got[prog]
                ro, _ = got[other]
                s = side[prog]
                E = s["E"]
                out = dict(pair=p["pair"], program=prog, analysis=s["analysis"],
                           what_differs_from_the_primary_analysis_of_Delphy=s.get("what_differs_from_the_primary_analysis_of_Delphy", ""),
                           differs_within_the_pair_in=p.get("differs_within_the_pair_in", ""), quantity=label, kind=kind, period=k,
                           status_of_the_period=(per.loc[k, "status"] if k else ""), source_table=s.get("estimates_name", ""),
                           source_version=s.get("estimates_version", ""), rule_3=NOT_READ, rule_fixed_before_values_seen=NOT_READ,
                           interim_values_displayed_to_the_track_before_the_rule=NOT_READ, state_of_the_table=state,
                           period_ends_on_the_same_date_in_both_programs=same_end,
                           summary_of_the_analysis_of_BEAST_X=p.get("summary", ""), pair_in_the_form_of_26_September=p.get("pair_in_the_form", ""), note_on_the_pair=p.get("note_on_the_pair", ""))
                ql = s.get("quantity_labels", {}).get(q, q)
                if r is None:
                    out.update(median="", meets_criteria_of_convergence=("not assessed: the table of estimates was not available" if E is None
                                                                         else "not assessed: the table of estimates holds no row of this quantity"),
                               source_row=f"quantity = {ql}")
                    vals.append(out)
                    continue
                m, lo, hi = num(r["median"]), num(r["hpd95_lower"]), num(r["hpd95_upper"])
                out.update(genomes=r["genomes"], last_collection_date=r["most_recent_collection_date"], chains_run=r["chains_run"],
                           chains_used=r["chains_retained"], states_used=r["states_used"], period_start=r["period_start"], period_end=r["period_end"],
                           end_used=r["end_used"], unit=r["unit"], ess_pooled=r["ess_pooled"], rhat_split=r["rhat_split"],
                           meets_criteria_of_convergence=r["meets_convergence_criterion"], source_row=f"quantity = {ql}" + (f"; end_used = {r['end_used']}" if k else ""),
                           ess_pooled_of_the_reproduction_number_series=r.get("ess_pooled_of_the_reproduction_number_series", ""),
                           rhat_split_of_the_reproduction_number_series=r.get("rhat_split_of_the_reproduction_number_series", ""),
                           meets_convergence_criterion_by_the_reproduction_number_series=r.get("meets_convergence_criterion_by_the_reproduction_number_series", ""),
                           states_with_R_set_to_0=r.get("states_with_R_set_to_0", ""), note_on_the_analysis=r.get("note_on_the_analysis", ""))
                if m is None:
                    out.update(median=r["median"], hpd_lo=r["hpd95_lower"], hpd_hi=r["hpd95_upper"])
                else:
                    out.update(median=m, hpd_lo=lo, hpd_hi=hi)
                if rd is not None:
                    out.update(median_date=rd["median"], hpd_early_date=rd["hpd95_lower"], hpd_late_date=rd["hpd95_upper"])
                if ro is not None:
                    mo, loo, hio = num(ro["median"]), num(ro["hpd95_lower"]), num(ro["hpd95_upper"])
                    if m is not None and loo is not None and hio is not None:
                        out["median_inside_interval_of_the_other_program"] = "yes" if loo <= m <= hio else "no"
                    if mo is not None and lo is not None and hi is not None:
                        out["median_of_the_other_program_inside_this_interval"] = "yes" if lo <= mo <= hi else "no"
                    b, d = (got["BEAST X"][0], got["Delphy"][0])
                    mb, md = num(b["median"]), num(d["median"])
                    dl, dh = num(d["hpd95_lower"]), num(d["hpd95_upper"])
                    if prog == "BEAST X" and mb is not None and md is not None:
                        out["difference_of_medians_BEAST_X_minus_Delphy"] = mb - md
                        if dl is not None and dh is not None and dh > dl:
                            out["difference_BEAST_X_minus_Delphy_as_share_of_the_width_of_the_Delphy_interval"] = (mb - md) / (dh - dl)
                vals.append(out)
        if p["model"] == "skygrid":
            for prog in ("Delphy", "BEAST X"):
                s = side[prog]
                I, E = s["I"], s["E"]
                if I is None:
                    traj.append(dict(pair=p["pair"], program=prog, analysis=s["analysis"],
                                     state_of_the_table="not complete: the table of intervals of this analysis was not available",
                                     summary_of_the_analysis_of_BEAST_X=p.get("summary", ""), pair_in_the_form_of_26_September=p.get("pair_in_the_form", ""), note_on_the_pair=p.get("note_on_the_pair", "")))
                    continue
                K = len(I)
                mids = [num(x) for x in I["mid_point_decimal"]]
                D = mids[0] - mids[1]
                for _, r in I.iterrows():
                    i = int(r["interval"])
                    newer = mids[0] + D / 2 - (i - 1) * D
                    older = newer - D if i < K else ""
                    traj.append(dict(pair=p["pair"], program=prog, analysis=s["analysis"], genomes=(E["genomes"].iloc[0] if E is not None else ""),
                                     last_collection_date=(E["most_recent_collection_date"].iloc[0] if E is not None else ""),
                                     chains_used=r["chains_retained"], states_used=r["states_used"], n_parameters=K, interval=i,
                                     older_edge_exact=older, newer_edge_exact=newer, lnNe_median=r["ln_Ne_tau_median"], lnNe_hpd_lo=r["ln_Ne_tau_hpd95_lower"],
                                     lnNe_hpd_hi=r["ln_Ne_tau_hpd95_upper"], Ne_tau_years_median=math.exp(float(r["ln_Ne_tau_median"])),
                                     Ne_tau_years_hpd_lo=math.exp(float(r["ln_Ne_tau_hpd95_lower"])), Ne_tau_years_hpd_hi=math.exp(float(r["ln_Ne_tau_hpd95_upper"])),
                                     state_of_the_table="complete", older_edge=r["older_edge"], newer_edge=r["newer_edge"], mid_point=r["mid_point"],
                                     ess_pooled=r["ess_pooled"], rhat_split=r["rhat_split"], meets_convergence_criterion=r["meets_convergence_criterion"],
                                     source_table=s.get("intervals_name", ""), source_version=s.get("intervals_version", ""),
                                     note_on_the_analysis=r.get("note_on_the_analysis", ""),
                                     summary_of_the_analysis_of_BEAST_X=p.get("summary", ""), pair_in_the_form_of_26_September=p.get("pair_in_the_form", ""), note_on_the_pair=p.get("note_on_the_pair", "")))
    V = pd.DataFrame(vals).reindex(columns=VALUE_COLUMNS).fillna("")
    V.to_csv(a.out_values, index=False)
    T = pd.DataFrame(traj).reindex(columns=TRAJ_COLUMNS).fillna("")
    T.to_csv(a.out_trajectories, index=False)
    print("pairs:", len(pairs), "| rows of the table of values:", len(V), "| rows of the table of trajectories:", len(T))
    print(V.groupby(["pair", "summary_of_the_analysis_of_BEAST_X", "program"]).size().to_string())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
