#!/usr/bin/env python3
"""combine_sensitivity_20261007.py (snapshot of 7 October 2026)

The analyses are those of the queue; the quantity P6_to_END is one more period wherever an analysis has it. The table of
the comparison of the removal sets is written by a script of its own, not by this one (the code for it is not in this
file, although its head below still lists it).

combine_sensitivity_20261007.py - joins the tables of the analyses that are derived into the common tables, in the
order of the schedule, and writes further tables from them:
  further_probabilities_sensitivity_20261007.csv   per analysis, from its posterior file: P(R > 1) and P(R < 1) of every period, and
                                                   P(R of period j < R of period i) for every two periods (share of the states)
  across_analyses_sensitivity_20261007.csv         per group of analyses and quantity: lowest and highest median with the analysis that
                                                   gives it, and how many 95 % intervals lie below 1, above 1 or hold 1 (numbers only)
  removal_sets_comparison_sensitivity_20261007.csv quantities (a) to (d) for every set and run,
                                                   v (largest absolute difference between the first and the second run of the same set),
                                                   the differences between every removal set and each of its controls; no outcome is stated
usage: python combine_sensitivity_20261007.py SCHEDULE.csv   (from the working directory)
"""
import os
import sys

import numpy as np
import pandas as pd

LABEL = '20261007'
D = 'delivery'
GROUP13 = ['allClasses', 'minusExcluded1046', 'retained1046', 'dupNone', 'dupKeepBest', 'drcOnly', 'upTo0831', 'upTo0817', 'gridWeekly32', 'gridCoarse10', 'gridCalendar20',
           'smooth15d', 'smooth60d']
P6E = 'P6_to_END'
QN = {P6E: 'P6_to_END (defined on 7 October 2026 after the collection dates of the new snapshot were seen and before any estimate)'}
PERIODS = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6', P6E]


def col_of(p, what='reproduction_number'):
    return f'{what}_{p}' + ('_own_end' if p in ('P3', 'P6') else '')


def main(schedule_csv):
    sched = pd.read_csv(schedule_csv, dtype=str, keep_default_na=False)
    mine = sched[sched['group'].isin(['sensitivity (the 13 analyses of Figure 2)', 'further analyses of the article', 'removal sets and their controls'])]
    in_queue = set(pd.read_csv(os.path.join('run', 'queue.csv'), dtype=str, keep_default_na=False)['analysis'])
    order = [a for a in dict.fromkeys(mine['analysis']) if a in in_queue]
    assert len(order) == len(in_queue) == 30, (len(order), len(in_queue))
    done = [a for a in order if os.path.exists(os.path.join(D, 'analyses', a, f'estimates_{a}_{LABEL}.csv'))]
    out = {}
    for kind in ('chains', 'estimates', 'intervals', 'trajectories', 'root_date_intervals', 'checks', 'halves', 'recomputation'):
        parts = []
        for a in done:
            p = os.path.join(D, 'analyses', a, f'{kind}_{a}_{LABEL}.csv')
            if os.path.exists(p):
                parts.append(pd.read_csv(p, dtype=str, keep_default_na=False))
        if parts:
            T = pd.concat(parts, ignore_index=True)
            name = {'checks': f'checks_by_analysis_sensitivity_{LABEL}.csv', 'halves': f'halves_of_the_chains_sensitivity_{LABEL}.csv',
                    'recomputation': f'recomputation_sensitivity_{LABEL}.csv'}.get(kind, f'{kind}_sensitivity_{LABEL}.csv')
            T.to_csv(os.path.join(D, name), index=False)
            out[kind] = T
    E = out.get('estimates')
    # ---- further probabilities ----
    rows = []
    for a in done:
        pf = os.path.join(D, 'analyses', a, f'posterior_{a}_{LABEL}.csv.gz')
        if not os.path.exists(pf):
            continue
        P = pd.read_csv(pf)
        e = E[E['analysis'] == a].iloc[0]
        base = dict(analysis=a, program=e['program'], setting=e['setting'], genome_set=e['genome_set'], genomes=e['genomes'], most_recent_collection_date=e['most_recent_collection_date'])
        per_a = [p for p in PERIODS if col_of(p) in P.columns]
        ev = [p for p in per_a if P[col_of(p)].notna().all()]
        for p in per_a:
            for q, f in ((f'P(R of {QN.get(p, p)} > 1)', lambda x: x > 1), (f'P(R of {QN.get(p, p)} < 1)', lambda x: x < 1)):
                rows.append(dict(base, quantity=q, value=(float(f(P[col_of(p)].values).mean()) if p in ev else 'not evaluated'), states_used=len(P), chains_retained=P['seed'].nunique()))
        for i in per_a:
            for j in per_a:
                if i == j:
                    continue
                rows.append(dict(base, quantity=f'P(R of {QN.get(j, j)} < R of {QN.get(i, i)})', value=(float((P[col_of(j)].values < P[col_of(i)].values).mean()) if (i in ev and j in ev) else 'not evaluated'),
                                 states_used=len(P), chains_retained=P['seed'].nunique()))
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(D, f'further_probabilities_sensitivity_{LABEL}.csv'), index=False)
    # ---- across analyses ----
    if E is not None:
        rows = []
        groups = {'the 13 sensitivity analyses of Figure 2': ['delphy_sens_' + x for x in GROUP13],
                  'cells of the coalescent prior': ['delphy_sens_cells1000', 'delphy_sens_cells16000', 'delphy_sens_cells400'],
                  'removal sets and their controls': [a for a in order if a.startswith('delphy_rem_')]}
        for gname, members in groups.items():
            have = [a for a in members if a in done]
            sub = E[E['analysis'].isin(have)]
            qs = [f'reproduction number, {QN.get(p, p)}' for p in PERIODS] + [f'growth rate, {QN.get(p, p)}' for p in PERIODS] + ['evolutionary rate', 'date of the most recent common ancestor (decimal year)', 'kappa',
                                                                                                         'P(R of P2 < R of P1)', 'P(R of P1 > 1)', 'P(R of P2 > 1)', 'P(R of P3 < 1)']
            for q in qs:
                s = sub[sub['quantity'] == q].copy()
                col = 'mean' if q.startswith('P(') else 'median'
                s['v'] = pd.to_numeric(s[col], errors='coerce')
                ok = s[s['v'].notna()]
                r = dict(group=gname, analyses_of_the_group=len(members), analyses_derived=len(have), quantity=q, value_compared=('probability' if col == 'mean' else 'median'),
                         analyses_with_a_value=len(ok), analyses_without_a_value=' '.join(s.loc[s['v'].isna(), 'analysis']))
                if len(ok):
                    lo, hi = ok.loc[ok['v'].idxmin()], ok.loc[ok['v'].idxmax()]
                    r.update(lowest=lo['v'], analysis_with_the_lowest=lo['analysis'], highest=hi['v'], analysis_with_the_highest=hi['analysis'])
                    if q.startswith('reproduction number'):
                        l_, u_ = pd.to_numeric(ok['hpd95_lower'], errors='coerce'), pd.to_numeric(ok['hpd95_upper'], errors='coerce')
                        r.update(intervals_that_lie_below_1=int((u_ < 1).sum()), analyses_whose_interval_lies_below_1=' '.join(ok.loc[u_ < 1, 'analysis']),
                                 intervals_that_lie_above_1=int((l_ > 1).sum()), analyses_whose_interval_lies_above_1=' '.join(ok.loc[l_ > 1, 'analysis']),
                                 intervals_that_hold_1=int(((l_ <= 1) & (u_ >= 1)).sum()))
                    r['analyses_that_do_not_meet_the_convergence_criterion'] = ' '.join(ok.loc[ok['meets_convergence_criterion'] == 'no', 'analysis'])
                rows.append(r)
        pd.DataFrame(rows).fillna('').to_csv(os.path.join(D, f'across_analyses_sensitivity_{LABEL}.csv'), index=False)
    return done, out


if __name__ == '__main__':
    done, _ = main(sys.argv[1])
    print(len(done), 'analyses in the combined tables')
