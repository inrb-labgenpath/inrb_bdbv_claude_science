#!/usr/bin/env python3
"""combine_sensitivity.py - joins the tables of the analyses that are derived into the common tables, in the
order of the schedule, and writes three further tables from them:
  further_probabilities_sensitivity_20261002.csv   per analysis, from its posterior file: P(R > 1) and P(R < 1) of every period, and
                                                   P(R of period j < R of period i) for every two periods (share of the states)
  across_analyses_sensitivity_20261002.csv         per group of analyses and quantity: lowest and highest median with the analysis that
                                                   gives it, and how many 95 % intervals lie below 1, above 1 or hold 1 (numbers only)
  removal_sets_comparison_sensitivity_20261002.csv quantities (a) to (d) for every set and run,
                                                   v (largest absolute difference between the first and the second run of the same set),
                                                   the differences between every removal set and each of its controls; no outcome is stated
usage: python combine_sensitivity.py SCHEDULE.csv   (from the working directory)
"""
import os
import sys

import numpy as np
import pandas as pd

LABEL = '20261002'
D = 'delivery'
GROUP13 = ['allClasses', 'minusExcluded1046', 'retained1046', 'dupNone', 'dupKeepBest', 'drcOnly', 'upTo0831', 'upTo0817', 'gridWeekly32', 'gridCoarse10', 'gridCalendar20',
           'smooth15d', 'smooth60d']
PERIODS = ['P1', 'P2', 'P3', 'P4', 'P5', 'P6']


def col_of(p, what='reproduction_number'):
    return f'{what}_{p}' + ('_own_end' if p in ('P3', 'P6') else '')


def main(schedule_csv):
    sched = pd.read_csv(schedule_csv, dtype=str, keep_default_na=False)
    mine = sched[sched['group'].isin(['sensitivity (the 13 analyses of Figure 2)', 'further analyses of the article', 'removal sets and their controls'])]
    order = list(dict.fromkeys(mine['analysis']))
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
        ev = [p for p in PERIODS if P[col_of(p)].notna().all()]
        for p in PERIODS:
            for q, f in ((f'P(R of {p} > 1)', lambda x: x > 1), (f'P(R of {p} < 1)', lambda x: x < 1)):
                rows.append(dict(base, quantity=q, value=(float(f(P[col_of(p)].values).mean()) if p in ev else 'not evaluated'), states_used=len(P), chains_retained=P['seed'].nunique()))
        for i in PERIODS:
            for j in PERIODS:
                if i == j:
                    continue
                rows.append(dict(base, quantity=f'P(R of {j} < R of {i})', value=(float((P[col_of(j)].values < P[col_of(i)].values).mean()) if (i in ev and j in ev) else 'not evaluated'),
                                 states_used=len(P), chains_retained=P['seed'].nunique()))
    if rows:
        pd.DataFrame(rows).to_csv(os.path.join(D, f'further_probabilities_sensitivity_{LABEL}.csv'), index=False)
    # ---- across analyses ----
    if E is not None:
        rows = []
        groups = {'the 13 sensitivity analyses of Figure 2': ['delphy_sens_' + x for x in GROUP13],
                  'cells of the coalescent prior and replicate': ['delphy_sens_cells1000', 'delphy_sens_cells16000', 'delphy_sens_cells400', 'delphy_sens_replicate'],
                  'removal sets and their controls': [a for a in order if a.startswith('delphy_rem_')]}
        for gname, members in groups.items():
            have = [a for a in members if a in done]
            sub = E[E['analysis'].isin(have)]
            qs = [f'reproduction number, {p}' for p in PERIODS] + [f'growth rate, {p}' for p in PERIODS] + ['evolutionary rate', 'date of the most recent common ancestor (decimal year)', 'kappa',
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
    # ---- removal sets: what the reading rule needs ----
    if E is not None:
        arms = {'B (thinned)': ['buniaThinned_rep1', 'buniaThinned_rep2', 'buniaThinned_rep3'], 'control A (of the thinned sets)': ['randomThinnedA_rep1', 'randomThinnedA_rep2', 'randomThinnedA_rep3'],
                'B (without recent Bunia)': ['withoutRecentBunia'], "control B (of 'without recent Bunia')": ['randomThinnedB_rep1', 'randomThinnedB_rep2', 'randomThinnedB_rep3']}
        Q = {'(a) median R of the last period (P3)': ('reproduction number, P3', 'median'),
             '(b) width of the 95 % HPD interval of R of the last period (P3)': ('reproduction number, P3', 'width'),
             '(c) median R of 1 August to the end (P6)': ('reproduction number, P6', 'median'),
             '(d) posterior SD / prior SD of the growth rate of 1 August to the end (P6)': ('growth rate, P6', 'ratio_of_the_posterior_to_the_prior_sd_of_the_growth_rate')}

        def val(a, q):
            s = E[(E['analysis'] == a) & (E['quantity'] == Q[q][0])]
            if not len(s):
                return np.nan
            s = s.iloc[0]
            if Q[q][1] == 'width':
                return pd.to_numeric(s['hpd95_upper'], errors='coerce') - pd.to_numeric(s['hpd95_lower'], errors='coerce')
            return pd.to_numeric(s[Q[q][1]], errors='coerce')
        rows, comp = [], []
        for arm, sets in arms.items():
            for s_ in sets:
                for run, a in (('first run', f'delphy_rem_{s_}'), ('second run', f'delphy_rem_{s_}_run2')):
                    if run == 'second run' and not arm.startswith('B'):
                        continue
                    e = E[E['analysis'] == a]
                    r = dict(arm=arm, set=s_, run=(run if arm.startswith('B') else 'only run'), analysis=a, derived='yes' if len(e) else 'no',
                             genomes=e['genomes'].iloc[0] if len(e) else '', chains_run=e['chains_run'].iloc[0] if len(e) else '', chains_retained=e['chains_retained'].iloc[0] if len(e) else '')
                    for q in Q:
                        r[q] = val(a, q) if len(e) else np.nan
                    rows.append(r)
        R = pd.DataFrame(rows)
        for q in Q:
            B = R[R['arm'].str.startswith('B')]
            v_by = {}
            for s_, g in B.groupby('set'):
                a1, a2 = g.loc[g['run'] == 'first run', q], g.loc[g['run'] == 'second run', q]
                if len(a1) and len(a2) and not (np.isnan(a1.iloc[0]) or np.isnan(a2.iloc[0])):
                    v_by[s_] = abs(float(a2.iloc[0]) - float(a1.iloc[0]))
            v = max(v_by.values()) if len(v_by) == 4 else np.nan
            for armB, armC in (('B (thinned)', 'control A (of the thinned sets)'), ('B (without recent Bunia)', "control B (of 'without recent Bunia')")):
                b2 = R[(R['arm'] == armB) & (R['run'] == 'second run')]
                b1 = R[(R['arm'] == armB) & (R['run'] == 'first run')]
                c = R[R['arm'] == armC]
                comp.append(dict(quantity=q, removal_sets=armB, controls=armC,
                                 values_of_the_removal_sets_second_run=' '.join(f'{x!r}' for x in b2[q]), values_of_the_removal_sets_first_run=' '.join(f'{x!r}' for x in b1[q]),
                                 values_of_the_controls=' '.join(f'{x!r}' for x in c[q]),
                                 lowest_of_the_removal_sets_second_run=b2[q].min(), highest_of_the_removal_sets_second_run=b2[q].max(),
                                 lowest_of_the_controls=c[q].min(), highest_of_the_controls=c[q].max(),
                                 absolute_difference_between_the_two_runs_by_set=' '.join(f'{k}: {x!r}' for k, x in v_by.items()),
                                 v_largest_absolute_difference_between_the_two_runs_of_a_set_of_arm_B=v, sets_of_arm_B_with_both_runs=len(v_by),
                                 lowest_of_the_controls_minus_v=c[q].min() - v, highest_of_the_controls_plus_v=c[q].max() + v,
                                 differences_removal_set_second_run_minus_control=' '.join(f'{sb}-{sc}: {float(xb) - float(xc)!r}' for sb, xb in zip(b2['set'], b2[q]) for sc, xc in zip(c['set'], c[q]))))
        R.to_csv(os.path.join(D, f'removal_sets_values_sensitivity_{LABEL}.csv'), index=False)
        pd.DataFrame(comp).to_csv(os.path.join(D, f'removal_sets_comparison_sensitivity_{LABEL}.csv'), index=False)
    return done, out


if __name__ == '__main__':
    done, _ = main(sys.argv[1])
    print(len(done), 'analyses in the combined tables')
