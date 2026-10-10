#!/usr/bin/env python3
"""primary_derived_quantities.py - quantities that the article derives from the primary analysis, computed for
the primary analysis of S (Delphy, smoothing fixed, set C00) with the functions of analysis.py and rt_lib.py.

Definitions are those of the units of the article (Table 1 with its caption P024; P089; P090; P022).
Rows are marked by the RULES that the article states, never by a date that it prints:
  pair with the highest median                    largest median of the reproduction number among the pairs of adjacent intervals
                                                  (a) among all pairs, (b) among the pairs that the function skygrid_tables does not
                                                  mark as lying before the lower bound of the date of the root (rule of
                                                  sensitivity_row); both marks are given
  interval with the highest median of Ne tau      largest median of the effective population size among the intervals; the two pairs at
                                                  its edges are the pairs in which this interval is the newer and the older interval
  pairs with an SD ratio of 0.90 or more          posterior SD of the growth rate / prior SD (prior_sd_contrast) >= 0.90
  period of a pair                                the period among P1, P2, P3 that holds the date on which the two intervals meet
  pairs of a period                               the pairs that enter the growth rate of the period with a weight that is not 0 (the
                                                  pairs whose mid-points span a part of the period)
usage: main(out_of_derive_sensitivity, periods_csv, pipeline, outdir)
"""
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

LABEL = '20261002'
ARTICLE = {  # unit of the article, label as it stands there with END in place of a date of the grid, name in analysis.py / rt_lib.py
    'P1': ('T1r1c0', '1 March - 15 May (before the declaration)', 'W1'), 'P2': ('T1r2c0', '15 May - 15 June', 'W2'),
    'P3': ('T1r3c0', '15 June - END (mid-point of the most recent Skygrid interval)', 'W3'), 'P4': ('T1r4c0', '15 June - 15 August', 'W3b'),
    'A1': ('T1r5c0', '1 April - 15 May', 'none in analysis.py (period of the article)'),
    'P5': ('T1r6c0', '15 June - 1 August', 'window 15Jun_to_1Aug of build_tables2.recent_extra'),
    'P6': ('T1r7c0', '1 August - END (mid-point of the most recent Skygrid interval)', 'window 1Aug_to_midpoint_of_most_recent_interval of build_tables2.recent_extra')}


def meets(ess, rhat):
    if ess != ess or rhat != rhat:
        return 'not assessed'
    return 'yes' if (ess >= 200 and rhat <= 1.05) else 'no'


def main(out, periods_csv, pipeline, outdir, analysis_label='delphy_fixed_C00'):
    sys.path.insert(0, pipeline)
    import analysis as A
    import rt_lib as L
    an, tabs = out['analysis'], out['tabs']
    g, gen = an.grid, an.gen
    K = g.K
    lp = an.lnNe()
    tau, fixed = A.precision_of(an)
    assert fixed
    per = pd.read_csv(periods_csv, dtype=str, keep_default_na=False).set_index('period')
    mid0 = float(g.mids[0])
    base = dict(analysis=analysis_label, program='Delphy 1.4.1', setting='Skygrid, smoothing fixed (program default)', genome_set='C00', genomes=an.n_tips,
                most_recent_collection_date=str(an.last_tip.date()))
    n_states = int(sum(len(x) for x in lp))
    tail = dict(chains_run=len(an.chains), chains_retained=len(an.kept), states_used=n_states)
    z = float(norm.ppf(0.975))
    # ---------------- periods of Table 1 ----------------
    defs = {}
    for p in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
        t0 = L.decimal_year(per.loc[p, 'start'])
        t1 = mid0 if per.loc[p, 'end'] == 'END' else L.decimal_year(per.loc[p, 'end'])
        defs[p] = (per.loc[p, 'start'], t0, t1, 'own' if per.loc[p, 'end'] == 'END' else 'fixed date')
    defs['A1'] = ('2026-04-01', L.decimal_year('2026-04-01'), L.decimal_year(per.loc['P1', 'end']), 'fixed date')
    rows, series = [], {}
    first_tip = an.first_tip
    for p in ('P1', 'P2', 'P3', 'P4', 'A1', 'P5', 'P6'):
        start, t0, t1, end_used = defs[p]
        c = L.window_contrast(g, t0, t1, 'interp')
        rs = [x @ c for x in lp]
        series[p] = rs
        r = A.r_summary(rs, gen, p)
        sR = L.summarise([gen.r2R(x) for x in rs], p)
        psd = L.prior_sd_contrast(c, tau)
        pr = A.prior_r_summary(c, tau, gen, K, 200000)
        med = r['growth_per_day_median']
        lo_sd, hi_sd = float(gen.r2R(-z * psd)), float(gen.r2R(z * psd))
        unit, lab, stored = ARTICLE[p]
        d_before = max(0, min((first_tip - pd.Timestamp(start)).days, int(round((t1 - t0) * L.YEAR_DAYS)))) if p in ('P1', 'A1') else ''
        rows.append(dict(base, period=(p if p != 'A1' else 'period 1 April - 15 May of the article'), unit_of_the_article=unit, label_of_the_article=lab, label_of_the_stored_script=stored,
                         period_start=start, period_end=str(L.dec2date(t1)), end_used=end_used, days=(t1 - t0) * L.YEAR_DAYS,
                         growth_rate_per_day_median=med, growth_rate_per_day_hpd95_lower=r['growth_per_day_hpd_lo'], growth_rate_per_day_hpd95_upper=r['growth_per_day_hpd_hi'],
                         R_median=r['R_median'], R_mean=sR['mean'], R_hpd95_lower=r['R_hpd_lo'], R_hpd95_upper=r['R_hpd_hi'], probability_that_R_exceeds_1=r['P_R_gt_1'],
                         doubling_time_days=(math.log(2) / med if med > 0 else ''), halving_time_days=(math.log(2) / -med if med < 0 else ''),
                         doubling_or_halving_time_at_the_lower_bound_of_the_growth_rate_days=math.log(2) / abs(r['growth_per_day_hpd_lo']) if r['growth_per_day_hpd_lo'] != 0 else '',
                         doubling_or_halving_time_at_the_upper_bound_of_the_growth_rate_days=math.log(2) / abs(r['growth_per_day_hpd_hi']) if r['growth_per_day_hpd_hi'] != 0 else '',
                         sign_of_the_growth_rate_at_its_bounds=('+' if r['growth_per_day_hpd_lo'] > 0 else '-') + ' / ' + ('+' if r['growth_per_day_hpd_hi'] > 0 else '-'),
                         prior_R_lower_from_the_prior_standard_deviation_of_the_growth_rate=lo_sd, prior_R_upper_from_the_prior_standard_deviation_of_the_growth_rate=hi_sd,
                         label_of_the_prior_range='from the prior standard deviation of the growth rate (rule of the caption of Table 1)',
                         prior_R_2_5_per_cent_quantile_simulated=pr['prior_R_q025'], prior_R_97_5_per_cent_quantile_simulated=pr['prior_R_q975'],
                         largest_difference_between_the_two_prior_ranges=max(abs(lo_sd - pr['prior_R_q025']), abs(hi_sd - pr['prior_R_q975'])),
                         prior_sd_of_the_growth_rate_per_day=psd / L.YEAR_DAYS, posterior_sd_of_the_growth_rate_per_day=r['growth_per_day_sd'],
                         sd_ratio=r['growth_per_day_sd'] / (psd / L.YEAR_DAYS), R_interval_lies_above_the_prior_range=('yes' if r['R_hpd_lo'] > hi_sd else 'no'),
                         R_interval_lies_below_the_prior_range=('yes' if r['R_hpd_hi'] < lo_sd else 'no'),
                         days_of_the_period_before_the_first_collection_date=d_before, first_collection_date=str(first_tip.date()),
                         states_with_R_set_to_0=r['n_samples_R_set_to_0'], ess_pooled=r['ess_pooled'], rhat_split=r['rhat_split'],
                         meets_convergence_criterion=meets(r['ess_pooled'], r['rhat_split']), ess_pooled_of_the_reproduction_number_series=sR['ess_pooled'],
                         rhat_split_of_the_reproduction_number_series=sR['rhat_split'], mid_point_of_the_most_recent_interval=str(L.dec2date(mid0)),
                         days_beyond_the_most_recent_mid_point=max(0.0, (t1 - mid0) * L.YEAR_DAYS), generation_time='gamma, mean 15.3 d, SD 9.3 d', **tail))
    T1 = pd.DataFrame(rows)
    # ---------------- intervals ----------------
    iv = tabs['intervals'].copy()
    imax = int(iv.loc[iv['Ne_tau_years_median'].idxmax(), 'interval'])
    I = pd.DataFrame(dict(interval=iv['interval'], older_edge=iv['older_edge'], newer_edge=iv['newer_edge'], mid_point=iv['mid_point'], genomes_collected=iv['n_genomes_collected'],
                          Ne_tau_years_median=iv['Ne_tau_years_median'], Ne_tau_years_hpd95_lower=iv['Ne_tau_years_hpd_lo'], Ne_tau_years_hpd95_upper=iv['Ne_tau_years_hpd_hi'],
                          ln_Ne_tau_median=iv['lnNe_median'], ess_pooled=iv['ess_pooled'], rhat_split=iv['rhat_split'],
                          meets_convergence_criterion=[meets(e, r_) for e, r_ in zip(iv['ess_pooled'], iv['rhat_split'])], status_by_the_stored_function=iv['status'],
                          highest_median_of_Ne_tau=['yes' if i == imax else 'no' for i in iv['interval']],
                          one_of_the_two_most_recent_intervals=['yes' if i <= 2 else 'no' for i in iv['interval']],
                          between_the_interval_with_the_highest_median_and_the_two_most_recent_intervals=['yes' if 2 < i < imax else 'no' for i in iv['interval']]))
    for k, v in base.items():
        I.insert(list(base).index(k), k, v)
    for k, v in tail.items():
        I[k] = v
    # ---------------- pairs ----------------
    rt = tabs['rt'].copy()
    pairs = [int(i) for i in rt['interval_pair']]
    Rp = {i: np.concatenate([gen.r2R(x @ L.pair_contrast(g, i)) for x in lp]) for i in pairs}
    Rper = {p: np.concatenate([gen.r2R(x) for x in series[p]]) for p in ('P1', 'P2', 'P3')}
    bounds = {p: (defs[p][1], defs[p][2]) for p in ('P1', 'P2', 'P3')}

    def period_of(i):
        tb = g.t_anchor - i * g.D
        for p, (a, b) in bounds.items():
            if a <= tb < b or (p == 'P3' and a <= tb <= b):
                return p
        return ''

    def pairs_of(p):
        c = L.window_contrast(g, bounds[p][0], bounds[p][1], 'interp')
        cum = np.cumsum(c)[:-1]                                 # weight of pair i (i = 1 .. K-1) in the growth rate of the period, up to the factor D
        return [i for i in pairs if abs(cum[i - 1]) > 1e-12]
    members = {p: pairs_of(p) for p in bounds}
    all_max = int(rt.loc[rt['R_median'].idxmax(), 'interval_pair'])
    ok = rt[rt['status'] == '']
    ok_max = int(ok.loc[ok['R_median'].idxmax(), 'interval_pair']) if len(ok) else all_max
    prow = []
    for _, r in rt.iterrows():
        i = int(r['interval_pair'])
        p = period_of(i)
        e = dict(base, pair=f'{i + 1} to {i}', newer_interval=i, older_interval=i + 1, date_on_which_the_two_intervals_meet=r['time_boundary'],
                 mid_point_of_the_older_interval=r['older_interval_mid'], mid_point_of_the_newer_interval=r['newer_interval_mid'],
                 R_median=r['R_median'], R_hpd95_lower=r['R_hpd_lo'], R_hpd95_upper=r['R_hpd_hi'], probability_that_R_exceeds_1=r['P_R_gt_1'],
                 growth_rate_per_day_median=r['growth_per_day_median'], growth_rate_per_day_hpd95_lower=r['growth_per_day_hpd_lo'], growth_rate_per_day_hpd95_upper=r['growth_per_day_hpd_hi'],
                 posterior_sd_of_the_growth_rate_per_day=r['growth_per_day_sd'], prior_sd_of_the_growth_rate_per_day=r['prior_sd_growth_per_day'],
                 sd_ratio=r['ratio_posterior_sd_to_prior_sd_growth'], sd_ratio_of_0_90_or_more=('yes' if r['ratio_posterior_sd_to_prior_sd_growth'] >= 0.90 else 'no'),
                 correlation_of_the_growth_rate_with_that_of_the_next_older_pair=r.get('corr_growth_with_next_older_pair', ''),
                 status_by_the_stored_function=r['status'], states_with_R_set_to_0=r['n_samples_R_set_to_0'], ess_pooled=r['ess_pooled'], rhat_split=r['rhat_split'],
                 meets_convergence_criterion=meets(r['ess_pooled'], r['rhat_split']),
                 highest_median_among_all_pairs=('yes' if i == all_max else 'no'),
                 highest_median_among_the_pairs_without_a_status_of_the_stored_function=('yes' if i == ok_max else 'no'),
                 newer_interval_is_the_interval_with_the_highest_median_of_Ne_tau=('yes' if i == imax else 'no'),
                 older_interval_is_the_interval_with_the_highest_median_of_Ne_tau=('yes' if i + 1 == imax else 'no'),
                 period_that_holds_the_date_on_which_the_intervals_meet=p,
                 periods_whose_growth_rate_the_pair_enters=' '.join(q for q in members if i in members[q]))
        for q in ('P1', 'P2', 'P3'):
            if i in members[q]:
                d = Rp[i] - Rper[q]
                lo, hi = L.hpd(d)
                am = np.argmax(np.column_stack([Rp[j] for j in members[q]]), axis=1)
                e.update({f'probability_that_R_of_the_pair_exceeds_R_of_{q}': float((Rp[i] > Rper[q]).mean()),
                          f'probability_that_the_pair_is_the_highest_of_the_pairs_of_{q}': float((am == members[q].index(i)).mean()),
                          f'pairs_of_{q}': ' '.join(f'{j + 1} to {j}' for j in members[q]), f'number_of_pairs_of_{q}': len(members[q]),
                          f'R_of_the_pair_minus_R_of_{q}_median': float(np.median(d)), f'R_of_the_pair_minus_R_of_{q}_hpd95_lower': lo,
                          f'R_of_the_pair_minus_R_of_{q}_hpd95_upper': hi})
        prow.append(dict(e, **tail))
    PT = pd.DataFrame(prow).fillna('')
    os.makedirs(outdir, exist_ok=True)
    T1.to_csv(os.path.join(outdir, f'primary_periods_of_table_1_sensitivity_{LABEL}.csv'), index=False)
    I.to_csv(os.path.join(outdir, f'primary_intervals_sensitivity_{LABEL}.csv'), index=False)
    PT.to_csv(os.path.join(outdir, f'primary_pairs_sensitivity_{LABEL}.csv'), index=False)
    return T1, I, PT
