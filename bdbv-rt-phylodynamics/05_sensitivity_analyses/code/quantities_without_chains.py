#!/usr/bin/env python3
"""quantities_without_chains.py - quantities that need no chain, from the
posterior file of the primary analysis of S (posterior_delphy_fixed_C00_20261002.csv.gz) with the functions of analysis.py and rt_lib.py.

  generation_times_sensitivity_20261002.csv   reproduction numbers of P1 to P6 under the generation time of the article and the four
                                              other generation times (list GT of build_tables2.py; functions GenTime.r2R,
                                              hpd, summarise)
  prior_alone_sensitivity_20261002.csv        range of the reproduction number under the smoothing prior alone for P1 to P6:
                                              (a) from the prior standard deviation of the growth rate (prior_sd_contrast;
                                              the growth rate under the prior is normal with mean 0, so that 95 % of it lies within
                                              1.959964 standard deviations; R by GenTime.r2R), (b) simulated 2.5 % and
                                              97.5 % quantiles of prior_r_summary (200,000 draws).  Both are given.
usage: python quantities_without_chains.py POSTERIOR.csv.gz ALIGNMENT.fasta PERIODS.csv PIPELINE_DIR OUTDIR
"""
import math
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

GT = [('mean 15.3 d, SD 9.3 d (reference)', 15.3, 9.3), ('mean 12 d, same coefficient of variation', 12.0, 12.0 * 9.3 / 15.3),
      ('mean 18 d, same coefficient of variation', 18.0, 18.0 * 9.3 / 15.3), ('mean 15.3 d, SD 7 d', 15.3, 7.0), ('mean 15.3 d, SD 12 d', 15.3, 12.0)]      # list GT of build_tables2.py (earlier analysis round)
ANALYSIS, SETTING, SET = 'delphy_fixed_C00', 'Skygrid, smoothing fixed (program default)', 'C00'


def meets(ess, rhat):
    if ess != ess or rhat != rhat:
        return 'not assessed'
    return 'yes' if (ess >= 200 and rhat <= 1.05) else 'no'


def main(posterior, fasta, periods_csv, pipeline, outdir):
    sys.path.insert(0, pipeline)
    import analysis as A
    import rt_lib as L
    os.makedirs(outdir, exist_ok=True)
    P = pd.read_csv(posterior)
    per = pd.read_csv(periods_csv, dtype=str, keep_default_na=False).set_index('period')
    dd = L.tip_dates_from_fasta(fasta)
    last_tip = dd.max()
    K = len([c for c in P.columns if c.startswith('ln_Ne_tau_')])
    g = L.Grid.from_cutoff(K, 0.8, last_tip)
    mid0 = float(g.mids[0])
    seeds = list(dict.fromkeys(P['seed']))
    tau = sorted(set(P['skygrid_precision']))
    assert len(tau) == 1, 'the precision is not fixed'
    tau = float(tau[0])
    base = dict(analysis=ANALYSIS, program='Delphy 1.4.1', setting=SETTING, genome_set=SET, genomes=len(dd), most_recent_collection_date=str(last_tip.date()))
    tail = dict(chains_run=24, chains_retained=len(seeds), states_used=len(P))
    rows, prior = [], []
    for p in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
        start = per.loc[p, 'start']
        t0 = L.decimal_year(start)
        t1 = mid0 if per.loc[p, 'end'] == 'END' else L.decimal_year(per.loc[p, 'end'])
        end = str(L.dec2date(t1)) if per.loc[p, 'end'] == 'END' else per.loc[p, 'end']
        end_used = 'own' if per.loc[p, 'end'] == 'END' else 'fixed date'
        col = f'growth_rate_per_day_{p}' + ('_own_end' if p in ('P3', 'P6') else '')
        rs = [P.loc[P['seed'] == s, col].values * L.YEAR_DAYS for s in seeds]          # growth rate per year of the log, per retained chain
        c = L.window_contrast(g, t0, t1, 'interp')
        lp = np.column_stack([P[f'ln_Ne_tau_{i}'].values for i in range(1, K + 1)])
        worst = float(np.max(np.abs(lp @ c - np.concatenate(rs))))
        assert worst < 1e-9, (p, worst)
        allr = np.concatenate(rs)
        sr = L.summarise(rs, p)
        for lab, m, sd in GT:
            gen = L.GenTime(m, sd)
            R = gen.r2R(allr)
            lo, hi = L.hpd(R)
            sR = L.summarise([gen.r2R(x) for x in rs], p)
            rows.append(dict(base, quantity=f'reproduction number, {p}', period_start=start, period_end=end, end_used=end_used, median=float(np.median(R)), mean=float(np.mean(R)),
                             hpd95_lower=lo, hpd95_upper=hi, ess_pooled=sr['ess_pooled'], rhat_split=sr['rhat_split'], meets_convergence_criterion=meets(sr['ess_pooled'], sr['rhat_split']),
                             **tail, unit='', generation_time=lab, generation_time_mean_days=m, generation_time_sd_days=round(sd, 3), generation_time_sd_days_as_used=sd,
                             states_with_R_set_to_0=gen.n_undefined(allr), probability_that_R_exceeds_1=float((allr > 0).mean()),
                             ess_pooled_of_the_reproduction_number_series=sR['ess_pooled'], rhat_split_of_the_reproduction_number_series=sR['rhat_split'],
                             meets_convergence_criterion_by_the_reproduction_number_series=meets(sR['ess_pooled'], sR['rhat_split']),
                             growth_rate_per_day_median=float(np.median(allr)) / L.YEAR_DAYS, period_length_days=(t1 - t0) * L.YEAR_DAYS,
                             mid_point_of_the_most_recent_interval=str(L.dec2date(mid0)), days_beyond_the_most_recent_mid_point=max(0.0, (t1 - mid0) * L.YEAR_DAYS)))
        gen = L.GenTime(15.3, 9.3)
        pr = A.prior_r_summary(c, tau, gen, K, 200000)
        sd_year = L.prior_sd_contrast(c, tau)
        z = float(norm.ppf(0.975))
        lo_sd, hi_sd = float(gen.r2R(-z * sd_year)), float(gen.r2R(z * sd_year))
        post_sd = float(np.std(allr)) / L.YEAR_DAYS
        R = gen.r2R(allr)
        plo, phi = L.hpd(R)
        d2 = lambda a_, b_: f'{a_:.2f}-{b_:.2f}'
        prior.append(dict(base, quantity=f'reproduction number, {p}', period_start=start, period_end=end, end_used=end_used,
                          skygrid_precision=tau, parameters_of_the_grid=K, interval_length_days=g.D * L.YEAR_DAYS,
                          prior_sd_of_the_growth_rate_per_day=sd_year / L.YEAR_DAYS, prior_sd_of_the_growth_rate_per_day_simulated=pr['prior_sd_growth_per_day_simulated'],
                          multiplier_of_the_standard_deviation=z,
                          prior_R_lower_from_the_standard_deviation=lo_sd, prior_R_upper_from_the_standard_deviation=hi_sd,
                          prior_R_2_5_per_cent_quantile_simulated=pr['prior_R_q025'], prior_R_97_5_per_cent_quantile_simulated=pr['prior_R_q975'],
                          prior_R_median_simulated=pr['prior_R_median'], prior_R_hpd95_lower_simulated=pr['prior_R_hpd_lo'], prior_R_hpd95_upper_simulated=pr['prior_R_hpd_hi'],
                          number_of_simulated_draws=pr['prior_n_sim'], seed_of_the_simulation='20260924 (default of the stored function prior_lnNe_samples)',
                          range_from_the_standard_deviation_two_decimals=d2(lo_sd, hi_sd), range_simulated_two_decimals=d2(pr['prior_R_q025'], pr['prior_R_q975']),
                          the_two_ranges_give_the_same_digits='yes' if d2(lo_sd, hi_sd) == d2(pr['prior_R_q025'], pr['prior_R_q975']) else 'no',
                          range_that_governs='from the standard deviation',
                          label_of_the_first_range='from the prior standard deviation of the growth rate (rule of the caption of Table 1)',
                          label_of_the_second_range='2.5 % and 97.5 % quantiles of simulated trajectories',
                          posterior_R_median=float(np.median(R)), posterior_R_hpd95_lower=plo, posterior_R_hpd95_upper=phi,
                          posterior_sd_of_the_growth_rate_per_day=post_sd, ratio_of_the_posterior_to_the_prior_sd_of_the_growth_rate=post_sd / (sd_year / L.YEAR_DAYS),
                          generation_time='mean 15.3 d, SD 9.3 d', period_length_days=(t1 - t0) * L.YEAR_DAYS, **tail,
                          rule='growth rate under the prior: normal, mean 0, standard deviation by the stored function prior_sd_contrast; range = R at -/+ 1.959964 standard deviations '
                               '(stored GenTime.r2R); simulated: stored prior_r_summary'))
    G = pd.DataFrame(rows)
    G.to_csv(os.path.join(outdir, 'generation_times_sensitivity_20261002.csv'), index=False)
    PR = pd.DataFrame(prior)
    PR.to_csv(os.path.join(outdir, 'prior_alone_sensitivity_20261002.csv'), index=False)
    return G, PR


if __name__ == '__main__':
    main(*sys.argv[1:6])
