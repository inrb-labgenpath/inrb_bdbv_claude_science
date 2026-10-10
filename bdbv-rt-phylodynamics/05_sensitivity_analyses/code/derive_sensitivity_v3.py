#!/usr/bin/env python3
"""derive_sensitivity_v3.py (round of October 2026)

The CODE is the same as that of derive_sensitivity_v2.py and derive_sensitivity.py (checked by program: the files have the same syntax
tree apart from this text); only this head text differs.

Retention, convergence and estimates of ONE analysis of Delphy, from its raw output, in the common form of the tables, in the order of
derive_rt.py of the earlier round: Chain (burn-in 30 %, root-date rule) -> chain that did not reach its last step -> density criterion
(apply_second_criterion, up to 50 sampled trees per chain, threshold 100 log units) -> tables.

COMPUTED BY THE FUNCTIONS OF analysis.py AND rt_lib.py:
  retention of every chain (Chain, the lines of derive_rt.py on chains that did not reach the last step, apply_second_criterion);
  growth rate and reproduction number of every period: median, 95 % HPD interval, pooled effective sample size, split R-hat, states with R
  set to 0, doubling time (window_contrast, Analysis.windows, r_summary, summarise, GenTime.r2R); evolutionary rate and kappa (summarise);
  date of the root: median and interval (tmrca_summary); ln(Ne tau) of every interval and the pairs of adjacent intervals
  (skygrid_tables); prior standard deviation of the growth rate (precision_of, prior_sd_contrast); ln(Ne tau) at a date (Grid.weights);
  every date and decimal year (decimal_year, dec2date, root_dates_from_height).
COMPUTED BY THIS SCRIPT ITSELF, because no function of analysis.py or rt_lib.py gives it (each is an arithmetic step on values of those
functions):
  1  the MEAN of the date of the root (numpy mean of the logged root heights; date = most recent collection date - rounded days);
     the means of the other quantities are the field 'mean' of summarise
  2  the four posterior PROBABILITIES: share of the states in which the condition holds (numpy mean of an indicator), with
     summarise applied to the indicator for its diagnostics; 'not assessed' where the indicator is constant in a chain
  3  Ne tau = exponential of the median and of the bounds of ln(Ne tau)
  4  the verdict 'meets_convergence_criterion' (effective sample size of at least 200 and R-hat of at most 1.05) from the two diagnostics
  5  period_length_days, days_beyond_the_most_recent_mid_point (differences of decimal years times 365)
  6  the columns on the tree of step 0 (difference at step 0, largest difference after step 0), read from the table of sampled trees that
     approximation_error returns
  7  the list of the shortest 95 % intervals of the root date in whole days (the functions return one of them)
  8  the dates of 'Ne tau at every 7th day': from the date of the mid-point of the most recent interval BACKWARD in steps of 7 days to the
     first such day on or after 1 January 2026. This is NOT the rule of the earlier example, which steps FORWARD from the date of the
     first mid-point after 1 January to the date of the most recent mid-point; on the grid of the example both rules give the same rows,
     on other grids the rows differ (on the grid of S this script gives one row more, at the start). The rule of this script is kept for
     every analysis; the difference is declared with the number of rows concerned.
  9  the ratio of the posterior to the prior standard deviation of the growth rate (quotient of two values of the functions)

CONSTANTS OF THE EARLIER ROUND in its scripts, and what is passed in their place:
  analysis.W_DATES            list of periods W1, W2, W3, W3b -> the periods P1, P2, P3, P4 of this round under the earlier keys,
                              and P5 added under the key 'P5' (fixed dates; evaluated by the rule of Analysis.windows)
  analysis.R0_START / R0_END  dates of P1 (set from the table of the periods)
  build_tables2.recent_extra  the period from 1 August to the mid-point of the most recent interval (P6) is computed by the lines
                              of that function (window_contrast, r_summary; 'not evaluated' where the end is not after the start)
  master.SFX, file names      the suffix 20261002 and the file names of this round
  derive_rt.py --steps        the number of steps of the command (1,500,000,000)
  most recent collection date read from the alignment of the analysis (tip_dates_from_fasta), never a constant
RULE OF Analysis.windows FOR PERIODS WITH FIXED DATES: a period is evaluated if its end lies on or before the most
recent collection date of the set and after its start; where the end lies beyond the mid-point of the most recent interval,
ln(Ne tau) at the end is the level of that interval; where the end lies after the most recent collection date the period is
not evaluated.
"""
import argparse
import csv
import gzip
import json
import math
import os
import re
import sys

import numpy as np
import pandas as pd

SUFFIX = '20261002'
BURNIN = 0.30
THRESHOLD = 100.0
MAX_TREES = 50
THIN_TO = 300000
EST_COLS = ['analysis', 'program', 'setting', 'genome_set', 'genomes', 'most_recent_collection_date', 'quantity', 'period_start', 'period_end', 'end_used',
            'median', 'mean', 'hpd95_lower', 'hpd95_upper', 'ess_pooled', 'rhat_split', 'meets_convergence_criterion', 'chains_run', 'chains_retained',
            'states_used', 'unit']
EST_EXTRA = ['states_with_R_set_to_0', 'ess_pooled_of_the_reproduction_number_series', 'rhat_split_of_the_reproduction_number_series',
             'meets_convergence_criterion_by_the_reproduction_number_series', 'period_length_days', 'mid_point_of_the_most_recent_interval',
             'days_beyond_the_most_recent_mid_point', 'note', 'posterior_sd_of_the_growth_rate_per_day', 'prior_sd_of_the_growth_rate_per_day',
             'ratio_of_the_posterior_to_the_prior_sd_of_the_growth_rate', 'skygrid_precision_fixed_at']
CHAIN_COLS = ['analysis', 'program', 'setting', 'genome_set', 'seed', 'command_line', 'started_utc', 'ended_utc', 'wall_time_s', 'exit_code', 'last_state_logged',
              'reached_last_step', 'started_again', 'root_date_rule_met', 'first_state_root_date_rule', 'density_criterion_met', 'n_sampled_trees_evaluated',
              'retained', 'why', 'n_states_after_burn_in', 'attempts', 'latest_root_date_logged', 'largest_difference_of_the_density_log_units',
              'end_of_the_last_attempt_by_the_ledger', 'density_difference_at_step_0', 'largest_density_difference_after_step_0', 'density_criterion_met_after_step_0',
              'discarded_for_the_tree_of_step_0_alone', 'first_tree_step_density_criterion', 'address_space_limit_kib', 'precision_is_estimated', 'smallest_precision',
              'largest_precision', 'entered_degenerate_state_by_the_definition_of_the_brief', 'files_cut_at_the_end', 'last_lines_of_standard_output_not_progress']
INT_COLS = ['analysis', 'program', 'setting', 'genome_set', 'interval', 'older_edge', 'newer_edge', 'mid_point', 'mid_point_decimal', 'genomes_collected',
            'ln_Ne_tau_median', 'ln_Ne_tau_hpd95_lower', 'ln_Ne_tau_hpd95_upper', 'ess_pooled', 'rhat_split', 'meets_convergence_criterion', 'status',
            'chains_retained', 'states_used', 'unit']
TRAJ_COLS = ['analysis', 'program', 'setting', 'genome_set', 'series', 'date', 'date_decimal', 'interval', 'older_date', 'newer_date', 'median', 'hpd95_lower',
             'hpd95_upper', 'unit', 'ess_pooled', 'rhat_split', 'note', 'meets_convergence_criterion', 'chains_retained', 'states_used']
NOT_EVALUATED_NO_LENGTH = 'not evaluated: no length'
NOT_EVALUATED_AFTER_LAST = 'not evaluated: the period ends after the most recent collection date of the set (rule of the stored function Analysis.windows)'
PROB_TEXT = "not applicable (probability: the value is in the column 'mean')"


def utc(t):
    import time
    return time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime(int(t)))


def read_text(path):
    try:
        return open(path, errors='replace').read().strip()
    except OSError:
        return ''


def whole_lines_copy(path, outdir, n_fields=None):
    """(file to read, note). A program that ends itself can leave a file that ends inside a line; a copy without the cut line is
    written to outdir, the file of the chain is not changed."""
    try:
        data = open(path, 'rb').read()
    except OSError:
        return path, 'file not found'
    whole = data.endswith(b'\n') or data.rstrip().endswith(b'End;')
    cut = len(data) > 0 and not whole
    body = data[:data.rfind(b'\n') + 1] if cut else data
    note = 'the file ends inside a line; the cut line was left out' if cut else ''
    if n_fields is not None:
        lines = body.split(b'\n')
        while len(lines) > 1 and lines[-1] == b'':
            lines.pop()
        if lines and not lines[-1].startswith(b'#') and lines[-1].count(b'\t') + 1 < n_fields:
            body = b'\n'.join(lines[:-1]) + b'\n'
            cut = True
            note = 'the last line of the file holds fewer fields than the header; it was left out'
    if not cut:
        return path, ''
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, os.path.basename(path))
    open(out, 'wb').write(body)
    return out, note


def last_state_of_log(path):
    last, nf = None, None
    try:
        with open(path, errors='replace') as fh:
            for line in fh:
                if line.startswith('#'):
                    continue
                if line.startswith('state\t'):
                    nf = len(line.rstrip('\n').rstrip('\t').split('\t'))
                    continue
                p = line.split('\t', 1)[0]
                if p.isdigit() and line.endswith('\n'):
                    last = int(p)
    except OSError:
        pass
    return last, nf


def not_progress_lines(path, n=3):
    out = []
    try:
        with open(path, 'rb') as fh:
            fh.seek(0, 2)
            size = fh.tell()
            fh.seek(max(0, size - 6000))
            tail = fh.read().decode(errors='replace').split('\n')
        for line in tail[1:]:
            if line.strip() and 'Step ' not in line and 'log_posterior' not in line:
                out.append(line.strip()[:200])
    except OSError:
        pass
    return ' | '.join(out[-n:])


def meets(ess, rhat):
    if ess != ess or rhat != rhat:
        return 'not assessed'
    return 'yes' if (ess >= 200 and rhat <= 1.05) else 'no'


def fmt(x):
    return x


def derive(analysis, setting, genome_set, seeds, steps, run_dir, fasta, pipeline, periods_csv, outdir, ledger_rows=None,
           num_parameters=20, cutoff=0.8, first_knot='', last_knot='', estimated=False, write=True, only_seeds=None, label_suffix=''):
    sys.path.insert(0, pipeline)
    import analysis as A
    import rt_lib as L
    os.makedirs(outdir, exist_ok=True)
    per = pd.read_csv(periods_csv, dtype=str, keep_default_na=False).set_index('period')
    # ---- constants of the earlier round replaced by the values of this round (see the head of the file) ----
    A.W_DATES = {'W1': (per.loc['P1', 'start'], per.loc['P1', 'end']), 'W2': (per.loc['P2', 'start'], per.loc['P2', 'end']),
                 'W3': (per.loc['P3', 'start'], None), 'W3b': (per.loc['P4', 'start'], per.loc['P4', 'end']),
                 'P5': (per.loc['P5', 'start'], per.loc['P5', 'end'])}
    A.R0_START, A.R0_END = per.loc['P1', 'start'], per.loc['P1', 'end']
    A.ERR_THRESHOLD = THRESHOLD
    P6_START = per.loc['P6', 'start']
    assert per.loc['P3', 'end'] == 'END' and per.loc['P6', 'end'] == 'END'
    KEY = {'P1': 'W1', 'P2': 'W2', 'P3': 'W3', 'P4': 'W3b', 'P5': 'P5'}
    dd = L.tip_dates_from_fasta(fasta)
    n_genomes, first_tip, last_tip = len(dd), dd.min(), dd.max()
    program = 'Delphy 1.4.1'
    base = dict(analysis=analysis + label_suffix, program=program, setting=setting, genome_set=genome_set)
    led = {}
    for r in (ledger_rows or []):
        if r['analysis'] == analysis:
            led.setdefault(r['seed'], []).append(r)
    san = os.path.join(outdir, 'copies_without_the_cut_line')
    chains, info_rows = [], {}
    use = [s for s in seeds if (only_seeds is None or s in only_seeds)]
    for seed in use:
        p = os.path.join(run_dir, f'{analysis}_s{seed}')
        last_txt, nf = last_state_of_log(p + '.log')
        f_log, n1 = whole_lines_copy(p + '.log', san, nf)
        f_trees, n2 = whole_lines_copy(p + '.trees', san)
        f_out, n3 = whole_lines_copy(p + '.stdout', san)
        cut_note = '; '.join(f'{k}: {v}' for k, v in (('log', n1), ('tree file', n2), ('standard output', n3)) if v)
        t0, t1 = read_text(p + '.start'), read_text(p + '.end')
        env = dict(x.split('=', 1) for x in read_text(p + '.env').split('\n') if '=' in x)
        lr = led.get(str(seed), [])
        row = dict(base, seed=seed, command_line=read_text(p + '.cmd'), started_utc=utc(t0) if t0 else '', ended_utc=utc(t1) if t1 else '',
                   wall_time_s=(int(t1) - int(t0)) if t0 and t1 else '', exit_code=read_text(p + '.exit'), last_state_logged=last_txt if last_txt is not None else '',
                   reached_last_step='yes' if last_txt == steps else 'no', attempts=read_text(p + '.attempt'),
                   started_again='yes' if (read_text(p + '.attempt') not in ('', '1')) else 'no',
                   end_of_the_last_attempt_by_the_ledger=(lr[-1]['outcome'] if lr else 'not assessed: no row of the ledger'),
                   address_space_limit_kib=env.get('address_space_limit_kib', ''), files_cut_at_the_end=cut_note,
                   last_lines_of_standard_output_not_progress=not_progress_lines(p + '.stdout'))
        try:
            c = A.Chain('delphy', f_log, fasta, int(seed), burnin=BURNIN, first_tip=first_tip, last_tip=last_tip, label=f'{analysis}_s{seed}')
            sp = int(np.diff(c.full['state'].values[:2])[0])
            if THIN_TO and sp < THIN_TO:                                   # lines of derive_rt.py of the earlier round
                keep = (c.full['state'] % THIN_TO == 0).values
                c.full = c.full[keep].reset_index(drop=True)
                c.root_dates_full = c.root_dates_full[keep]
                c.collapsed_mask = c.collapsed_mask[keep]
                c.n_logged = len(c.full)
                c.nb = int(c.n_logged * BURNIN)
                c.post = c.full.iloc[c.nb:].reset_index(drop=True)
            c.info['steps'] = steps
            c.info['complete'] = bool(c.info['last_state'] >= steps)
            if not c.info['complete'] and c.retained:
                c.retained = False
                c.why = f"incomplete (stopped at step {c.info['last_state']} of {steps})"
            c.f_out, c.f_trees = f_out, f_trees
            chains.append(c)
        except Exception as exc:
            row.update(retained='no', why=f'not assessed: the log could not be read ({str(exc)[:120]})', root_date_rule_met='not assessed',
                       density_criterion_met='not assessed', n_sampled_trees_evaluated=0, n_states_after_burn_in=0)
        info_rows[seed] = row
    grid = (L.Grid.from_knots(num_parameters, first_knot, last_knot, last_tip) if first_knot else L.Grid.from_cutoff(num_parameters, cutoff, last_tip))
    an = A.Analysis(analysis, 'delphy', 'Skygrid', os.path.basename(fasta), fasta, chains, grid=grid, gen=L.GenTime(15.3, 9.3))
    an.progress = {c.seed: c.f_out for c in chains}
    an.trees = {c.seed: c.f_trees for c in chains}
    # ---- density criterion: the function, chain by chain so that a file that cannot be read stops one chain only ----
    for c in chains:
        sub = A.Analysis(analysis, 'delphy', 'Skygrid', os.path.basename(fasta), fasta, [c], grid=grid, gen=an.gen)
        sub.progress, sub.trees = {c.seed: c.f_out}, {c.seed: c.f_trees}
        try:
            A.apply_second_criterion(sub, THRESHOLD, MAX_TREES)
            c.density_assessed = True
        except Exception as exc:
            c.density_assessed = False
            c.density_error = str(exc)[:160]
            if c.retained:
                c.retained = False
                c.why = f'density criterion not assessed ({c.density_error}); the chain is left out and reported'
    for c in chains:
        row = info_rows[str(c.seed)] if str(c.seed) in info_rows else info_rows[c.seed]
        E = getattr(c, 'err', None)
        d0 = dmax_after = np.nan
        met_after = 'not assessed'
        if c.density_assessed and E is not None and len(E):
            e0 = E[E['state'] == 0]
            d0 = float(e0['err'].iloc[0]) if len(e0) else np.nan
            ea = E[E['state'] > 0]
            dmax_after = float(ea['err'].max()) if len(ea) else np.nan
            met_after = ('yes' if (len(ea) and dmax_after > THRESHOLD) else 'no')
        flagged = c.info.get('flagged_by_density_criterion')
        prec_lo, prec_hi = c.info.get('skygrid_precision_min'), c.info.get('skygrid_precision_max')
        row.update(root_date_rule_met='yes' if c.info['collapsed'] else 'no',
                   first_state_root_date_rule=(int(c.info['first_collapsed_state']) if c.info['collapsed'] else ''),
                   density_criterion_met=(('yes' if flagged else 'no') if c.density_assessed else 'not assessed'),
                   n_sampled_trees_evaluated=c.info.get('n_trees_checked', 0), retained='yes' if c.retained else 'no', why=c.why,
                   n_states_after_burn_in=len(c.post), latest_root_date_logged=c.info['latest_root_date'],
                   largest_difference_of_the_density_log_units=c.info.get('approx_error_max', np.nan),
                   density_difference_at_step_0=d0, largest_density_difference_after_step_0=dmax_after, density_criterion_met_after_step_0=met_after,
                   discarded_for_the_tree_of_step_0_alone=('yes' if (c.density_assessed and flagged and not c.info['collapsed'] and c.info['complete'] and met_after == 'no') else 'no'),
                   first_tree_step_density_criterion=(int(c.info['first_tree_state_error_above_threshold']) if (c.density_assessed and flagged) else ''),
                   precision_is_estimated='yes' if estimated else 'no', smallest_precision=prec_lo, largest_precision=prec_hi,
                   entered_degenerate_state_by_the_definition_of_the_brief=('yes' if c.info['collapsed'] else 'no'))
    CH = pd.DataFrame([info_rows[s] for s in use])
    for col in CHAIN_COLS:
        if col not in CH.columns:
            CH[col] = ''
    CH = CH[CHAIN_COLS]
    out = dict(chains=CH, analysis=an)
    kept = an.kept
    n_run, n_kept = len(use), len(kept)
    n_states = int(sum(len(c.post) for c in kept))
    common = dict(base, genomes=n_genomes, most_recent_collection_date=str(last_tip.date()))
    tail = dict(chains_run=n_run, chains_retained=n_kept, states_used=n_states)
    mid0 = float(grid.mids[0])
    mid0_date = str(L.dec2date(mid0))
    t_last = L.decimal_year(last_tip)
    est = []

    def period_cols(t0, t1):
        return dict(period_length_days=(t1 - t0) * L.YEAR_DAYS, mid_point_of_the_most_recent_interval=mid0_date,
                    days_beyond_the_most_recent_mid_point=max(0.0, (t1 - mid0) * L.YEAR_DAYS))

    # ---- periods: definitions by the functions (also where no chain is retained) ----
    W = an.windows('interp', 'midpoint')
    PD = {}
    for P in ('P1', 'P2', 'P3', 'P4', 'P5'):
        w = W[KEY[P]]
        state = 'evaluated'
        if not w['t1'] > w['t0']:
            state = NOT_EVALUATED_NO_LENGTH
        elif not w['evaluable']:
            state = NOT_EVALUATED_AFTER_LAST
        PD[P] = dict(t0=w['t0'], t1=w['t1'], c=w['c'], state=state, note=w['note'], start=per.loc[P, 'start'],
                     end=(str(L.dec2date(w['t1'])) if P == 'P3' else per.loc[P, 'end']), end_used=('own' if P == 'P3' else 'fixed date'))
    t_b = L.decimal_year(P6_START)                       # lines of build_tables2.recent_extra
    if mid0 <= t_b + 1e-9:
        PD['P6'] = dict(t0=t_b, t1=mid0, c=None, state=NOT_EVALUATED_NO_LENGTH, note='the mid-point of the most recent interval lies on or before the start of the period',
                        start=P6_START, end=mid0_date, end_used='own')
    else:
        PD['P6'] = dict(t0=t_b, t1=mid0, c=L.window_contrast(grid, t_b, mid0, 'interp'), state='evaluated', note='', start=P6_START, end=mid0_date, end_used='own')
    out['periods'] = PD
    if n_kept == 0:
        for P in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
            d = PD[P]
            for q, unit in ((f'growth rate, {P}', 'per day'), (f'reproduction number, {P}', '')):
                txt = d['state'] if d['state'] != 'evaluated' else 'not assessed: no chain retained'
                est.append(dict(common, quantity=q, period_start=d['start'], period_end=d['end'], end_used=d['end_used'], median=txt, mean=txt, hpd95_lower=txt,
                                hpd95_upper=txt, ess_pooled=txt, rhat_split=txt, meets_convergence_criterion=txt, unit=unit, note=d['note'], **tail, **period_cols(d['t0'], d['t1'])))
        for q in ('doubling time, P1', 'evolutionary rate', 'date of the most recent common ancestor (date)', 'date of the most recent common ancestor (decimal year)', 'kappa',
                  'P(R of P2 < R of P1)', 'P(R of P1 > 1)', 'P(R of P2 > 1)', 'P(R of P3 < 1)'):
            txt = 'not assessed: no chain retained'
            est.append(dict(common, quantity=q, median=txt, mean=txt, hpd95_lower=txt, hpd95_upper=txt, ess_pooled=txt, rhat_split=txt, meets_convergence_criterion=txt, **tail))
        E_ = pd.DataFrame(est)
        for col in EST_COLS + EST_EXTRA:
            if col not in E_.columns:
                E_[col] = ''
        out.update(estimates=E_[EST_COLS + EST_EXTRA].fillna(''), intervals=pd.DataFrame(columns=INT_COLS), trajectories=pd.DataFrame(columns=TRAJ_COLS), posterior=None,
                   root_date_intervals=None)
        return out
    lp = an.lnNe()
    gen = an.gen
    tau_fixed, fixed_precision = A.precision_of(an)
    RS, RR = {}, {}
    for P in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
        d = PD[P]
        pc_ = period_cols(d['t0'], d['t1'])
        if d['state'] != 'evaluated':
            for q, unit in ((f'growth rate, {P}', 'per day'), (f'reproduction number, {P}', '')):
                est.append(dict(common, quantity=q, period_start=d['start'], period_end=d['end'], end_used=d['end_used'], median=d['state'], mean=d['state'],
                                hpd95_lower=d['state'], hpd95_upper=d['state'], ess_pooled=d['state'], rhat_split=d['state'], meets_convergence_criterion=d['state'],
                                unit=unit, note=d['note'], **tail, **pc_))
            if P == 'P1':
                est.append(dict(common, quantity='doubling time, P1', period_start=d['start'], period_end=d['end'], end_used=d['end_used'], median=d['state'], mean=d['state'],
                                hpd95_lower=d['state'], hpd95_upper=d['state'], ess_pooled=d['state'], rhat_split=d['state'], meets_convergence_criterion=d['state'],
                                unit='days', note=d['note'], **tail, **pc_))
            continue
        rs = [x @ d['c'] for x in lp]                                      # growth rate per year of the log, per retained chain
        RS[P] = rs
        r = A.r_summary(rs, gen, P)
        sg = L.summarise(rs, P, 1.0 / L.YEAR_DAYS)
        Rs = [gen.r2R(x) for x in rs]
        RR[P] = Rs
        sR = L.summarise(Rs, P)
        ok = meets(r['ess_pooled'], r['rhat_split'])
        if fixed_precision:                                                # functions precision_of and prior_sd_contrast
            psd = L.prior_sd_contrast(d['c'], tau_fixed) / L.YEAR_DAYS
            pc_ = dict(pc_, posterior_sd_of_the_growth_rate_per_day=r['growth_per_day_sd'], prior_sd_of_the_growth_rate_per_day=psd,
                       ratio_of_the_posterior_to_the_prior_sd_of_the_growth_rate=r['growth_per_day_sd'] / psd, skygrid_precision_fixed_at=tau_fixed)
        est.append(dict(common, quantity=f'growth rate, {P}', period_start=d['start'], period_end=d['end'], end_used=d['end_used'], median=r['growth_per_day_median'],
                        mean=sg['mean'], hpd95_lower=r['growth_per_day_hpd_lo'], hpd95_upper=r['growth_per_day_hpd_hi'], ess_pooled=r['ess_pooled'],
                        rhat_split=r['rhat_split'], meets_convergence_criterion=ok, unit='per day', note=d['note'], **tail, **pc_))
        est.append(dict(common, quantity=f'reproduction number, {P}', period_start=d['start'], period_end=d['end'], end_used=d['end_used'], median=r['R_median'],
                        mean=sR['mean'], hpd95_lower=r['R_hpd_lo'], hpd95_upper=r['R_hpd_hi'], ess_pooled=r['ess_pooled'], rhat_split=r['rhat_split'],
                        meets_convergence_criterion=ok, unit='', states_with_R_set_to_0=r['n_samples_R_set_to_0'],
                        ess_pooled_of_the_reproduction_number_series=sR['ess_pooled'], rhat_split_of_the_reproduction_number_series=sR['rhat_split'],
                        meets_convergence_criterion_by_the_reproduction_number_series=meets(sR['ess_pooled'], sR['rhat_split']), note=d['note'], **tail, **pc_))
        if P == 'P1':
            est.append(dict(common, quantity='doubling time, P1', period_start=d['start'], period_end=d['end'], end_used=d['end_used'], median=r['doubling_time_d'],
                            mean='not assessed: the stored function gives no mean (the doubling time is ln 2 / median growth rate)',
                            hpd95_lower=r['doubling_time_d_lo'], hpd95_upper=r['doubling_time_d_hi'], ess_pooled=r['ess_pooled'], rhat_split=r['rhat_split'],
                            meets_convergence_criterion=ok, unit='days', note=d['note'], **tail, **pc_))
    s = L.summarise(an.series('clock.rate'), 'clock rate')
    est.append(dict(common, quantity='evolutionary rate', median=s['median'], mean=s['mean'], hpd95_lower=s['hpd_lo'], hpd95_upper=s['hpd_hi'], ess_pooled=s['ess_pooled'],
                    rhat_split=s['rhat_split'], meets_convergence_criterion=meets(s['ess_pooled'], s['rhat_split']), unit='substitutions per site and year', **tail))
    tm = A.tmrca_summary(an)
    rh_all = np.concatenate(an.root_height())
    mean_days = float(np.mean(rh_all * L.YEAR_DAYS))
    okt = meets(tm['ess_pooled'], tm['rhat_split'])
    est.append(dict(common, quantity='date of the most recent common ancestor (date)', median=tm['tmrca_median'],
                    mean=str((an.last_tip - pd.Timedelta(days=int(round(mean_days)))).date()), hpd95_lower=tm['tmrca_hpd_early'], hpd95_upper=tm['tmrca_hpd_late'],
                    ess_pooled=tm['ess_pooled'], rhat_split=tm['rhat_split'], meets_convergence_criterion=okt, unit='date (day resolution, from the root height of the log)', **tail))
    est.append(dict(common, quantity='date of the most recent common ancestor (decimal year)', median=tm['tmrca_decimal_median'], mean=an.t_last - mean_days / L.YEAR_DAYS,
                    hpd95_lower=tm['tmrca_decimal_early'], hpd95_upper=tm['tmrca_decimal_late'], ess_pooled=tm['ess_pooled'], rhat_split=tm['rhat_split'],
                    meets_convergence_criterion=okt, unit='decimal year (year of 365 days)', **tail))
    s = L.summarise(an.series('kappa'), 'kappa')
    est.append(dict(common, quantity='kappa', median=s['median'], mean=s['mean'], hpd95_lower=s['hpd_lo'], hpd95_upper=s['hpd_hi'], ess_pooled=s['ess_pooled'],
                    rhat_split=s['rhat_split'], meets_convergence_criterion=meets(s['ess_pooled'], s['rhat_split']), unit='', **tail))
    if estimated:
        s = L.summarise(an.series('skygrid.precision'), 'precision')
        est.append(dict(common, quantity='precision of the Skygrid', median=s['median'], mean=s['mean'], hpd95_lower=s['hpd_lo'], hpd95_upper=s['hpd_hi'],
                        ess_pooled=s['ess_pooled'], rhat_split=s['rhat_split'], meets_convergence_criterion=meets(s['ess_pooled'], s['rhat_split']), unit='', **tail))

    def prob(q, ind, P, both=None):
        d = PD[P]
        if ind is None:
            txt = (d['state'] if d['state'] != 'evaluated' else (PD[both]['state'] if both else 'not assessed'))
            est.append(dict(common, quantity=q, period_start='' if both else d['start'], period_end='' if both else d['end'], end_used=d['end_used'], median=txt, mean=txt,
                            hpd95_lower=txt, hpd95_upper=txt, ess_pooled=txt, rhat_split=txt, meets_convergence_criterion=txt, unit='probability', **tail))
            return
        const = any(np.std(x) == 0 for x in ind)
        s_ = L.summarise(ind, q)
        pc2 = {} if both else period_cols(d['t0'], d['t1'])
        est.append(dict(common, quantity=q, period_start='' if both else d['start'], period_end='' if both else d['end'], end_used=d['end_used'], median=PROB_TEXT,
                        mean=float(np.mean(np.concatenate(ind))), hpd95_lower=PROB_TEXT, hpd95_upper=PROB_TEXT,
                        ess_pooled='' if const else s_['ess_pooled'], rhat_split='' if const else s_['rhat_split'],
                        meets_convergence_criterion=('not assessed: the indicator is constant in at least one chain' if const else meets(s_['ess_pooled'], s_['rhat_split'])),
                        unit='probability', **tail, **pc2))
    prob('P(R of P2 < R of P1)', [(b < a_).astype(float) for a_, b in zip(RR['P1'], RR['P2'])] if ('P1' in RR and 'P2' in RR) else None, 'P1', both='P2')
    prob('P(R of P1 > 1)', [(x > 1).astype(float) for x in RR['P1']] if 'P1' in RR else None, 'P1')
    prob('P(R of P2 > 1)', [(x > 1).astype(float) for x in RR['P2']] if 'P2' in RR else None, 'P2')
    prob('P(R of P3 < 1)', [(x < 1).astype(float) for x in RR['P3']] if 'P3' in RR else None, 'P3')
    E_ = pd.DataFrame(est)
    for col in EST_COLS + EST_EXTRA:
        if col not in E_.columns:
            E_[col] = ''
    out['estimates'] = E_[EST_COLS + EST_EXTRA].fillna('')
    # ---- intervals and trajectories (skygrid_tables) ----
    tabs = A.skygrid_tables(an)
    K = grid.K
    iv = tabs['intervals']
    IV = pd.DataFrame(dict(interval=iv['interval'], older_edge=iv['older_edge'], newer_edge=iv['newer_edge'], mid_point=iv['mid_point'],
                           mid_point_decimal=[grid.mids[i - 1] if i < K else '' for i in iv['interval']], genomes_collected=iv['n_genomes_collected'],
                           ln_Ne_tau_median=iv['lnNe_median'], ln_Ne_tau_hpd95_lower=iv['lnNe_hpd_lo'], ln_Ne_tau_hpd95_upper=iv['lnNe_hpd_hi'],
                           ess_pooled=iv['ess_pooled'], rhat_split=iv['rhat_split'],
                           meets_convergence_criterion=[meets(e, r_) for e, r_ in zip(iv['ess_pooled'], iv['rhat_split'])], status=iv['status']))
    for k, v in base.items():
        IV[k] = v
    IV['chains_retained'], IV['states_used'], IV['unit'] = n_kept, n_states, 'ln(years)'
    out['intervals'] = IV[INT_COLS]
    tr = []
    open_note = f'interval {K} is open towards the past (everything before the ' + ('first knot date' if first_knot else 'cut-off') + '); it has no mid-point'
    for _, r in iv.iterrows():
        i = int(r['interval'])
        tr.append(dict(base, series='Ne tau at the mid-point of an interval', date=r['mid_point'], date_decimal=(grid.mids[i - 1] if i < K else ''), interval=i,
                       older_date='', newer_date='', median=math.exp(r['lnNe_median']), hpd95_lower=math.exp(r['lnNe_hpd_lo']), hpd95_upper=math.exp(r['lnNe_hpd_hi']),
                       unit='years', ess_pooled=r['ess_pooled'], rhat_split=r['rhat_split'], note=(open_note if i == K else r['status']),
                       meets_convergence_criterion=meets(r['ess_pooled'], r['rhat_split'])))
    end_day = pd.Timestamp(L.dec2date(mid0))
    days = []
    d_ = end_day
    while d_ >= pd.Timestamp('2026-01-01'):
        days.append(d_)
        d_ -= pd.Timedelta(days=7)
    for d_ in sorted(days):
        t = L.decimal_year(d_)
        w = grid.weights(t, 'interp')
        s_ = L.summarise([x @ w for x in lp], 'lnNe')
        tr.append(dict(base, series='Ne tau at every 7th day', date=str(d_.date()), date_decimal=t, interval=grid.interval_of(t), older_date='', newer_date='',
                       median=math.exp(s_['median']), hpd95_lower=math.exp(s_['hpd_lo']), hpd95_upper=math.exp(s_['hpd_hi']), unit='years', ess_pooled=s_['ess_pooled'],
                       rhat_split=s_['rhat_split'], note='ln(Ne tau) read by linear interpolation between interval mid-points',
                       meets_convergence_criterion=meets(s_['ess_pooled'], s_['rhat_split'])))
    rt = tabs['rt']
    for _, r in rt.iterrows():
        i = int(r['interval_pair'])
        t_b_ = grid.t_anchor - i * grid.D
        cmn = dict(base, date=r['time_boundary'], date_decimal=t_b_, interval=f'{i + 1} to {i}', older_date=r['older_interval_mid'], newer_date=r['newer_interval_mid'],
                   ess_pooled=r['ess_pooled'], rhat_split=r['rhat_split'], note=r['status'], meets_convergence_criterion=meets(r['ess_pooled'], r['rhat_split']))
        tr.append(dict(cmn, series='reproduction number between successive mid-points', median=r['R_median'], hpd95_lower=r['R_hpd_lo'], hpd95_upper=r['R_hpd_hi'], unit=''))
        tr.append(dict(cmn, series='growth rate between successive mid-points', median=r['growth_per_day_median'], hpd95_lower=r['growth_per_day_hpd_lo'],
                       hpd95_upper=r['growth_per_day_hpd_hi'], unit='per day'))
    TR = pd.DataFrame(tr)
    TR['chains_retained'], TR['states_used'] = n_kept, n_states
    out['trajectories'] = TR[TRAJ_COLS]
    out['tabs'] = tabs
    # ---- posterior file ----
    frames = []
    for c, x in zip(kept, lp):
        f = pd.DataFrame(dict(seed=c.seed, state=c.post['state'].values, evolutionary_rate=c.post['clock.rate'].values, root_height_years=c.post[c.rh_col].values,
                              root_date=[str(v) for v in L.root_dates_from_height(c.post[c.rh_col].values, last_tip)], kappa=c.post['kappa'].values))
        for i, col in enumerate(c.lp_cols, 1):
            f[f'ln_Ne_tau_{i}'] = c.post[col].values
        f['skygrid_precision'] = c.post['skygrid.precision'].values
        for P in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
            nm = P + ('_own_end' if P in ('P3', 'P6') else '')
            f[f'growth_rate_per_day_{nm}'] = (x @ PD[P]['c']) / L.YEAR_DAYS if PD[P]['state'] == 'evaluated' else np.nan
        for P in ('P1', 'P2', 'P3', 'P4', 'P5', 'P6'):
            nm = P + ('_own_end' if P in ('P3', 'P6') else '')
            f[f'reproduction_number_{nm}'] = gen.r2R(x @ PD[P]['c']) if PD[P]['state'] == 'evaluated' else np.nan
        frames.append(f)
    out['posterior'] = pd.concat(frames, ignore_index=True)
    # ---- shortest 95 % intervals of the root date in whole days ----
    days_all = np.sort(np.round(rh_all * L.YEAR_DAYS).astype(int))
    n = len(days_all)
    m = int(np.ceil(0.95 * n))
    widths = days_all[m - 1:] - days_all[:n - m + 1]
    wmin = int(widths.min())
    js = np.where(widths == wmin)[0]
    cand = sorted(set((int(days_all[j]), int(days_all[j + m - 1])) for j in js))
    f_date = lambda dys: str((an.last_tip - pd.Timedelta(days=int(dys))).date())
    out['root_date_intervals'] = dict(base, genomes=n_genomes, most_recent_collection_date=str(last_tip.date()), states_used=n, states_inside_an_interval_at_least=m,
                                      width_of_the_shortest_interval_days=wmin, number_of_shortest_intervals_in_whole_days=len(cand),
                                      earliest_of_them=f'{f_date(max(c_[1] for c_ in cand))} to {f_date(max(c_[0] for c_ in cand if c_[1] == max(x[1] for x in cand)))}',
                                      latest_of_them=f'{f_date(min(c_[1] for c_ in cand))} to {f_date(min(c_[0] for c_ in cand if c_[1] == min(x[1] for x in cand)))}',
                                      interval_that_the_stored_function_returns=f"{tm['tmrca_hpd_early']} to {tm['tmrca_hpd_late']}",
                                      chains_retained=n_kept)
    return out


def write_tables(out, analysis, outdir):
    os.makedirs(outdir, exist_ok=True)
    out['chains'].to_csv(os.path.join(outdir, f'chains_{analysis}_{SUFFIX}.csv'), index=False)
    out['estimates'].to_csv(os.path.join(outdir, f'estimates_{analysis}_{SUFFIX}.csv'), index=False)
    out['intervals'].to_csv(os.path.join(outdir, f'intervals_{analysis}_{SUFFIX}.csv'), index=False)
    out['trajectories'].to_csv(os.path.join(outdir, f'trajectories_{analysis}_{SUFFIX}.csv'), index=False)
    if out.get('posterior') is not None:
        out['posterior'].to_csv(os.path.join(outdir, f'posterior_{analysis}_{SUFFIX}.csv.gz'), index=False, compression={'method': 'gzip', 'mtime': 0})
    if out.get('root_date_intervals') is not None:
        pd.DataFrame([out['root_date_intervals']]).to_csv(os.path.join(outdir, f'root_date_intervals_{analysis}_{SUFFIX}.csv'), index=False)
