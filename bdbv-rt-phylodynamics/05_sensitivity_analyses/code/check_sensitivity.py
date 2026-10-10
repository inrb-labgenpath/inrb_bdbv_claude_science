#!/usr/bin/env python3
"""check_sensitivity.py - checks of one analysis that do NOT call the functions of analysis.py and rt_lib.py:
  recompute()   the estimates of the table of estimates, recomputed from the posterior file (median, mean, shortest 95 % interval,
                probabilities, pooled effective sample size, split R-hat) with code written for this check; a second interval of
                exactly the same width is reported as a tie
  halves()      the two halves of the retained chains (in the order of their seeds: the first floor(n/2) chains against the rest):
                difference of the medians in per cent of the width of the 95 % interval of the whole analysis
The definitions are those of analysis.py and rt_lib.py (year of 365 days; generation time in years = days / 365.25; ln(Ne tau) read by
linear interpolation between interval mid-points, beyond the most recent mid-point the level of the most recent interval)."""
import math

import numpy as np
import pandas as pd

MEAN_D, SD_D = 15.3, 9.3


def dec_year(date):
    ts = pd.Timestamp(date)
    y0 = pd.Timestamp(ts.year, 1, 1)
    return ts.year + (ts - y0).days / (pd.Timestamp(ts.year + 1, 1, 1) - y0).days


def shortest_interval(x, cred=0.95):
    """(lower, upper, number of intervals of the same shortest width, lowest lower bound among them, highest lower bound among them)"""
    v = np.sort(np.asarray(x, float))
    n = v.size
    m = int(math.ceil(cred * n))
    best, where = None, []
    for j in range(0, n - m + 1):
        w = v[j + m - 1] - v[j]
        if best is None or w < best:
            best, where = w, [j]
        elif w == best:
            where.append(j)
    lows = sorted(set(float(v[j]) for j in where))
    return float(v[where[0]]), float(v[where[0] + m - 1]), len(lows), lows[0], lows[-1], float(best)


def r_to_R(r_per_day, mean_d=MEAN_D, sd_d=SD_D):
    tg, sg = mean_d / 365.25, sd_d / 365.25
    shape, rate = (tg / sg) ** 2, tg / sg ** 2
    b = 1.0 + (np.asarray(r_per_day, float) * 365.0) / rate
    out = np.zeros_like(b)
    pos = b > 0
    out[pos] = b[pos] ** shape
    return out


def split_rhat(chains):
    n = min(len(c) for c in chains)
    h = n // 2
    parts = []
    for c in chains:
        c = np.asarray(c[:n], float)
        parts += [c[:h], c[h:2 * h]]
    means = np.array([p.mean() for p in parts])
    vars_ = np.array([p.var(ddof=1) for p in parts])
    W = vars_.mean()
    if W == 0:
        return float('nan')
    B = h * means.var(ddof=1)
    return float(math.sqrt(((h - 1) / h * W + B / h) / W))


def pooled_ess(chains):
    """multi-chain effective sample size (within- and between-chain variance, Geyer's initial monotone sequence of sums of pairs)"""
    n = min(len(c) for c in chains)
    m = len(chains)
    if n < 8:
        return float('nan')
    X = np.array([np.asarray(c[:n], float) for c in chains])
    ac = np.zeros((m, n))
    for k in range(m):
        d = X[k] - X[k].mean()
        f = np.fft.rfft(d, 2 * n)
        ac[k] = np.fft.irfft(f * np.conj(f))[:n] / n
    cv = ac[:, 0] * n / (n - 1.0)
    W = cv.mean()
    if W == 0:
        return float('nan')
    vp = W * (n - 1.0) / n + (X.mean(axis=1).var(ddof=1) if m > 1 else 0.0)
    rho = 1.0 - (W - ac.mean(axis=0)) / vp
    rho[0] = 1.0
    tau, prev, t = -1.0, float('inf'), 0
    while t + 1 < n:
        pair = rho[t] + rho[t + 1]
        if pair < 0:
            break
        pair = min(pair, prev)
        prev = pair
        tau += 2.0 * pair
        t += 2
    tau = max(tau, 1.0 / math.log10(max(m * n, 10)))
    return float(m * n / tau)


def ln_ne_at(P, t, mids, K):
    """ln(Ne tau) at decimal year t for every state: linear interpolation between the mid-points; rows of P = states"""
    cols = np.column_stack([P[f'ln_Ne_tau_{i}'].values for i in range(1, K + 1)])
    order = np.argsort(mids)
    xs = np.asarray(mids)[order]
    out = np.empty(cols.shape[0])
    if t >= xs[-1]:
        return cols[:, order[-1]].copy()
    if t <= xs[0]:
        return cols[:, order[0]].copy()
    j = int(np.searchsorted(xs, t, side='right')) - 1
    f = (t - xs[j]) / (xs[j + 1] - xs[j])
    out = (1 - f) * cols[:, order[j]] + f * cols[:, order[j + 1]]
    return out


def close(a, b, rel=1e-9):
    try:
        a, b = float(a), float(b)
    except (TypeError, ValueError):
        return str(a) == str(b)
    if a != a and b != b:
        return True
    return abs(a - b) <= rel * max(1.0, abs(a), abs(b))


def recompute(posterior, estimates, last_tip, anchor_decimal, D, K, periods):
    """periods: dict P -> (start date, end as decimal year or date, evaluated?).  Returns one row per compared cell."""
    P = posterior
    seeds = list(dict.fromkeys(P['seed']))
    by = {s: P[P['seed'] == s] for s in seeds}
    mids = anchor_decimal - (np.arange(1, K + 1) - 0.5) * D
    E = estimates.set_index('quantity')
    rows = []

    def put(q, what, mine, stored, tie=''):
        ok = close(mine, stored)
        rows.append(dict(quantity=q, value=what, recomputed=mine, in_the_table_of_estimates=stored,
                         agrees=('yes' if ok else ('tie: another interval of exactly the same width' if tie else 'no')), note=tie))

    def series_checks(q, chains, scale=1.0, bounds=True, diag=True, diag_cols=('ess_pooled', 'rhat_split')):
        allx = np.concatenate(chains) * scale
        put(q, 'median', float(np.median(allx)), E.loc[q, 'median'])
        put(q, 'mean', float(np.mean(allx)), E.loc[q, 'mean'])
        if bounds:
            lo, hi, n_same, low_first, low_last, width = shortest_interval(allx)
            s_lo, s_hi = float(E.loc[q, 'hpd95_lower']), float(E.loc[q, 'hpd95_upper'])
            tie = ''
            if not (close(lo, s_lo) and close(hi, s_hi)) and close(hi - lo, s_hi - s_lo, 1e-12):
                tie = f'{n_same} intervals of the width {width!r}'
            put(q, 'lower bound of the 95 % interval', lo, s_lo, tie)
            put(q, 'upper bound of the 95 % interval', hi, s_hi, tie)
        if diag:
            put(q, diag_cols[0], pooled_ess(chains), E.loc[q, diag_cols[0]])
            put(q, diag_cols[1], split_rhat(chains), E.loc[q, diag_cols[1]])

    R = {}
    for p, (start, end_dec, evaluated) in periods.items():
        if not evaluated:
            continue
        t0 = dec_year(start)
        g = {s: (ln_ne_at(by[s], end_dec, mids, K) - ln_ne_at(by[s], t0, mids, K)) / (end_dec - t0) / 365.0 for s in seeds}
        col = f'growth_rate_per_day_{p}' + ('_own_end' if p in ('P3', 'P6') else '')
        worst = max(float(np.max(np.abs(g[s] - by[s][col].values))) for s in seeds)
        rows.append(dict(quantity=f'growth rate, {p}', value='largest difference between the growth rate recomputed from ln(Ne tau) and the column of the posterior file',
                         recomputed=worst, in_the_table_of_estimates=0.0, agrees='yes' if worst < 1e-9 else 'no', note=''))
        series_checks(f'growth rate, {p}', [g[s] for s in seeds])
        R[p] = {s: r_to_R(g[s]) for s in seeds}
        q = f'reproduction number, {p}'
        allR = np.concatenate([R[p][s] for s in seeds])
        put(q, 'median', float(np.median(allR)), E.loc[q, 'median'])
        put(q, 'mean', float(np.mean(allR)), E.loc[q, 'mean'])
        lo, hi, n_same, _, _, width = shortest_interval(allR)
        s_lo, s_hi = float(E.loc[q, 'hpd95_lower']), float(E.loc[q, 'hpd95_upper'])
        tie = f'{n_same} intervals of the width {width!r}' if (not (close(lo, s_lo) and close(hi, s_hi)) and close(hi - lo, s_hi - s_lo, 1e-12)) else ''
        put(q, 'lower bound of the 95 % interval', lo, s_lo, tie)
        put(q, 'upper bound of the 95 % interval', hi, s_hi, tie)
        put(q, 'states_with_R_set_to_0', int((allR == 0).sum()), E.loc[q, 'states_with_R_set_to_0'])
        put(q, 'ess_pooled_of_the_reproduction_number_series', pooled_ess([R[p][s] for s in seeds]), E.loc[q, 'ess_pooled_of_the_reproduction_number_series'])
        put(q, 'rhat_split_of_the_reproduction_number_series', split_rhat([R[p][s] for s in seeds]), E.loc[q, 'rhat_split_of_the_reproduction_number_series'])
        if p == 'P1':
            allg = np.concatenate([g[s] for s in seeds])
            med = float(np.median(allg))
            glo, ghi, *_ = shortest_interval(allg)
            put('doubling time, P1', 'median', math.log(2) / med if med > 0 else float('nan'), E.loc['doubling time, P1', 'median'])
            put('doubling time, P1', 'lower bound of the 95 % interval', math.log(2) / ghi if glo > 0 else float('nan'), E.loc['doubling time, P1', 'hpd95_lower'])
            put('doubling time, P1', 'upper bound of the 95 % interval', math.log(2) / glo if glo > 0 else float('nan'), E.loc['doubling time, P1', 'hpd95_upper'])
    series_checks('evolutionary rate', [by[s]['evolutionary_rate'].values for s in seeds])
    series_checks('kappa', [by[s]['kappa'].values for s in seeds])
    # root date
    days = {s: by[s]['root_height_years'].values * 365.0 for s in seeds}
    alld = np.concatenate([days[s] for s in seeds])
    f = lambda d: str((pd.Timestamp(last_tip) - pd.Timedelta(days=int(round(d)))).date())
    lo, hi, n_same, low_first, low_last, width = shortest_interval(alld)
    qd, qy = 'date of the most recent common ancestor (date)', 'date of the most recent common ancestor (decimal year)'
    put(qd, 'median', f(float(np.median(alld))), E.loc[qd, 'median'])
    put(qd, 'mean', f(float(np.mean(alld))), E.loc[qd, 'mean'])
    tie = ''
    if not (f(hi) == E.loc[qd, 'hpd95_lower'] and f(lo) == E.loc[qd, 'hpd95_upper']):
        w_st = (pd.Timestamp(E.loc[qd, 'hpd95_upper']) - pd.Timestamp(E.loc[qd, 'hpd95_lower'])).days
        if abs(w_st - (hi - lo)) < 0.5:
            tie = f'{n_same} intervals of the width {width!r} days'
    put(qd, 'lower bound of the 95 % interval (earlier date)', f(hi), E.loc[qd, 'hpd95_lower'], tie)
    put(qd, 'upper bound of the 95 % interval (later date)', f(lo), E.loc[qd, 'hpd95_upper'], tie)
    t_last = dec_year(last_tip)
    put(qy, 'median', t_last - float(np.median(alld)) / 365.0, E.loc[qy, 'median'])
    put(qy, 'mean', t_last - float(np.mean(alld)) / 365.0, E.loc[qy, 'mean'])
    put(qy, 'lower bound of the 95 % interval', t_last - hi / 365.0, E.loc[qy, 'hpd95_lower'], tie)
    put(qy, 'upper bound of the 95 % interval', t_last - lo / 365.0, E.loc[qy, 'hpd95_upper'], tie)
    put(qd, 'ess_pooled', pooled_ess([by[s]['root_height_years'].values for s in seeds]), E.loc[qd, 'ess_pooled'])
    put(qd, 'rhat_split', split_rhat([by[s]['root_height_years'].values for s in seeds]), E.loc[qd, 'rhat_split'])
    # probabilities
    pr = []
    if 'P1' in R and 'P2' in R:
        pr.append(('P(R of P2 < R of P1)', {s: (R['P2'][s] < R['P1'][s]).astype(float) for s in seeds}))
    if 'P1' in R:
        pr.append(('P(R of P1 > 1)', {s: (R['P1'][s] > 1).astype(float) for s in seeds}))
    if 'P2' in R:
        pr.append(('P(R of P2 > 1)', {s: (R['P2'][s] > 1).astype(float) for s in seeds}))
    if 'P3' in R:
        pr.append(('P(R of P3 < 1)', {s: (R['P3'][s] < 1).astype(float) for s in seeds}))
    for q, ind in pr:
        put(q, 'probability (column mean)', float(np.mean(np.concatenate([ind[s] for s in seeds]))), E.loc[q, 'mean'])
        if E.loc[q, 'ess_pooled'] not in ('', None) and str(E.loc[q, 'ess_pooled']) != 'nan':
            put(q, 'ess_pooled', pooled_ess([ind[s] for s in seeds]), E.loc[q, 'ess_pooled'])
            put(q, 'rhat_split', split_rhat([ind[s] for s in seeds]), E.loc[q, 'rhat_split'])
        else:
            const = any(np.std(ind[s]) == 0 for s in seeds)
            rows.append(dict(quantity=q, value='the indicator is constant in at least one chain', recomputed='yes' if const else 'no', in_the_table_of_estimates='yes',
                             agrees='yes' if const else 'no', note=''))
    return pd.DataFrame(rows)


def halves(posterior, estimates, last_tip):
    P = posterior
    seeds = sorted(dict.fromkeys(P['seed']))
    n = len(seeds)
    E = estimates.set_index('quantity')
    rows = []
    if n < 2:
        for q in ('reproduction number, P1', 'reproduction number, P2', 'reproduction number, P3', 'evolutionary rate', 'date of the most recent common ancestor (decimal year)'):
            rows.append(dict(quantity=q, chains_in_the_first_half=n, chains_in_the_second_half=0, difference_of_the_medians_in_per_cent_of_the_width_of_the_95_interval='not assessed: fewer than 2 retained chains'))
        return pd.DataFrame(rows)
    h1, h2 = seeds[:n // 2], seeds[n // 2:]
    a, b = P[P['seed'].isin(h1)], P[P['seed'].isin(h2)]
    spec = [('reproduction number, P1', 'reproduction_number_P1', 1.0), ('reproduction number, P2', 'reproduction_number_P2', 1.0),
            ('reproduction number, P3', 'reproduction_number_P3_own_end', 1.0), ('evolutionary rate', 'evolutionary_rate', 1.0),
            ('date of the most recent common ancestor (decimal year)', 'root_height_years', -1.0)]
    for q, col, sign in spec:
        try:
            width = float(E.loc[q, 'hpd95_upper']) - float(E.loc[q, 'hpd95_lower'])
            m1, m2 = float(np.median(a[col].values)), float(np.median(b[col].values))
            if np.isnan(m1) or np.isnan(m2):
                raise ValueError('the period is not evaluated')
            diff = sign * (m2 - m1)
            rows.append(dict(quantity=q, chains_in_the_first_half=len(h1), chains_in_the_second_half=len(h2), seeds_of_the_first_half=' '.join(map(str, h1)),
                             seeds_of_the_second_half=' '.join(map(str, h2)), median_of_the_first_half=(m1 if sign > 0 else dec_year(last_tip) - m1),
                             median_of_the_second_half=(m2 if sign > 0 else dec_year(last_tip) - m2), difference_second_minus_first=diff,
                             width_of_the_95_interval_of_the_whole_analysis=width,
                             difference_of_the_medians_in_per_cent_of_the_width_of_the_95_interval=100.0 * diff / width if width > 0 else float('nan'),
                             unit=('years' if sign < 0 else str(E.loc[q, 'unit']))))
        except (KeyError, ValueError) as exc:
            rows.append(dict(quantity=q, chains_in_the_first_half=len(h1), chains_in_the_second_half=len(h2),
                             difference_of_the_medians_in_per_cent_of_the_width_of_the_95_interval=f'not assessed: {exc}'))
    return pd.DataFrame(rows)
