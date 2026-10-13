#!/usr/bin/env python3
"""rt_lib.py - functions added to derive_rt.py for the 24 Sep 2026 round.

Conventions (unchanged from derive_rt.py of the pipeline unless stated)
  time               decimal years as written in the chain logs: a year has 365 days (tip dates of 2025/2026)
  grid               K parameters; interval i (1 = most recent) spans [t_anchor - i*D, t_anchor - (i-1)*D] for
                     i = 1..K-1; interval K is everything older.  Grid anchored at the last tip: t_anchor = last tip,
                     D = cutoff / (K - 1).  Grid fixed by calendar dates: t_anchor = last knot date,
                     D = (last knot - first knot) / (K - 1).
  Skygrid type       staircase (ln Ne constant inside an interval); checked from the log column skygrid.isloglinear
  growth rate        per year of the log; R = (1 + r sd^2/mu)^(mu^2/sd^2), generation time in years = days / 365.25
                     (pipeline convention).  For 1 + r sd^2/mu <= 0 the relation has the limit R = 0; such samples
                     are counted and set to 0 (the pipeline dropped them).
  window growth      W-type quantities: (ln Ne(t1) - ln Ne(t0)) / (t1 - t0) with ln Ne(t) read from the trajectory:
                       'interp'  linear interpolation between interval mid-points (primary reading); beyond the
                                 mid-point of the most recent interval the level of that interval is used
                       'stair'   level of the interval that contains the date
                       'extrap'  as 'interp', but linear extrapolation from the two most recent intervals beyond
                                 the mid-point of the most recent interval
  ESS                per chain: estimator of the pipeline (initial positive sequence).  Pooled: multi-chain estimator
                     (within- and between-chain variance, Geyer initial monotone sequence; Gelman et al., BDA3 /
                     Vehtari et al. 2021, without rank normalisation); the sum of the per-chain values is given too.
  split-R-hat        estimator of the pipeline (chains split in halves); rank-normalised split-R-hat in addition.
"""
import math
import re

import numpy as np
import pandas as pd

YEAR_DAYS = 365.0


# ----------------------------------------------------------------------------------------------------
# dates
# ----------------------------------------------------------------------------------------------------
def decimal_year(ts):
    ts = pd.Timestamp(ts)
    start = pd.Timestamp(ts.year, 1, 1)
    end = pd.Timestamp(ts.year + 1, 1, 1)
    return ts.year + (ts - start).days / (end - start).days


def dec2ts(d):
    y = int(math.floor(d))
    start = pd.Timestamp(y, 1, 1)
    nd = (pd.Timestamp(y + 1, 1, 1) - start).days
    return (start + pd.Timedelta(days=(d - y) * nd)).round("s")      # round: a date given as year + day/365 must map back to that day


def dec2date(d):
    return dec2ts(d).date()


def root_dates_from_height(root_height_years, last_tip):
    """day-resolution root date: last tip - round(rootHeight * 365) days (array of numpy datetime64[D])"""
    days = np.round(np.asarray(root_height_years, float) * YEAR_DAYS).astype(int)
    return np.datetime64(pd.Timestamp(last_tip).date()) - days.astype("timedelta64[D]")


# ----------------------------------------------------------------------------------------------------
# summaries
# ----------------------------------------------------------------------------------------------------
def hpd(x, cred=0.95):
    x = np.sort(np.asarray(x, float))
    x = x[~np.isnan(x)]
    n = len(x)
    if n == 0:
        return float("nan"), float("nan")
    m = int(np.ceil(cred * n))
    widths = x[m - 1:] - x[:n - m + 1]
    j = widths.argmin()
    return float(x[j]), float(x[j + m - 1])


def _acf_pipeline(x):
    n = len(x)
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n]
    return ac / (np.arange(n, 0, -1) * x.var())


def ess_single(x):
    """ESS estimator of the pipeline (derive_rt.ess)"""
    x = np.asarray(x, float)
    n = len(x)
    x = x - x.mean()
    if n < 4 or x.std() == 0:
        return float("nan")
    acf = _acf_pipeline(x)
    s = 0.0
    for k in range(1, n // 2):
        pair = acf[2 * k - 1] + acf[2 * k] if 2 * k < n else acf[2 * k - 1]
        if pair < 0:
            break
        s += pair
    return n / (1 + 2 * s)


def _autocov(x):
    n = len(x)
    x = x - x.mean()
    f = np.fft.rfft(x, 2 * n)
    return np.fft.irfft(f * np.conj(f))[:n] / n          # biased autocovariance, lag 0 = variance (ddof 0)


def ess_multichain(chains):
    """multi-chain ESS (BDA3 / Stan): chains = list of 1-d arrays, cut to the shortest length"""
    n = min(len(c) for c in chains)
    X = np.array([np.asarray(c[:n], float) for c in chains])
    m = X.shape[0]
    if n < 8:
        return float("nan")
    acov = np.array([_autocov(x) for x in X])
    chain_var = acov[:, 0] * n / (n - 1.0)
    W = chain_var.mean()
    if W == 0:
        return float("nan")
    var_plus = W * (n - 1.0) / n
    if m > 1:
        var_plus += X.mean(axis=1).var(ddof=1)
    rho = 1.0 - (W - acov.mean(axis=0)) / var_plus
    rho[0] = 1.0
    # Geyer: sums of adjacent pairs, initial positive then monotone
    tau = -1.0
    prev = np.inf
    t = 0
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


def rhat_split(chains):
    """split-R-hat of the pipeline (derive_rt.rhat)"""
    m = min(len(c) for c in chains)
    X = np.array([np.asarray(c[:m]) for c in chains], float)
    half = m // 2
    X = np.concatenate([X[:, :half], X[:, half:2 * half]], axis=0)
    M, N = X.shape
    W = X.var(axis=1, ddof=1).mean()
    B = N * X.mean(axis=1).var(ddof=1)
    if W == 0:
        return float("nan")
    return float(np.sqrt(((N - 1) / N * W + B / N) / W))


def rhat_rank(chains):
    """rank-normalised split-R-hat (Vehtari et al. 2021), maximum of the bulk and the folded version"""
    from scipy.stats import norm, rankdata
    m = min(len(c) for c in chains)
    X = np.array([np.asarray(c[:m]) for c in chains], float)
    half = m // 2
    X = np.concatenate([X[:, :half], X[:, half:2 * half]], axis=0)

    def z(v):
        r = rankdata(v.ravel(), method="average").reshape(v.shape)
        return norm.ppf((r - 0.375) / (v.size + 0.25))

    def rh(v):
        N = v.shape[1]
        W = v.var(axis=1, ddof=1).mean()
        B = N * v.mean(axis=1).var(ddof=1)
        return float(np.sqrt(((N - 1) / N * W + B / N) / W)) if W > 0 else float("nan")

    return max(rh(z(X)), rh(z(np.abs(X - np.median(X)))))


def summarise(chains, name=None, scale=1.0):
    """chains: list of 1-d arrays (post-burn-in samples of one quantity, one array per retained chain)"""
    allx = np.concatenate(chains) * scale
    lo, hi = hpd(allx)
    per = [ess_single(c) for c in chains]
    out = dict(quantity=name, median=float(np.nanmedian(allx)), hpd_lo=lo, hpd_hi=hi, mean=float(np.nanmean(allx)),
               sd=float(np.nanstd(allx)), n_samples=int(len(allx)), n_chains=len(chains),
               ess_pooled=ess_multichain(chains), ess_sum_of_chains=float(np.nansum(per)),
               ess_min_chain=float(np.nanmin(per)) if len(per) else float("nan"),
               rhat_split=rhat_split(chains) if len(chains) > 0 else float("nan"),
               rhat_rank=rhat_rank(chains) if len(chains) > 0 else float("nan"))
    out["mcse_median_approx"] = float(1.2533 * out["sd"] / math.sqrt(out["ess_pooled"])) if out["ess_pooled"] == out["ess_pooled"] and out["ess_pooled"] > 0 else float("nan")
    return out


# ----------------------------------------------------------------------------------------------------
# generation time
# ----------------------------------------------------------------------------------------------------
class GenTime:
    def __init__(self, mean_d=15.3, sd_d=9.3):
        self.mean_d, self.sd_d = float(mean_d), float(sd_d)
        tg, sg = mean_d / 365.25, sd_d / 365.25
        self.A = (tg / sg) ** 2       # shape
        self.B = tg / sg ** 2         # rate per year

    def r2R(self, r):
        r = np.asarray(r, float)
        base = 1.0 + r / self.B
        with np.errstate(invalid="ignore"):
            return np.where(base > 0, np.power(np.maximum(base, 0.0), self.A), 0.0)

    def n_undefined(self, r):
        return int((1.0 + np.asarray(r, float) / self.B <= 0).sum())


# ----------------------------------------------------------------------------------------------------
# grid
# ----------------------------------------------------------------------------------------------------
class Grid:
    def __init__(self, K, t_anchor, D, last_tip, anchored="last tip"):
        self.K = int(K)
        self.t_anchor = float(t_anchor)
        self.D = float(D)
        self.t_last_tip = float(decimal_year(last_tip))
        self.last_tip = pd.Timestamp(last_tip)
        self.anchored = anchored
        self.mids = self.t_anchor - (np.arange(1, self.K + 1) - 0.5) * self.D

    @classmethod
    def from_cutoff(cls, K, cutoff, last_tip):
        return cls(K, decimal_year(last_tip), float(cutoff) / (int(K) - 1), last_tip, "last tip")

    @classmethod
    def from_knots(cls, K, first_knot, last_knot, last_tip):
        t1, t0 = decimal_year(last_knot), decimal_year(first_knot)
        return cls(K, t1, (t1 - t0) / (int(K) - 1), last_tip, f"calendar knots {first_knot} .. {last_knot}")

    def interval_of(self, t):
        i = int(np.floor((self.t_anchor - t) / self.D + 1e-12)) + 1
        return min(max(i, 1), self.K)

    def edges(self, i):
        """(older edge, newer edge) of interval i in decimal years; the older edge of interval K is -inf"""
        newer = self.t_anchor - (i - 1) * self.D
        older = self.t_anchor - i * self.D if i < self.K else -np.inf
        return older, newer

    def weights(self, t, mode="interp"):
        """vector w with ln Ne(t) = w . lnNe (lnNe[0] = interval 1)"""
        w = np.zeros(self.K)
        if mode == "stair":
            w[self.interval_of(t) - 1] = 1.0
            return w
        m = self.mids
        if t >= m[0]:
            if mode == "extrap":
                x = (t - m[1]) / (m[0] - m[1])          # > 1 beyond the most recent mid-point
                w[0], w[1] = x, 1.0 - x
            else:
                w[0] = 1.0
            return w
        if t <= m[self.K - 1]:
            w[self.K - 1] = 1.0
            return w
        j = int(np.floor((m[0] - t) / self.D))            # t lies between mids[j] (newer) and mids[j+1] (older)
        j = min(max(j, 0), self.K - 2)
        x = (m[j] - t) / self.D
        w[j], w[j + 1] = 1.0 - x, x
        return w

    def table(self):
        rows = []
        for i in range(1, self.K + 1):
            older, newer = self.edges(i)
            rows.append(dict(interval=i,
                             older_edge=str(dec2date(older)) if i < self.K else "open (before the cutoff)",
                             newer_edge=str(dec2date(newer)), mid_point=str(dec2date(self.mids[i - 1])) if i < self.K else "",
                             older_edge_decimal=older if i < self.K else np.nan, newer_edge_decimal=newer,
                             width_days=self.D * YEAR_DAYS if i < self.K else np.nan))
        return pd.DataFrame(rows)


def window_contrast(grid, t0, t1, mode="interp"):
    """weights c with growth rate [per year] = c . lnNe for the window t0 -> t1 (decimal years)"""
    return (grid.weights(t1, mode) - grid.weights(t0, mode)) / (t1 - t0)


def pair_contrast(grid, i):
    """adjacent intervals: older interval i+1 -> newer interval i"""
    c = np.zeros(grid.K)
    c[i - 1], c[i] = 1.0 / grid.D, -1.0 / grid.D
    return c


def prior_lnNe_samples(K, tau, n=200000, seed=20260924):
    """GMRF prior alone: first-order random walk, increments N(0, 1/tau); the level of interval K is set to 0
    (the prior is flat in the overall level, contrasts do not depend on it)"""
    rng = np.random.default_rng(seed)
    inc = rng.normal(0.0, 1.0 / math.sqrt(tau), size=(n, K - 1))      # inc[:, j] = lnNe_{j+1} - lnNe_{j+2} (newer - older)
    g = np.zeros((n, K))
    g[:, :K - 1] = np.cumsum(inc[:, ::-1], axis=1)[:, ::-1]
    return g


def prior_sd_contrast(c, tau):
    """exact prior SD of c . lnNe for a contrast (sum of c = 0)"""
    cum = np.cumsum(c)[:-1]
    return float(math.sqrt((cum ** 2).sum() / tau))


# ----------------------------------------------------------------------------------------------------
# logs
# ----------------------------------------------------------------------------------------------------
def load_log(path):
    df = pd.read_csv(path, sep="\t", comment="#")
    return df.dropna(axis=1, how="all")


def lp_columns(df):
    cols = [c for c in df.columns if c.startswith("skygrid.logPopSize")]
    return sorted(cols, key=lambda c: int(c.replace("skygrid.logPopSize", "") or 1))


def tip_dates_from_fasta(path):
    d = []
    with open_any(path) as fh:
        for line in fh:
            if line.startswith(">"):
                d.append(line[1:].strip().split("|")[-1])
    dd = pd.to_datetime(pd.Series(d), format="ISO8601", errors="coerce")
    assert dd.notna().all(), "FASTA headers must end with |YYYY-MM-DD"
    return dd


_PROGRESS = re.compile(r"Step (\d+), log_posterior = ([-\d.e+]+), log_G = ([-\d.e+]+), log_coal = ([-\d.e+]+), "
                       r"log_other_priors = ([-\d.e+]+), num_muts = (\d+), T = ([-\d.e+]+), t_MRCA = (\d{4}-\d{2}-\d{2})")


def open_any(path, mode="rt"):
    import gzip
    return gzip.open(path, mode, errors="replace") if str(path).endswith(".gz") else open(path, mode.replace("t", ""), errors="replace")


def parse_delphy_stdout(path):
    """progress lines of a Delphy stdout file, or the reduced table written by compact_outputs.py (*.progress.tsv.gz)"""
    if str(path).endswith((".tsv", ".tsv.gz")):
        return pd.read_csv(path, sep="\t")
    rows = []
    with open_any(path) as fh:
        for line in fh:
            m = _PROGRESS.search(line)
            if m:
                g = m.groups()
                rows.append((int(g[0]), float(g[1]), float(g[2]), float(g[3]), float(g[4]), int(g[5]), float(g[6]), g[7]))
    return pd.DataFrame(rows, columns=["state", "log_posterior", "log_G", "log_coal", "log_other_priors", "num_muts",
                                       "T", "t_MRCA"])


# ----------------------------------------------------------------------------------------------------
# trees (Delphy / BEAST NEXUS with numbered taxa)
# ----------------------------------------------------------------------------------------------------
_ANN = re.compile(r"\[&[^\]]*\]")


def read_nexus_trees(path, want_states=None, first_only=False):
    labels, taxl, trees = {}, [], []
    in_tr = in_tax = False
    with open_any(path) as fh:
        for line in fh:
            ls = line.strip()
            low = ls.lower()
            if low.startswith("taxlabels"):
                in_tax = True
                continue
            if in_tax:
                if ls.startswith(";") or ls == "":
                    in_tax = False
                else:
                    taxl.append(ls.rstrip(";").strip("'"))
                    if ls.endswith(";"):
                        in_tax = False
                continue
            if low.startswith("translate"):
                in_tr = True
                continue
            if in_tr:
                if ls.startswith(";"):
                    in_tr = False
                    continue
                m = re.match(r"(\d+)\s+(\S+?),?;?$", ls)
                if m:
                    labels[int(m.group(1))] = m.group(2).strip("'").rstrip(",")
                if ls.endswith(";"):
                    in_tr = False
                continue
            if low.startswith("tree "):
                m = re.match(r"tree\s+(\S+)\s*(\[[^\]]*\])?\s*=\s*(\[&[RU]\])?\s*(.*);\s*$", ls)
                name = m.group(1)
                ms = re.match(r"STATE_(\d+)", name)
                st = int(ms.group(1)) if ms else None
                if want_states is not None and st not in want_states:
                    continue
                trees.append((st if st is not None else name, _ANN.sub("", m.group(4))))
                if first_only:
                    break
    if not labels:
        labels = {i + 1: l for i, l in enumerate(taxl)}
    return labels, trees


def parse_newick(newick):
    """returns nodes = list of [children, label, branch length]; children are created before parents; root = last"""
    s = newick
    i, n = 0, len(s)
    nodes, stack = [], []
    root = None
    while i < n:
        ch = s[i]
        if ch == "(":
            stack.append([])
            i += 1
        elif ch == ",":
            i += 1
        elif ch == ")":
            kids = stack.pop()
            i += 1
            j = i
            while j < n and s[j] not in ":,()":
                j += 1
            i = j
            bl = 0.0
            if i < n and s[i] == ":":
                j = i + 1
                while j < n and s[j] not in ",()":
                    j += 1
                bl = float(s[i + 1:j])
                i = j
            nodes.append([kids, None, bl])
            if stack:
                stack[-1].append(len(nodes) - 1)
            else:
                root = len(nodes) - 1
        elif ch in " \t;":
            i += 1
        else:
            j = i
            while j < n and s[j] not in ":,()":
                j += 1
            lab = s[i:j]
            i = j
            bl = 0.0
            if i < n and s[i] == ":":
                j = i + 1
                while j < n and s[j] not in ",()":
                    j += 1
                bl = float(s[i + 1:j])
                i = j
            nodes.append([[], lab, bl])
            stack[-1].append(len(nodes) - 1)
    return nodes, root


def node_times(newick, labels, dec_of_label, tol=5e-4):
    nodes, root = parse_newick(newick)
    t = [None] * len(nodes)
    tips, ints = [], []
    for k, (kids, lab, bl) in enumerate(nodes):
        if not kids:
            key = labels[int(lab)] if lab.isdigit() else lab.strip("'")
            t[k] = dec_of_label[key]
            tips.append(t[k])
        else:
            cand = [t[c] - nodes[c][2] for c in kids]
            assert max(cand) - min(cand) < tol, f"inconsistent node time (spread {max(cand) - min(cand):.2e} yr)"
            t[k] = float(np.mean(cand))
            ints.extend([t[k]] * (len(kids) - 1))          # a node with m children counts as m-1 coalescent events
    return np.array(tips), np.array(ints)


def skygrid_tree_terms(tips, ints, grid):
    """per interval: number of coalescent events and exposure integral of C(k,2) dt (pair-years)"""
    K = grid.K
    ev = np.concatenate([np.column_stack([tips, np.ones(len(tips))]), np.column_stack([ints, -np.ones(len(ints))])])
    bounds = grid.t_anchor - grid.D * np.arange(1, K)
    ev = np.concatenate([ev, np.column_stack([bounds, np.zeros(K - 1)])])
    ev = ev[np.argsort(-ev[:, 0], kind="stable")]
    k = 0
    tprev = ev[0, 0]
    n_i, A_i = np.zeros(K), np.zeros(K)
    for tt, typ in ev:
        dt = tprev - tt
        if dt > 0 and k >= 2:
            A_i[grid.interval_of(0.5 * (tprev + tt)) - 1] += k * (k - 1) / 2.0 * dt
        if typ == 1:
            k += 1
        elif typ == -1:
            n_i[grid.interval_of(tt) - 1] += 1
            k -= 1
        tprev = tt
    return n_i, A_i


def exact_skygrid_logdensity(n_i, A_i, lnN, time_unit_days=1.0):
    """exact log density of the coalescent times under a staircase N(t); time unit of the density: days by default
    (Delphy's log_coal), years for time_unit_days = 365"""
    lnN = np.asarray(lnN, float)
    return float(-(n_i * lnN).sum() - (A_i / np.exp(lnN)).sum() - n_i.sum() * math.log(YEAR_DAYS / time_unit_days))


# ----------------------------------------------------------------------------------------------------
# Delphy prior components (verified against the log columns of Delphy 1.4.1, see degenerate_state report)
#   skygrid (log column)   = log_coal + (n_tips - 1) ln(365) + (K-1)/2 ln(tau / 2 pi) - tau/2 sum (lnNe_i - lnNe_i+1)^2
#   log_other_priors       = - tau/2 sum (lnNe_i - lnNe_i+1)^2 + barrier + priors of the other parameters (+ constant)
#   barrier                = - sum_i ( min(lnNe_i - ln(loc), 0) / ln(1 - scale) )^2 ,  loc = 1/365 yr, scale = 0.30 by default
# ----------------------------------------------------------------------------------------------------
def gmrf_term(lnNe, tau):
    lnNe = np.asarray(lnNe, float)
    return -0.5 * tau * ((lnNe[:, :-1] - lnNe[:, 1:]) ** 2).sum(axis=1)


def barrier_term(lnNe, loc_years=1.0 / 365.0, scale=0.30):
    x = np.minimum(np.asarray(lnNe, float) - math.log(loc_years), 0.0) / math.log(1.0 - scale)
    return -(x ** 2).sum(axis=1)
