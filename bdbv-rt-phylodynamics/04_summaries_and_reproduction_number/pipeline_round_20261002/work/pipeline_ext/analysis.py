#!/usr/bin/env python3
"""analysis.py - chain handling and derived quantities (24 Sep 2026 round); uses rt_lib.py.

Pre-specified windows (fixed before any 24 Sep estimate was computed)
  W1   1 Mar 2026 -> 15 May 2026                      early-phase R0 (primary definition)
  W2   15 May 2026 -> 15 Jun 2026
  W3   15 Jun 2026 -> mid-point of the most recent interval (last date at which the trajectory is defined
       without extrapolation)
  W3b  15 Jun 2026 -> 15 Aug 2026
  reading of ln Ne at a date: linear interpolation between interval mid-points ('interp').
"""
import math
import os

import numpy as np
import pandas as pd

import rt_lib as L

W_DATES = {"W1": ("2026-03-01", "2026-05-15"), "W2": ("2026-05-15", "2026-06-15"),
           "W3": ("2026-06-15", None), "W3b": ("2026-06-15", "2026-08-15")}
R0_START, R0_END = "2026-03-01", "2026-05-15"


# ----------------------------------------------------------------------------------------------------
class Chain:
    def __init__(self, engine, log, fasta, seed, stdout=None, burnin=0.30, collapse_days=1.0, label=None,
                 extra_logs=None, first_tip=None, last_tip=None):
        self.engine, self.log, self.fasta, self.seed, self.stdout = engine, log, fasta, seed, stdout
        self.label = label or os.path.splitext(os.path.basename(log))[0]
        full = L.load_log(log)
        if extra_logs:                                   # resumed BEAST segments: states continue
            for p in extra_logs:
                nxt = L.load_log(p)
                nxt = nxt[nxt["state"] > full["state"].iloc[-1]]
                full = pd.concat([full, nxt], ignore_index=True)
        self.full = full
        if first_tip is None:
            dd = L.tip_dates_from_fasta(fasta)
            first_tip, last_tip = dd.min(), dd.max()
        self.first_tip, self.last_tip = pd.Timestamp(first_tip), pd.Timestamp(last_tip)
        self.rh_col = "rootHeight" if engine == "delphy" else "treeModel.rootHeight"
        rd = L.root_dates_from_height(full[self.rh_col].values, self.last_tip)
        thr = np.datetime64((self.first_tip - pd.Timedelta(days=collapse_days)).date())
        coll = rd >= thr
        self.root_dates_full = rd
        self.collapsed_mask = coll
        self.n_logged = len(full)
        self.nb = int(self.n_logged * burnin)
        self.post = full.iloc[self.nb:].reset_index(drop=True)
        self.lp_cols = L.lp_columns(full)
        self.is_sky = len(self.lp_cols) > 0
        self.info = dict(chain=self.label, engine=engine, seed=seed, n_states_logged=self.n_logged,
                         last_state=int(full["state"].iloc[-1]), n_post_burnin=len(self.post),
                         collapse_threshold=str(thr), collapsed=bool(coll.any()),
                         first_collapsed_state=(int(full["state"].values[np.argmax(coll)]) if coll.any() else np.nan),
                         n_collapsed_states=int(coll.sum()),
                         n_regular_states_after_first_collapse=(int((~coll[np.argmax(coll):]).sum()) if coll.any() else np.nan),
                         latest_root_date=str(rd.max()))
        if self.is_sky:
            lp = full[self.lp_cols].values
            self.info.update(min_lnNe_any_interval=float(lp.min()),
                             interval_of_min_lnNe=int(np.unravel_index(lp.argmin(), lp.shape)[1] + 1),
                             skygrid_precision_min=float(full["skygrid.precision"].min()),
                             skygrid_precision_max=float(full["skygrid.precision"].max()))
            if "skygrid.isloglinear" in full.columns:
                self.info["skygrid_type"] = "log-linear" if full["skygrid.isloglinear"].max() > 0 else "staircase"
        self.retained = not self.info["collapsed"]
        self.why = "regular" if self.retained else "entered the degenerate state (root date >= first tip - 1 d)"


class Analysis:
    def __init__(self, name, engine, model, dataset, fasta, chains, grid=None, gen=None, spec=""):
        self.name, self.engine, self.model, self.dataset, self.fasta = name, engine, model, dataset, fasta
        self.chains = chains
        self.grid = grid
        self.gen = gen or L.GenTime()
        self.spec = spec
        dd = L.tip_dates_from_fasta(fasta)
        self.tip_dates = dd
        self.n_tips, self.first_tip, self.last_tip = len(dd), dd.min(), dd.max()
        self.t_last = L.decimal_year(self.last_tip)

    @property
    def kept(self):
        return [c for c in self.chains if c.retained]

    # ------------------------------------------------------------------ series per retained chain
    def series(self, col, chains=None):
        return [c.post[col].values.astype(float) for c in (chains or self.kept)]

    def lnNe(self, chains=None):
        return [c.post[c.lp_cols].values.astype(float) for c in (chains or self.kept)]

    def root_height(self, chains=None):
        return [c.post[c.rh_col].values.astype(float) for c in (chains or self.kept)]

    def contrast_series(self, cvec, chains=None):
        return [lp @ cvec for lp in self.lnNe(chains)]

    def w3_end(self):
        return float(self.grid.mids[0])

    def windows(self, mode="interp", w3_variant="midpoint"):
        """dict name -> (t0, t1, contrast vector, label)"""
        g = self.grid
        out = {}
        for k, (a, b) in W_DATES.items():
            t0 = L.decimal_year(a)
            if k == "W3":
                if w3_variant == "midpoint":
                    t1 = self.w3_end()
                    c = L.window_contrast(g, t0, t1, mode)
                    lab = f"{a} -> {L.dec2date(t1)} (mid-point of the most recent interval)"
                elif w3_variant == "constant_to_last_tip":
                    t1 = self.t_last
                    c = L.window_contrast(g, t0, t1, "interp" if mode != "stair" else "stair")
                    lab = f"{a} -> {self.last_tip.date()} (last tip; level of the most recent interval held constant)"
                else:
                    t1 = self.t_last
                    c = L.window_contrast(g, t0, t1, "extrap")
                    lab = f"{a} -> {self.last_tip.date()} (last tip; linear extrapolation from the two most recent intervals)"
            else:
                t1 = L.decimal_year(b)
                c = L.window_contrast(g, t0, t1, mode)
                lab = f"{a} -> {b}"
            ok = (t1 <= self.t_last + 1e-9) and (t1 > t0)
            note = ""
            if t1 > g.mids[0] + 1e-9 and k != "W3":
                note = "end date later than the mid-point of the most recent interval: level of that interval used"
            if not ok:
                note = "window ends after the last tip: not evaluated"
            out[k] = dict(t0=t0, t1=t1, c=c, label=lab, days=(t1 - t0) * L.YEAR_DAYS, evaluable=ok, note=note)
        return out


# ----------------------------------------------------------------------------------------------------
def r_summary(chains_r, gen, name, extra=None):
    """growth rate (per year) -> summary of r per day and of R"""
    allr = np.concatenate(chains_r)
    R = gen.r2R(allr)
    lo, hi = L.hpd(R)
    rlo, rhi = L.hpd(allr / L.YEAR_DAYS)
    s = L.summarise(chains_r, name)
    med_r = float(np.median(allr)) / L.YEAR_DAYS
    out = dict(quantity=name, R_median=float(np.median(R)), R_hpd_lo=lo, R_hpd_hi=hi, P_R_gt_1=float((allr > 0).mean()),
               R_sd=float(np.std(R)), growth_per_day_median=med_r, growth_per_day_hpd_lo=rlo, growth_per_day_hpd_hi=rhi,
               growth_per_day_sd=float(np.std(allr)) / L.YEAR_DAYS,
               doubling_time_d=(math.log(2) / med_r if med_r > 0 else np.nan),
               doubling_time_d_lo=(math.log(2) / rhi if rlo > 0 else np.nan),
               doubling_time_d_hi=(math.log(2) / rlo if rlo > 0 else np.nan),
               halving_time_d=(math.log(2) / -med_r if med_r < 0 else np.nan),
               n_samples_R_set_to_0=gen.n_undefined(allr), n_samples=int(len(allr)), n_chains=len(chains_r),
               ess_pooled=s["ess_pooled"], ess_sum_of_chains=s["ess_sum_of_chains"], ess_min_chain=s["ess_min_chain"],
               rhat_split=s["rhat_split"], rhat_rank=s["rhat_rank"])
    if extra:
        out.update(extra)
    return out


def prior_r_summary(cvec, tau, gen, K, nsim=200000):
    g = L.prior_lnNe_samples(K, tau, nsim)
    r = g @ cvec
    R = gen.r2R(r)
    lo, hi = L.hpd(R)
    q = np.quantile(R, [0.025, 0.975])
    return dict(prior_sd_growth_per_day=L.prior_sd_contrast(cvec, tau) / L.YEAR_DAYS,
                prior_sd_growth_per_day_simulated=float(np.std(r)) / L.YEAR_DAYS,
                prior_R_median=float(np.median(R)), prior_R_hpd_lo=lo, prior_R_hpd_hi=hi,
                prior_R_q025=float(q[0]), prior_R_q975=float(q[1]), prior_P_R_gt_1=float((r > 0).mean()),
                prior_n_sim=nsim)


def flag_ratio(x):
    if x != x:
        return ""
    if x >= 0.9:
        return "prior-like (posterior not narrower than the smoothing prior)"
    if x >= 0.7:
        return "weakly constrained by data"
    return "constrained by data"


def tmrca_summary(an, chains=None):
    rh = an.root_height(chains)
    allrh = np.concatenate(rh)
    days = allrh * L.YEAR_DAYS
    med = float(np.median(days))
    lo, hi = L.hpd(days)                       # lo = smallest height = latest date
    s = L.summarise(rh, "root height (yr)")
    f = lambda d: str((an.last_tip - pd.Timedelta(days=int(round(d)))).date())
    return dict(tmrca_median=f(med), tmrca_hpd_early=f(hi), tmrca_hpd_late=f(lo), root_height_days_median=med,
                root_height_days_hpd_lo=lo, root_height_days_hpd_hi=hi, root_height_days_sd=float(np.std(days)),
                tmrca_decimal_median=an.t_last - med / L.YEAR_DAYS, tmrca_decimal_early=an.t_last - hi / L.YEAR_DAYS,
                tmrca_decimal_late=an.t_last - lo / L.YEAR_DAYS,
                ess_pooled=s["ess_pooled"], ess_sum_of_chains=s["ess_sum_of_chains"], ess_min_chain=s["ess_min_chain"],
                rhat_split=s["rhat_split"], rhat_rank=s["rhat_rank"], n_chains=len(rh), n_samples=int(len(allrh)))


def precision_of(an):
    c = an.kept[0] if an.kept else an.chains[0]
    lo, hi = c.info.get("skygrid_precision_min"), c.info.get("skygrid_precision_max")
    return (lo if lo == hi else np.nan), (lo == hi)


# ----------------------------------------------------------------------------------------------------
def skygrid_tables(an, mcc_counts=None, tree_counts=None, n_report=None):
    """returns dict of DataFrames: rt (pairs), intervals, windows (prior vs posterior), r0defs, corr"""
    g, gen = an.grid, an.gen
    K = g.K
    tau, fixed = precision_of(an)
    tm = tmrca_summary(an)
    t_lo = tm["tmrca_decimal_early"]           # lower (earliest) HPD bound of the tMRCA
    lp = an.lnNe()
    n_report = n_report or (K - 1)
    # ---------------- intervals
    tips_dec = np.array([L.decimal_year(d) for d in an.tip_dates])
    irows = []
    for i in range(1, K + 1):
        older, newer = g.edges(i)
        s = L.summarise([x[:, i - 1] for x in lp], f"lnNe_{i}")
        Ne = np.exp(np.concatenate([x[:, i - 1] for x in lp]))
        nlo, nhi = L.hpd(np.log(Ne))
        n_g = int(((tips_dec > older) & (tips_dec <= newer + 1e-12)).sum()) if i > 1 else int((tips_dec > older).sum())
        if newer <= t_lo:
            dom = "prior-dominated (interval older than the lower 95 % HPD bound of the tMRCA)"
        elif older < t_lo:
            dom = "partly older than the lower 95 % HPD bound of the tMRCA"
        else:
            dom = ""
        row = dict(interval=i, older_edge=str(L.dec2date(older)) if i < K else "open", newer_edge=str(L.dec2date(newer)),
                   mid_point=str(L.dec2date(g.mids[i - 1])) if i < K else "", width_days=g.D * L.YEAR_DAYS if i < K else np.nan,
                   lnNe_median=s["median"], lnNe_hpd_lo=s["hpd_lo"], lnNe_hpd_hi=s["hpd_hi"], lnNe_sd=s["sd"],
                   Ne_tau_years_median=float(np.median(Ne)), Ne_tau_years_hpd_lo=float(np.exp(nlo)), Ne_tau_years_hpd_hi=float(np.exp(nhi)),
                   Ne_tau_days_median=float(np.median(Ne)) * L.YEAR_DAYS,
                   n_genomes_collected=n_g, status=dom, ess_pooled=s["ess_pooled"], ess_sum_of_chains=s["ess_sum_of_chains"],
                   ess_min_chain=s["ess_min_chain"], rhat_split=s["rhat_split"], rhat_rank=s["rhat_rank"])
        if mcc_counts is not None:
            row["n_coalescent_events_MCC_tree"] = int(mcc_counts[i - 1])
        if tree_counts is not None:
            row["n_coalescent_events_posterior_mean"] = float(tree_counts["mean"][i - 1])
            row["n_coalescent_events_posterior_q025"] = float(tree_counts["q025"][i - 1])
            row["n_coalescent_events_posterior_q975"] = float(tree_counts["q975"][i - 1])
        irows.append(row)
    intervals = pd.DataFrame(irows)
    # ---------------- pairs
    prows, rser = [], {}
    for i in range(1, n_report + 1):
        c = L.pair_contrast(g, i)
        rs = [x @ c for x in lp]
        rser[i] = rs
        older_i = i + 1
        _, newer_edge_old = g.edges(older_i)
        older_edge_old, _ = g.edges(older_i)
        if newer_edge_old <= t_lo:
            dom = "prior-dominated (older interval lies before the lower 95 % HPD bound of the tMRCA)"
        elif older_edge_old < t_lo:
            dom = "older interval partly before the lower 95 % HPD bound of the tMRCA"
        else:
            dom = ""
        row = r_summary(rs, gen, f"pair {i}", dict(interval_pair=i, newer_interval=i, older_interval=i + 1,
                                                  time_boundary=str(L.dec2date(g.t_anchor - i * g.D)),
                                                  newer_interval_mid=str(L.dec2date(g.mids[i - 1])),
                                                  older_interval_mid=str(L.dec2date(g.mids[i])) if i + 1 < K else "",
                                                  status=dom))
        if fixed:
            pr = prior_r_summary(c, tau, gen, K, 100000)
            row.update(pr)
            row["skygrid_precision"] = tau
            row["ratio_posterior_sd_to_prior_sd_growth"] = row["growth_per_day_sd"] / pr["prior_sd_growth_per_day"]
            row["data_constraint"] = flag_ratio(row["ratio_posterior_sd_to_prior_sd_growth"])
        prows.append(row)
    rt = pd.DataFrame(prows)
    # correlation of growth rates of adjacent pairs
    crow = []
    for i in range(1, n_report):
        a = np.concatenate(rser[i]); b = np.concatenate(rser[i + 1])
        crow.append(dict(pair_newer=i, pair_older=i + 1, boundary_newer=str(L.dec2date(g.t_anchor - i * g.D)),
                         boundary_older=str(L.dec2date(g.t_anchor - (i + 1) * g.D)),
                         posterior_correlation_growth_rates=float(np.corrcoef(a, b)[0, 1]),
                         prior_correlation_growth_rates=0.0,
                         R_median_newer=float(rt.loc[i - 1, "R_median"]), R_median_older=float(rt.loc[i, "R_median"])))
    corr = pd.DataFrame(crow)
    if len(corr):
        rt["corr_growth_with_next_older_pair"] = list(corr["posterior_correlation_growth_rates"]) + [np.nan] * (len(rt) - len(corr))
    # ---------------- windows: prior vs posterior
    wrows = []
    variants = [("interp", "midpoint", "primary reading: ln Ne interpolated between interval mid-points"),
                ("stair", "midpoint", "sensitivity: level of the interval that contains the date"),
                ("interp", "constant_to_last_tip", "sensitivity (W3 only): level held constant to the last tip"),
                ("interp", "extrapolate_to_last_tip", "sensitivity (W3 only): linear extrapolation to the last tip")]
    wser = {}
    for mode, w3v, vlab in variants:
        W = an.windows(mode, w3v)
        for k, w in W.items():
            if w3v != "midpoint" and k != "W3":
                continue
            if not w["evaluable"]:
                wrows.append(dict(window=k, reading=vlab, definition=w["label"], window_days=w["days"], note=w["note"]))
                continue
            rs = [x @ w["c"] for x in lp]
            key = f"{k}|{mode}|{w3v}"
            wser[key] = rs
            row = r_summary(rs, gen, k, dict(window=k, reading=vlab, definition=w["label"], window_days=w["days"],
                                             start_date=str(L.dec2date(w["t0"])), end_date=str(L.dec2date(w["t1"])),
                                             note=w["note"], is_primary_definition=bool(mode == "interp" and w3v == "midpoint")))
            if fixed:
                pr = prior_r_summary(w["c"], tau, gen, K, 200000)
                row.update(pr)
                row["skygrid_precision"] = tau
                row["ratio_posterior_sd_to_prior_sd_growth"] = row["growth_per_day_sd"] / pr["prior_sd_growth_per_day"]
                row["data_constraint"] = flag_ratio(row["ratio_posterior_sd_to_prior_sd_growth"])
            wrows.append(row)
    windows = pd.DataFrame(wrows)
    # ---------------- R0 definitions
    t0, t1 = L.decimal_year(R0_START), L.decimal_year(R0_END)
    i_old, i_new = g.interval_of(t0), g.interval_of(t1)
    drows = []

    def add(defn, label, cvec, days, extra=None):
        rs = [x @ cvec for x in lp]
        e = dict(definition=defn, description=label, window_days=days)
        if extra:
            e.update(extra)
        row = r_summary(rs, gen, defn, e)
        if fixed:
            psd = L.prior_sd_contrast(cvec, tau) / L.YEAR_DAYS
            row["prior_sd_growth_per_day"] = psd
            row["ratio_posterior_sd_to_prior_sd_growth"] = row["growth_per_day_sd"] / psd
        drows.append(row)

    def anchors(io, inw):
        c = np.zeros(K)
        dist = g.mids[inw - 1] - g.mids[io - 1]
        c[inw - 1], c[io - 1] = 1.0 / dist, -1.0 / dist
        return c, dist * L.YEAR_DAYS

    add("W1 (primary)", f"growth of ln Ne between {R0_START} and {R0_END}, ln Ne interpolated between interval mid-points",
        L.window_contrast(g, t0, t1, "interp"), (t1 - t0) * L.YEAR_DAYS, dict(is_primary_definition=True))
    add("W1 staircase reading", f"level of the interval containing {R0_END} minus level of the interval containing {R0_START}, "
        "divided by the elapsed time between the two dates", L.window_contrast(g, t0, t1, "stair"), (t1 - t0) * L.YEAR_DAYS,
        dict(interval_old=i_old, interval_new=i_new))
    c, d = anchors(i_old, i_new)
    add("(i) intervals containing the dates", f"intervals {i_old} and {i_new}, divided by the distance between their mid-points "
        f"({L.dec2date(g.mids[i_old - 1])} to {L.dec2date(g.mids[i_new - 1])}; definition used so far)", c, d,
        dict(interval_old=i_old, interval_new=i_new))
    for sh, lab in ((1, "one interval older"), (-1, "one interval newer")):
        io, inw = i_old + sh, i_new + sh
        if 1 <= inw < io <= K - 1:
            c, d = anchors(io, inw)
            add(f"(ii) both anchor intervals shifted, {lab}", f"intervals {io} and {inw} "
                f"({L.dec2date(g.mids[io - 1])} to {L.dec2date(g.mids[inw - 1])})", c, d, dict(interval_old=io, interval_new=inw))
    for sh, lab in ((1, "older anchor one interval older"), (-1, "older anchor one interval newer")):
        io = i_old + sh
        if i_new < io <= K - 1:
            c, d = anchors(io, i_new)
            add(f"(ii-b) {lab}", f"intervals {io} and {i_new} ({L.dec2date(g.mids[io - 1])} to {L.dec2date(g.mids[i_new - 1])})",
                c, d, dict(interval_old=io, interval_new=i_new))
    inside = [i for i in range(1, K) if t0 - 1e-9 <= g.mids[i - 1] <= t1 + 1e-9]
    if len(inside) >= 2:
        tt = g.mids[[i - 1 for i in inside]]
        c = np.zeros(K)
        c[[i - 1 for i in inside]] = (tt - tt.mean()) / ((tt - tt.mean()) ** 2).sum()
        add("(iii) least-squares slope", f"slope of a least-squares line through ln Ne at the mid-points of intervals "
            f"{min(inside)}-{max(inside)} ({L.dec2date(tt.min())} to {L.dec2date(tt.max())}), per posterior sample", c,
            (tt.max() - tt.min()) * L.YEAR_DAYS, dict(interval_old=max(inside), interval_new=min(inside)))
    r0defs = pd.DataFrame(drows)
    return dict(rt=rt, intervals=intervals, windows=windows, r0defs=r0defs, corr=corr, tmrca=tm, wser=wser, rser=rser,
                precision=tau, precision_fixed=fixed)


def headline(an, tabs=None):
    """one row per headline quantity"""
    rows = []
    base = dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, n_genomes=an.n_tips,
                last_tip=str(an.last_tip.date()), chains_run=len(an.chains), chains_retained=len(an.kept))

    def put(q, med, lo, hi, s, unit, extra=None):
        r = dict(base)
        r.update(quantity=q, median=med, hpd_lo=lo, hpd_hi=hi, unit=unit, ess_pooled=s.get("ess_pooled"),
                 ess_sum_of_chains=s.get("ess_sum_of_chains"), ess_min_chain=s.get("ess_min_chain"),
                 rhat_split=s.get("rhat_split"), rhat_rank=s.get("rhat_rank"), n_chains=s.get("n_chains"),
                 n_samples=s.get("n_samples"))
        if extra:
            r.update(extra)
        rows.append(r)

    if not an.kept:
        return pd.DataFrame(rows)
    s = L.summarise(an.series("clock.rate"), "clock rate", 1e3)
    put("clock rate", s["median"], s["hpd_lo"], s["hpd_hi"], s, "1e-3 substitutions/site/year")
    tm = tmrca_summary(an)
    put("tMRCA", tm["tmrca_median"], tm["tmrca_hpd_early"], tm["tmrca_hpd_late"], tm, "date (day resolution, from root height)",
        dict(root_height_days_median=tm["root_height_days_median"]))
    if an.model.startswith("exp"):
        g = an.series("exponential.growthRate")
        rs = r_summary(g, an.gen, "exponential growth")
        put("growth rate (whole period, exponential model)", rs["growth_per_day_median"], rs["growth_per_day_hpd_lo"],
            rs["growth_per_day_hpd_hi"], rs, "per day")
        put("doubling time (whole period, exponential model)", rs["doubling_time_d"], rs["doubling_time_d_lo"],
            rs["doubling_time_d_hi"], rs, "days")
        put("R (whole period, exponential model)", rs["R_median"], rs["R_hpd_lo"], rs["R_hpd_hi"], rs, "",
            dict(P_R_gt_1=rs["P_R_gt_1"]))
    elif tabs is not None:
        W = tabs["windows"]
        Wp = W[W.get("is_primary_definition") == True]
        for _, w in Wp.iterrows():
            nm = {"W1": "R0, early phase (W1)", "W2": "R, W2", "W3": "R, W3", "W3b": "R, W3b"}[w["window"]]
            put(f"{nm}: {w['definition']}", w["R_median"], w["R_hpd_lo"], w["R_hpd_hi"], w.to_dict(), "",
                dict(P_R_gt_1=w["P_R_gt_1"], growth_per_day_median=w["growth_per_day_median"], doubling_time_d=w["doubling_time_d"],
                     ratio_posterior_sd_to_prior_sd_growth=w.get("ratio_posterior_sd_to_prior_sd_growth"),
                     data_constraint=w.get("data_constraint")))
        rt = tabs["rt"]
        pk = rt.loc[rt["R_median"].idxmax()]
        put(f"peak R(t): pair {int(pk['interval_pair'])}, boundary {pk['time_boundary']}", pk["R_median"], pk["R_hpd_lo"],
            pk["R_hpd_hi"], pk.to_dict(), "", dict(P_R_gt_1=pk["P_R_gt_1"],
                                                   ratio_posterior_sd_to_prior_sd_growth=pk.get("ratio_posterior_sd_to_prior_sd_growth"),
                                                   data_constraint=pk.get("data_constraint"), status=pk.get("status")))
        ok = rt[rt["status"] == ""]
        if len(ok):
            pk = ok.loc[ok["R_median"].idxmax()]
            put(f"peak R(t) among pairs not flagged as prior-dominated: pair {int(pk['interval_pair'])}, boundary {pk['time_boundary']}",
                pk["R_median"], pk["R_hpd_lo"], pk["R_hpd_hi"], pk.to_dict(), "",
                dict(P_R_gt_1=pk["P_R_gt_1"], ratio_posterior_sd_to_prior_sd_growth=pk.get("ratio_posterior_sd_to_prior_sd_growth"),
                     data_constraint=pk.get("data_constraint")))
        la = rt.iloc[0]
        put(f"latest-interval R: pair 1, boundary {la['time_boundary']}", la["R_median"], la["R_hpd_lo"], la["R_hpd_hi"],
            la.to_dict(), "", dict(P_R_gt_1=la["P_R_gt_1"],
                                   ratio_posterior_sd_to_prior_sd_growth=la.get("ratio_posterior_sd_to_prior_sd_growth"),
                                   data_constraint=la.get("data_constraint")))
        if not tabs["precision_fixed"]:
            s = L.summarise(an.series("skygrid.precision"), "Skygrid precision")
            put("Skygrid precision (estimated)", s["median"], s["hpd_lo"], s["hpd_hi"], s, "")
    return pd.DataFrame(rows)


def stays_below(rt, col, thr, below=True):
    """boundary date of the oldest pair from which (towards the present) the column stays below thr"""
    v = rt.sort_values("interval_pair")[col].values          # pair 1 = most recent
    ok = v < thr
    first = None
    for i in range(len(v)):                                    # i = 0 is the most recent pair
        if ok[:i + 1].all():
            first = i
        else:
            break
    if first is None:
        return None
    return rt.sort_values("interval_pair").iloc[first]


def per_chain_diagnostics(an, tabs=None, params=None):
    """per-chain ESS (pipeline estimator) and pooled ESS / split-R-hat for every reported parameter"""
    rows = []
    kept = an.kept
    ser = {}
    ser["clock.rate"] = an.series("clock.rate", an.chains)
    ser["root height (tMRCA)"] = an.root_height(an.chains)
    if an.chains[0].is_sky:
        lps = an.lnNe(an.chains)
        K = lps[0].shape[1]
        for i in range(1, K + 1):
            ser[f"lnNe interval {i}"] = [x[:, i - 1] for x in lps]
        if an.grid is not None:
            for i in range(1, K):
                c = L.pair_contrast(an.grid, i)
                ser[f"growth rate pair {i}"] = [x @ c for x in lps]
            for mode, w3v in (("interp", "midpoint"),):
                for k, w in an.windows(mode, w3v).items():
                    if w["evaluable"]:
                        ser[f"growth rate {k}"] = [x @ w["c"] for x in lps]
        if "skygrid.precision" in an.chains[0].post.columns and an.chains[0].info["skygrid_precision_min"] != an.chains[0].info["skygrid_precision_max"]:
            ser["skygrid.precision"] = an.series("skygrid.precision", an.chains)
    else:
        ser["exponential.growthRate"] = an.series("exponential.growthRate", an.chains)
        if "exponential.popSize" in an.chains[0].post.columns:
            ser["ln exponential.popSize"] = [np.log(x) for x in an.series("exponential.popSize", an.chains)]
    # sampler-level quantities: log-posterior, tree likelihood, number of mutations on the tree (Delphy), tree length
    for p in ("kappa", "treeLength", "posterior_for_Delphy", "prior_for_Delphy", "likelihood_really_logG", "numMuts", "skygrid",
              "joint", "prior", "likelihood", "coalescent"):
        if p in an.chains[0].post.columns:
            v = an.series(p, an.chains)
            if all(np.nanstd(x) > 0 for x in v):
                ser[p] = v
    keep_idx = [j for j, c in enumerate(an.chains) if c.retained]
    for name, lst in ser.items():
        for j, c in enumerate(an.chains):
            x = lst[j]
            rows.append(dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, parameter=name,
                             scope="chain", chain=c.label, seed=c.seed, retained=c.retained, n_samples=len(x),
                             ess=L.ess_single(x), median=float(np.median(x)), rhat_split=np.nan, rhat_rank=np.nan))
        if keep_idx:
            k = [lst[j] for j in keep_idx]
            s = L.summarise(k, name)
            ok = (s["ess_pooled"] >= 200) and (s["rhat_split"] <= 1.05)
            rows.append(dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, parameter=name,
                             scope="pooled retained chains", chain=f"{len(k)} chains", seed="", retained=True,
                             n_samples=s["n_samples"], ess=s["ess_pooled"], ess_sum_of_chains=s["ess_sum_of_chains"],
                             median=s["median"], rhat_split=s["rhat_split"], rhat_rank=s["rhat_rank"],
                             criterion_ess200_rhat105_met=bool(ok)))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------
# degenerate state: description of the two states from the chain logs, progress lines and sampled trees
# ----------------------------------------------------------------------------------------------------
def degenerate_rows(an, progress_paths, trees_paths, barrier_loc_years=1.0 / 365.0, barrier_scale=0.30, n_trees=20,
                    transient_steps=100e6, slide_margin_steps=30e6, settle_steps=100e6, spec_label=""):
    """one row per chain and part ('regular' / 'degenerate'); progress_paths, trees_paths: dict seed -> path"""
    g = an.grid
    tau, fixed = precision_of(an)
    labels_fa = [l[1:].strip() for l in L.open_any(an.fasta) if l.startswith(">")]
    dec = {l: L.decimal_year(l.split("|")[-1]) for l in labels_fa}
    rows = []
    for c in an.chains:
        full = c.full.copy()
        prog = L.parse_delphy_stdout(progress_paths[c.seed])
        a = full.merge(prog[["state", "log_posterior", "log_G", "log_coal", "log_other_priors"]], on="state", how="inner")
        idx = full["state"].isin(a["state"]).values
        lp = a[c.lp_cols].values.astype(float)
        a["gmrf"] = L.gmrf_term(lp, tau) if fixed else np.nan
        a["barrier"] = L.barrier_term(lp, barrier_loc_years, barrier_scale)
        a["other"] = a["log_other_priors"] - a["gmrf"] - a["barrier"]
        a["root_date"] = c.root_dates_full[idx]
        a["coll"] = c.collapsed_mask[idx]
        a["minlp"] = lp.min(axis=1)
        a["argminlp"] = lp.argmin(axis=1) + 1
        last = float(a["state"].iloc[-1])
        parts = {}
        if c.info["collapsed"]:
            fc = float(c.info["first_collapsed_state"])
            parts["regular (before entry)"] = a[(a["state"] >= transient_steps) & (a["state"] <= fc - slide_margin_steps)]
            deg = a[a["state"] >= fc + settle_steps]
            parts["degenerate"] = deg if len(deg) >= 10 else a[a["state"] >= fc]
        else:
            parts["regular (after 30 % burn-in)"] = a[a["state"] >= c.full["state"].iloc[c.nb]]
        trees_all = None
        for pname, d in parts.items():
            row = dict(analysis=an.name, specification=spec_label or an.spec, engine=an.engine, dataset=an.dataset, chain=c.label,
                       seed=c.seed, steps=int(last), chain_entered_degenerate_state=bool(c.info["collapsed"]),
                       first_step_in_degenerate_state=c.info["first_collapsed_state"], part=pname, n_logged_states=len(d))
            if len(d) == 0:
                rows.append(row)
                continue
            rd = pd.to_datetime(d["root_date"].values)
            med_rd = pd.Timestamp(np.sort(rd.values)[len(rd) // 2])
            mode_int = int(pd.Series(d["argminlp"]).mode().iloc[0])
            row.update(first_step=int(d["state"].iloc[0]), last_step=int(d["state"].iloc[-1]),
                       root_date_median=str(med_rd.date()), root_date_earliest=str(rd.min().date()), root_date_latest=str(rd.max().date()),
                       clock_rate_median_x1e3=float(d["clock.rate"].median() * 1e3),
                       log_posterior_median=float(d["log_posterior"].median()),
                       tree_log_likelihood_median=float(d["log_G"].median()),
                       coalescent_log_density_delphy_median=float(d["log_coal"].median()),
                       smoothing_term_median=float(d["gmrf"].median()), low_population_barrier_term_median=float(d["barrier"].median()),
                       other_priors_median=float(d["other"].median()),
                       lowest_lnNe_median=float(d["minlp"].median()), lowest_lnNe_min=float(d["minlp"].min()),
                       interval_with_lowest_lnNe_most_frequent=mode_int,
                       lowest_Ne_tau_days_median=float(np.exp(d["minlp"].median()) * L.YEAR_DAYS),
                       n_mutations_median=float(d["numMuts"].median()) if "numMuts" in d else np.nan,
                       tree_length_years_median=float(d["treeLength"].median()) if "treeLength" in d else np.nan)
            tp = trees_paths.get(c.seed)
            if tp and os.path.exists(tp) and n_trees > 0:
                if trees_all is None:
                    labels, trees_all = L.read_nexus_trees(tp)
                    tree_by_state = dict(trees_all)
                st = [s for s in d["state"].values if s in tree_by_state]
                if len(st) > n_trees:
                    st = [st[i] for i in np.linspace(0, len(st) - 1, n_trees).astype(int)]
                ex, de, n9 = [], [], []
                dd = d.set_index("state")
                for s_ in st:
                    tips, ints = L.node_times(tree_by_state[s_], labels, dec)
                    n_i, A_i = L.skygrid_tree_terms(tips, ints, g)
                    lnN = dd.loc[s_, c.lp_cols].values.astype(float)
                    ex.append(L.exact_skygrid_logdensity(n_i, A_i, lnN))
                    de.append(float(dd.loc[s_, "log_coal"]))
                    n9.append(n_i[mode_int - 1])
                if ex:
                    err = np.array(de) - np.array(ex)
                    row.update(n_trees_checked=len(ex), coalescent_log_density_exact_median=float(np.median(ex)),
                               delphy_minus_exact_median=float(np.median(err)), delphy_minus_exact_min=float(err.min()),
                               delphy_minus_exact_max=float(err.max()),
                               coalescent_events_in_lowest_interval_median=float(np.median(n9)))
            rows.append(row)
    return pd.DataFrame(rows)


def collapse_counts(an, spec_label=""):
    ch = an.chains
    fc = [c.info["first_collapsed_state"] for c in ch if c.info["collapsed"]]
    steps = [c.info["last_state"] for c in ch]
    at_risk = sum(min(c.info["first_collapsed_state"], c.info["last_state"]) if c.info["collapsed"] else c.info["last_state"] for c in ch)
    return dict(analysis=an.name, specification=spec_label or an.spec, engine=an.engine, model=an.model, dataset=an.dataset,
                n_genomes=an.n_tips, chains_run=len(ch), chains_entered_degenerate_state=len(fc), chains_retained=len(an.kept),
                steps_per_chain=int(np.median(steps)), seeds=" ".join(str(c.seed) for c in ch),
                seeds_entered=" ".join(str(c.seed) for c in ch if c.info["collapsed"]),
                first_entry_steps_millions=" ".join(f"{x / 1e6:.1f}" for x in sorted(fc)),
                steps_at_risk_millions=at_risk / 1e6,
                entries_per_1000M_steps_at_risk=(len(fc) / (at_risk / 1e9)) if at_risk > 0 else np.nan,
                chains_returning_to_regular_state=sum(1 for c in ch if c.info["collapsed"] and c.info["n_regular_states_after_first_collapse"] > 50))


# ----------------------------------------------------------------------------------------------------
# second criterion for the degenerate state (added 24 Sep): Delphy's coalescent log-density against the exact
# log-density of the sampled tree under the logged population sizes
# ----------------------------------------------------------------------------------------------------
ERR_THRESHOLD = 100.0     # log units


def approximation_error(an, chain, progress_path, trees_path, max_trees=60):
    labels_fa = [l[1:].strip() for l in L.open_any(an.fasta) if l.startswith(">")]
    dec = {l: L.decimal_year(l.split("|")[-1]) for l in labels_fa}
    labels, trees = L.read_nexus_trees(trees_path)
    if len(trees) > max_trees:
        idx = np.unique(np.linspace(0, len(trees) - 1, max_trees).astype(int))
        trees = [trees[i] for i in idx]
    prog = L.parse_delphy_stdout(progress_path).set_index("state")
    full = chain.full.set_index("state")
    rows = []
    for st, nw in trees:
        if st not in prog.index or st not in full.index:
            continue
        tips, ints = L.node_times(nw, labels, dec)
        n_i, A_i = L.skygrid_tree_terms(tips, ints, an.grid)
        lnN = full.loc[st, chain.lp_cols].values.astype(float)
        ex = L.exact_skygrid_logdensity(n_i, A_i, lnN)
        rows.append(dict(state=int(st), coal_delphy=float(prog.loc[st, "log_coal"]), coal_exact=ex,
                         err=float(prog.loc[st, "log_coal"]) - ex, root_tree=float(ints.min()),
                         n_events_max_interval=int(n_i.max()), interval_max_events=int(n_i.argmax() + 1)))
    return pd.DataFrame(rows)


def apply_second_criterion(an, threshold=ERR_THRESHOLD, max_trees=60):
    """flags chains whose sampled trees show Delphy's coalescent log-density above the exact value by more than
    'threshold'; updates chain.retained / chain.why; returns a table (one row per chain)"""
    out = []
    for c in an.chains:
        E = approximation_error(an, c, an.progress[c.seed], an.trees[c.seed], max_trees)
        c.err = E
        bad = E[E["err"] > threshold]
        post = E[E["state"] >= c.full["state"].iloc[c.nb]]
        c.info.update(n_trees_checked=len(E), approx_error_max=float(E["err"].max()) if len(E) else np.nan,
                      approx_error_median_post_burnin=float(post["err"].median()) if len(post) else np.nan,
                      approx_error_max_post_burnin=float(post["err"].max()) if len(post) else np.nan,
                      first_tree_state_error_above_threshold=(int(bad["state"].iloc[0]) if len(bad) else np.nan),
                      flagged_by_density_criterion=bool(len(bad) > 0))
        if len(bad) and c.retained:
            c.retained = False
            c.why = (f"Delphy's coalescent log-density exceeds the exact value by more than {threshold:g} log units "
                     f"(first tree at step {int(bad['state'].iloc[0])}); root-date rule not met")
        elif len(bad):
            c.why += f"; density criterion also met (first tree at step {int(bad['state'].iloc[0])})"
        out.append(dict(analysis=an.name, chain=c.label, seed=c.seed, root_date_rule=c.info["collapsed"],
                        density_criterion=c.info["flagged_by_density_criterion"], retained=c.retained,
                        first_step_root_date_rule=c.info["first_collapsed_state"],
                        first_tree_step_density_criterion=c.info["first_tree_state_error_above_threshold"],
                        approx_error_median_post_burnin=c.info["approx_error_median_post_burnin"],
                        approx_error_max=c.info["approx_error_max"], why=c.why))
    return pd.DataFrame(out)


# ----------------------------------------------------------------------------------------------------
# pooled MCC tree (delphy_mcc removes the first 30 % of the trees of its input file by count)
# ----------------------------------------------------------------------------------------------------
def pooled_trees_file(an, out_path, burnin=0.30, pad_tree_from=None):
    """writes one NEXUS file: [padding trees][post-burn-in trees of all retained chains]; the padding (copies of
    the first post-burn-in tree) is what delphy_mcc removes as its built-in 30 % burn-in.
    Returns (n_pooled, n_padding)."""
    header = None
    pooled = []
    for c in an.kept:
        tp = an.trees[c.seed]
        hdr, lines = [], []
        with L.open_any(tp) as fh:
            for line in fh:
                if line.startswith("tree "):
                    lines.append(line.rstrip("\n"))
                elif not lines:
                    hdr.append(line.rstrip("\n"))
        if header is None:
            header = [h for h in hdr if not h.startswith("#") or h.startswith("#NEXUS")]
        n = len(lines)
        keep = lines[int(n * burnin):]
        for k, ln in enumerate(keep):
            pooled.append(ln)
    N = len(pooled)
    B = int(math.floor(3 * N / 7.0))
    while int(math.floor(0.3 * (B + N))) < B:
        B -= 1
    while int(math.floor(0.3 * (B + N))) > B:
        B += 1
    assert int(math.floor(0.3 * (B + N))) == B
    import re as _re
    with open(out_path, "w") as fh:
        fh.write("\n".join(header) + "\n")
        k = 0
        for ln in [pooled[0]] * B + pooled:
            fh.write(_re.sub(r"^tree \S+", f"tree STATE_{k}", ln) + "\n")
            k += 1
        fh.write("End;\n")
    return N, B


def coalescent_counts(an, mcc_path=None, max_trees_per_chain=40):
    """coalescent events per interval: MCC tree and posterior summary over sampled post-burn-in trees"""
    labels_fa = [l[1:].strip() for l in L.open_any(an.fasta) if l.startswith(">")]
    dec = {l: L.decimal_year(l.split("|")[-1]) for l in labels_fa}
    out = {}
    if mcc_path:
        labels, trees = L.read_nexus_trees(mcc_path)
        tips, ints = L.node_times(trees[0][1], labels, dec, tol=5e-3)
        n_i, A_i = L.skygrid_tree_terms(tips, ints, an.grid)
        out["mcc"] = n_i
        out["mcc_root"] = float(ints.min())
    cnt = []
    for c in an.kept:
        labels, trees = L.read_nexus_trees(an.trees[c.seed])
        trees = trees[int(len(trees) * 0.30):]
        if len(trees) > max_trees_per_chain:
            idx = np.unique(np.linspace(0, len(trees) - 1, max_trees_per_chain).astype(int))
            trees = [trees[i] for i in idx]
        for st, nw in trees:
            tips, ints = L.node_times(nw, labels, dec)
            n_i, _ = L.skygrid_tree_terms(tips, ints, an.grid)
            cnt.append(n_i)
    cnt = np.array(cnt)
    out["posterior"] = dict(mean=cnt.mean(axis=0), q025=np.quantile(cnt, 0.025, axis=0), q975=np.quantile(cnt, 0.975, axis=0),
                            n_trees=len(cnt))
    return out


def exponential_table(an, label=None):
    g = an.series("exponential.growthRate")
    rs = r_summary(g, an.gen, "exponential growth")
    tm = tmrca_summary(an)
    s = L.summarise(an.series("clock.rate"), "clock", 1e3)
    row = dict(analysis=label or an.name, engine=an.engine, dataset=an.dataset, n_genomes=an.n_tips, first_tip=str(an.first_tip.date()),
               last_tip=str(an.last_tip.date()), chains_run=len(an.chains), chains_retained=len(an.kept),
               growth_rate_per_year_median=rs["growth_per_day_median"] * L.YEAR_DAYS, growth_rate_per_day_median=rs["growth_per_day_median"],
               growth_rate_per_day_hpd_lo=rs["growth_per_day_hpd_lo"], growth_rate_per_day_hpd_hi=rs["growth_per_day_hpd_hi"],
               doubling_time_d_median=rs["doubling_time_d"], doubling_time_d_hpd_lo=rs["doubling_time_d_lo"],
               doubling_time_d_hpd_hi=rs["doubling_time_d_hi"], R_median=rs["R_median"], R_hpd_lo=rs["R_hpd_lo"], R_hpd_hi=rs["R_hpd_hi"],
               P_R_gt_1=rs["P_R_gt_1"], growth_ess_pooled=rs["ess_pooled"], growth_rhat_split=rs["rhat_split"],
               tmrca_median=tm["tmrca_median"], tmrca_hpd_early=tm["tmrca_hpd_early"], tmrca_hpd_late=tm["tmrca_hpd_late"],
               root_height_ess_pooled=tm["ess_pooled"], root_height_rhat_split=tm["rhat_split"],
               clock_rate_x1e3_median=s["median"], clock_rate_x1e3_hpd_lo=s["hpd_lo"], clock_rate_x1e3_hpd_hi=s["hpd_hi"],
               clock_ess_pooled=s["ess_pooled"], clock_rhat_split=s["rhat_split"],
               generation_time_mean_d=an.gen.mean_d, generation_time_sd_d=an.gen.sd_d)
    if "exponential.popSize" in an.kept[0].post.columns:
        p = L.summarise([np.log(x) for x in an.series("exponential.popSize")], "ln N0")
        row.update(ln_final_popsize_median=p["median"], ln_final_popsize_ess_pooled=p["ess_pooled"], ln_final_popsize_rhat_split=p["rhat_split"])
    return row
