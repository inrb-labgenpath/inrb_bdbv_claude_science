#!/usr/bin/env python
"""beast_analysis.py - BEAST X cross-check of the Delphy analyses (24 Sep 2026 round).

usage (from the workspace root, conda environment with BEAST X for the tree step):
    python work/pipeline_ext/beast_analysis.py --in <directory with one sub-directory per chain> --out work/deliverables [--test]

Models (XML files in work/beast/xml, written by beast_make_xml.py):
  A  Skygrid, 20 parameters, cutoff 0.8 yr, precision fixed at 4.06335 (as in the Delphy analyses), HKY, strict clock
  B  as A, precision estimated (prior Gamma(shape 0.001, scale 1000))
  C  'INRB-like': Skygrid with 32 parameters, last transition point at 224 d (cutoff 0.613699 yr), precision estimated, GTR without rate
     heterogeneity, strict clock; primary genome set and INRB-like retained set.  Settings reconstructed from the wording of a public post;
     agreement is indicative, not a replication.
  D  exponential growth, HKY, strict clock
Burn-in 30 % (12 M of 40 M states), fixed when the XML files were written.  Chains that did not reach the requested length are not used.
--test: accepts incomplete chains and writes to the given directory with the prefix TEST_ (code test only).
"""
import argparse, glob, gzip, json, math, os, re, shutil, subprocess, sys, time
import numpy as np
import pandas as pd
sys.path.insert(0, "work/pipeline_ext")
import rt_lib as L, analysis as A

SFX = "20260924"
CHAIN_LENGTH = 40_000_000
BURNIN = 0.30
LABEL_DEFAULT = ("conditional on the chain not having entered the degenerate state within 1,500 M steps; "
                 "retained chains are a selected minority of the chains run")
MODELS = {
    "A_skyfix_primary": dict(xml="beastA_skyfix_primary_20260924", short="A", dataset="primary", kind="skygrid",
                             label="A: Skygrid (20 parameters, cutoff 0.8 yr), precision fixed at 4.06335, HKY, strict clock",
                             delphy="Skygrid | primary | 8000 cells", delphy_alt="Skygrid | primary | 1000 cells", like_for_like="same model and grid as the Delphy primary analysis"),
    "B_skyest_primary": dict(xml="beastB_skyest_primary_20260924", short="B", dataset="primary", kind="skygrid",
                             label="B: Skygrid (20 parameters, cutoff 0.8 yr), precision estimated (Gamma(0.001, 1000) prior), HKY, strict clock",
                             delphy="Skygrid | primary | 8000 cells", delphy_alt="Skygrid | primary | 1000 cells",
                             like_for_like="same grid as the Delphy primary analysis; precision estimated in BEAST X, fixed at 4.06335 in Delphy"),
    "C_inrblike_primary": dict(xml="beastC_inrblike_primary_20260924", short="C", dataset="primary", kind="skygrid",
                               label="C (INRB-like): Skygrid (32 parameters, cutoff 0.613699 yr), precision estimated, GTR, strict clock; primary genome set",
                               delphy="Skygrid | gridWeekly32 | 8000 cells", delphy_alt="Skygrid | gridWeekly32 | 1000 cells",
                               like_for_like="same grid as the Delphy weekly-grid analysis; BEAST X: GTR and estimated precision, Delphy: HKY and fixed precision"),
    "C_inrblike_inrbRetained": dict(xml="beastC_inrblike_inrbRetained_20260924", short="C-inrbRetained", dataset="inrbRetained", kind="skygrid",
                                    label="C (INRB-like): Skygrid (32 parameters, cutoff 0.613699 yr), precision estimated, GTR, strict clock; INRB-like retained set",
                                    delphy="Skygrid | inrbRetained | 8000 cells", delphy_alt="Skygrid | inrbRetained | 1000 cells",
                                    like_for_like="same genome set; grids differ (BEAST X 32 parameters of 7.2 d, Delphy 20 parameters of 15.4 d): windows and "
                                                  "scalar parameters are compared, interval pairs are not"),
    "D_exp_primary": dict(xml="beastD_exp_primary_20260924", short="D", dataset="primary", kind="exponential",
                          label="D: exponential growth, HKY, strict clock", delphy="exponential | primary | 8000 cells",
                          delphy_alt="exponential | primary | 1000 cells", like_for_like="same model as the Delphy exponential analysis"),
}
ALN = {"primary": "inputs/aln/primary_20260924.fasta", "inrbRetained": "inputs/aln/inrbRetained_20260924.fasta"}


def xml_params(m):
    return json.load(open(f"work/beast/xml/beast_xml_parameters_{MODELS[m]['xml']}.json"))


# ----------------------------------------------------------------------------------------------------- chains
def load_model(m, indir, test=False):
    info = MODELS[m]
    P = xml_params(m)
    fasta = ALN[info["dataset"]]
    dd = L.tip_dates_from_fasta(fasta)
    chains = []
    for d in sorted(glob.glob(os.path.join(indir, m + "_s*"))):
        seed = int(d.rsplit("_s", 1)[1])
        logs = sorted(glob.glob(os.path.join(d, info["xml"] + ".log*")))
        if not logs:
            continue
        c = A.Chain("beast", logs[0], fasta, seed, burnin=BURNIN, first_tip=dd.min(), last_tip=dd.max(), label=os.path.basename(d))
        c.dir = d
        c.info["steps"] = CHAIN_LENGTH
        c.info["complete"] = bool(c.info["last_state"] >= CHAIN_LENGTH)
        spacing = np.diff(c.full["state"].values)
        assert len(set(spacing)) == 1 and spacing[0] == int(P["log_every"]), (d, set(spacing))
        if not c.info["complete"]:
            if test:
                c.why += f" [TEST: incomplete chain accepted, state {c.info['last_state']}]"
            else:
                c.retained = False
                c.why = f"incomplete (state {c.info['last_state']} of {CHAIN_LENGTH})"
        chains.append(c)
    grid = None
    if info["kind"] == "skygrid":
        grid = L.Grid.from_cutoff(int(P["num_parameters"]), float(P["cutoff"]), dd.max())
    an = A.Analysis("beastx_" + m, "beast", info["label"], info["dataset"], fasta, chains, grid=grid, spec=info["xml"] + ".xml")
    an.key, an.P = m, P
    if test:                                     # equal length is needed for the pooled diagnostics
        n = min(len(c.post) for c in an.kept)
        for c in an.kept:
            c.post = c.post.iloc[:n].reset_index(drop=True)
    return an


def tree_span_minimum(an):
    """per chain: lowest ln Ne among the intervals that overlap the tree, over all logged states after the first 1 M states"""
    rows = []
    g = an.grid
    for c in an.chains:
        f = c.full[c.full["state"] >= 1_000_000]
        lp = f[c.lp_cols].values
        rh = f[c.rh_col].values
        iroot = np.minimum(np.floor(rh / g.D + 1e-12).astype(int) + 1, g.K)          # interval that contains the root
        mask = np.arange(1, g.K + 1)[None, :] <= iroot[:, None]
        v = np.where(mask, lp, np.inf)
        mn = v.min(axis=1)
        imin = v.argmin(axis=1) + 1
        j = int(mn.argmin())
        low = mn < math.log(1 / 365.0)
        rd = L.root_dates_from_height(rh, an.last_tip)
        # the same, leaving out the interval that contains the root (in which two lineages remain)
        v2 = np.where(np.arange(1, g.K + 1)[None, :] < iroot[:, None], lp, np.inf)
        mn2 = v2.min(axis=1)
        rows.append(dict(chain=c.label, seed=c.seed, lowest_lnNe_within_tree_span=float(mn.min()), interval_of_lowest=int(imin[j]),
                         state_of_lowest=int(f["state"].values[j]), root_date_at_state_of_lowest=str(rd[j]),
                         fraction_of_states_with_Ne_tau_below_1_day=float(low.mean()),
                         fraction_of_states_with_Ne_tau_below_7_days=float((mn < math.log(7 / 365.0)).mean()),
                         fraction_of_those_states_in_which_the_lowest_interval_contains_the_root=(float((imin[low] == iroot[low]).mean()) if low.any() else np.nan),
                         latest_root_date_among_those_states=(str(rd[low].max()) if low.any() else ""),
                         median_root_date_among_those_states=(str(np.sort(rd[low])[low.sum() // 2]) if low.any() else ""),
                         median_root_date_other_states=(str(np.sort(rd[~low])[(~low).sum() // 2]) if (~low).any() else ""),
                         lowest_lnNe_excluding_the_interval_with_the_root=float(mn2[np.isfinite(mn2)].min()) if np.isfinite(mn2).any() else np.nan,
                         fraction_of_states_with_Ne_tau_below_1_day_excluding_the_interval_with_the_root=float((mn2 < math.log(1 / 365.0)).mean()),
                         latest_root_date=str(rd.max()), first_tip=str(an.first_tip.date())))
    return pd.DataFrame(rows)


def retention_row_beast(an):
    ch = an.chains
    return dict(analysis=an.name, engine="BEAST X 10.5.0", model=an.model, dataset=an.dataset, n_genomes=an.n_tips, first_tip=str(an.first_tip.date()),
                last_tip=str(an.last_tip.date()), chains_run=len(ch), states_per_chain=CHAIN_LENGTH, burnin_states=int(BURNIN * CHAIN_LENGTH),
                discarded_root_date_rule=sum(1 for c in ch if c.info["collapsed"]),
                discarded_density_criterion_only=0,
                density_criterion_note="not applicable: BEAST X evaluates the coalescent density without a time grid (checked against an independent implementation, "
                                       "density_check_three_way table)",
                discarded_incomplete=sum(1 for c in ch if not c.info["complete"] and not c.info["collapsed"]),
                chains_retained=len(an.kept), seeds_run=" ".join(str(c.seed) for c in ch), seeds_retained=" ".join(str(c.seed) for c in an.kept),
                states_reached_min=min(c.info["last_state"] for c in ch), states_reached_max=max(c.info["last_state"] for c in ch),
                latest_root_date_any_chain=max(c.info["latest_root_date"] for c in ch))


def extra_diagnostics(an):
    """ESS and split-R-hat of the logged parameters of the substitution model that the shared diagnostics function does not cover"""
    rows = []
    cols = [c for c in an.chains[0].post.columns if c.startswith("gtr.rates.") or c.startswith("frequencies") or c in ("kappa", "alpha", "meanRate")]
    have = set(A.per_chain_diagnostics(an)["parameter"])
    keep_idx = [j for j, c in enumerate(an.chains) if c.retained]
    for name in cols:
        if name in have:
            continue
        lst = an.series(name, an.chains)
        if not all(np.nanstd(x) > 0 for x in lst):
            continue
        for j, c in enumerate(an.chains):
            x = lst[j]
            rows.append(dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, parameter=name, scope="chain", chain=c.label, seed=c.seed,
                             retained=c.retained, n_samples=len(x), ess=L.ess_single(x), median=float(np.median(x)), rhat_split=np.nan, rhat_rank=np.nan))
        if keep_idx:
            k = [lst[j] for j in keep_idx]
            s = L.summarise(k, name)
            rows.append(dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, parameter=name, scope="pooled retained chains",
                             chain=f"{len(k)} chains", seed="", retained=True, n_samples=s["n_samples"], ess=s["ess_pooled"], ess_sum_of_chains=s["ess_sum_of_chains"],
                             median=s["median"], rhat_split=s["rhat_split"], rhat_rank=s["rhat_rank"],
                             criterion_ess200_rhat105_met=bool((s["ess_pooled"] >= 200) and (s["rhat_split"] <= 1.05))))
    return pd.DataFrame(rows)


# ----------------------------------------------------------------------------------------------------- quantities from per-chain series
class Series:
    """per-chain post-burn-in series of one analysis (either engine)"""
    def __init__(self, name, engine, lnNe, clock, rh, grid, last_tip, growth=None):
        self.name, self.engine, self.lnNe, self.clock, self.rh, self.grid, self.growth = name, engine, lnNe, clock, rh, grid, growth
        self.last_tip = pd.Timestamp(last_tip)
        self.t_last = L.decimal_year(last_tip)
        self.gen = L.GenTime()


def series_from_analysis(an):
    k = an.kept
    return Series(an.name, "BEAST X", an.lnNe() if an.grid is not None else None, an.series("clock.rate"), an.root_height(), an.grid, an.last_tip,
                  growth=(an.series("exponential.growthRate") if an.grid is None else None))


_PS = {}


def series_from_delphy(label, samples_path):
    if "df" not in _PS:
        _PS["df"] = pd.read_csv(samples_path, low_memory=False)
    d = _PS["df"][_PS["df"]["analysis"] == label]
    assert len(d), f"no Delphy samples for {label}"
    chains = [x for _, x in d.groupby("seed", sort=True)]
    n = min(len(x) for x in chains)
    assert all(len(x) == n for x in chains)
    clock = [x["clock_rate"].values for x in chains]; rh = [x["root_height_years"].values for x in chains]
    # last tip and grid are taken from the alignment and the specification of the analysis, not from the rounded columns of the sample file
    sets = d["set"].unique()
    assert len(sets) == 1, sets
    last_tip = L.tip_dates_from_fasta(f"inputs/aln/{sets[0]}_{SFX}.fasta").max()
    assert abs(L.decimal_year(last_tip) - float((d["root_age_decimal_year"] + d["root_height_years"]).median())) < 2e-4, (label, last_tip)
    if label.startswith("exponential"):
        return Series(label, "Delphy", None, clock, rh, None, last_tip, growth=[x["growth_rate_per_year_exponential"].values for x in chains])
    K = int(d["n_parameters"].iloc[0]); D = float(d["interval_days"].iloc[0]) / L.YEAR_DAYS
    cut = [c for c in (0.8, 0.613699) if abs(D * (K - 1) - c) < 1e-5]
    assert len(cut) == 1, (label, K, D * (K - 1))
    grid = L.Grid.from_cutoff(K, cut[0], last_tip)
    assert abs(grid.t_anchor - float(d["grid_anchor_decimal_year"].iloc[0])) < 1e-4
    lp = [x[[f"lnNe_{i:02d}" for i in range(1, K + 1)]].values.astype(float) for x in chains]
    return Series(label, "Delphy", lp, clock, rh, grid, last_tip)


def _last_tip_of(d):
    t = (d["root_age_decimal_year"] + d["root_height_years"]).median()
    return L.dec2date(t)


def quantities(S):
    """dict quantity -> dict(median, lo, hi, ess, rhat, mcse, kind, P_gt_1)"""
    out = {}

    def put(q, s, kind, scale_mcse=None, extra=None):
        mc = s.get("mcse_median_approx")
        if mc is None:
            mc = 1.2533 * s["R_sd"] / math.sqrt(s["ess_pooled"]) if s.get("ess_pooled", 0) and s["ess_pooled"] > 0 else np.nan
        r = dict(median=s.get("median", s.get("R_median")), lo=s.get("hpd_lo", s.get("R_hpd_lo")), hi=s.get("hpd_hi", s.get("R_hpd_hi")),
                 ess=s["ess_pooled"], rhat=s["rhat_split"], mcse=mc, kind=kind, P_gt_1=s.get("P_R_gt_1", np.nan), n_chains=s["n_chains"], n_samples=s["n_samples"])
        if extra:
            r.update(extra)
        out[q] = r

    put("clock rate (1e-3 subst/site/yr)", L.summarise(S.clock, "clock", 1e3), "scalar")
    put("root height (days before the last tip)", L.summarise(S.rh, "rh", L.YEAR_DAYS), "scalar")
    if S.grid is None:
        rs = A.r_summary(S.growth, S.gen, "exp")
        put("exponential model: R of the whole period", rs, "R")
        g = L.summarise(S.growth, "g", 1.0 / L.YEAR_DAYS)
        put("exponential model: growth rate (per day)", g, "scalar")
        return out
    g = S.grid
    for k, (a, b) in A.W_DATES.items():
        t0 = L.decimal_year(a)
        t1 = float(g.mids[0]) if k == "W3" else L.decimal_year(b)
        if t1 > S.t_last + 1e-9 or t1 <= t0:
            continue
        c = L.window_contrast(g, t0, t1, "interp")
        lab = f"R {k}: {a} -> " + (f"{L.dec2date(t1)} (mid-point of the most recent interval)" if k == "W3" else b)
        put(lab, A.r_summary([x @ c for x in S.lnNe], S.gen, k), "R", extra=dict(window=k, window_days=(t1 - t0) * L.YEAR_DAYS))
    # two sub-windows that were NOT pre-specified (exploratory, defined after the Delphy estimates were seen)
    t_a, t_b, t_c = L.decimal_year("2026-06-15"), L.decimal_year("2026-08-01"), float(g.mids[0])
    for key, t0, t1, lab in (("X1", t_a, t_b, "R exploratory: 2026-06-15 -> 2026-08-01"),
                             ("X2", t_b, t_c, f"R exploratory: 2026-08-01 -> {L.dec2date(t_c)} (mid-point of the most recent interval)")):
        if (t1 - t0) * L.YEAR_DAYS < MIN_WINDOW_DAYS or t1 > t_c + 1e-9:
            continue
        c = L.window_contrast(g, t0, t1, "interp")
        put(lab, A.r_summary([x @ c for x in S.lnNe], S.gen, key), "R", extra=dict(window=key, window_days=(t1 - t0) * L.YEAR_DAYS, exploratory=True))
    # the same periods with the end date of the primary analysis, for grids whose most recent interval is shorter (mid-point later than that date)
    if t_c > T_END_PRIMARY + 1.0 / L.YEAR_DAYS:
        for key, t0, lab in (("W3p", t_a, f"R W3 with the end date of the primary analysis: 2026-06-15 -> {L.dec2date(T_END_PRIMARY)}"),
                             ("X2p", t_b, f"R exploratory with the end date of the primary analysis: 2026-08-01 -> {L.dec2date(T_END_PRIMARY)}")):
            c = L.window_contrast(g, t0, T_END_PRIMARY, "interp")
            put(lab, A.r_summary([x @ c for x in S.lnNe], S.gen, key), "R",
                extra=dict(window=key, window_days=(T_END_PRIMARY - t0) * L.YEAR_DAYS, exploratory=key.startswith("X")))
    for i in range(1, g.K):
        c = L.pair_contrast(g, i)
        put(f"R pair {i:02d} (boundary {L.dec2date(g.t_anchor - i * g.D)})", A.r_summary([x @ c for x in S.lnNe], S.gen, f"pair {i}"), "R",
            extra=dict(pair=i, boundary=str(L.dec2date(g.t_anchor - i * g.D))))
    return out


EXPLORATORY_LABEL = "exploratory, defined after the estimates were seen"
MIN_WINDOW_DAYS = 14.0                           # shorter periods are not evaluated (genome sets that end early)
# mid-point of the most recent interval of the primary grid (20 parameters, cutoff 0.8 yr, last tip 2026-09-05)
T_END_PRIMARY = float(L.Grid.from_cutoff(20, 0.8, pd.Timestamp("2026-09-05")).mids[0])


def post_rows(analysis, engine, model_text, dataset, chains_run, chains_retained, q, last_tip, precision=None, note="", n_genomes=None):
    """rows of the table for the post: one row per reported parameter of one analysis"""
    rows = []
    base = dict(analysis=analysis, engine=engine, model=model_text, genome_set=dataset, n_genomes=n_genomes, last_collection_date=str(pd.Timestamp(last_tip).date()),
                chains_run=chains_run, chains_retained=chains_retained, note=note)

    def add(par, unit, med, lo, hi, ess, rhat, label="", extra=None):
        r = dict(base); r.update(parameter=par, unit=unit, median=med, hpd95_lo=lo, hpd95_hi=hi, ess_pooled=ess, rhat_split=rhat,
                                 criterion_ess200_rhat105_met=bool(ess >= 200 and rhat <= 1.05) if (ess == ess and rhat == rhat) else "", label=label)
        if extra:
            r.update(extra)
        rows.append(r)
    for k, v in q.items():
        if k.startswith("clock rate"):
            add("clock rate", "1e-3 substitutions/site/year", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"])
        elif k.startswith("root height"):
            lt = pd.Timestamp(last_tip)
            add("tMRCA", "date (from the root height; interval bounds: earlier date first)", str((lt - pd.Timedelta(days=round(v["median"]))).date()), str((lt - pd.Timedelta(days=round(v["hi"]))).date()),
                str((lt - pd.Timedelta(days=round(v["lo"]))).date()), v["ess"], v["rhat"], extra=dict(root_height_days_median=v["median"], root_height_days_hpd95_lo=v["lo"], root_height_days_hpd95_hi=v["hi"]))
        elif k.startswith("R W3 with the end date of the primary analysis"):
            add("R W3, end date of the primary analysis", "", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"],
                extra=dict(period=k.split(": ", 1)[1], P_R_gt_1=v.get("P_gt_1", np.nan)))
        elif k.startswith("R W1") or k.startswith("R W2") or k.startswith("R W3:") or k.startswith("R W3 "):
            add(k.split(":")[0], "", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"], extra=dict(period=k.split(": ", 1)[1] if ": " in k else "", P_R_gt_1=v.get("P_gt_1", np.nan)))
        elif k.startswith("R exploratory with the end date of the primary analysis"):
            add("R exploratory 1 Aug to the end date of the primary analysis", "", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"], label=EXPLORATORY_LABEL,
                extra=dict(period=k.split(": ", 1)[1], P_R_gt_1=v.get("P_gt_1", np.nan)))
        elif k.startswith("R exploratory"):
            add("R exploratory " + ("15 Jun to 1 Aug" if "2026-06-15" in k else "1 Aug to the mid-point of the most recent interval"), "", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"], label=EXPLORATORY_LABEL,
                extra=dict(period=k.split(": ", 1)[1], P_R_gt_1=v.get("P_gt_1", np.nan)))
        elif k.startswith("exponential model: R"):
            add("R of the whole period (exponential model)", "", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"], extra=dict(P_R_gt_1=v.get("P_gt_1", np.nan)))
        elif k.startswith("exponential model: growth"):
            add("growth rate (exponential model)", "per day", v["median"], v["lo"], v["hi"], v["ess"], v["rhat"],
                extra=dict(doubling_time_days_at_median=(math.log(2) / v["median"] if v["median"] > 0 else np.nan)))
    is_sky = any(k.startswith("R W1") for k in q)
    if is_sky and not any(k.startswith("R exploratory: 2026-08-01") for k in q):
        add("R exploratory 1 Aug to the mid-point of the most recent interval", "", np.nan, np.nan, np.nan, np.nan, np.nan, label=EXPLORATORY_LABEL,
            extra=dict(period=f"not evaluated: the genome set ends on {pd.Timestamp(last_tip).date()}, the period would be shorter than {MIN_WINDOW_DAYS:.0f} days"))
    if precision is not None:
        add("Skygrid precision", "", precision["median"], precision["hpd_lo"], precision["hpd_hi"], precision["ess_pooled"], precision["rhat_split"])
    return rows


def compare(qb, qd, model_key, delphy_label, what):
    rows = []
    for q, b in qb.items():
        if q not in qd:
            continue
        d = qd[q]
        diff = b["median"] - d["median"]
        se = math.sqrt(b["mcse"] ** 2 + d["mcse"] ** 2) if (b["mcse"] == b["mcse"] and d["mcse"] == d["mcse"]) else np.nan
        ok_b = bool(b["ess"] >= 200 and b["rhat"] <= 1.05); ok_d = bool(d["ess"] >= 200 and d["rhat"] <= 1.05)
        rows.append(dict(beast_model=MODELS[model_key]["label"], delphy_analysis=delphy_label, comparison=what, quantity=q,
                         delphy_median=d["median"], delphy_hpd_lo=d["lo"], delphy_hpd_hi=d["hi"], delphy_ess_pooled=d["ess"], delphy_rhat_split=d["rhat"],
                         delphy_chains=d["n_chains"], beast_median=b["median"], beast_hpd_lo=b["lo"], beast_hpd_hi=b["hi"], beast_ess_pooled=b["ess"],
                         beast_rhat_split=b["rhat"], beast_chains=b["n_chains"], difference_beast_minus_delphy=diff,
                         difference_relative_to_delphy_hpd_width=diff / (d["hi"] - d["lo"]) if d["hi"] > d["lo"] else np.nan,
                         z_difference_over_monte_carlo_error=diff / se if se and se == se and se > 0 else np.nan,
                         beast_median_inside_delphy_hpd=bool(d["lo"] <= b["median"] <= d["hi"]), delphy_median_inside_beast_hpd=bool(b["lo"] <= d["median"] <= b["hi"]),
                         hpd_overlap_fraction_of_narrower=max(0.0, min(b["hi"], d["hi"]) - max(b["lo"], d["lo"])) / min(b["hi"] - b["lo"], d["hi"] - d["lo"])
                         if min(b["hi"] - b["lo"], d["hi"] - d["lo"]) > 0 else np.nan,
                         beast_P_R_gt_1=b["P_gt_1"], delphy_P_R_gt_1=d["P_gt_1"], beast_criteria_met=ok_b, delphy_criteria_met=ok_d))
    return rows


# ----------------------------------------------------------------------------------------------------- trees
_NODE_ANN = re.compile(r"\[&([^\]]*)\]")


def read_annotated_tree(path):
    """first tree of a NEXUS file: (dict number -> label, newick string with annotations)"""
    op = gzip.open if str(path).endswith(".gz") else open
    labels, in_tr, newick = {}, False, None
    with op(path, "rt") as fh:
        for line in fh:
            s = line.strip()
            if s.lower().startswith("translate"):
                in_tr = True
                continue
            if in_tr:
                if s.startswith(";"):
                    in_tr = False
                    continue
                m = re.match(r"(\d+)\s+'?([^',;]+)'?,?;?$", s)
                if m:
                    labels[m.group(1)] = m.group(2)
                if s.endswith(";"):
                    in_tr = False
                continue
            if s.lower().startswith("tree "):
                newick = s[s.index("=") + 1:].strip()
                newick = re.sub(r"^\[&[RU]\]\s*", "", newick)
                break
    return labels, newick


def clades_with_support(path, acc_index):
    """dict bitmask -> posterior support for every internal node (except the root) of the first tree in the file"""
    labels, nw = read_annotated_tree(path)
    i, n = 0, len(nw)
    stack, out = [], {}
    cur = None

    def read_annotation(j):
        if j < n and nw[j] == "[":
            k = nw.index("]", j)
            return nw[j + 1:k], k + 1
        return None, j

    def skip_branch(j):
        if j < n and nw[j] == ":":
            j += 1
            _, j = read_annotation(j)
            m = re.match(r"[-+0-9.eE]+", nw[j:])
            j += m.end() if m else 0
        return j

    while i < n:
        ch = nw[i]
        if ch == "(":
            stack.append(0)
            i += 1
        elif ch == ",":
            i += 1
        elif ch == ")":
            mask = stack.pop()
            i += 1
            m = re.match(r"[^\[\]():,;]*", nw[i:])          # optional internal label
            i += m.end()
            ann, i = read_annotation(i)
            i = skip_branch(i)
            post = np.nan
            if ann:
                mm = re.search(r"(?:^|,|&)posterior=([-+0-9.eE]+)", ann)
                if mm:
                    post = float(mm.group(1))
            if stack:
                out[mask] = post
                stack[-1] |= mask
            else:
                root_mask = mask
        elif ch == ";":
            break
        else:
            m = re.match(r"'[^']*'|[^\[\]():,;]+", nw[i:])
            tok = m.group(0)
            i += m.end()
            _, i = read_annotation(i)
            i = skip_branch(i)
            lab = labels.get(tok, tok).strip("'")
            acc = lab.split("|")[0].split(".")[0]
            stack[-1] |= 1 << acc_index[acc]
    return out, root_mask


def compare_clades(delphy_tree, beast_tree, fasta, thr=0.9):
    accs = sorted(l[1:].strip().split("|")[0].split(".")[0] for l in open(fasta) if l.startswith(">"))
    assert len(set(accs)) == len(accs)
    idx = {a: i for i, a in enumerate(accs)}
    cd, rd = clades_with_support(delphy_tree, idx)
    cb, rb = clades_with_support(beast_tree, idx)
    assert rd == rb == (1 << len(accs)) - 1, "the two trees do not contain the same genomes"
    hd = {m: p for m, p in cd.items() if p == p and p >= thr}
    hb = {m: p for m, p in cb.items() if p == p and p >= thr}

    def conflicts(a, b):
        x = a & b
        return x != 0 and x != a and x != b

    rows = []
    n_conf_d = n_conf_b = 0
    for src, H, other_all, other_high in (("Delphy", hd, cb, hb), ("BEAST X", hb, cd, hd)):
        for m, p in H.items():
            conf = [o for o in other_high if conflicts(m, o)]
            conf_any = [o for o in other_all if conflicts(m, o)]
            rows.append(dict(clade_supported_in=src, clade_size=bin(m).count("1"), support=p, present_in_other_summary_tree=m in other_all,
                             support_in_other_summary_tree=other_all.get(m, np.nan),
                             conflicts_with_a_clade_of_support_ge_threshold_in_other=len(conf) > 0, n_conflicting_clades_ge_threshold=len(conf),
                             highest_support_among_conflicting_clades_of_other_tree=max([other_all[o] for o in conf_any], default=np.nan),
                             first_accessions=" ".join([a for a in accs if m >> idx[a] & 1][:3])))
            if conf:
                if src == "Delphy":
                    n_conf_d += 1
                else:
                    n_conf_b += 1
    T = pd.DataFrame(rows)
    shared = sum(1 for m in hd if m in hb)
    summ = dict(threshold=thr, n_genomes=len(accs), internal_nodes_delphy_tree=len(cd), internal_nodes_beast_tree=len(cb),
                clades_ge_threshold_delphy=len(hd), clades_ge_threshold_beast=len(hb), clades_ge_threshold_in_both=shared,
                delphy_clades_ge_threshold_present_in_beast_tree=sum(1 for m in hd if m in cb),
                beast_clades_ge_threshold_present_in_delphy_tree=sum(1 for m in hb if m in cd),
                delphy_clades_ge_threshold_in_conflict_with_beast_clade_ge_threshold=n_conf_d,
                beast_clades_ge_threshold_in_conflict_with_delphy_clade_ge_threshold=n_conf_b,
                clades_of_size_ge_2_in_both_trees=sum(1 for m in cd if m in cb))
    return summ, T


def beast_summary_tree(an, out_nexus, workdir, heights="median"):
    os.makedirs(workdir, exist_ok=True)
    files, n_trees = [], 0
    for c in an.kept:
        src = sorted(glob.glob(os.path.join(c.dir, MODELS[an.key]["xml"] + ".trees*")))
        assert src, c.dir
        dst = os.path.join(workdir, c.label + ".trees")
        if src[0].endswith(".gz"):
            with gzip.open(src[0], "rb") as fi, open(dst, "wb") as fo:
                shutil.copyfileobj(fi, fo)
        else:
            shutil.copy(src[0], dst)
        files.append(dst)
    burn = int(BURNIN * CHAIN_LENGTH)
    last = min(c.info["last_state"] for c in an.kept)
    if last < CHAIN_LENGTH:                      # only reached with --test (incomplete chains are not retained otherwise)
        burn = int(BURNIN * last)
    comb = os.path.join(workdir, an.key + "_combined.trees")
    p1 = subprocess.run(["logcombiner", "-trees", "-burnin", str(burn)] + files + [comb], capture_output=True, text=True)
    open(os.path.join(workdir, "logcombiner.stdout"), "w").write(p1.stdout[-5000:] + p1.stderr[-3000:])
    assert os.path.exists(comb) and os.path.getsize(comb) > 0, p1.stderr[-500:]
    n_trees = sum(1 for l in open(comb) if l.startswith("tree "))
    tmp = os.path.join(workdir, an.key + ".mcc.nexus")
    p2 = subprocess.run(["treeannotator", "-type", "mcc", "-heights", heights, comb, tmp], capture_output=True, text=True)
    open(os.path.join(workdir, "treeannotator.stdout"), "w").write(p2.stdout[-5000:] + p2.stderr[-3000:])
    assert os.path.exists(tmp) and os.path.getsize(tmp) > 0, p2.stderr[-500:]
    txt = open(tmp).read()
    hdr = (f"#NEXUS\n[Maximum clade credibility (MCC) tree, TreeAnnotator 10.5.0 (-type mcc -heights {heights}). Input: {n_trees} trees pooled with LogCombiner 10.5.0 "
           f"from {len(an.kept)} BEAST X 10.5.0 chains of model {MODELS[an.key]['label']} (seeds {' '.join(str(c.seed) for c in an.kept)}), "
           f"states {burn} to {last} of each chain, one tree every {an.P['tree_every']} states.]\n")
    open(out_nexus, "w").write(hdr + txt.replace("#NEXUS", "", 1).lstrip("\n"))
    return n_trees


def samples_frame_beast(an):
    fr = []
    g = an.grid
    for c in an.kept:
        p = c.post
        d = pd.DataFrame({"engine": "beastx", "model": an.model, "set": an.dataset, "analysis": "BEAST X | " + an.key, "coalescent_prior_cells": np.nan,
                          "seed": c.seed, "state": p["state"].values, "clock_rate": p["clock.rate"].values, "root_height_years": p[c.rh_col].values})
        d["root_age_decimal_year"] = an.t_last - d["root_height_years"]
        d["root_date"] = L.root_dates_from_height(d["root_height_years"].values, an.last_tip).astype(str)
        for col in ("joint", "prior", "likelihood", "kappa", "treeLength"):
            if col in p.columns:
                d[col] = p[col].values
        if c.is_sky:
            lp = p[c.lp_cols].values
            for i in range(lp.shape[1]):
                d[f"lnNe_{i + 1:02d}"] = lp[:, i]
            d["skygrid_precision"] = p["skygrid.precision"].values
            d["n_parameters"] = g.K; d["interval_days"] = g.D * 365; d["grid_anchor_decimal_year"] = g.t_anchor
            for k, w in an.windows("interp", "midpoint").items():
                if w["evaluable"]:
                    d[f"growth_rate_per_year_{k}"] = lp @ w["c"]
        else:
            d["growth_rate_per_year_exponential"] = p["exponential.growthRate"].values
            d["exponential_final_popsize_years"] = p["exponential.popSize"].values
        d["sample_spacing_states"] = int(an.P["log_every"])
        fr.append(d)
    return pd.concat(fr, ignore_index=True) if fr else pd.DataFrame()


def manifest_rows(an):
    rows = []
    for c in an.chains:
        cmd, started, ended = "", np.nan, np.nan
        p = os.path.join(c.dir, "command.txt")
        if os.path.exists(p):
            cmd = " ".join(open(p).read().split())
            cmd = re.sub(r"\S*/xml/", "", cmd)
        logs = sorted(glob.glob(os.path.join(c.dir, MODELS[an.key]["xml"] + ".log*")))
        hdr = [l for l in L.open_any(logs[0]) if l.startswith("#")][:4]
        opts = next((h[2:].strip() for h in hdr if h.startswith("# -")), "")
        opts = re.sub(r"\S*/xml/", "", opts)
        wall = np.nan
        so = os.path.join(c.dir, "stdout.txt")
        if os.path.exists(so):
            t = open(so, errors="replace").read()
            m = re.search(r"([0-9.]+) (hours|minutes|seconds)\s*$", t.strip().splitlines()[-1]) if t.strip() else None
            m2 = re.findall(r"^\s*([0-9.]+) (hours|minutes|seconds)\s*$", t, flags=re.M)
            if m2:
                v, u = m2[-1]
                wall = float(v) * {"hours": 3600, "minutes": 60, "seconds": 1}[u]
        rows.append(dict(engine="BEAST X 10.5.0 (BEAGLE 4.0.1)", model=an.model, set=an.dataset, analysis=an.key, seed=c.seed, steps_requested=CHAIN_LENGTH,
                         steps_reached=c.info["last_state"], threads=(re.search(r"-threads (\d+)", opts).group(1) if re.search(r"-threads (\d+)", opts) else ""),
                         coalescent_prior_cells="", low_population_barrier="", double_half_time_d="",
                         n_parameters=(an.P["num_parameters"] or ""), grid=(f"cutoff {an.P['cutoff']} yr" if an.P["cutoff"] else ""),
                         log_every=an.P["log_every"], tree_every=an.P["tree_every"], where_run="Modal", remote_job=os.environ.get("BEAST_JOB_ID", ""),
                         command_line=("beast " + opts) if opts else cmd, wall_time_s=wall, output_found=True,
                         entered_degenerate_state_root_date_rule=c.info["collapsed"], first_step_root_date_rule=c.info["first_collapsed_state"],
                         density_criterion_met="not applicable", retained="yes" if c.retained else "no", why=c.why, xml_sha256=an.P["sha256"]))
    return rows


# ----------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--out", default="work/deliverables")
    ap.add_argument("--delphy-tables", default="work/deliverables")
    ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    t0 = time.time()
    os.makedirs(a.out, exist_ok=True)
    pre = "TEST_" if a.test else ""
    AN, TABS, QB = {}, {}, {}
    for m in MODELS:
        an = load_model(m, a.indir, a.test)
        AN[m] = an
        print(m, "chains", len(an.chains), "retained", len(an.kept), "states", min(c.info["last_state"] for c in an.chains), "-", max(c.info["last_state"] for c in an.chains),
              round(time.time() - t0), "s", flush=True)
    # ---- retention and degenerate-state check
    ret = pd.DataFrame([retention_row_beast(an) for an in AN.values()])
    span = []
    for m, an in AN.items():
        if an.grid is not None:
            s = tree_span_minimum(an); s.insert(0, "analysis", an.name); span.append(s)
    span = pd.concat(span, ignore_index=True)
    agg = span.groupby("analysis").agg(lowest_lnNe_within_tree_span_any_chain=("lowest_lnNe_within_tree_span", "min"),
                                       largest_fraction_of_states_with_Ne_tau_below_1_day=("fraction_of_states_with_Ne_tau_below_1_day", "max"),
                                       largest_fraction_of_states_with_Ne_tau_below_7_days=("fraction_of_states_with_Ne_tau_below_7_days", "max"),
                                       lowest_lnNe_excluding_the_interval_with_the_root_any_chain=("lowest_lnNe_excluding_the_interval_with_the_root", "min"),
                                       largest_fraction_below_1_day_excluding_the_interval_with_the_root=(
                                           "fraction_of_states_with_Ne_tau_below_1_day_excluding_the_interval_with_the_root", "max")).reset_index()
    ret = ret.merge(agg, on="analysis", how="left")
    ret.to_csv(f"{a.out}/{pre}beast_retention_{SFX}.csv", index=False)
    span.to_csv(f"{a.out}/{pre}beast_degenerate_state_check_{SFX}.csv", index=False)
    # ---- tables per model
    head, diag, rts, ivs, wins, r0s = [], [], [], [], [], []
    for m, an in AN.items():
        if not an.kept:
            continue
        if an.grid is not None:
            tb = A.skygrid_tables(an)
            TABS[m] = tb
            h = A.headline(an, tb)
            for nm, lst in (("rt", rts), ("intervals", ivs), ("windows", wins), ("r0defs", r0s)):
                x = tb[nm].copy(); x.insert(0, "beast_model", MODELS[m]["label"]); x.insert(1, "dataset", an.dataset)
                x.insert(2, "chains_run", len(an.chains)); x.insert(3, "chains_retained", len(an.kept)); lst.append(x)
        else:
            h = A.headline(an)
        h.insert(4, "setting", "BEAST X 10.5.0, 8 chains x 40 M states, 30 % burn-in")
        h.insert(5, "coalescent_prior_cells", np.nan)
        rr = ret[ret["analysis"] == an.name].iloc[0]
        h["discarded_root_date_rule"] = rr["discarded_root_date_rule"]; h["discarded_density_criterion_only"] = 0
        h["label"] = ""
        head.append(h)
        d = pd.concat([A.per_chain_diagnostics(an), extra_diagnostics(an)], ignore_index=True); d.insert(0, "analysis_label", "BEAST X | " + m); diag.append(d)
        QB[m] = quantities(series_from_analysis(an))
        print(m, "tables", round(time.time() - t0), "s", flush=True)
    HEADB = pd.concat(head, ignore_index=True)
    DIAGB = pd.concat(diag, ignore_index=True)
    DIAGB.to_csv(f"{a.out}/{pre}chain_diagnostics_beast_{SFX}.csv", index=False)
    pd.concat(rts, ignore_index=True).to_csv(f"{a.out}/{pre}Rt_skygrid_beast_{SFX}.csv", index=False)
    pd.concat(ivs, ignore_index=True).to_csv(f"{a.out}/{pre}skygrid_ne_intervals_beast_{SFX}.csv", index=False)
    pd.concat(wins, ignore_index=True).to_csv(f"{a.out}/{pre}Rt_windows_beast_{SFX}.csv", index=False)
    pd.concat(r0s, ignore_index=True).to_csv(f"{a.out}/{pre}R0_definitions_beast_{SFX}.csv", index=False)
    # ---- headline summary of both engines
    hd = pd.read_csv(f"{a.delphy_tables}/phylo_summary_delphy_{SFX}.csv")
    ALLH = pd.concat([hd, HEADB], ignore_index=True)
    ALLH["engine"] = ALLH["engine"].replace({"delphy": "Delphy 1.4.1", "beast": "BEAST X 10.5.0"})
    ALLH.to_csv(f"{a.out}/{pre}phylo_summary_{SFX}.csv", index=False)
    # ---- comparison with Delphy
    samples = f"{a.delphy_tables}/posterior_samples_delphy_{SFX}.csv.gz"
    comp, QD = [], {}
    for m in QB:
        for lab, what in ((MODELS[m]["delphy"], MODELS[m]["like_for_like"]), (MODELS[m]["delphy_alt"], MODELS[m]["like_for_like"])):
            if lab not in QD:
                QD[lab] = quantities(series_from_delphy(lab, samples))
            qd = QD[lab]
            qb = QB[m]
            if m == "C_inrblike_inrbRetained":
                qb = {k: v for k, v in qb.items() if not k.startswith("R pair")}
                # windows are compared by name; W3 ends at the mid-point of the most recent interval of each grid
                qb = {(k if not k.startswith("R W3:") else "R W3 (end = mid-point of the most recent interval of each grid)"): v for k, v in qb.items()}
                qd = {(k if not k.startswith("R W3:") else "R W3 (end = mid-point of the most recent interval of each grid)"): v for k, v in qd.items()
                      if not k.startswith("R pair")}
            comp += compare(qb, qd, m, lab, what)
    COMP = pd.DataFrame(comp)
    is_pair = COMP["quantity"].str.startswith("R pair")
    COMP[~is_pair].to_csv(f"{a.out}/{pre}beast_vs_delphy_{SFX}.csv", index=False)
    COMP[is_pair].to_csv(f"{a.out}/{pre}beast_vs_delphy_Rt_{SFX}.csv", index=False)
    # ---- table for the post: one row per analysis and reported parameter
    RETD = pd.read_csv(f"{a.delphy_tables}/retention_by_analysis_{SFX}.csv").set_index("analysis")
    POST = []
    DMOD = {"Skygrid | primary | 8000 cells": ("Skygrid, 20 parameters, cutoff 0.8 yr, precision fixed at 4.06335, HKY, strict clock; 8,000 coalescent-prior cells (primary analysis)", "primary"),
            "Skygrid | gridWeekly32 | 8000 cells": ("Skygrid, 32 parameters, cutoff 0.613699 yr, precision fixed (same double-half time of 30 d), HKY, strict clock; 8,000 cells", "primary"),
            "Skygrid | inrbRetained | 8000 cells": ("Skygrid, 20 parameters, cutoff 0.8 yr, precision fixed at 4.06335, HKY, strict clock; 8,000 cells", "inrbRetained"),
            "exponential | primary | 8000 cells": ("exponential growth, HKY, strict clock; 8,000 cells", "primary")}
    for lab, (txt, ds) in DMOD.items():
        if lab in QD:
            sd_ = series_from_delphy(lab, samples)
            POST += post_rows("Delphy | " + lab, "Delphy 1.4.1", txt, ds, int(RETD.loc[lab, "chains_run"]), int(RETD.loc[lab, "chains_retained"]), QD[lab], sd_.last_tip,
                              note="posterior sample of the delivered analysis", n_genomes=sum(1 for l_ in L.open_any(f"inputs/aln/{ds}_{SFX}.fasta") if l_.startswith(">")))
    SHORT = {"A_skyfix_primary": "precision fixed as in Delphy", "B_skyest_primary": "precision estimated", "C_inrblike_primary": "settings of the previous post, primary genome set",
             "C_inrblike_inrbRetained": "settings of the previous post, genome set retained by the INRB-like rules", "D_exp_primary": "exponential growth"}
    for m, an in AN.items():
        if m not in QB:
            continue
        prec = L.summarise(an.series("skygrid.precision"), "precision") if (an.grid is not None and not TABS[m]["precision_fixed"]) else None
        POST += post_rows("BEAST X | " + m + " (" + SHORT[m] + ")", "BEAST X 10.5.0", MODELS[m]["label"], an.dataset, len(an.chains), len(an.kept), QB[m], an.last_tip, precision=prec,
                          note=("INRB-like: settings reconstructed from the wording of a public post; agreement is indicative, not a replication" if m.startswith("C_") else ""),
                          n_genomes=an.n_tips)
    POST = pd.DataFrame(POST)
    n_ok = 0
    for lab, an_name, pairs in (("Skygrid | primary | 8000 cells", "delphy_skygrid_primary_cells8000",
                                 (("clock rate", "clock rate"), ("tMRCA", "tMRCA"), ("R W1", "R0, early phase (W1)"), ("R W2", "R, W2:"), ("R W3", "R, W3:"))),
                                ("exponential | primary | 8000 cells", "delphy_exponential_primary_cells8000",
                                 (("clock rate", "clock rate"), ("tMRCA", "tMRCA"), ("R of the whole period (exponential model)", "R (whole period"),
                                  ("growth rate (exponential model)", "growth rate (whole period")))):
        for par, qn in pairs:
            p_ = POST[(POST["analysis"] == "Delphy | " + lab) & (POST["parameter"] == par)]
            h_ = hd[(hd["analysis"] == an_name) & hd["quantity"].str.startswith(qn)]
            assert len(p_) == 1 and len(h_) == 1, (lab, par, len(p_), len(h_))
            for c1, c2 in (("median", "median"), ("hpd95_lo", "hpd_lo"), ("hpd95_hi", "hpd_hi")):
                v1, v2 = p_[c1].iloc[0], h_[c2].iloc[0]
                same = (str(v1) == str(v2)) if par == "tMRCA" else (abs(float(v1) - float(v2)) <= 1e-9 * max(1.0, abs(float(v2))))
                assert same, ("Delphy row of the table for the post differs from the delivered summary", lab, par, c1, v1, v2)
                n_ok += 1
    print("Delphy rows of the table for the post equal to the delivered summary table:", n_ok, "values", flush=True)
    SENS = pd.read_csv(f"{a.delphy_tables}/sensitivity_{SFX}.csv", low_memory=False).set_index("analysis")
    X1, X2 = "exploratory_window_15Jun_to_1Aug", "exploratory_window_1Aug_to_midpoint_of_most_recent_interval"
    n_ok2, newly = 0, []
    for lab in ("primary | 8000 cells", "gridWeekly32 | 8000 cells", "inrbRetained | 8000 cells"):
        r_ = SENS.loc[lab]
        p_ = POST[POST["analysis"] == "Delphy | Skygrid | " + lab].set_index("parameter")
        mp = {"clock rate": ("clock_rate_x1e3_median", "clock_rate_x1e3_hpd_lo", "clock_rate_x1e3_hpd_hi"), "tMRCA": ("tmrca_median", "tmrca_hpd_early", "tmrca_hpd_late"),
              "R W1": ("R0_W1_median", "R0_W1_hpd_lo", "R0_W1_hpd_hi"), "R W2": ("R_W2_median", "R_W2_hpd_lo", "R_W2_hpd_hi"), "R W3": ("R_W3_median", "R_W3_hpd_lo", "R_W3_hpd_hi"),
              "R exploratory 15 Jun to 1 Aug": (X1 + "_R_median", X1 + "_R_hpd_lo", X1 + "_R_hpd_hi"),
              "R exploratory 1 Aug to the mid-point of the most recent interval": (X2 + "_R_median", X2 + "_R_hpd_lo", X2 + "_R_hpd_hi")}
        for par, cols in mp.items():
            for c1, c2 in zip(("median", "hpd95_lo", "hpd95_hi"), cols):
                v1, v2 = p_.loc[par, c1], r_[c2]
                if pd.isna(v2):                  # not in the Delphy table (exploratory periods were computed for the primary analysis only)
                    if not pd.isna(v1):
                        newly.append((lab, par))
                    continue
                same = (str(v1) == str(v2)) if par == "tMRCA" else (abs(float(v1) - float(v2)) <= 1e-9 * max(1.0, abs(float(v2))))
                assert same, ("Delphy row of the table for the post differs from the delivered sensitivity table", lab, par, c1, v1, v2)
                n_ok2 += 1
    newly = sorted(set(newly))
    print("Delphy rows of the table for the post equal to the delivered sensitivity table:", n_ok2, "values; computed for this table from the delivered posterior samples:", newly, flush=True)
    for lab, par in newly:
        k_ = (POST["analysis"] == "Delphy | Skygrid | " + lab) & (POST["parameter"] == par)
        POST.loc[k_, "note"] = "computed for this table from the posterior sample of the delivered analysis (not in the delivered tables)"
    k_ = POST["analysis"].str.startswith("Delphy | Skygrid | gridWeekly32") & POST["parameter"].str.contains("end date of the primary analysis")
    POST.loc[k_, "note"] = "computed for this table from the posterior sample of the delivered analysis (not in the delivered tables)"
    first = ["analysis", "engine", "model", "genome_set", "n_genomes", "last_collection_date", "chains_run", "chains_retained", "parameter", "period", "unit", "median", "hpd95_lo", "hpd95_hi", "ess_pooled", "rhat_split", "criterion_ess200_rhat105_met", "label"]
    POST = POST[[c for c in first if c in POST.columns] + [c for c in POST.columns if c not in first]]
    POST.to_csv(f"{a.out}/{pre}post_table_{SFX}.csv", index=False)
    print("post table rows", len(POST), "| analyses", POST["analysis"].nunique(), flush=True)
    # ---- precision (models with estimated precision)
    prow = []
    for m, an in AN.items():
        if an.grid is not None and an.kept and not TABS[m]["precision_fixed"]:
            s = L.summarise(an.series("skygrid.precision"), "precision")
            dh = [L.YEAR_DAYS * (math.log(2) ** 2) * an.grid.D * x for x in (s["median"], s["hpd_lo"], s["hpd_hi"])]
            prow.append(dict(beast_model=MODELS[m]["label"], dataset=an.dataset, precision_median=s["median"], precision_hpd_lo=s["hpd_lo"], precision_hpd_hi=s["hpd_hi"],
                             ess_pooled=s["ess_pooled"], rhat_split=s["rhat_split"], interval_days=an.grid.D * L.YEAR_DAYS,
                             double_half_time_d_at_median=dh[0], double_half_time_d_at_hpd_lo=dh[1], double_half_time_d_at_hpd_hi=dh[2],
                             sd_of_change_of_lnNe_per_interval_at_median=1 / math.sqrt(s["median"]),
                             delphy_fixed_precision=4.06335, delphy_double_half_time_d=30.0,
                             note="double-half time = ln(2)^2 x interval x precision (relation used for the Delphy setting: 30 d with 15.37-d intervals = 4.06335)"))
    pd.DataFrame(prow).to_csv(f"{a.out}/{pre}beast_skygrid_precision_{SFX}.csv", index=False)
    # ---- posterior samples, manifest
    PSB = pd.concat([samples_frame_beast(an) for an in AN.values() if an.kept], ignore_index=True)
    PSB.to_csv(f"{a.out}/{pre}posterior_samples_beast_{SFX}.csv.gz", index=False, float_format="%.8g", compression="gzip")
    MANB = pd.DataFrame([r for an in AN.values() for r in manifest_rows(an)])
    MANB.to_csv(f"{a.out}/{pre}run_manifest_beast_{SFX}.csv", index=False)
    # ---- summary tree of model A and clade comparison
    try:
        anA = AN["A_skyfix_primary"]
        out_tree = f"{a.out}/{pre}beast_skygrid_mcc_{SFX}.nexus"
        n_trees = beast_summary_tree(anA, out_tree, f"work/beast/{pre}summary_tree")
        summ, CT = compare_clades(f"{a.delphy_tables}/delphy_skygrid_mcc_{SFX}.nexus", out_tree, anA.fasta, 0.9)
        summ.update(beast_trees=n_trees, beast_chains=len(anA.kept), delphy_tree=f"delphy_skygrid_mcc_{SFX}.nexus", beast_tree=f"beast_skygrid_mcc_{SFX}.nexus",
                    note="clades of the two summary trees only; taxa matched by accession without version")
        pd.DataFrame([summ]).to_csv(f"{a.out}/{pre}beast_vs_delphy_clades_summary_{SFX}.csv", index=False)
        CT.to_csv(f"{a.out}/{pre}beast_vs_delphy_clades_{SFX}.csv", index=False)
        print("clades:", summ, flush=True)
    except Exception as e:                       # the tree step needs the BEAST X programs on the PATH
        print("summary tree step failed:", repr(e), flush=True)
        if not a.test:
            raise
    pd.set_option("display.width", 250); pd.set_option("display.max_columns", 40); pd.set_option("display.max_colwidth", 50)
    print(ret[["analysis", "chains_run", "discarded_root_date_rule", "discarded_incomplete", "chains_retained", "states_reached_min", "states_reached_max",
               "lowest_lnNe_within_tree_span_any_chain", "largest_fraction_of_states_with_Ne_tau_below_1_day"]].to_string(index=False))
    pooled = DIAGB[DIAGB["scope"] == "pooled retained chains"]
    print("pooled diagnostics rows", len(pooled), "not meeting the criteria", int((pooled["criterion_ess200_rhat105_met"] == False).sum()))
    x = COMP[~is_pair][["beast_model", "delphy_analysis", "quantity", "delphy_median", "beast_median", "beast_hpd_lo", "beast_hpd_hi", "beast_ess_pooled", "beast_rhat_split",
                        "z_difference_over_monte_carlo_error", "beast_median_inside_delphy_hpd"]].copy()
    x["beast_model"] = x["beast_model"].str.slice(0, 16); x["quantity"] = x["quantity"].str.slice(0, 40); x["delphy_analysis"] = x["delphy_analysis"].str.slice(10, 40)
    print(x.round(3).to_string(index=False))
    print("done", round(time.time() - t0), "s")


if __name__ == "__main__":
    main()
