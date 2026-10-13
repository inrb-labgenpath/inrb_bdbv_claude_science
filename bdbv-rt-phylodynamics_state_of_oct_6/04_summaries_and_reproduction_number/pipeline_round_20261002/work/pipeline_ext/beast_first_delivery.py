#!/usr/bin/env python
"""beast_first_delivery.py - first set of files of the BEAST X cross-check of the 24 Sep round.

usage (workspace root, conda environment `beast`):
    python work/pipeline_ext/beast_first_delivery.py --in <directory with one sub-directory per chain> --out <dir>
           --tables work/report_inputs/tables --versions <json: Delphy table -> sha256> [--mcc <BEAST MCC tree of model A>] [--test]

Writes (prefix TEST_ with --test, which accepts incomplete chains: code test only, the values are not results)
  beast_estimates_<SFX>.csv                     one row per BEAST X analysis and quantity
  beast_ne_intervals_<SFX>.csv                  ln(Ne tau) per interval and analysis, columns of the Delphy table of intervals
  beast_posterior_<analysis>_<SFX>.csv.gz       every logged state after burn-in of every retained chain, all columns of the BEAST log
  delphy_counterparts_estimates_<SFX>.csv       the same quantities of the Delphy counterparts, copied from the Delphy tables (text of the cells), with table and row
  delphy_counterparts_ne_intervals_<SFX>.csv    intervals of the Delphy counterparts (primary analysis: rows of the Delphy table; other two: computed from the Delphy posterior sample, marked)
  beast_vs_delphy_reading_<SFX>.csv             the reading rule (column reading_rule_version) applied to the rows of the two tables of estimates
  beast_first_delivery_facts_<SFX>.json         counts used by the note

Definitions are those of the Delphy tables: the functions of rt_lib.py / analysis.py / beast_analysis.py are called, nothing is re-implemented:
  growth rate of a period = (ln Ne tau(t1) - ln Ne tau(t0)) / (t1 - t0), ln Ne tau(t) read by linear interpolation between the mid-points of the intervals
  (constant beyond the mid-point of the most recent interval); R = (1 + r sd^2 / mu)^(mu^2 / sd^2), mu 15.3 d, sd 9.3 d; year of 365 d;
  decimal date = year + (day of the year - 1) / 365; W3 and the second exploratory period end at the mid-point of the most recent interval of the analysis.
Retention of chains (rules fixed before the first BEAST X run): a chain is used if it reached the requested length and if no logged state
  has a root date on or after the day before the first collection date; burn-in = first 30 % of the logged states.
"""
import argparse, glob, gzip, hashlib, json, math, os, re, sys
import numpy as np
import pandas as pd
sys.path.insert(0, "work/pipeline_ext")
import rt_lib as L, analysis as A
import beast_analysis as B

SFX = B.SFX
RULE_VERSION = os.environ.get("RULE_VERSION_ID", "program_comparison_reading_rule_20260926.md")          # the reading rule applied (file program_comparison_reading_rule_20260926.md)
RULE_NUMBER = os.environ.get("RULE_VERSION_NUMBER", "12")
EXPL = B.EXPLORATORY_LABEL
SHORT = {"A_skyfix_primary": "A", "B_skyest_primary": "B", "C_inrblike_primary": "C, primary set", "C_inrblike_inrbRetained": "C, INRB-like retained set", "D_exp_primary": "D"}
# Delphy counterpart of each BEAST X analysis: (label in the sensitivity table / table of the exponential model, label in the posterior sample file, what differs apart from the program)
COUNTERPART = {
    "A_skyfix_primary": dict(table="sensitivity", key="primary | 8000 cells", samples="Skygrid | primary | 8000 cells", name="Delphy primary analysis",
                             differs="priors of the parameters other than the population sizes; Delphy evaluates the coalescent prior on a time grid (8,000 cells) and has its low-population barrier switched on"),
    "B_skyest_primary": dict(table="sensitivity", key="primary | 8000 cells", samples="Skygrid | primary | 8000 cells", name="Delphy primary analysis",
                             differs="smoothing: precision estimated in BEAST X (prior Gamma(0.001, 1000)), fixed at 4.06335 in Delphy; otherwise as for A"),
    "C_inrblike_primary": dict(table="sensitivity", key="gridWeekly32 | 8000 cells", samples="Skygrid | gridWeekly32 | 8000 cells", name="Delphy, weekly grid (gridWeekly32)",
                               differs="substitution model (GTR in BEAST X, HKY in Delphy); smoothing (precision estimated in BEAST X, fixed at 8.64226 in Delphy); otherwise as for A"),
    "C_inrblike_inrbRetained": dict(table="sensitivity", key="inrbRetained | 8000 cells", samples="Skygrid | inrbRetained | 8000 cells", name="Delphy, INRB-like retained set (inrbRetained)",
                                    differs="grid (32 parameters, cut-off 0.613699 yr in BEAST X; 20 parameters, cut-off 0.8 yr in Delphy), therefore other end date of W3; substitution model (GTR / HKY); "
                                            "smoothing (precision estimated / fixed at 4.06335); otherwise as for A"),
    "D_exp_primary": dict(table="early_window_exponential", key="primary | 8000 cells", samples="exponential | primary | 8000 cells", name="Delphy, exponential model",
                          differs="prior of the growth rate (Laplace, mean 0, scale 100 per year in BEAST X; Laplace, mean 0.001, scale 30.70 per year in Delphy); prior of the population size (1/x in BEAST X); "
                                  "Delphy bounds the population size below (1 day); otherwise as for A"),
}
DIFFERS_FROM_PRIMARY = {
    "A_skyfix_primary": "program only (same grid, same fixed precision, HKY, same genomes); see the note for the priors of the other parameters",
    "B_skyest_primary": "smoothing estimated",
    "C_inrblike_primary": "other grid (32 parameters, cut-off 0.613699 yr); smoothing estimated; other substitution model (GTR)",
    "C_inrblike_inrbRetained": "other genome set (515 genomes, last collection date 2026-08-09); other grid; smoothing estimated; other substitution model (GTR)",
    "D_exp_primary": "other tree prior (exponential growth)",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


# ------------------------------------------------------------------------------------------------ BEAST X rows
def left_out_text(an):
    out = [f"seed {c.seed}: {c.why}" for c in an.chains if not c.retained]
    return "; ".join(out) if out else "no chain left out"


def base_cols(an, P):
    kept = an.kept
    n_used = int(sum(len(c.post) for c in kept))
    return dict(analysis=an.key, analysis_short=SHORT[an.key], engine="BEAST X 10.5.0", model=an.model, genome_set=an.dataset, n_genomes=an.n_tips,
                first_collection_date=str(an.first_tip.date()), last_collection_date=str(an.last_tip.date()),
                chains_run=len(an.chains), chains_retained=len(kept), chains_left_out=len(an.chains) - len(kept), criterion_chain_left_out=left_out_text(an),
                seeds_retained=" ".join(str(c.seed) for c in kept), states_requested_per_chain=B.CHAIN_LENGTH,
                burnin_states=int(kept[0].post["state"].iloc[0]) if kept else np.nan, log_every_states=int(P["log_every"]),
                n_states_used=n_used, n_states_used_per_chain=" ".join(str(len(c.post)) for c in kept),
                skygrid_parameters=(an.grid.K if an.grid is not None else np.nan), skygrid_cutoff_years=(float(P["cutoff"]) if an.grid is not None else np.nan),
                skygrid_interval_days=(an.grid.D * L.YEAR_DAYS if an.grid is not None else np.nan),
                grid_anchor_last_collection_date_decimal=(an.grid.t_anchor if an.grid is not None else np.nan),
                generation_time_mean_d=an.gen.mean_d, generation_time_sd_d=an.gen.sd_d, year_length_days=L.YEAR_DAYS)


def crit(ess, rhat):
    return bool(ess >= 200 and rhat <= 1.05) if (ess == ess and rhat == rhat) else ""


WINDOW_DEF = ("growth rate = (ln Ne tau(end) - ln Ne tau(start)) / (end - start), ln Ne tau read by linear interpolation between the mid-points of the intervals; "
              "R = (1 + r sd^2 / mu)^(mu^2 / sd^2) with r per year of the decimal dates (365 d) and the generation time (mean 15.3 d, SD 9.3 d) converted into years with 365.25 d, "
              "as in the delivered tables (README of the pipeline, section 3.5); growth rate per day = growth rate per year / 365")


def windows_of(S):
    """periods of one Skygrid analysis (either program): list of dicts; 'rs' = summary of analysis.r_summary or None if the period is not evaluated"""
    g = S.grid
    t_end = float(g.mids[0])
    defs = [("W1", "R, 1 March - 15 May (W1)", "2026-03-01", "2026-05-15", "fixed before the estimates"),
            ("W2", "R, 15 May - 15 June (W2)", "2026-05-15", "2026-06-15", "fixed before the estimates"),
            ("W3", "R, 15 June - mid-point of the most recent interval (W3)", "2026-06-15", None, "fixed before the estimates"),
            ("W3b", "R, 15 June - 15 August (W3b)", "2026-06-15", "2026-08-15", "fixed before the estimates"),
            ("X1", "R, 15 June - 1 August", "2026-06-15", "2026-08-01", EXPL),
            ("X2", "R, 1 August - mid-point of the most recent interval", "2026-08-01", None, EXPL)]
    if t_end > B.T_END_PRIMARY + 1.0 / L.YEAR_DAYS:
        e = str(L.dec2date(B.T_END_PRIMARY))
        defs += [("W3p", "R, 15 June - end date of the Delphy primary analysis", "2026-06-15", e, "additional row: W3 with the end date of the Delphy primary analysis"),
                 ("X2p", "R, 1 August - end date of the Delphy primary analysis", "2026-08-01", e, EXPL + "; additional row with the end date of the Delphy primary analysis")]
    out = []
    for key, name, a_, b_, status in defs:
        t0 = L.decimal_year(a_)
        t1 = t_end if b_ is None else (B.T_END_PRIMARY if key in ("W3p", "X2p") else L.decimal_year(b_))
        days = (t1 - t0) * L.YEAR_DAYS
        w = dict(key=key, name=name, status=status, start=a_, end=str(L.dec2date(t1)), days=round(days, 4), end_is_midpoint=bool(b_ is None), note="", rs=None, shown=None)
        if t1 > S.t_last + 1e-9:
            w["note"] = f"not evaluated: the period ends after the last collection date ({S.last_tip.date()})"
        elif days < 1.0:
            w["note"] = (f"not evaluated: the period has no length (the mid-point of the most recent interval is {L.dec2ts(t1)})")
        else:
            if b_ is not None and key not in ("W3p", "X2p") and t1 > t_end + 1e-9:
                w["note"] = "the end date is later than the mid-point of the most recent interval: the level of that interval is used from the mid-point on"
            c = L.window_contrast(g, t0, t1, "interp")
            w["rs"] = A.r_summary([x @ c for x in S.lnNe], S.gen, key)
            if key.startswith("X") and days < B.MIN_WINDOW_DAYS:
                w["shown"] = False
                w["note"] = (f"period of {days:.1f} days, shorter than {B.MIN_WINDOW_DAYS:.0f} days: not evaluated in the table for the post (limit written into beast_analysis.py on 25 Sep, "
                             "before any chain was complete); given here because the row was asked for; the value is the growth between the two most recent intervals")
        out.append(w)
    return out


def skygrid_rows(an, P):
    base = base_cols(an, P)
    S = B.series_from_analysis(an)
    Q = B.quantities(S)                                     # the dictionary the table for the post is built from: used as a check below
    rows = []
    s = L.summarise(S.clock, "clock", 1e3)
    rows.append(dict(base, quantity="evolutionary rate", quantity_key="rate", unit="1e-3 substitutions/site/year", status="", median=s["median"], hpd95_lo=s["hpd_lo"], hpd95_hi=s["hpd_hi"],
                     ess_pooled=s["ess_pooled"], rhat_split=s["rhat_split"], criterion_ess200_rhat105_met=crit(s["ess_pooled"], s["rhat_split"]),
                     definition="log column clock.rate x 1000"))
    assert abs(s["median"] - Q["clock rate (1e-3 subst/site/yr)"]["median"]) < 1e-12
    tm = A.tmrca_summary(an)
    rows.append(dict(base, quantity="tMRCA", quantity_key="tMRCA", unit="date; hpd95_lo = earlier date", status="", median=tm["tmrca_median"], hpd95_lo=tm["tmrca_hpd_early"], hpd95_hi=tm["tmrca_hpd_late"],
                     root_height_days_median=tm["root_height_days_median"], root_height_days_hpd95_lo=tm["root_height_days_hpd_lo"], root_height_days_hpd95_hi=tm["root_height_days_hpd_hi"],
                     ess_pooled=tm["ess_pooled"], rhat_split=tm["rhat_split"], criterion_ess200_rhat105_met=crit(tm["ess_pooled"], tm["rhat_split"]),
                     definition="last collection date - round(root height x 365) days; median and 95 % HPD of the root height (log column treeModel.rootHeight, years)"))
    for w in windows_of(S):
        r = dict(base, quantity=w["name"], quantity_key=w["key"], unit="", status=w["status"], period_start=w["start"], period_end=w["end"], period_days=w["days"],
                 period_end_is_midpoint_of_most_recent_interval=w["end_is_midpoint"], definition=WINDOW_DEF, note=w["note"])
        rs = w["rs"]
        if rs is not None:
            r.update(median=rs["R_median"], hpd95_lo=rs["R_hpd_lo"], hpd95_hi=rs["R_hpd_hi"], P_R_gt_1=rs["P_R_gt_1"], growth_rate_per_day_median=rs["growth_per_day_median"],
                     growth_rate_per_day_hpd95_lo=rs["growth_per_day_hpd_lo"], growth_rate_per_day_hpd95_hi=rs["growth_per_day_hpd_hi"], n_samples_R_set_to_0=rs["n_samples_R_set_to_0"],
                     ess_pooled=rs["ess_pooled"], rhat_split=rs["rhat_split"], criterion_ess200_rhat105_met=crit(rs["ess_pooled"], rs["rhat_split"]))
            if w["shown"] is False:
                r["shown_in_the_table_for_the_post"] = False
            qk = [k for k, v in Q.items() if v.get("window") == w["key"]]                   # check against beast_analysis.quantities
            if qk:
                assert abs(Q[qk[0]]["median"] - rs["R_median"]) < 1e-12 and abs(Q[qk[0]]["lo"] - rs["R_hpd_lo"]) < 1e-12, (an.key, w["key"])
        rows.append(r)
    if "skygrid.precision" in an.kept[0].post.columns and an.kept[0].post["skygrid.precision"].nunique() > 1:
        s = L.summarise(an.series("skygrid.precision"), "precision")
        rows.append(dict(base, quantity="Skygrid precision (estimated)", quantity_key="precision", unit="", status="", median=s["median"], hpd95_lo=s["hpd_lo"], hpd95_hi=s["hpd_hi"],
                         ess_pooled=s["ess_pooled"], rhat_split=s["rhat_split"], criterion_ess200_rhat105_met=crit(s["ess_pooled"], s["rhat_split"]), definition="log column skygrid.precision"))
    return rows


def exponential_rows(an, P):
    base = base_cols(an, P)
    et = A.exponential_table(an)
    rows = []
    rows.append(dict(base, quantity="evolutionary rate", quantity_key="rate", unit="1e-3 substitutions/site/year", status="", median=et["clock_rate_x1e3_median"], hpd95_lo=et["clock_rate_x1e3_hpd_lo"],
                     hpd95_hi=et["clock_rate_x1e3_hpd_hi"], ess_pooled=et["clock_ess_pooled"], rhat_split=et["clock_rhat_split"],
                     criterion_ess200_rhat105_met=crit(et["clock_ess_pooled"], et["clock_rhat_split"]), definition="log column clock.rate x 1000"))
    tm = A.tmrca_summary(an)
    rows.append(dict(base, quantity="tMRCA", quantity_key="tMRCA", unit="date; hpd95_lo = earlier date", status="", median=tm["tmrca_median"], hpd95_lo=tm["tmrca_hpd_early"], hpd95_hi=tm["tmrca_hpd_late"],
                     root_height_days_median=tm["root_height_days_median"], root_height_days_hpd95_lo=tm["root_height_days_hpd_lo"], root_height_days_hpd95_hi=tm["root_height_days_hpd_hi"],
                     ess_pooled=tm["ess_pooled"], rhat_split=tm["rhat_split"], criterion_ess200_rhat105_met=crit(tm["ess_pooled"], tm["rhat_split"]),
                     definition="last collection date - round(root height x 365) days; median and 95 % HPD of the root height (log column treeModel.rootHeight, years)"))
    ge, gr = et["growth_ess_pooled"], et["growth_rhat_split"]
    rows.append(dict(base, quantity="growth rate (exponential model)", quantity_key="growth", unit="per day", status="", median=et["growth_rate_per_day_median"], hpd95_lo=et["growth_rate_per_day_hpd_lo"],
                     hpd95_hi=et["growth_rate_per_day_hpd_hi"], ess_pooled=ge, rhat_split=gr, criterion_ess200_rhat105_met=crit(ge, gr), definition="log column exponential.growthRate (per year) / 365"))
    rows.append(dict(base, quantity="doubling time (exponential model)", quantity_key="doubling", unit="days", status="", median=et["doubling_time_d_median"], hpd95_lo=et["doubling_time_d_hpd_lo"],
                     hpd95_hi=et["doubling_time_d_hpd_hi"], ess_pooled=ge, rhat_split=gr, criterion_ess200_rhat105_met=crit(ge, gr),
                     definition="ln 2 / growth rate: at the median of the growth rate and at the two bounds of its 95 % HPD interval (as in the Delphy table); ESS and R-hat are those of the growth rate"))
    rows.append(dict(base, quantity="R of the whole period (exponential model)", quantity_key="R_exp", unit="", status="", median=et["R_median"], hpd95_lo=et["R_hpd_lo"], hpd95_hi=et["R_hpd_hi"],
                     P_R_gt_1=et["P_R_gt_1"], growth_rate_per_day_median=et["growth_rate_per_day_median"], growth_rate_per_day_hpd95_lo=et["growth_rate_per_day_hpd_lo"],
                     growth_rate_per_day_hpd95_hi=et["growth_rate_per_day_hpd_hi"], ess_pooled=ge, rhat_split=gr, criterion_ess200_rhat105_met=crit(ge, gr),
                     definition="R from the growth rate with the gamma generation time; ESS and R-hat are those of the growth rate"))
    return rows


HPD_METHOD = ("shortest interval that holds m = ceil(0.95 n) of the n states used: the states are sorted, and of all runs of m consecutive values the one with the smallest difference between its last and its first "
              "value is taken (function hpd of the pipeline)")


def dec2ts_exact(d):
    """decimal year -> time stamp (years of the calendar; microseconds), not rounded to the second"""
    y = int(math.floor(d))
    start = pd.Timestamp(y, 1, 1)
    nd = (pd.Timestamp(y + 1, 1, 1) - start).days
    return (start + pd.Timedelta(days=(d - y) * nd)).round("ms")


def add_grid_columns(tb, grid, n_states):
    """number of parameters, edges and mid-point as date and time, how the HPD interval is computed"""
    tb = tb.copy()
    iso = lambda t: t.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3]
    oe, ne, mp, od, nd_ = [], [], [], [], []
    for i in tb["interval"].astype(int):
        older, newer = grid.edges(i)
        oe.append(iso(dec2ts_exact(older)) if i < grid.K else "open"); ne.append(iso(dec2ts_exact(newer))); mp.append(iso(dec2ts_exact(grid.mids[i - 1])) if i < grid.K else "")
        od.append(older if i < grid.K else np.nan); nd_.append(newer)
    tb["interval_number_1_is_the_most_recent"] = tb["interval"].astype(int)
    tb["n_parameters"] = int(grid.K)
    tb["older_edge_date_time"] = oe; tb["newer_edge_date_time"] = ne; tb["mid_point_date_time"] = mp
    tb["older_edge_decimal_year"] = od; tb["newer_edge_decimal_year"] = nd_
    tb["date_time_convention"] = "decimal year of the pipeline (year + day of the year / days of the year) converted back to date and time; no time zone is implied: collection dates are days"
    tb["n_states_used_for_the_interval"] = int(n_states)
    m = int(np.ceil(0.95 * n_states))
    tb["hpd95_method"] = HPD_METHOD
    tb["hpd95_m_states_inside"] = m
    tb["hpd95_same_run_as_with_k_equal_floor_of_0.95_n_steps"] = bool(m - 1 == int(np.floor(0.95 * n_states)))
    return tb


def intervals_frame(an, mcc_counts=None, tree_counts=None):
    tb = A.skygrid_tables(an, mcc_counts=mcc_counts, tree_counts=tree_counts)["intervals"].copy()
    for c_ in ("n_coalescent_events_MCC_tree", "n_coalescent_events_posterior_mean", "n_coalescent_events_posterior_q025", "n_coalescent_events_posterior_q975"):
        if c_ not in tb.columns:
            tb[c_] = np.nan
    tb.insert(0, "analysis", an.key); tb.insert(1, "engine", "BEAST X 10.5.0"); tb.insert(2, "model", an.model); tb.insert(3, "genome_set", an.dataset)
    tb.insert(4, "coalescent_prior_cells", ""); tb.insert(5, "setting", f"BEAST X 10.5.0, {len(an.chains)} chains x {B.CHAIN_LENGTH // 1_000_000} M states, 30 % burn-in")
    tb.insert(6, "chains_run", len(an.chains)); tb.insert(7, "chains_retained", len(an.kept)); tb.insert(8, "n_states_used", int(sum(len(c.post) for c in an.kept)))
    return add_grid_columns(tb, an.grid, int(sum(len(c.post) for c in an.kept)))


def tree_counts_of(an, max_trees_per_chain=40):
    """posterior number of coalescent events per interval from the sampled trees after burn-in (function of the Delphy table)"""
    an.trees = {}
    for c in an.kept:
        src = sorted(glob.glob(os.path.join(c.dir, B.MODELS[an.key]["xml"] + ".trees*")))
        if not src:
            return None
        an.trees[c.seed] = src[0]
    return A.coalescent_counts(an, mcc_path=None, max_trees_per_chain=max_trees_per_chain)["posterior"]


def posterior_frame(an):
    fr = []
    for c in an.kept:
        d = c.post.copy()
        d.insert(0, "analysis", an.key); d.insert(1, "chain", c.label); d.insert(2, "seed", c.seed)
        fr.append(d)
    return pd.concat(fr, ignore_index=True)


# ------------------------------------------------------------------------------------------------ Delphy counterparts from the Delphy tables
class Delivered:
    def __init__(self, tables, versions):
        self.dir, self.V, self.cache = tables, versions, {}

    def table(self, name):
        fn = f"{name}_{SFX}.csv"
        if fn not in self.cache:
            p = os.path.join(self.dir, fn)
            assert fn in self.V, f"no version recorded for {fn}"
            assert sha256(p) == self.V[fn]["sha256"], f"{fn}: the local copy is not the delivered version"
            self.cache[fn] = pd.read_csv(p, dtype=str, keep_default_na=False)
        return fn, self.cache[fn]

    def has(self, name, col, key, extra=None):
        fn, d = self.table(name)
        m = d[col] == key
        for k_, v_ in (extra or {}).items():
            m &= d[k_] == v_
        return int(m.sum()) == 1

    def row(self, name, col, key, extra=None):
        fn, d = self.table(name)
        m = d[col] == key
        if extra:
            for k_, v_ in extra.items():
                m &= d[k_] == v_
        idx = list(d.index[m])
        assert len(idx) == 1, (fn, key, extra, len(idx))
        return fn, self.V[fn]["version_id"], idx[0] + 1, d.loc[idx[0]]


SENS = {"rate": ("evolutionary rate", "1e-3 substitutions/site/year", "clock_rate_x1e3_median", "clock_rate_x1e3_hpd_lo", "clock_rate_x1e3_hpd_hi", None, "clock_rate_ess_pooled", "clock_rate_rhat_split"),
        "tMRCA": ("tMRCA", "date; hpd95_lo = earlier date", "tmrca_median", "tmrca_hpd_early", "tmrca_hpd_late", None, "root_age_ess_pooled", "root_age_rhat_split"),
        "W1": ("R, 1 March - 15 May (W1)", "", "R0_W1_median", "R0_W1_hpd_lo", "R0_W1_hpd_hi", "R0_W1_P_gt_1", "R0_W1_ess_pooled", "R0_W1_rhat_split"),
        "W2": ("R, 15 May - 15 June (W2)", "", "R_W2_median", "R_W2_hpd_lo", "R_W2_hpd_hi", "R_W2_P_gt_1", "R_W2_ess_pooled", "R_W2_rhat_split"),
        "W3": ("R, 15 June - mid-point of the most recent interval (W3)", "", "R_W3_median", "R_W3_hpd_lo", "R_W3_hpd_hi", "R_W3_P_gt_1", "R_W3_ess_pooled", "R_W3_rhat_split"),
        "W3b": ("R, 15 June - 15 August (W3b)", "", "R_W3b_median", "R_W3b_hpd_lo", "R_W3b_hpd_hi", "R_W3b_P_gt_1", "R_W3b_ess_pooled", "R_W3b_rhat_split"),
        "X1": ("R, 15 June - 1 August", "", *[f"exploratory_window_15Jun_to_1Aug_{x}" for x in ("R_median", "R_hpd_lo", "R_hpd_hi", "P_R_gt_1", "ess_pooled", "rhat_split")]),
        "X2": ("R, 1 August - mid-point of the most recent interval", "", *[f"exploratory_window_1Aug_to_midpoint_of_most_recent_interval_{x}" for x in ("R_median", "R_hpd_lo", "R_hpd_hi", "P_R_gt_1", "ess_pooled", "rhat_split")])}
EXPO = {"rate": ("evolutionary rate", "1e-3 substitutions/site/year", "clock_rate_x1e3_median", "clock_rate_x1e3_hpd_lo", "clock_rate_x1e3_hpd_hi", None, "clock_ess_pooled", "clock_rhat_split"),
        "tMRCA": ("tMRCA", "date; hpd95_lo = earlier date", "tmrca_median", "tmrca_hpd_early", "tmrca_hpd_late", None, "root_height_ess_pooled", "root_height_rhat_split"),
        "growth": ("growth rate (exponential model)", "per day", "growth_rate_per_day_median", "growth_rate_per_day_hpd_lo", "growth_rate_per_day_hpd_hi", None, "growth_ess_pooled", "growth_rhat_split"),
        "doubling": ("doubling time (exponential model)", "days", "doubling_time_d_median", "doubling_time_d_hpd_lo", "doubling_time_d_hpd_hi", None, "growth_ess_pooled", "growth_rhat_split"),
        "R_exp": ("R of the whole period (exponential model)", "", "R_median", "R_hpd_lo", "R_hpd_hi", "P_R_gt_1", "growth_ess_pooled", "growth_rhat_split")}
PERIODS = {"W1": ("2026-03-01", "2026-05-15"), "W2": ("2026-05-15", "2026-06-15"), "W3b": ("2026-06-15", "2026-08-15"), "X1": ("2026-06-15", "2026-08-01")}
DIAG_PAR = {"rate": "clock.rate", "tMRCA": "root height (tMRCA)", "W1": "growth rate W1", "W2": "growth rate W2", "W3": "growth rate W3", "W3b": "growth rate W3b",
            "growth": "exponential.growthRate", "doubling": "exponential.growthRate", "R_exp": "exponential.growthRate"}


DELIVERED_ROW = "cell of a delivered table (text of the cell)"
COMPUTED_ROW = ("computed on 26 Sep 2026 from the delivered posterior sample with the functions of the delivered tables; NOT a value of a delivered table "
                "(the delivered table gives the exploratory periods for the primary analysis only)")


def delphy_rows(DV, samples_path, samples_version):
    rows, checks = [], []
    done = set()
    for m, cp in COUNTERPART.items():
        if (cp["table"], cp["key"]) in done:
            continue
        done.add((cp["table"], cp["key"]))
        fn, vid, rno, r = DV.row(cp["table"], "analysis", cp["key"])
        used_by = "; ".join(SHORT[k] for k, v in COUNTERPART.items() if (v["table"], v["key"]) == (cp["table"], cp["key"]))
        mp = SENS if cp["table"] == "sensitivity" else EXPO
        W = None
        if cp["table"] == "sensitivity":
            W = {w["key"]: w for w in windows_of(B.series_from_delphy(cp["samples"], samples_path))}
            # check of the functions: the periods of this analysis that ARE in the Delphy table, computed from the posterior sample
            n_c, worst = 0, 0.0
            for key in ("W1", "W2", "W3", "W3b", "X1", "X2"):
                name, unit, cm, cl, ch, cpv, ce, cr = mp[key]
                if r[cm] == "" or W[key]["rs"] is None:
                    continue
                rs = W[key]["rs"]
                for col, mine in ((cm, rs["R_median"]), (cl, rs["R_hpd_lo"]), (ch, rs["R_hpd_hi"]), (cpv, rs["P_R_gt_1"]), (ce, rs["ess_pooled"]), (cr, rs["rhat_split"])):
                    n_c += 1
                    worst = max(worst, abs(float(r[col]) - mine) / max(abs(float(r[col])), 1e-300))
            assert worst < 1e-9, (cp["key"], worst)
            checks.append(dict(delphy_analysis=cp["samples"], cells_of_the_delivered_table_recomputed=n_c, largest_relative_difference=worst))
        for key, (name, unit, cm, cl, ch, cpv, ce, cr) in mp.items():
            out = dict(delphy_analysis=cp["samples"], delphy_analysis_name=cp["name"], counterpart_of_beast_analysis=used_by, engine="Delphy 1.4.1", genome_set=r["dataset"], n_genomes=r["n_genomes"],
                       last_collection_date=r["last_tip"], coalescent_prior_cells=r["coalescent_prior_cells"], chains_run=r["chains_run"], chains_retained=r["chains_retained"],
                       quantity=name, quantity_key=key, unit=unit, status=(EXPL if key.startswith("X") else "fixed before the estimates" if key.startswith("W") else ""),
                       median=r[cm], hpd95_lo=r[cl], hpd95_hi=r[ch], P_R_gt_1=(r[cpv] if cpv else ""), ess_pooled=r[ce], rhat_split=r[cr],
                       source=DELIVERED_ROW, source_table=fn, source_version_id=vid, source_row=rno, source_row_key=f"analysis = {cp['key']}",
                       source_columns=" ".join(x for x in (cm, cl, ch, cpv, ce, cr) if x), note="")
            if cp["table"] == "sensitivity":
                out.update(skygrid_precision=r["skygrid_precision"], grid=r["grid"], chains_discarded_root_date_rule=r["discarded_root_date_rule"],
                           chains_discarded_density_criterion_only=r["discarded_density_criterion_only"], chains_discarded_incomplete=r["discarded_incomplete"])
                if key == "W3":
                    out.update(period_text_of_the_table=r["R_W3_window"])
                elif key == "X1":
                    out.update(period_text_of_the_table=r["exploratory_window_15Jun_to_1Aug_definition"])
                elif key == "X2":
                    out.update(period_text_of_the_table=r["exploratory_window_1Aug_to_midpoint_of_most_recent_interval_definition"])
                if key in W:
                    out.update(period_start=W[key]["start"], period_end=W[key]["end"], period_days=W[key]["days"], period_end_is_midpoint_of_most_recent_interval=W[key]["end_is_midpoint"])
                if key in W and r[cm] == "":
                    w = W[key]
                    if w["rs"] is None:
                        out.update(source="", source_table="", source_version_id="", source_row="", source_columns="", note="cell empty in the delivered table; " + w["note"])
                    else:
                        rs = w["rs"]
                        out.update(median=rs["R_median"], hpd95_lo=rs["R_hpd_lo"], hpd95_hi=rs["R_hpd_hi"], P_R_gt_1=rs["P_R_gt_1"], ess_pooled=rs["ess_pooled"], rhat_split=rs["rhat_split"],
                                   source=COMPUTED_ROW, source_table=os.path.basename(samples_path), source_version_id=samples_version, source_row="",
                                   source_row_key=f"analysis = {cp['samples']}", source_columns="lnNe_01 ... ; seed; state",
                                   note=("cell empty in the delivered table. " + w["note"]).strip())
            # number of states used: pooled row of the Delphy table of chain diagnostics
            # (periods without a row of their own in that table: row of the evolutionary rate, same states)
            par = DIAG_PAR.get(key, "clock.rate")
            if not DV.has("chain_diagnostics", "analysis_label", cp["samples"], extra=dict(scope="pooled retained chains", parameter=par)):
                par = "clock.rate"
            f2, v2, rno2, r2 = DV.row("chain_diagnostics", "analysis_label", cp["samples"], extra=dict(scope="pooled retained chains", parameter=par))
            out.update(n_states_used=r2["n_samples"], n_states_source=f"{f2}, version {v2}, row {rno2} (analysis_label = {cp['samples']}, parameter = {par}, pooled)")
            rows.append(out)
    return pd.DataFrame(rows), checks


def intervals_from_series(S, set_name):
    """interval rows of a Delphy analysis from its posterior sample: the statements of analysis.skygrid_tables for the intervals"""
    g = S.grid
    dd = L.tip_dates_from_fasta(f"inputs/aln/{set_name}_{SFX}.fasta")
    tips_dec = np.array([L.decimal_year(t) for t in dd])
    tm_lo = S.t_last - L.hpd(np.concatenate(S.rh) * L.YEAR_DAYS)[1] / L.YEAR_DAYS           # lower (earliest) HPD bound of the tMRCA
    rows = []
    for i in range(1, g.K + 1):
        older, newer = g.edges(i)
        s = L.summarise([v[:, i - 1] for v in S.lnNe], f"lnNe_{i}")
        Ne = np.exp(np.concatenate([v[:, i - 1] for v in S.lnNe]))
        nlo, nhi = L.hpd(np.log(Ne))
        n_g = int(((tips_dec > older) & (tips_dec <= newer + 1e-12)).sum()) if i > 1 else int((tips_dec > older).sum())
        dom = ("prior-dominated (interval older than the lower 95 % HPD bound of the tMRCA)" if newer <= tm_lo else
               "partly older than the lower 95 % HPD bound of the tMRCA" if older < tm_lo else "")
        rows.append(dict(coalescent_prior_cells="8000", setting="8,000 cells", chains_run=len(S.lnNe), chains_retained=len(S.lnNe), interval=i,
                         older_edge=str(L.dec2date(older)) if i < g.K else "open", newer_edge=str(L.dec2date(newer)), mid_point=str(L.dec2date(g.mids[i - 1])) if i < g.K else "",
                         width_days=g.D * L.YEAR_DAYS if i < g.K else np.nan, lnNe_median=s["median"], lnNe_hpd_lo=s["hpd_lo"], lnNe_hpd_hi=s["hpd_hi"], lnNe_sd=s["sd"],
                         Ne_tau_years_median=float(np.median(Ne)), Ne_tau_years_hpd_lo=float(np.exp(nlo)), Ne_tau_years_hpd_hi=float(np.exp(nhi)),
                         Ne_tau_days_median=float(np.median(Ne)) * L.YEAR_DAYS, n_genomes_collected=n_g, status=dom, ess_pooled=s["ess_pooled"],
                         ess_sum_of_chains=s["ess_sum_of_chains"], ess_min_chain=s["ess_min_chain"], rhat_split=s["rhat_split"], rhat_rank=s["rhat_rank"]))
    return pd.DataFrame(rows)


CHECK_COLS = ["older_edge", "newer_edge", "mid_point", "width_days", "lnNe_median", "lnNe_hpd_lo", "lnNe_hpd_hi", "lnNe_sd", "Ne_tau_years_median", "Ne_tau_years_hpd_lo", "Ne_tau_years_hpd_hi",
              "Ne_tau_days_median", "n_genomes_collected", "status", "ess_pooled", "ess_sum_of_chains", "ess_min_chain", "rhat_split", "rhat_rank"]


def delphy_intervals(DV, samples_path, samples_version):
    out = []
    # primary analysis: rows of the Delphy table
    fn, d = DV.table("skygrid_ne_intervals")
    x = d[d["coalescent_prior_cells"] == "8000"].copy()
    assert len(x) == 20
    # check of the function on the primary analysis: the rows computed from the Delphy posterior sample against the rows of the Delphy table
    mine = intervals_from_series(B.series_from_delphy("Skygrid | primary | 8000 cells", samples_path), "primary")
    n_cells, n_equal_text, worst = 0, 0, 0.0
    for (_, rd), (_, rm) in zip(x.iterrows(), mine.iterrows()):
        assert int(rd["interval"]) == int(rm["interval"])
        for c_ in CHECK_COLS:
            n_cells += 1
            a_, b_ = rd[c_], rm[c_]
            try:
                fa, fb = float(a_), float(b_)
                if fa != fa and fb != fb:
                    n_equal_text += 1
                    continue
                rel = abs(fa - fb) / max(abs(fa), 1e-300)
                worst = max(worst, rel)
                n_equal_text += int(fa == fb)
            except (TypeError, ValueError):
                b_ = "" if (isinstance(b_, float) and b_ != b_) else str(b_)
                assert str(a_) == b_, (c_, a_, b_)
                n_equal_text += 1
    assert worst < 1e-9, worst
    check = dict(delphy_analysis="Skygrid | primary | 8000 cells", rows=20, cells_compared=n_cells, cells_equal=n_equal_text, largest_relative_difference=worst)
    Sp = B.series_from_delphy("Skygrid | primary | 8000 cells", samples_path)
    x = add_grid_columns(x, Sp.grid, int(sum(len(v) for v in Sp.lnNe)))
    x.insert(0, "delphy_analysis", "Skygrid | primary | 8000 cells"); x.insert(1, "source", "row of a delivered table (text of the cells); the columns from interval_number_1_is_the_most_recent on are added here")
    x.insert(2, "source_table", fn); x.insert(3, "source_version_id", DV.V[fn]["version_id"]); x.insert(4, "source_row", [i + 1 for i in x.index])
    x["check"] = ""
    out.append(x)
    for key, lab, st in (("gridWeekly32 | 8000 cells", "Skygrid | gridWeekly32 | 8000 cells", "primary"), ("inrbRetained | 8000 cells", "Skygrid | inrbRetained | 8000 cells", "inrbRetained")):
        Sy = B.series_from_delphy(lab, samples_path)
        y = add_grid_columns(intervals_from_series(Sy, st), Sy.grid, int(sum(len(v) for v in Sy.lnNe)))
        y.insert(0, "delphy_analysis", lab)
        y.insert(1, "source", "computed on 26 Sep 2026 from the delivered posterior sample with the statements of the delivered table; NOT a row of a delivered table "
                              "(the delivered table of intervals holds the primary genome set on the primary grid only)")
        y.insert(2, "source_table", os.path.basename(samples_path)); y.insert(3, "source_version_id", samples_version); y.insert(4, "source_row", "")
        y["check"] = (f"the same statements applied to the posterior sample of the primary analysis give the 20 delivered rows: {n_cells} cells compared, {n_equal_text} equal, "
                      f"largest relative difference {worst:.1e}")
        out.append(y)
    return pd.concat(out, ignore_index=True), [check]


# ------------------------------------------------------------------------------------------------ reading rule
def num(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return np.nan


SEEN_COORD_TMRCA = "no: criterion written after interim tMRCA values of analysis C (incomplete chains) had been seen"
SEEN_COORD_OTHER = "no unmasked value found in the record"
SEEN_COORD_NONE = "not stated for this quantity (quantity outside rule 2)"
RECORD = "beast_values_seen_record_20260926.md"
SEEN_TRACK = {
    "tMRCA": f"yes: root height of incomplete chains of the hedge (item b of {RECORD}, 25 Sep 02:52 UTC) and tMRCA of incomplete chains of the second run (item c3, 25 Sep 18:45 UTC); root height of the first run (item a1)",
    "rate": f"yes: evolutionary rate of incomplete chains of the first run (item a1 of {RECORD}, 25 Sep 02:26 UTC) and of the hedge (item b, 25 Sep 02:52 UTC)",
    "W": f"yes: R of this period from incomplete chains of the hedge (item b of {RECORD}, 25 Sep 02:52 UTC)",
    "X": f"no unmasked value of this period was displayed to the track according to its record ({RECORD}); two archived test texts hold such values (section 3 of the record); R of W1, W2, W3, W3b of incomplete chains had been displayed (item b)",
    "D": f"yes: growth rate of incomplete chains of the first run (item a1 of {RECORD}) and growth rate and R of the hedge (item b, 25 Sep 02:52 UTC)",
}
PAIR_LABEL = {("A_skyfix_primary", "Delphy primary analysis"): "same model; differs in the program",
              ("B_skyest_primary", "Delphy primary analysis"): "differs in smoothing (estimated)",
              ("C_inrblike_primary", "Delphy primary analysis"): "differs in grid, smoothing and substitution model; last period ends on another date",
              ("C_inrblike_primary", "Delphy counterpart"): "differs in substitution model and smoothing",
              ("C_inrblike_inrbRetained", "Delphy primary analysis"): "differs in genome set, grid, smoothing and substitution model; last period ends on another date",
              ("C_inrblike_inrbRetained", "Delphy counterpart"): "differs in grid, smoothing and substitution model; last period ends on another date",
              ("D_exp_primary", "Delphy primary analysis"): "differs in the tree prior (exponential growth)",
              ("D_exp_primary", "Delphy counterpart"): "same model; differs in the program"}


DESIGN_PAIR = {("A_skyfix_primary", "Delphy primary analysis"): "A - primary analysis (same model)",
               ("B_skyest_primary", "Delphy primary analysis"): "B - primary analysis (differs in the program and in the smoothing, which B estimates)",
               ("C_inrblike_primary", "Delphy counterpart"): "C on 555 genomes - analysis with the weekly grid of 32 parameters (differs in the program, in the substitution model and in the smoothing, which C estimates)",
               ("C_inrblike_inrbRetained", "Delphy counterpart"): "C on the 515 genomes retained in the previous post - analysis of these genomes (differs in the program, the grid, the substitution model and the smoothing; "
                                                                  "the last period ends on another date)",
               ("D_exp_primary", "Delphy counterpart"): "D - exponential model (same model)"}
DESIGN_QUANT = {"skygrid": ("W1", "W2", "W3", "rate", "tMRCA"), "exp": ("growth", "doubling", "R_exp", "rate", "tMRCA")}
DESIGN_TABLE_ONLY = ("X1", "X2")
PRIORS_NOTE = ("the priors of the parameters other than the population sizes are those of each program (evolutionary rate: CTMC scale prior in BEAST X, uniform in Delphy; kappa and base frequencies: "
               "log-normal(1, 1.25) and uniform in BEAST X, defaults of Delphy)")
PRIORS_NOTE_D = ("the priors are those of each program (growth rate: Laplace, location 0, scale 100 per year in BEAST X; Laplace, location 0.001, scale 30.701135 per year in Delphy, its default; population size at the last "
                 "collection date: 1/x in both programs (default of Delphy); Delphy replaces a population size N(t) below 0.0027397 yr (1 day; default of --v0-pop-min-pop) by this value in the coalescent density, BEAST X has "
                 "no such floor; evolutionary rate, kappa, base frequencies as for A)")


def reading_rows(E, DLP):
    """'agree' if the median of each analysis lies inside the 95 % HPD interval of the other"""
    rows = []
    E = E.set_index(["analysis", "quantity_key"]).sort_index(); DL = DLP.set_index(["delphy_analysis", "quantity_key"]).sort_index()
    PRIM = "Skygrid | primary | 8000 cells"
    for m in COUNTERPART:
        pairs = [(PRIM, "Delphy primary analysis")]
        if COUNTERPART[m]["samples"] != PRIM:
            pairs.append((COUNTERPART[m]["samples"], "Delphy counterpart"))
        for ref, what in pairs:
            if m == "D_exp_primary":
                keys = ("rate", "tMRCA") if what == "Delphy primary analysis" else ("rate", "tMRCA", "growth", "doubling", "R_exp")
            else:
                keys = ("W1", "W2", "W3", "rate", "tMRCA") + (("X1", "X2", "W3b") if True else ())
            same_genomes = E.loc[(m, "rate"), "genome_set"] == "primary"
            if what == "Delphy primary analysis" and (m, "W3p") in E.index:
                keys = keys + ("W3p", "X2p")                        # BEAST X period with the end date of the Delphy primary analysis against W3 / second exploratory period of Delphy
            for k in keys:
                kd = {"W3p": "W3", "X2p": "X2"}.get(k, k)
                if (ref, kd) not in DL.index or (m, k) not in E.index:
                    continue
                d, b = DL.loc[(ref, kd)], E.loc[(m, k)]
                if str(d["median"]) in ("", "nan") or str(b["median"]) in ("", "nan"):
                    continue
                if k == "tMRCA":
                    f = lambda s_: float((pd.Timestamp(s_) - pd.Timestamp("2026-01-01")).days)
                    dm, dl, dh = f(d["median"]), f(d["hpd95_lo"]), f(d["hpd95_hi"]); bm, bl, bh = f(b["median"]), f(b["hpd95_lo"]), f(b["hpd95_hi"])
                    unit = "difference and width in days"
                else:
                    dm, dl, dh = num(d["median"]), num(d["hpd95_lo"]), num(d["hpd95_hi"]); bm, bl, bh = num(b["median"]), num(b["hpd95_lo"]), num(b["hpd95_hi"])
                    unit = b["unit"]
                if any(v != v for v in (dm, dl, dh, bm, bl, bh)):
                    continue
                in1, in2 = bool(dl <= bm <= dh), bool(bl <= dm <= bh)
                applies = bool(what == "Delphy primary analysis" and same_genomes and k in ("W1", "W2", "W3", "rate", "tMRCA"))
                why = []
                if what != "Delphy primary analysis":
                    why.append("the Delphy analysis is not the primary analysis")
                if not same_genomes:
                    why.append("other genome set")
                if k in ("W3p", "X2p"):
                    why.append("additional row: period of the BEAST X analysis with the end date of the Delphy primary analysis")
                elif k not in ("W1", "W2", "W3", "rate", "tMRCA"):
                    why.append("quantity not named in rule 2")
                # how rule and design are applied in the post (rule file, version 9, entries of 26 Sep 02:31 and 02:39 UTC): reading within the pairs of C.1; the column of Table S1 is filled for R of the three periods,
                # the evolutionary rate and, in D, the growth rate alone (doubling time and R are transformations of it) - not for exploratory periods, not for the tMRCA, not where the periods of the two programs end on
                # different dates
                in_pair = (m, what) in DESIGN_PAIR
                is_period = k in ("W1", "W2", "W3", "W3b", "X1", "X2", "W3p", "X2p")
                same_end = (str(b.get("period_end", "")) == str(d.get("period_end", ""))) if is_period else None
                carries = ("growth", "rate") if m == "D_exp_primary" else ("W1", "W2", "W3", "rate")
                no_fill = []
                if not in_pair:
                    no_fill.append("not a pair of section C.1")
                if k in ("X1", "X2", "X2p"):
                    no_fill.append("exploratory period (shown, not read)")
                if k == "tMRCA":
                    no_fill.append("tMRCA: the criterion was written after interim values had been seen; estimates and difference in days are given without a reading")
                if k in ("W3b", "W3p"):
                    no_fill.append("additional period of the track, not a quantity of the design")
                if is_period and same_end is False:
                    no_fill.append("the periods of the two programs end on different dates")
                if k in ("doubling", "R_exp"):
                    no_fill.append("transformation of the growth rate: the reading is given for the growth rate alone (entry of 02:39 UTC)")
                if k not in carries and not no_fill:
                    no_fill.append("quantity does not carry the reading in the post")
                filled = bool(in_pair and k in carries and not no_fill)
                rows.append(dict(beast_analysis=m, beast_analysis_short=SHORT[m], compared_with=what, delphy_analysis=ref, quantity=b["quantity"], quantity_key=k, unit=unit,
                                 reading_within_a_pair_of_C1=in_pair, periods_of_the_two_programs_end_on_the_same_date=("" if same_end is None else same_end),
                                 reading_given_in_the_post_by_the_entries_of_0231_and_0239=filled, why_the_column_is_not_filled="; ".join(no_fill),
                                 delphy_median=d["median"], delphy_hpd95_lo=d["hpd95_lo"], delphy_hpd95_hi=d["hpd95_hi"], delphy_period_end=d.get("period_end", ""),
                                 beast_median=b["median"], beast_hpd95_lo=b["hpd95_lo"], beast_hpd95_hi=b["hpd95_hi"], beast_period_end=b.get("period_end", ""),
                                 beast_median_inside_delphy_hpd=in1, delphy_median_inside_beast_hpd=in2, reading=("agree" if in1 and in2 else "differ"),
                                 difference_beast_minus_delphy=bm - dm, difference_as_share_of_the_width_of_the_delphy_hpd=((bm - dm) / (dh - dl) if dh > dl else np.nan),
                                 rule_3_of_the_reading_rule_applies=applies,
                                 rule_fixed_before_values_seen=(SEEN_COORD_TMRCA if k == "tMRCA" else SEEN_COORD_OTHER if (k == "rate" or k[0] in "WX" or k == "R_exp") else SEEN_COORD_NONE),
                                 displayed_to_the_track_before_the_rule_was_written=(SEEN_TRACK["tMRCA"] if k == "tMRCA" else SEEN_TRACK["rate"] if k == "rate" else
                                                                                     SEEN_TRACK["W"] if k in ("W1", "W2", "W3", "W3b") else SEEN_TRACK["X"] if k in ("X1", "X2", "W3p", "X2p") else SEEN_TRACK["D"]),
                                 pair_label=PAIR_LABEL[(m, what)],
                                 pair_of_the_design_C1=((m, what) in DESIGN_PAIR),
                                 pair_of_the_design_C1_text=DESIGN_PAIR.get((m, what), ""),
                                 quantity_of_the_design_C2=bool((m, what) in DESIGN_PAIR and (k in DESIGN_QUANT["exp" if m == "D_exp_primary" else "skygrid"] or k in DESIGN_TABLE_ONLY)),
                                 place_in_the_design=("" if (m, what) not in DESIGN_PAIR else "figure and table" if k in DESIGN_QUANT["exp" if m == "D_exp_primary" else "skygrid"] else
                                                      "table only (exploratory period)" if k in DESIGN_TABLE_ONLY else "not a quantity of the design"),
                                 remark_of_the_track_on_the_label_of_the_pair=("" if (m, what) not in DESIGN_PAIR else PRIORS_NOTE_D if m == "D_exp_primary" else PRIORS_NOTE),
                                 remark=("" if applies else "outside rule 3 (" + "; ".join(why) + "): same computation, shown for information"),
                                 differs_apart_from_the_program=(DIFFERS_FROM_PRIMARY[m] if what == "Delphy primary analysis" else COUNTERPART[m]["differs"]),
                                 reading_rule_version=RULE_VERSION, reading_rule_version_number=RULE_NUMBER, delphy_value_source=d["source"], delphy_source_table=d["source_table"], delphy_source_version_id=d["source_version_id"],
                                 delphy_source_row=d["source_row"]))
    return pd.DataFrame(rows)


# ------------------------------------------------------------------------------------------------ main
XML_DIR = os.environ.get("BEAST_XML_DIR", "work/beast/xml_regen")
PROGRAM_LINES = [("program_version_line", r"^\s*BEAST v[\w.\-]+"), ("beagle_version_line", r"Using BEAGLE library"), ("beagle_resource_line", r"Using BEAGLE resource"),
                 ("beagle_threads_line", r"Using \d+ threads? for CPU"), ("likelihood_threads_line", r"Likelihood computation is using a pool of"), ("development_parsers_line", r"development parsers"),
                 ("seed_line", r"^Random number seed:"), ("command_line_as_printed", r"^# -threads ")]


def program_output(indir, AN):
    """what the program printed (standard output of every chain), the settings of the chain in the XML file, and the states of the log: one row per chain. Lines are quoted with their number."""
    recorded = {}
    rp = os.path.join(indir, "xml_sha256.txt")
    if os.path.exists(rp):
        for l in open(rp):
            if l.strip():
                h, n = l.split()
                recorded[os.path.basename(n)] = h
    rows = []
    for m, an in AN.items():
        xml = os.path.join(XML_DIR, f"beast{m}_{SFX}.xml")
        xs = sha256(xml) if os.path.exists(xml) else ""
        xt = open(xml).read() if os.path.exists(xml) else ""
        g = lambda rx: (int(re.search(rx, xt).group(1)) if re.search(rx, xt) else None)
        requested, log_every, tree_every = g(r'<mcmc [^>]*chainLength="(\d+)"'), g(r'<log id="fileLog" logEvery="(\d+)"'), g(r'<logTree [^>]*logEvery="(\d+)"')
        for c in an.chains:
            d = os.path.dirname(c.log)
            so = os.path.join(d, "stdout.txt") if d else None
            r = dict(analysis=m, analysis_short=SHORT[m], seed=int(c.seed), directory_of_the_chain=os.path.basename(d) if d else "", file_read=(os.path.basename(d) + "/stdout.txt") if d else "")
            lines = open(so, errors="replace").read().split("\n") if so and os.path.exists(so) else []
            r["lines_of_the_file"] = len(lines)
            for key, rx in PROGRAM_LINES:
                hit = [(i + 1, l) for i, l in enumerate(lines) if re.search(rx, l)]
                r[key] = hit[0][1].strip() if hit else ""
                r[key + "_number"] = hit[0][0] if hit else ""
                r[key + "_times_printed"] = len(hit)
            mv = re.search(r"BEAST v([\w.\-]+)", r["program_version_line"]); r["program_version_printed"] = mv.group(1).rstrip(",") if mv else ""
            mb = re.search(r"BEAGLE library v([\w.\-]+(?: \([^)]*\))?)", r["beagle_version_line"]); r["beagle_version_printed"] = mb.group(1) if mb else ""
            mt = re.search(r"Using (\d+) threads? for CPU", r["beagle_threads_line"]); r["beagle_threads_per_instance_printed"] = int(mt.group(1)) if mt else ""
            ml = re.search(r"pool of (\d+) threads", r["likelihood_threads_line"]); r["likelihood_thread_pool_printed"] = int(ml.group(1)) if ml else ""
            r["beagle_instances_printed"] = r["beagle_resource_line_times_printed"]
            cm = os.path.join(d, "command.txt") if d else None
            r["command_line_of_the_run"] = open(cm).read().strip() if cm and os.path.exists(cm) else ""
            mo = re.search(r"-threads (\d+)", r["command_line_of_the_run"]); r["threads_option_of_the_command_line"] = int(mo.group(1)) if mo else ""
            r["xml_file"] = os.path.basename(xml); r["xml_sha256_of_the_copy_read"] = xs; r["xml_sha256_recorded_by_the_run"] = recorded.get(os.path.basename(xml), "")
            r["xml_copy_read_is_the_file_of_the_run"] = bool(xs) and xs == recorded.get(os.path.basename(xml), "")
            r["states_requested_in_the_xml"] = requested; r["logging_interval_of_the_parameters_in_the_xml"] = log_every; r["logging_interval_of_the_trees_in_the_xml"] = tree_every
            st = c.full["state"].to_numpy()
            r["states_reached_last_state_of_the_log"] = int(st[-1]); r["logged_states_in_the_log"] = int(len(st))
            r["logging_interval_observed_in_the_log"] = "; ".join(str(int(x)) for x in sorted(set(np.diff(st).tolist())))
            r["chain_complete"] = bool(st[-1] >= B.CHAIN_LENGTH)
            r["burn_in_fraction"] = B.BURNIN
            r["burn_in_logged_states_not_used"] = int(c.nb)
            r["burn_in_last_state_not_used"] = int(st[c.nb - 1]) if c.nb > 0 else ""
            r["burn_in_states"] = int(st[c.nb]) if c.nb < len(st) else ""          # states of the chain before the first state used
            kept = c in an.kept
            r["chain_retained"] = bool(kept)
            r["first_state_used"] = int(c.post["state"].iloc[0]) if kept and len(c.post) else ""
            r["last_state_used"] = int(c.post["state"].iloc[-1]) if kept and len(c.post) else ""
            r["logged_states_used_after_burn_in"] = int(len(c.post)) if kept else 0
            rows.append(r)
    return pd.DataFrame(rows)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--in", dest="indir", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--tables", default="work/report_inputs/tables")
    ap.add_argument("--versions", required=True)
    ap.add_argument("--mcc", default=None)
    ap.add_argument("--test", action="store_true")
    a = ap.parse_args()
    os.makedirs(a.out, exist_ok=True)
    pre = "TEST_" if a.test else ""
    V = json.load(open(a.versions))
    DV = Delivered(a.tables, V)
    AN = {m: B.load_model(m, a.indir, a.test) for m in B.MODELS}
    assert sum(len(an.chains) for an in AN.values()) == 40, {m: len(an.chains) for m, an in AN.items()}
    facts = dict(test=a.test, analyses={})
    est, ivs = [], []
    for m, an in AN.items():
        P = an.P
        st = [c.info["last_state"] for c in an.chains]
        facts["analyses"][m] = dict(chains_run=len(an.chains), chains_retained=len(an.kept), chains_complete=int(sum(s_ >= B.CHAIN_LENGTH for s_ in st)), states_min=int(min(st)), states_max=int(max(st)),
                                    left_out=left_out_text(an), root_date_rule=int(sum(c.info["collapsed"] for c in an.chains)), latest_root_date=max(c.info["latest_root_date"] for c in an.chains),
                                    n_states_used=int(sum(len(c.post) for c in an.kept)), first_state_used=int(min(c.post["state"].iloc[0] for c in an.kept)) if an.kept else None,
                                    last_state_used=int(max(c.post["state"].iloc[-1] for c in an.kept)) if an.kept else None, log_every=int(P["log_every"]), tree_every=int(P["tree_every"]),
                                    n_log_columns=int(an.chains[0].full.shape[1]))
        if not an.kept:
            est.append(dict(analysis=m, analysis_short=SHORT[m], engine="BEAST X 10.5.0", model=an.model, genome_set=an.dataset, chains_run=len(an.chains), chains_retained=0,
                            criterion_chain_left_out=left_out_text(an), quantity="no estimate: no chain retained"))
            continue
        if an.grid is not None:
            est += skygrid_rows(an, P)
            mc = None
            if a.mcc and m == "A_skyfix_primary" and os.path.exists(a.mcc):
                labels_fa = [l[1:].strip() for l in L.open_any(an.fasta) if l.startswith(">")]
                dec = {l: L.decimal_year(l.split("|")[-1]) for l in labels_fa}
                labels, trees = L.read_nexus_trees(a.mcc)
                tips, ints = L.node_times(trees[0][1], labels, dec, tol=5e-3)
                mc, _ = L.skygrid_tree_terms(tips, ints, an.grid)
            try:
                tc = tree_counts_of(an)
            except Exception as e:                       # the table is written without these columns rather than not at all
                tc = None
                facts["analyses"][m]["tree_counts_error"] = repr(e)[:300]
            facts["analyses"][m]["trees_in_coalescent_counts"] = (int(tc["n_trees"]) if tc else 0)
            iv = intervals_frame(an, mcc_counts=mc, tree_counts=tc)
            iv["coalescent_events_note"] = (f"posterior columns: {tc['n_trees']} sampled trees after burn-in (at most 40 per chain, evenly spaced)" if tc else "sampled trees not evaluated") + \
                                           ("; MCC column: summary tree of model A" if mc is not None else "; MCC column empty: a summary tree is built for model A only" if m != "A_skyfix_primary" else
                                            "; MCC column empty: summary tree not given")
            ivs.append(iv)
        else:
            est += exponential_rows(an, P)
        pf = posterior_frame(an)
        assert len(pf) == facts["analyses"][m]["n_states_used"]
        pth = f"{a.out}/{pre}beast_posterior_{m}_{SFX}.csv.gz"
        pf.to_csv(pth, index=False, compression=dict(method="gzip", mtime=0))
        facts["analyses"][m].update(posterior_file=os.path.basename(pth), posterior_rows=int(len(pf)), posterior_columns=int(pf.shape[1]))
        print(m, "rows of the posterior sample", len(pf), "columns", pf.shape[1], flush=True)
    E = pd.DataFrame(est)
    CRIT_ON = {"rate": "the evolutionary rate (log column clock.rate)", "tMRCA": "the height of the root (log column treeModel.rootHeight)", "precision": "the precision (log column skygrid.precision)",
               "growth": "the growth rate (log column exponential.growthRate)", "doubling": "the growth rate (log column exponential.growthRate), of which the doubling time is a transformation",
               "R_exp": "the growth rate (log column exponential.growthRate), of which R is a transformation"}
    if "quantity_key" in E.columns:
        E["criteria_of_convergence_computed_on"] = [CRIT_ON.get(k, "the GROWTH RATE of the period (difference of ln Ne tau between the two ends of the period, interpolated between the mid-points of the intervals, "
                                                                  "divided by its length), of which R is a transformation") if isinstance(k, str) else "" for k in E["quantity_key"]]
        E["definition"] = [(str(d_) + "; pooled ESS and split R-hat computed on " + c_ + " (as in the delivered Delphy tables)") if isinstance(d_, str) and d_ and "pooled ESS and split R-hat computed on" not in d_ else d_
                           for d_, c_ in zip(E["definition"], E["criteria_of_convergence_computed_on"])]
        E["pooled_ess_function"] = "rt_lib.ess_multichain: multi-chain estimator over the retained chains, chains not split (source quoted in the note and in pooled_ess_and_rhat_source_20260926.md)"
        E["split_rhat_function"] = "rt_lib.rhat_split: each retained chain split into two halves"
    first = ["analysis", "analysis_short", "engine", "model", "genome_set", "n_genomes", "first_collection_date", "last_collection_date", "quantity", "quantity_key", "status", "unit",
             "period_start", "period_end", "period_days", "period_end_is_midpoint_of_most_recent_interval", "median", "hpd95_lo", "hpd95_hi", "P_R_gt_1",
             "growth_rate_per_day_median", "growth_rate_per_day_hpd95_lo", "growth_rate_per_day_hpd95_hi", "root_height_days_median", "root_height_days_hpd95_lo", "root_height_days_hpd95_hi",
             "chains_run", "chains_retained", "chains_left_out", "criterion_chain_left_out", "ess_pooled", "rhat_split", "criterion_ess200_rhat105_met", "n_states_used"]
    E = E[[c for c in first if c in E.columns] + [c for c in E.columns if c not in first]]
    E.to_csv(f"{a.out}/{pre}beast_estimates_{SFX}.csv", index=False)
    IV = pd.concat(ivs, ignore_index=True) if ivs else pd.DataFrame()
    facts["hpd95_method"] = HPD_METHOD
    facts["hpd95_same_as_floor_convention"] = (sorted(set(bool(x) for x in IV["hpd95_same_run_as_with_k_equal_floor_of_0.95_n_steps"])) if len(IV) else [])
    IV.to_csv(f"{a.out}/{pre}beast_ne_intervals_{SFX}.csv", index=False)
    sp = os.path.join(a.tables, f"posterior_samples_delphy_{SFX}.csv.gz")
    assert sha256(sp) == V[os.path.basename(sp)]["sha256"], "posterior sample of Delphy: the local copy is not the delivered version"
    DLP, dchecks = delphy_rows(DV, sp, V[os.path.basename(sp)]["version_id"])
    DLP.to_csv(f"{a.out}/{pre}delphy_counterparts_estimates_{SFX}.csv", index=False)
    DIV, checks = delphy_intervals(DV, sp, V[os.path.basename(sp)]["version_id"])
    DIV.to_csv(f"{a.out}/{pre}delphy_counterparts_ne_intervals_{SFX}.csv", index=False)
    R = reading_rows(E, DLP)
    R.to_csv(f"{a.out}/{pre}beast_vs_delphy_reading_{SFX}.csv", index=False)
    Rr = R[R["rule_3_of_the_reading_rule_applies"]]
    Rd = R[R["pair_of_the_design_C1"] & R["quantity_of_the_design_C2"]]
    Rp = R[R["reading_given_in_the_post_by_the_entries_of_0231_and_0239"]]
    facts.update(estimates_rows=int(len(E)), intervals_rows=int(len(IV)), delphy_rows=int(len(DLP)), delphy_interval_rows=int(len(DIV)), delphy_interval_checks=checks, delphy_estimate_checks=dchecks,
                 delphy_rows_delivered=int((DLP['source'] == DELIVERED_ROW).sum()), delphy_rows_computed=int((DLP['source'] == COMPUTED_ROW).sum()),
                 reading_rows=int(len(R)), design_rows=int(len(Rd)), design_rows_figure=int((Rd['place_in_the_design'] == 'figure and table').sum()), design_pairs=int(Rd['pair_of_the_design_C1_text'].nunique()),
                 post_rows=int(len(Rp)), post_agree=int((Rp["reading"] == "agree").sum()), post_differ=int((Rp["reading"] == "differ").sum()),
                 post_differ_rows=[dict(analysis=r_["beast_analysis_short"], quantity_key=r_["quantity_key"]) for _, r_ in Rp[Rp["reading"] == "differ"].iterrows()],
                 post_rows_per_analysis={k_: int(v_) for k_, v_ in Rp.groupby("beast_analysis_short", sort=False).size().items()},
                 reading_rows_rule3=int(len(Rr)), reading_agree_rule3=int((Rr["reading"] == "agree").sum()), reading_differ_rule3=int((Rr["reading"] == "differ").sum()),
                 rows_not_meeting_ess_rhat=int((E["criterion_ess200_rhat105_met"].astype(str) == "False").sum()),
                 rows_with_ess_rhat=int(E["criterion_ess200_rhat105_met"].astype(str).isin(["True", "False"]).sum()))
    bad = E[E["criterion_ess200_rhat105_met"].astype(str) == "False"]
    facts["rows_not_meeting"] = [dict(analysis=r_["analysis_short"], quantity_key=r_["quantity_key"], ess_pooled=float(r_["ess_pooled"]), rhat_split=float(r_["rhat_split"])) for _, r_ in bad.iterrows()]
    facts["ess_min"] = float(pd.to_numeric(E["ess_pooled"], errors="coerce").min()); facts["rhat_max"] = float(pd.to_numeric(E["rhat_split"], errors="coerce").max())
    # source of the functions of the criteria of convergence, quoted with line numbers
    import inspect
    src_file = inspect.getsourcefile(L)
    src_lines = open(src_file).read().split("\n")
    QUOTE = []
    for fn_ in (L._autocov, L.ess_multichain, L.rhat_split, L.hpd, L._acf_pipeline, L.ess_single):
        ls_, n0_ = inspect.getsourcelines(fn_)
        QUOTE.append(dict(function=fn_.__name__, first_line=n0_, last_line=n0_ + len(ls_) - 1, source="".join(ls_).rstrip("\n")))
    facts["criteria_source"] = dict(file="pipeline_ext/rt_lib.py", sha256_of_the_file=sha256(src_file), lines_of_the_file=len(src_lines), functions=QUOTE)
    with open(f"{a.out}/{pre}pooled_ess_and_rhat_source_20260926.md", "w") as fh:
        fh.write("# Source of the functions that compute the criteria of convergence and the HPD interval (phylodynamics track)\n\n")
        fh.write(f"File `pipeline_ext/rt_lib.py` of the archive of scripts (SHA-256 of the file {sha256(src_file)}; {len(src_lines)} lines). The same functions computed the delivered Delphy tables and compute the BEAST X tables.\n\n")
        fh.write("- pooled ESS of the tables (`ess_pooled`): `ess_multichain`, with `_autocov`; chains are NOT split; chains are cut to the length of the shortest.\n- split R-hat (`rhat_split`): `rhat_split`; each chain is split into two halves.\n"
                 "- ESS of a single chain (`ess_min_chain`, `ess_sum_of_chains`): `ess_single`, with `_acf_pipeline` (the estimator of the pipeline supplied with the task).\n- 95 % HPD interval: `hpd`.\n"
                 "- For a reproduction number both criteria are computed on the growth rate of the period, of which R is a transformation (R itself is set to 0 where the transformation is undefined).\n\n")
        for q_ in QUOTE:
            fh.write(f"## `{q_['function']}` (lines {q_['first_line']} to {q_['last_line']})\n\n```python\n{q_['source']}\n```\n\n")
    PO = program_output(a.indir, AN)
    PO.to_csv(f"{a.out}/{pre}beast_program_output_{SFX}.csv", index=False)
    uniq = lambda c: sorted(set(str(x) for x in PO[c]))
    ex = PO.iloc[0]
    facts["program"] = dict(chains=int(len(PO)), file_quoted=ex["file_read"],
                            quoted={k: dict(line_number=(int(ex[k + "_number"]) if ex[k + "_number"] != "" else None), line=ex[k]) for k, _ in PROGRAM_LINES},
                            program_versions_printed=uniq("program_version_printed"), beagle_versions_printed=uniq("beagle_version_printed"), beagle_instances_printed=uniq("beagle_instances_printed"),
                            chains_with_the_same_version_line=int((PO["program_version_line"] == ex["program_version_line"]).sum()),
                            chains_with_the_same_beagle_line=int((PO["beagle_version_line"] == ex["beagle_version_line"]).sum()),
                            chains_with_development_parsers_line=int((PO["development_parsers_line"] != "").sum()),
                            xml_copies_equal_to_the_files_of_the_run=int(PO["xml_copy_read_is_the_file_of_the_run"].sum()),
                            per_analysis={m: dict(threads_option=sorted(set(int(x) for x in g_["threads_option_of_the_command_line"] if x != "")),
                                                  beagle_threads_printed=sorted(set(int(x) for x in g_["beagle_threads_per_instance_printed"] if x != "")),
                                                  likelihood_pool_printed=sorted(set(int(x) for x in g_["likelihood_thread_pool_printed"] if x != "")),
                                                  states_requested=sorted(set(int(x) for x in g_["states_requested_in_the_xml"])),
                                                  states_reached=[int(x) for x in g_["states_reached_last_state_of_the_log"]],
                                                  seeds=[int(x) for x in g_["seed"]],
                                                  logging_interval_xml=sorted(set(int(x) for x in g_["logging_interval_of_the_parameters_in_the_xml"])),
                                                  logging_interval_observed=sorted(set(g_["logging_interval_observed_in_the_log"])),
                                                  burn_in_states=sorted(set(int(x) for x in g_["burn_in_states"] if x != "")),
                                                  burn_in_logged_states_not_used=sorted(set(int(x) for x in g_["burn_in_logged_states_not_used"])),
                                                  logged_states_in_the_log=sorted(set(int(x) for x in g_["logged_states_in_the_log"])),
                                                  chains_complete=int(g_["chain_complete"].sum()), chains_retained=int(g_["chain_retained"].sum()),
                                                  logged_states_used_per_chain=[int(x) for x in g_["logged_states_used_after_burn_in"]],
                                                  logged_states_used=int(g_["logged_states_used_after_burn_in"].sum()),
                                                  first_state_used=sorted(set(int(x) for x in g_["first_state_used"] if x != "")),
                                                  last_state_used=sorted(set(int(x) for x in g_["last_state_used"] if x != "")))
                                          for m, g_ in PO.groupby("analysis", sort=False)})
    for m in facts["program"]["per_analysis"]:
        assert facts["program"]["per_analysis"][m]["logged_states_used"] == facts["analyses"][m]["n_states_used"], m
    json.dump(facts, open(f"{a.out}/{pre}beast_first_delivery_facts_{SFX}.json", "w"), indent=1)
    print("estimates", E.shape, "| intervals", IV.shape, "| Delphy rows", DLP.shape, "| Delphy intervals", DIV.shape, "| reading", R.shape, flush=True)


if __name__ == "__main__":
    main()
