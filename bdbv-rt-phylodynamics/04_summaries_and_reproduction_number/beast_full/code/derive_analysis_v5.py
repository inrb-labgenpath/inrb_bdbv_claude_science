#!/usr/bin/env python3
"""derive_analysis_v5.py - tables of ONE analysis in the common form.

Differences from the earlier versions of this script (no estimate is changed; validate_derive_v4.sh and validate_derive_v5.sh
compare, cell for cell as text, the tables of the versions without the new columns):
  - a chain of the configuration can carry its own number of states (key steps_of_this_chain): a chain that is CUT at the
    saved state from which it was continued is complete at that state; without the key the chain has the number of states
    of the analysis.  reached_last_step of such a chain reads 'not applicable: the chain is cut at its saved state
    (N of M states)'.
  - the table of the chains has the further columns saved_state_from_which_the_chain_was_continued, logged_states_read,
    logged_states_discarded_as_burn_in (the number that the class analysis.Chain discards: the whole part of 0.30 times the
    logged states that were read) and start_end_wall_time_exit_code_and_command_line_are_those_of, which says to which
    attempt the columns started_utc, ended_utc, wall_time_s, exit_code and command_line belong (for a chain that was
    continued they are those of the CONTINUATION, which ran only the states after the saved state).
  - every table (chains, estimates, intervals, trajectories, root_date_intervals) has a further column note_on_the_analysis
    with the text of the configuration (empty where the configuration gives none): it marks an analysis whose chains were
    continued.
  - the note of a row of a doubling time and the note of its period are joined; the note of the series 'Ne tau at every
    7th day' holds the rule of its dates where the configuration gives it; exp() of a value too large for a floating-point
    number gives inf instead of ending the script; the table of the shortest intervals of the date of the root is written
    for Delphy only.

Medians, 95 % intervals, pooled ESS, split R-hat, the reproduction number and the rules for chains come from the
functions of the analysis round of 24 September (work/pipeline_ext/analysis.py, rt_lib.py):
  chain, burn-in, root-date rule    analysis.Chain (first 30 % of the logged states discarded)
  density criterion                 analysis.apply_second_criterion (as derive_rt.py calls it: Delphy with the Skygrid
                                    only, up to 50 sampled trees per chain, threshold 100 log units)
  growth rate of a period           rt_lib.window_contrast(grid, start, end, 'interp') applied to ln(Ne tau)
  reproduction number               rt_lib.GenTime(mean, SD).r2R through analysis.r_summary
  median, mean, 95 % interval,      rt_lib.summarise (rt_lib.hpd, rt_lib.ess_multichain, rt_lib.rhat_split)
  pooled ESS, split R-hat
  date of the root                  analysis.tmrca_summary
  intervals of the Skygrid          analysis.skygrid_tables
Computed by THIS script, because the functions of analysis.py and rt_lib.py return no such value:
  mean of the growth rate per day   numpy mean of the pooled series of the retained chains, divided by rt_lib.YEAR_DAYS
                                    (analysis.r_summary returns no mean)
  mean date of the root             numpy mean of the pooled root heights in days, rounded to a whole day and subtracted from
                                    the most recent collection date (analysis.tmrca_summary returns no mean)
  probabilities                     mean of the indicator series over the pooled states (P(R of P2 < R of P1), P(R of P3 < 1));
                                    P(R > 1) of a period is the value P_R_gt_1 of analysis.r_summary; pooled ESS and split
                                    R-hat of an indicator series by rt_lib.ess_multichain and rt_lib.rhat_split
  shortest intervals of the date    every interval of whole days of the smallest width that holds 95 % of the pooled states
  of the root in whole days         (table root_date_intervals_*); the interval of the table of estimates is the one of
                                    analysis.tmrca_summary (rt_lib.hpd)
  series 'Ne tau at every 7th day'  the dates are set by this script; the value at a date is rt_lib.Grid.weights(date,
                                    'interp') applied to ln(Ne tau), summarised by rt_lib.summarise, and exp() of the results
  Ne tau of the tables              exp() of the median and of the bounds of ln(Ne tau) that the functions return
Everything else is passed to the functions of analysis.py and rt_lib.py, and their results are written into the columns
of the common form.

Values of this round passed where a script of the earlier round holds a constant of that round:
  analysis.W_DATES (periods W1, W2, W3, W3b)       the periods are read from the table of the periods of this round
                                                    (P1 to P6); P5 and P6 are not in that dictionary
  beast_analysis.MODELS, ALN, SFX, CHAIN_LENGTH     names of files, alignment, suffix and chain length are arguments
  beast_analysis.T_END_PRIMARY                      not used: P3 and P6 end at the mid-point of the most recent
                                                    interval of the set and grid of the analysis
  beast_analysis.MIN_WINDOW_DAYS (14 days)          not used: no smallest length is set for P3 or P6
  beast_analysis.xml_params (work/beast/xml)        number of parameters and cut-off are read from the parameter record
                                                    of the input file of this round
Rule for a period with fixed dates (it is the rule of analysis.Analysis.windows, which wrote the table of 24 September):
  the period is evaluated if its end is not later than the most recent collection date of the set and later than its
  start; where the end lies beyond the mid-point of the most recent interval, the function uses the level of the
  most recent interval from that mid-point on.  A period that ends at END (P3, P6) is evaluated if END lies after its start.

Usage: derive_analysis_v5.py --config analysis.json --out DIR
"""
import argparse
import datetime
import glob
import json
import math
import os
import sys

import numpy as np
import pandas as pd

SUFFIX = "20261002"
ESS_MIN, RHAT_MAX = 200.0, 1.05
NA_PROB = "not applicable (probability: the value is in the column 'mean')"
NO_MEAN_DT = "not assessed: the stored function gives no mean (the doubling time is ln 2 / median growth rate)"

CHAIN_COLUMNS = ["analysis", "program", "setting", "genome_set", "seed", "command_line", "started_utc", "ended_utc", "wall_time_s",
                 "exit_code", "last_state_logged", "reached_last_step", "started_again", "root_date_rule_met",
                 "first_state_root_date_rule", "density_criterion_met", "n_sampled_trees_evaluated", "retained", "why",
                 "n_states_after_burn_in", "attempts", "latest_root_date_logged", "largest_difference_of_the_density_log_units",
                 "end_of_the_last_attempt_by_the_ledger", "density_difference_at_step_0", "largest_density_difference_after_step_0",
                 "density_criterion_met_after_step_0", "continued_from_a_saved_state", "files_of_the_chain",
                 "saved_state_from_which_the_chain_was_continued", "logged_states_read", "logged_states_discarded_as_burn_in", "note_on_the_analysis",
                 "start_end_wall_time_exit_code_and_command_line_are_those_of"]
EST_COLUMNS = ["analysis", "program", "setting", "genome_set", "genomes", "most_recent_collection_date", "quantity", "period_start",
               "period_end", "end_used", "median", "mean", "hpd95_lower", "hpd95_upper", "ess_pooled", "rhat_split",
               "meets_convergence_criterion", "chains_run", "chains_retained", "states_used", "unit",
               "period_length_days", "mid_point_of_the_most_recent_interval", "days_beyond_the_most_recent_mid_point",
               "states_with_R_set_to_0", "ess_pooled_of_the_reproduction_number_series", "rhat_split_of_the_reproduction_number_series",
               "meets_convergence_criterion_by_the_reproduction_number_series", "note", "note_on_the_analysis"]
INT_COLUMNS = ["analysis", "program", "setting", "genome_set", "interval", "older_edge", "newer_edge", "mid_point", "mid_point_decimal",
               "genomes_collected", "ln_Ne_tau_median", "ln_Ne_tau_hpd95_lower", "ln_Ne_tau_hpd95_upper", "ess_pooled", "rhat_split",
               "meets_convergence_criterion", "status", "chains_retained", "states_used", "unit", "note_on_the_analysis"]
TRAJ_COLUMNS = ["analysis", "program", "setting", "genome_set", "series", "date", "date_decimal", "interval", "older_date", "newer_date",
                "median", "hpd95_lower", "hpd95_upper", "unit", "ess_pooled", "rhat_split", "note", "meets_convergence_criterion",
                "chains_retained", "states_used", "note_on_the_analysis"]


def exp_(x):
    """exp(x); a value too large for a floating-point number gives inf instead of ending the script (seen on made-up logs of a test only)"""
    try:
        return math.exp(x)
    except OverflowError:
        return float("inf")


def meets(ess, rhat):
    if ess != ess or rhat != rhat:
        return "not assessed: the diagnostics have no value"
    return "yes" if (ess >= ESS_MIN and rhat <= RHAT_MAX) else "no"


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--config", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = json.load(open(a.config))
    sys.path.insert(0, cfg["pipeline"])
    import analysis as A
    import rt_lib as L
    os.makedirs(a.out, exist_ok=True)
    an_name, engine, model = cfg["analysis"], cfg["engine"], cfg["model"]
    base = dict(analysis=an_name, program=cfg["program"], setting=cfg["setting"], genome_set=cfg["genome_set"])
    gen = L.GenTime(float(cfg["generation_time_mean_days"]), float(cfg["generation_time_sd_days"]))
    dd = L.tip_dates_from_fasta(cfg["alignment"])
    steps = int(cfg["steps"])
    NOTE = str(cfg.get("note_on_the_analysis", "") or "")          # change C

    # ------------------------------------------------------------------------------------------------ chains
    chains, meta = [], {}
    for ch in cfg["chains"]:
        if not ch.get("log") or not os.path.exists(ch["log"]):
            meta[int(ch["seed"])] = dict(ch, no_log=True)
            continue
        c = A.Chain("delphy" if engine == "delphy" else "beast", ch["log"], cfg["alignment"], int(ch["seed"]), burnin=0.30,
                    first_tip=dd.min(), last_tip=dd.max(), label=ch["chain"], extra_logs=ch.get("extra_logs") or None)
        steps_c = int(ch.get("steps_of_this_chain") or steps)          # change A
        c.info["steps"] = steps_c
        c.info["complete"] = bool(c.info["last_state"] >= steps_c)
        if not c.info["complete"] and c.retained:
            c.retained = False
            c.why = (f"incomplete (stopped at step {c.info['last_state']} of {steps_c})" if engine == "delphy"
                     else f"incomplete (state {c.info['last_state']} of {steps_c})")
        elif not c.info["complete"]:
            c.why += f"; incomplete (stopped at step {c.info['last_state']} of {steps_c})"
        chains.append(c)
        meta[int(ch["seed"])] = ch
    grid = None
    if model == "skygrid":
        grid = L.Grid.from_cutoff(int(cfg["num_parameters"]), float(cfg["cutoff"]), dd.max())
    an = A.Analysis(an_name, engine, "Skygrid" if grid is not None else "exponential growth", os.path.basename(cfg["alignment"]),
                    cfg["alignment"], chains, grid=grid, gen=gen)
    density_evaluated = False
    if engine == "delphy" and grid is not None:
        an.progress = {c.seed: meta[c.seed]["stdout"] for c in chains}
        an.trees = {c.seed: meta[c.seed]["trees"] for c in chains}
        A.apply_second_criterion(an, 100.0, 50)
        density_evaluated = True
    if engine == "delphy":
        dens_na = ("not evaluated: the stored pipeline evaluates the density criterion for chains of the Skygrid only "
                   "(derive_rt.py; the exact density of rt_lib.py is that of a staircase of population sizes)")
    else:
        dens_na = "not applicable"

    rows = []
    for ch in cfg["chains"]:
        seed = int(ch["seed"])
        c = next((x for x in chains if x.seed == seed), None)
        r = dict(base, seed=seed, command_line=ch.get("command_line", ""), started_utc=ch.get("started_utc", ""), ended_utc=ch.get("ended_utc", ""),
                 wall_time_s=ch.get("wall_time_s", ""), exit_code=ch.get("exit_code", ""), attempts=ch.get("attempts", ""),
                 started_again=ch.get("started_again", ""), end_of_the_last_attempt_by_the_ledger=ch.get("end_of_the_last_attempt_by_the_ledger", ""),
                 continued_from_a_saved_state=ch.get("continued_from_a_saved_state", "no"), files_of_the_chain=ch.get("files_of_the_chain", ""),
                 saved_state_from_which_the_chain_was_continued=ch.get("saved_state_from_which_the_chain_was_continued", ""),          # change B
                 logged_states_read=(c.n_logged if c is not None else 0), logged_states_discarded_as_burn_in=(c.nb if c is not None else 0))
        r["start_end_wall_time_exit_code_and_command_line_are_those_of"] = ch.get("start_end_wall_time_exit_code_and_command_line_are_those_of", "")          # change D
        if c is None:
            r.update(last_state_logged="", reached_last_step="no", root_date_rule_met="not assessed: the chain has no log",
                     first_state_root_date_rule="", density_criterion_met=dens_na if not density_evaluated else "not assessed: the chain has no log",
                     n_sampled_trees_evaluated="", retained="no", why="no log", n_states_after_burn_in=0, latest_root_date_logged="",
                     largest_difference_of_the_density_log_units="", density_difference_at_step_0="", largest_density_difference_after_step_0="",
                     density_criterion_met_after_step_0="")
            rows.append(r)
            continue
        cut_here = bool(ch.get("steps_of_this_chain")) and int(ch["steps_of_this_chain"]) < steps          # change E
        r.update(last_state_logged=c.info["last_state"],
                 reached_last_step=(f"not applicable: the chain is cut at its saved state ({c.info['last_state']} of {steps} states)" if cut_here
                                    else ("yes" if c.info["complete"] else "no")),
                 root_date_rule_met="yes" if c.info["collapsed"] else "no",
                 first_state_root_date_rule=("" if not c.info["collapsed"] else int(c.info["first_collapsed_state"])),
                 retained="yes" if c.retained else "no", why=c.why, n_states_after_burn_in=len(c.post),
                 latest_root_date_logged=c.info["latest_root_date"])
        if density_evaluated:
            E = c.err
            e0 = E[E["state"] == 0]
            later = E[E["state"] > 0]
            r.update(density_criterion_met="yes" if c.info["flagged_by_density_criterion"] else "no", n_sampled_trees_evaluated=c.info["n_trees_checked"],
                     largest_difference_of_the_density_log_units=c.info["approx_error_max"],
                     density_difference_at_step_0=(float(e0["err"].iloc[0]) if len(e0) else "not assessed: the tree of step 0 is not among the sampled trees that were evaluated"),
                     largest_density_difference_after_step_0=(float(later["err"].max()) if len(later) else ""),
                     density_criterion_met_after_step_0=("yes" if (len(later) and float(later["err"].max()) > 100.0) else "no"))
        else:
            r.update(density_criterion_met=dens_na, n_sampled_trees_evaluated=dens_na if engine == "delphy" else "not applicable",
                     largest_difference_of_the_density_log_units=dens_na if engine == "delphy" else "not applicable",
                     density_difference_at_step_0=dens_na if engine == "delphy" else "not applicable",
                     largest_density_difference_after_step_0=dens_na if engine == "delphy" else "not applicable",
                     density_criterion_met_after_step_0=dens_na if engine == "delphy" else "not applicable")
        rows.append(r)
    CH = pd.DataFrame(rows).assign(note_on_the_analysis=NOTE)[CHAIN_COLUMNS]
    CH.to_csv(os.path.join(a.out, f"chains_{an_name}_{SUFFIX}.csv"), index=False)

    kept = an.kept
    common = dict(base, genomes=an.n_tips, most_recent_collection_date=str(an.last_tip.date()), chains_run=len(cfg["chains"]),
                  chains_retained=len(kept), states_used=int(sum(len(c.post) for c in kept)))
    est = []

    def put(quantity, unit, median="", mean="", lo="", hi="", ess="", rhat="", crit=None, **extra):
        r = dict(common, quantity=quantity, unit=unit, median=median, mean=mean, hpd95_lower=lo, hpd95_upper=hi, ess_pooled=ess, rhat_split=rhat,
                 meets_convergence_criterion=(crit if crit is not None else (meets(ess, rhat) if ess != "" else "")),
                 period_start="", period_end="", end_used="", period_length_days="", mid_point_of_the_most_recent_interval="",
                 days_beyond_the_most_recent_mid_point="", states_with_R_set_to_0="", ess_pooled_of_the_reproduction_number_series="",
                 rhat_split_of_the_reproduction_number_series="", meets_convergence_criterion_by_the_reproduction_number_series="", note="")
        r.update(extra)
        est.append(r)

    if not kept:
        put("no chain is retained: no estimate", "", crit="not assessed: no chain is retained")
        pd.DataFrame(est).assign(note_on_the_analysis=NOTE)[EST_COLUMNS].to_csv(os.path.join(a.out, f"estimates_{an_name}_{SUFFIX}.csv"), index=False)
        print(an_name, "no chain retained")
        return 0

    def growth_rows(rs, label_g, label_R, label_dt, pinfo, doubling):
        """rs: growth rate per year, one array per retained chain (convention of analysis.py and rt_lib.py)"""
        s = A.r_summary(rs, gen, label_g)
        allr = np.concatenate(rs)
        Rser = [gen.r2R(x) for x in rs]
        sR = L.summarise(Rser, label_R)
        put(label_g, "per day", s["growth_per_day_median"], float(np.mean(allr / L.YEAR_DAYS)), s["growth_per_day_hpd_lo"], s["growth_per_day_hpd_hi"],
            s["ess_pooled"], s["rhat_split"], **pinfo)
        put(label_R, "", s["R_median"], sR["mean"], s["R_hpd_lo"], s["R_hpd_hi"], s["ess_pooled"], s["rhat_split"],
            states_with_R_set_to_0=s["n_samples_R_set_to_0"], ess_pooled_of_the_reproduction_number_series=sR["ess_pooled"],
            rhat_split_of_the_reproduction_number_series=sR["rhat_split"],
            meets_convergence_criterion_by_the_reproduction_number_series=meets(sR["ess_pooled"], sR["rhat_split"]), **pinfo)
        if doubling:
            pin = dict(pinfo)
            notes = [pin.pop("note", ""), ("" if s["doubling_time_d"] == s["doubling_time_d"] else
                                           "the stored function gives a doubling time only where the median growth rate is above 0; "
                                           "it gives bounds only where the lower bound of the interval of the growth rate is above 0")]
            put(label_dt, "days", s["doubling_time_d"], NO_MEAN_DT, s["doubling_time_d_lo"], s["doubling_time_d_hi"], s["ess_pooled"], s["rhat_split"],
                note="; ".join(x for x in notes if x), **pin)
        return Rser, s

    def probability(quantity, ind, pinfo, value=None):
        """ind: indicator series (0/1), one array per retained chain"""
        p = float(np.concatenate(ind).mean()) if value is None else value
        if any(float(np.std(x)) == 0.0 for x in ind):
            put(quantity, "probability", NA_PROB, p, NA_PROB, NA_PROB, "", "", crit="not assessed: the indicator is constant in at least one chain", **pinfo)
        else:
            e, rh = L.ess_multichain(ind), L.rhat_split(ind)
            put(quantity, "probability", NA_PROB, p, NA_PROB, NA_PROB, e, rh, **pinfo)

    post_cols = {}
    periods = pd.read_csv(cfg["periods"], dtype=str, keep_default_na=False)
    periods = periods[periods["period"].isin(["P1", "P2", "P3", "P4", "P5", "P6"])]
    Rs = {}
    if grid is not None:
        lp = an.lnNe()
        mid0 = float(grid.mids[0])
        mid0_date = str(L.dec2date(mid0))
        for _, p in periods.iterrows():
            k, start, end = p["period"], p["start"], p["end"]
            t0 = L.decimal_year(start)
            own = (end == "END")
            t1 = mid0 if own else L.decimal_year(end)
            pinfo = dict(period_start=start, period_end=(mid0_date if own else end), end_used=("own" if own else "fixed date"),
                         period_length_days=(t1 - t0) * L.YEAR_DAYS, mid_point_of_the_most_recent_interval=mid0_date,
                         days_beyond_the_most_recent_mid_point=max(0.0, (t1 - mid0) * L.YEAR_DAYS))
            why_not = ""
            if t1 <= t0:
                why_not = "not evaluated: no length"
            elif (not own) and t1 > an.t_last + 1e-9:
                why_not = ("not evaluated: the period ends after the most recent collection date of the set "
                           "(rule of the stored function analysis.Analysis.windows)")
            if why_not:
                pinfo["period_length_days"] = (pinfo["period_length_days"] if t1 > t0 else 0.0)
                for q, u in ((f"growth rate, {k}", "per day"), (f"reproduction number, {k}", "")) + (((f"doubling time, {k}", "days"),) if k == "P1" else ()):
                    put(q, u, why_not, why_not, why_not, why_not, "", "", crit=why_not, **pinfo)
                continue
            if (not own) and t1 > mid0 + 1e-9:
                pinfo["note"] = "end date later than the mid-point of the most recent interval: level of that interval used"
            c = L.window_contrast(grid, t0, t1, "interp")
            rs = [x @ c for x in lp]
            sfx = k + ("_own_end" if own else "")
            Rser, s = growth_rows(rs, f"growth rate, {k}", f"reproduction number, {k}", f"doubling time, {k}", pinfo, doubling=(k == "P1"))
            Rs[k] = (rs, Rser, pinfo, s)
            post_cols[f"growth_rate_per_day_{sfx}"] = [x / L.YEAR_DAYS for x in rs]
            post_cols[f"reproduction_number_{sfx}"] = Rser
    else:
        rs = an.series("exponential.growthRate")
        pinfo = dict(period_start="not applicable", period_end="not applicable", end_used="not applicable",
                     period_length_days="not applicable", mid_point_of_the_most_recent_interval="not applicable",
                     days_beyond_the_most_recent_mid_point="not applicable",
                     note="constant exponential growth: one growth rate for the whole tree")
        Rser, s = growth_rows(rs, "growth rate, constant exponential growth", "reproduction number, constant exponential growth",
                              "doubling time, constant exponential growth", pinfo, doubling=True)
        Rs["exp"] = (rs, Rser, pinfo, s)
        post_cols["growth_rate_per_day_constant_exponential_growth"] = [x / L.YEAR_DAYS for x in rs]
        post_cols["reproduction_number_constant_exponential_growth"] = Rser
        na = "not applicable: under constant exponential growth the growth rate is one number for the whole tree"
        for _, p in periods.iterrows():
            for q, u in ((f"growth rate, {p['period']}", "per day"), (f"reproduction number, {p['period']}", "")) + (
                    ((f"doubling time, {p['period']}", "days"),) if p["period"] == "P1" else ()):
                put(q, u, na, na, na, na, "", "", crit=na, period_start=p["start"], period_end=("" if p["end"] == "END" else p["end"]),
                    end_used=("own" if p["end"] == "END" else "fixed date"))

    # ---- evolutionary rate, date of the root, kappa, precision
    s = L.summarise(an.series("clock.rate"), "clock rate")
    put("evolutionary rate", "substitutions per site and year", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"])
    tm = A.tmrca_summary(an)
    rh_all = np.concatenate(an.root_height())
    mean_days = float(np.mean(rh_all * L.YEAR_DAYS))
    put("date of the most recent common ancestor (date)", "date (day resolution, from the root height of the log)", tm["tmrca_median"],
        str((an.last_tip - pd.Timedelta(days=int(round(mean_days)))).date()), tm["tmrca_hpd_early"], tm["tmrca_hpd_late"], tm["ess_pooled"], tm["rhat_split"])
    put("date of the most recent common ancestor (decimal year)", "decimal year (year of 365 days)", tm["tmrca_decimal_median"],
        an.t_last - mean_days / L.YEAR_DAYS, tm["tmrca_decimal_early"], tm["tmrca_decimal_late"], tm["ess_pooled"], tm["rhat_split"])
    if "kappa" in kept[0].post.columns:
        s = L.summarise(an.series("kappa"), "kappa")
        put("kappa", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"])
    else:
        na = "not applicable: the substitution model of the setting has no kappa (GTR: six relative rates)"
        put("kappa", "", na, na, na, na, "", "", crit=na)
    if grid is not None:
        tau, fixed = A.precision_of(an)
        if fixed:
            na = "not applicable: the precision is fixed"
            put("precision of the Skygrid", "", na, na, na, na, "", "", crit=na, note=f"fixed at {tau}")
        else:
            s = L.summarise(an.series("skygrid.precision"), "precision")
            put("precision of the Skygrid", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"])
    # ---- probabilities
    if grid is not None:
        def pi(k):
            return {x: Rs[k][2][x] for x in ("period_start", "period_end", "end_used", "period_length_days", "mid_point_of_the_most_recent_interval",
                                             "days_beyond_the_most_recent_mid_point")}
        if "P1" in Rs and "P2" in Rs:
            probability("P(R of P2 < R of P1)", [(b < c_).astype(float) for b, c_ in zip(Rs["P2"][1], Rs["P1"][1])], dict(end_used="fixed date"))
        else:
            put("P(R of P2 < R of P1)", "probability", crit="not evaluated: one of the two periods is not evaluated")
        for k, q, f in (("P1", "P(R of P1 > 1)", lambda r: (r > 0)), ("P2", "P(R of P2 > 1)", lambda r: (r > 0)), ("P3", "P(R of P3 < 1)", None)):
            if k not in Rs:
                put(q, "probability", crit="not evaluated: the period is not evaluated")
                continue
            if f is not None:
                probability(q, [f(x).astype(float) for x in Rs[k][0]], pi(k), value=Rs[k][3]["P_R_gt_1"])
            else:
                probability(q, [(x < 1.0).astype(float) for x in Rs[k][1]], pi(k))
    else:
        probability("P(R > 1), constant exponential growth", [(x > 0).astype(float) for x in Rs["exp"][0]],
                    dict(note="constant exponential growth: one growth rate for the whole tree"), value=Rs["exp"][3]["P_R_gt_1"])
        na = "not applicable: under constant exponential growth the growth rate is one number for the whole tree"
        for q in ("P(R of P2 < R of P1)", "P(R of P1 > 1)", "P(R of P2 > 1)", "P(R of P3 < 1)"):
            put(q, "probability", na, na, na, na, "", "", crit=na)
    E = pd.DataFrame(est).assign(note_on_the_analysis=NOTE)[EST_COLUMNS]
    E.to_csv(os.path.join(a.out, f"estimates_{an_name}_{SUFFIX}.csv"), index=False)

    # ------------------------------------------------------------------------------------- intervals, trajectories
    if grid is not None:
        tabs = A.skygrid_tables(an)
        I = tabs["intervals"]
        K = grid.K
        irows = []
        for _, r in I.iterrows():
            i = int(r["interval"])
            irows.append(dict(base, interval=i, older_edge=r["older_edge"], newer_edge=r["newer_edge"], mid_point=r["mid_point"],
                              mid_point_decimal=(float(grid.mids[i - 1]) if i < K else ""), genomes_collected=int(r["n_genomes_collected"]),
                              ln_Ne_tau_median=r["lnNe_median"], ln_Ne_tau_hpd95_lower=r["lnNe_hpd_lo"], ln_Ne_tau_hpd95_upper=r["lnNe_hpd_hi"],
                              ess_pooled=r["ess_pooled"], rhat_split=r["rhat_split"], meets_convergence_criterion=meets(r["ess_pooled"], r["rhat_split"]),
                              status=r["status"], chains_retained=len(kept), states_used=common["states_used"], unit="ln(years)"))
        pd.DataFrame(irows).assign(note_on_the_analysis=NOTE)[INT_COLUMNS].to_csv(os.path.join(a.out, f"intervals_{an_name}_{SUFFIX}.csv"), index=False)
        tr = []
        tb = dict(base, chains_retained=len(kept), states_used=common["states_used"])
        for _, r in I.iterrows():
            i = int(r["interval"])
            tr.append(dict(tb, series="Ne tau at the mid-point of an interval", date=(r["mid_point"] if i < K else ""),
                           date_decimal=(float(grid.mids[i - 1]) if i < K else ""), interval=i, older_date="", newer_date="",
                           median=exp_(r["lnNe_median"]), hpd95_lower=exp_(r["lnNe_hpd_lo"]), hpd95_upper=exp_(r["lnNe_hpd_hi"]), unit="years",
                           ess_pooled=r["ess_pooled"], rhat_split=r["rhat_split"],
                           note=(r["status"] if i < K else f"interval {K} is open towards the past (everything before the cut-off); it has no mid-point"),
                           meets_convergence_criterion=meets(r["ess_pooled"], r["rhat_split"])))
        # every 7th day: from the date of the first interval mid-point after 1 January 2026 (intervals 1 to K-1) to the date of the
        # mid-point of the most recent interval; the first and the last date are included; Ne tau at the START of the day
        first_mid = min(float(m) for m in grid.mids[:K - 1] if m > L.decimal_year(cfg["first_date_not_before"]))
        d, last_day = pd.Timestamp(L.dec2date(first_mid)), pd.Timestamp(L.dec2date(mid0))
        days = []
        while d < last_day:
            days.append(d)
            d += pd.Timedelta(days=7)
        days.append(last_day)
        for d in days:
            t = L.decimal_year(str(d.date()))
            w = grid.weights(t, "interp")
            s = L.summarise([x @ w for x in lp], "lnNe")
            tr.append(dict(tb, series="Ne tau at every 7th day", date=str(d.date()), date_decimal=t, interval=grid.interval_of(t), older_date="", newer_date="",
                           median=exp_(s["median"]), hpd95_lower=exp_(s["hpd_lo"]), hpd95_upper=exp_(s["hpd_hi"]), unit="years",
                           ess_pooled=s["ess_pooled"], rhat_split=s["rhat_split"], note="ln(Ne tau) read by linear interpolation between interval mid-points"
                           + (("; " + cfg["rule_of_the_dates_of_the_series"]) if cfg.get("rule_of_the_dates_of_the_series") else ""),
                           meets_convergence_criterion=meets(s["ess_pooled"], s["rhat_split"])))
        RT = tabs["rt"]
        for series, unit, cm, cl, chh in (("reproduction number between successive mid-points", "", "R_median", "R_hpd_lo", "R_hpd_hi"),
                                          ("growth rate between successive mid-points", "per day", "growth_per_day_median", "growth_per_day_hpd_lo", "growth_per_day_hpd_hi")):
            for _, r in RT.iterrows():
                i = int(r["interval_pair"])
                tr.append(dict(tb, series=series, date=r["time_boundary"], date_decimal=float(grid.t_anchor - i * grid.D), interval=f"{i + 1} to {i}",
                               older_date=r["older_interval_mid"], newer_date=r["newer_interval_mid"], median=r[cm], hpd95_lower=r[cl], hpd95_upper=r[chh],
                               unit=unit, ess_pooled=r["ess_pooled"], rhat_split=r["rhat_split"], note=r["status"],
                               meets_convergence_criterion=meets(r["ess_pooled"], r["rhat_split"])))
        pd.DataFrame(tr).assign(note_on_the_analysis=NOTE)[TRAJ_COLUMNS].to_csv(os.path.join(a.out, f"trajectories_{an_name}_{SUFFIX}.csv"), index=False)

    # -------------------------------------------------------------------------------------------------- posterior
    fr = []
    for j, c in enumerate(kept):
        p = c.post
        d = pd.DataFrame({"seed": c.seed, "state": p["state"].values, "evolutionary_rate": p["clock.rate"].values,
                          "root_height_years": p[c.rh_col].values})
        d["root_date"] = L.root_dates_from_height(p[c.rh_col].values, an.last_tip).astype(str)
        d["kappa"] = p["kappa"].values if "kappa" in p.columns else ""
        if grid is not None:
            lpj = p[c.lp_cols].values
            for i in range(lpj.shape[1]):
                d[f"ln_Ne_tau_{i + 1}"] = lpj[:, i]
            d["skygrid_precision"] = p["skygrid.precision"].values
        else:
            d["exponential_growth_rate_per_year"] = p["exponential.growthRate"].values
            d["exponential_population_size_at_the_most_recent_date_years"] = p["exponential.popSize"].values
        for col in [x for x in p.columns if x.startswith("gtr.rates")]:
            d[col] = p[col].values
        for name, lst in post_cols.items():
            d[name] = lst[j]
        fr.append(d)
    PS = pd.concat(fr, ignore_index=True)
    PS.to_csv(os.path.join(a.out, f"posterior_{an_name}_{SUFFIX}.csv.gz"), index=False, compression={"method": "gzip", "mtime": 0})

    # ---------------------------------------------------------------- root date: shortest intervals in whole days (Delphy)
    # (this list is written for every analysis of Delphy: Delphy logs root heights of whole days; BEAST X does not)
    if engine == "delphy":
        days = np.sort(np.round(rh_all * L.YEAR_DAYS).astype(int))
        n = len(days)
        m = int(np.ceil(0.95 * n))
        widths = days[m - 1:] - days[:n - m + 1]
        w = int(widths.min())
        js = np.where(widths == w)[0]
        ivs = sorted({(int(days[j]), int(days[j + m - 1])) for j in js})
        f = lambda h: str((an.last_tip - pd.Timedelta(days=int(h))).date())
        whole = bool(np.all(np.abs(rh_all * L.YEAR_DAYS - np.round(rh_all * L.YEAR_DAYS)) < 0.05))
        pd.DataFrame([dict(base, states_used=n, states_in_a_95_per_cent_interval=m, every_logged_root_height_is_a_whole_number_of_days="yes" if whole else "no",
                           width_of_the_shortest_interval_days=w, number_of_shortest_intervals=len(ivs),
                           earliest_shortest_interval=f"{f(ivs[-1][1])} to {f(ivs[-1][0])}", latest_shortest_interval=f"{f(ivs[0][1])} to {f(ivs[0][0])}",
                           interval_of_the_table_of_estimates=f"{tm['tmrca_hpd_early']} to {tm['tmrca_hpd_late']}",
                           all_shortest_intervals="; ".join(f"{f(hi)} to {f(lo)}" for lo, hi in reversed(ivs)))]).assign(note_on_the_analysis=NOTE).to_csv(
            os.path.join(a.out, f"root_date_intervals_{an_name}_{SUFFIX}.csv"), index=False)
    open(os.path.join(a.out, f"derivation_{an_name}.json"), "w").write(json.dumps(dict(
        analysis=an_name, derived_utc=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        chains_run=len(cfg["chains"]), chains_with_a_log=len(chains), chains_retained=len(kept), states_used=common["states_used"],
        python=sys.version.split()[0], numpy=np.__version__, pandas=pd.__version__), indent=1) + "\n")
    print(an_name, "chains", len(cfg["chains"]), "retained", len(kept), "states used", common["states_used"], "rows of estimates", len(E))
    return 0


if __name__ == "__main__":
    sys.exit(main())
