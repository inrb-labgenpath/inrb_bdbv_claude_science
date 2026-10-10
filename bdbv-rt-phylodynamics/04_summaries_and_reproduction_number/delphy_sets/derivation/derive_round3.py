#!/usr/bin/env python3
"""derive_round3.py - retention, convergence and estimates of ONE analysis of the Delphy sets (round of October 2026).

usage: python derive_round3.py --analysis delphy_fixed_C95 --run-root RUN_ROOT --pipeline TREE/work/pipeline_ext --outdir OUT

Compared with the first version (derive_round.py, kept unchanged), every number is computed in the same way; what is
added: one more column of the table of chains (end of the last attempt, read from the ledger), and a table that
compares what this script computes in its own code with the functions of analysis.py and rt_lib.py where such a
function exists (own_code_against_stored_<analysis>_20261002.csv).

A. Computed by the functions of analysis.py and rt_lib.py (work/pipeline_ext, imported unchanged)
  chains, burn-in, root-date rule            analysis.Chain
  density criterion                          analysis.apply_second_criterion (called chain by chain, see D)
  grid of the Skygrid                        rt_lib.Grid.from_cutoff
  weights of ln(Ne tau) at a date            rt_lib.Grid.weights / rt_lib.window_contrast, mode 'interp'
  growth rate -> reproduction number         rt_lib.GenTime.r2R
  growth rate and reproduction number of a period: median, HPD bounds, ESS, split R-hat, doubling time
                                             analysis.r_summary
  evolutionary rate, kappa, precision: median, mean, HPD bounds, ESS, split R-hat
                                             rt_lib.summarise
  date of the most recent common ancestor: median, HPD bounds, ESS, split R-hat
                                             analysis.tmrca_summary
  intervals of ln(Ne tau), reproduction number between successive mid-points
                                             analysis.skygrid_tables (tables 'intervals' and 'rt')
  root date of a state                       rt_lib.root_dates_from_height
  HPD interval, ESS and split R-hat wherever they appear in B
                                             rt_lib.hpd, rt_lib.ess_multichain, rt_lib.rhat_split

B. Computed by THIS script in its own code, because analysis.py and rt_lib.py hold no function for it
  1. the MEAN of the growth rate per day and of the reproduction number of a period: arithmetic mean of the pooled
     states (numpy.mean); analysis.r_summary gives no mean
  2. the MEAN date of the most recent common ancestor: most recent collection date minus the arithmetic mean of the
     root height (numpy.mean), as a date by the rounding of analysis.tmrca_summary and as decimal year
  3. the four POSTERIOR PROBABILITIES P(R of P2 < R of P1), P(R of P1 > 1), P(R of P2 > 1), P(R of P3 < 1): share of
     the pooled states for which the relation holds (mean of an indicator); their ESS and split R-hat are those of
     the indicator, by the functions of A.  analysis.r_summary gives P(R > 1) as the share of states with a growth
     rate above 0; the table of C compares the two
  4. Ne tau at every 7th day: weights by rt_lib.Grid.weights; median = numpy.median of exp(ln Ne tau), HPD bounds
     = exp of rt_lib.hpd of ln(Ne tau), which is how analysis.skygrid_tables treats an interval
  5. the degenerate state: running median of 50 logged states of the precision
     (pandas rolling(50, min_periods=50).median(), the form used by the run script of the estimated smoothing)
  6. the last state of a log, read as text; whether a file ends inside a line
  7. the periods P5 and P6 and the common end of P3 and P6: the DATES are passed by this script; the growth rate and
     the reproduction number of these periods are computed by the functions of A

C. own_code_against_stored_<analysis>_20261002.csv: for B.1 the mean by rt_lib.summarise of the same series, for B.3
   the P(R > 1) of analysis.r_summary; differences are reported, no value of the tables is replaced.

D. Values of this round that are passed where a script of the earlier round holds a constant or a default written for
   another snapshot: number of steps (derive_rt.py --steps), suffix of the file names, label of the analysis, the
   alignment (from which the functions read the first and the most recent collection date).  P1 to P4 are the
   dictionary analysis.W_DATES (W1, W2, W3, W3b); the script stops if the periods file of this round does not give
   the same dates.  analysis.apply_second_criterion loops over the chains of an analysis and stops at the first
   chain whose files it cannot read; it is therefore called once per chain (the list of chains of the analysis
   object is set to that chain for the call), so that a chain with unreadable files is reported as 'not assessed' and
   does not hide the others.  The function itself is unchanged.

This file holds the code of derive_round2.py unchanged; only the header and the names of the scripts that it calls
differ.
"""
import argparse
import csv
import json
import math
import os
import re
import sys
import traceback

import numpy as np
import pandas as pd

SUFFIX = "20261002"
STEPS = 1500000000
BURNIN = 0.30
K_PARAM, CUTOFF = 20, 0.8
GEN_MEAN, GEN_SD = 15.3, 9.3
DENSITY_THRESHOLD, MAX_TREES = 100.0, 50
ESS_MIN, RHAT_MAX = 200.0, 1.05
PRECISION_LIMIT, PRECISION_WINDOW = 0.05, 50
PROGRAM = "Delphy 1.4.1"
NA = "not assessed"
STORED_NAME = {"P1": "W1", "P2": "W2", "P3": "W3", "P4": "W3b"}


def crit(ess, rhat):
    if ess is None or rhat is None or ess != ess or rhat != rhat:
        return NA
    return "yes" if (ess >= ESS_MIN and rhat <= RHAT_MAX) else "no"


def utc_of(epoch_text):
    import time
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime(int(epoch_text)))


def read_text(p, default=""):
    try:
        return open(p).read().strip()
    except OSError:
        return default


def last_state_of_log(path):
    """last state of the log (first column of the last complete line), by reading the file as text"""
    last, n = None, 0
    try:
        with open(path, errors="replace") as fh:
            for line in fh:
                if not line.endswith("\n"):
                    break
                f = line.split("\t", 1)[0]
                if f.isdigit():
                    last = int(f); n += 1
    except OSError:
        return None, 0
    return last, n


def whole_lines_copy(path, outdir, min_fields=None):
    """(path to read, note).  If the file ends inside a line (a program that ends itself does not write its buffers
    out), a copy without the cut line is written to outdir and returned; the file of the chain is not changed."""
    try:
        with open(path, "rb") as fh:
            data = fh.read()
    except OSError:
        return path, "file not found"
    whole = data.endswith(b"\n") or data.rstrip().endswith(b"End;")     # a complete tree file of Delphy ends with 'End;' without a line end
    cut = len(data) > 0 and not whole
    body = data[:data.rfind(b"\n") + 1] if cut else data
    note = "the file ends inside a line; the cut line was left out" if cut else ""
    if min_fields is not None:
        lines = body.split(b"\n")
        while len(lines) > 1 and lines[-1] == b"":
            lines.pop()
        if lines and not lines[-1].startswith(b"#") and lines[-1].count(b"\t") + 1 < min_fields:
            body = b"\n".join(lines[:-1]) + b"\n"; cut = True
            note = "the last line of the file holds fewer fields than the header; it was left out"
    if not cut:
        return path, ""
    os.makedirs(outdir, exist_ok=True)
    out = os.path.join(outdir, os.path.basename(path))
    with open(out, "wb") as fh:
        fh.write(body)
    return out, note


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--analysis", required=True)
    ap.add_argument("--run-root", required=True)
    ap.add_argument("--pipeline", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--only-seeds", default="", help="TEST of the code only: restrict to these seeds; the output is marked")
    a = ap.parse_args()
    sys.path.insert(0, os.path.abspath(a.pipeline))
    import rt_lib as L
    import analysis as A
    os.makedirs(a.outdir, exist_ok=True)
    rr = a.run_root
    sched = [r for r in csv.DictReader(open(f"{rr}/aln/chain_schedule_round_20261002.csv", newline=""))
             if r["track"] == "Delphy sets" and r["analysis"] == a.analysis]
    assert sched, f"no chain of {a.analysis} in the schedule"
    assert all(int(r["steps"]) == STEPS and float(r["burn_in_share"]) == BURNIN for r in sched)
    test_only = bool(a.only_seeds)
    if test_only:
        keep = set(a.only_seeds.split(","))
        sched = [r for r in sched if r["seed"] in keep]
    gset, setting = sched[0]["genome_set"], sched[0]["setting"]
    estimated = a.analysis.startswith("delphy_estimated_")
    fasta = f"{rr}/aln/{gset}_{SUFFIX}.fasta"
    dd = L.tip_dates_from_fasta(fasta)
    A.ERR_THRESHOLD = DENSITY_THRESHOLD
    ledger = {}
    if os.path.exists(f"{rr}/ledger.csv"):
        for r in csv.DictReader(open(f"{rr}/ledger.csv", newline="")):
            if r["analysis"] == a.analysis:
                ledger.setdefault(r["seed"], []).append(r)

    # ------------------------------------------------------------------ chains, as derive_rt.py builds them
    chains, crow = [], {}
    for r in sched:
        seed = r["seed"]; p = f"{rr}/runs/{a.analysis}/{a.analysis}_s{seed}"
        att = ledger.get(seed, [])
        outside = [x for x in att if x["outcome"] == "ended_from_outside"]
        last_att = max(att, key=lambda x: int(x["attempt"])) if att else None
        if last_att is None:
            end_last = f"{NA}: the chain is not in the ledger"
        elif last_att["outcome"] == "ended_from_outside":
            end_last = "ended from outside" + (", not started again" if "NOT started again" in last_att["note"] else ", to be started again")
        else:
            end_last = {"completed": "completed", "ended_by_program": "ended by the program", "needs_attention": "needs attention",
                        "running": "running"}.get(last_att["outcome"], last_att["outcome"])
        san = os.path.join(a.outdir, "copies_without_the_cut_line")
        nf = None
        try:
            with open(p + ".log", errors="replace") as fh:
                for line in fh:
                    if line.startswith("state\t"):
                        nf = len(line.rstrip("\n").rstrip("\t").split("\t")); break
        except OSError:
            pass
        f_log, n1 = whole_lines_copy(p + ".log", san, nf)
        f_trees, n2 = whole_lines_copy(p + ".trees", san)
        f_out, n3 = whole_lines_copy(p + ".stdout", san)
        cut_note = "; ".join(f"{k}: {v}" for k, v in (("log", n1), ("trees", n2), ("standard output", n3)) if v)
        last_txt, n_txt = last_state_of_log(f_log)
        st, en = read_text(p + ".start"), read_text(p + ".end")
        row = dict(analysis=a.analysis, program=PROGRAM, setting=setting, genome_set=gset, seed=int(seed),
                   command_line=read_text(p + ".cmd"),
                   started_utc=utc_of(st) if st else "", ended_utc=utc_of(en) if en else "",
                   wall_time_s=(int(en) - int(st)) if st and en else "", exit_code=read_text(p + ".exit", NA),
                   last_state_logged=last_txt if last_txt is not None else NA,
                   reached_last_step="yes" if last_txt == STEPS else "no",
                   started_again=("no" if not outside else
                                  f"yes ({len(outside)} time(s)): " + " | ".join(f"attempt {x['attempt']} ended from outside, exit code "
                                                                                  f"{x['exit_code'] or 'none left'}; {x['note']}" for x in outside)),
                   root_date_rule_met=NA, first_state_root_date_rule=NA, density_criterion_met=NA,
                   n_sampled_trees_evaluated=NA, retained="no", why="", n_states_after_burn_in=NA,
                   address_space_limit_kib=dict(l.split("=", 1) for l in read_text(p + ".env").splitlines() if "=" in l).get("address_space_limit_kib", NA),
                   attempts=len(att) if att else NA,
                   precision_running_median_below_limit=NA, first_state_precision_running_median_below_limit=NA,
                   smallest_running_median_of_precision=NA, smallest_precision=NA, latest_root_date_logged=NA,
                   last_lines_of_standard_output_not_progress="", files_cut_at_the_end=cut_note,
                   end_of_the_last_attempt_by_the_ledger=end_last)
        tail = read_text(p + ".stdout")[-6000:]
        notp = [l.strip() for l in tail.splitlines() if l.strip() and not re.match(r"^[0-9.]+ M steps/s, Step ", l.strip())]
        if row["exit_code"] not in ("0", NA):
            row["last_lines_of_standard_output_not_progress"] = " | ".join(notp[-3:])[:500]
        crow[seed] = row
        if n_txt < 2:
            row["why"] = f"{NA}: the log holds {n_txt} complete state(s); the stored functions need at least 2"
            continue
        try:
            c = A.Chain("delphy", f_log, fasta, int(seed), burnin=BURNIN, first_tip=dd.min(), last_tip=dd.max(),
                        label=f"{a.analysis}_s{seed}")
            assert c.info["last_state"] == last_txt, (c.info["last_state"], last_txt)
        except Exception as e:                                   # a log that the reader of analysis.py cannot read
            row["why"] = f"{NA}: the stored reader could not read the log ({type(e).__name__}: {str(e)[:200]})"
            continue
        sp = int(np.diff(c.full["state"].values[:2])[0])
        assert sp == 300000, f"log spacing {sp}"
        c.info["steps"] = STEPS
        c.info["complete"] = bool(c.info["last_state"] >= STEPS)
        if not c.info["complete"] and c.retained:
            c.retained = False
            c.why = f"incomplete (stopped at step {c.info['last_state']} of {STEPS})"
        c.paths = p; c.f_trees = f_trees; c.f_out = f_out
        chains.append(c)
    grid = L.Grid.from_cutoff(K_PARAM, CUTOFF, dd.max())
    an = A.Analysis(a.analysis, "delphy", "Skygrid", os.path.basename(fasta), fasta, chains, grid=grid, gen=L.GenTime(GEN_MEAN, GEN_SD))
    an.progress = {c.seed: c.f_out for c in chains}
    an.trees = {c.seed: c.f_trees for c in chains}
    all_chains = list(chains)
    for c in all_chains:                                          # function of analysis.py, one chain per call
        an.chains = [c]
        try:
            A.apply_second_criterion(an, DENSITY_THRESHOLD, MAX_TREES)
            c.density_assessed = True
        except Exception as e:
            c.density_assessed = False
            c.density_error = f"{type(e).__name__}: {str(e)[:200]}"
            print(f"ATTENTION {a.analysis} seed {c.seed}: the density criterion could not be evaluated ({c.density_error})")
            if c.retained:
                # a chain whose density criterion cannot be evaluated is not retained silently: it is reported and left out
                c.retained = False
                c.why = f"density criterion {NA} ({c.density_error}); chain left out and reported"
                print(f"ATTENTION {a.analysis} seed {c.seed}: the chain met no rule and is left out because one criterion is {NA}")
    an.chains = all_chains
    for c in all_chains:
        row = crow[str(c.seed)]
        f = c.full
        row.update(root_date_rule_met="yes" if c.info["collapsed"] else "no",
                   first_state_root_date_rule=("" if not c.info["collapsed"] else int(c.info["first_collapsed_state"])),
                   density_criterion_met=(("yes" if c.info["flagged_by_density_criterion"] else "no") if c.density_assessed
                                          else f"{NA}: {c.density_error}"),
                   n_sampled_trees_evaluated=(int(c.info["n_trees_checked"]) if c.density_assessed else NA),
                   retained="yes" if c.retained else "no", why=c.why, n_states_after_burn_in=int(c.info["n_post_burnin"]),
                   latest_root_date_logged=c.info["latest_root_date"], smallest_precision=float(f["skygrid.precision"].min()))
        rm = f["skygrid.precision"].rolling(PRECISION_WINDOW, min_periods=PRECISION_WINDOW).median()
        row.update(precision_running_median_below_limit="yes" if (rm < PRECISION_LIMIT).any() else "no",
                   first_state_precision_running_median_below_limit=(int(f["state"][rm < PRECISION_LIMIT].iloc[0]) if (rm < PRECISION_LIMIT).any() else ""),
                   smallest_running_median_of_precision=(float(rm.min()) if rm.notna().any() else NA))
        if c.density_assessed and len(c.err):
            row["largest_difference_of_the_density_log_units"] = float(c.err["err"].max())
    ch = pd.DataFrame([crow[r["seed"]] for r in sched])
    ch["entered_degenerate_state_by_the_definition_of_the_brief"] = [
        (NA if x.root_date_rule_met == NA else ("yes" if (x.root_date_rule_met == "yes" or x.precision_running_median_below_limit == "yes") else "no"))
        for x in ch.itertuples()]
    ch.to_csv(f"{a.outdir}/chains_{a.analysis}_{SUFFIX}.csv", index=False)

    base = dict(analysis=a.analysis, program=PROGRAM, setting=setting, genome_set=gset, genomes=an.n_tips,
                most_recent_collection_date=str(an.last_tip.date()))
    tail_cols = dict(chains_run=len(sched), chains_retained=len(an.kept))
    meta = dict(analysis=a.analysis, test_run_on_a_subset_of_chains=test_only, genome_set=gset, genomes=an.n_tips,
                first_collection_date=str(an.first_tip.date()), most_recent_collection_date=str(an.last_tip.date()),
                grid_anchor_decimal=grid.t_anchor, interval_width_days=grid.D * L.YEAR_DAYS, own_end_decimal=float(grid.mids[0]),
                own_end_date=str(L.dec2date(grid.mids[0])), chains_run=len(sched), chains_retained=len(an.kept),
                chains_read_by_the_stored_reader=len(all_chains), smoothing_estimated=estimated,
                driver="derive_round3.py", stored_modules=dict(analysis=os.path.abspath(A.__file__), rt_lib=os.path.abspath(L.__file__)))

    # ------------------------------------------------------------------ periods of this round
    per = list(csv.DictReader(open(f"{rr}/aln/periods_round_20261002.csv", newline="")))
    ends = {r["genome_set"]: r for r in csv.DictReader(open(f"{rr}/aln/periods_ends_by_set_20261002.csv", newline=""))}
    alt = {r["set"]: r for r in csv.DictReader(open(f"{rr}/aln/genome_set_alignments_20261002.csv", newline=""))}
    for p in per:                                                 # P1 to P4 must be the windows of analysis.W_DATES
        if p["period"] in STORED_NAME:
            w0, w1 = A.W_DATES[STORED_NAME[p["period"]]]
            assert p["start"] == w0 and (p["end"] == w1 or (p["end"] == "END" and w1 is None)), p
    own_end = float(grid.mids[0])
    e = ends[gset]
    assert f"{own_end:.6f}" == f"{float(e['mid_point_of_the_most_recent_interval_decimal']):.6f}", (own_end, e)
    assert str(L.dec2date(own_end)) == e["mid_point_of_the_most_recent_interval"] and str(an.last_tip.date()) == e["most_recent_collection_date"]
    cset = e["set_of_the_common_end"]
    common_end = float(L.Grid.from_cutoff(K_PARAM, CUTOFF, pd.Timestamp(alt[cset]["most_recent_collection_date"])).mids[0])
    assert f"{common_end:.6f}" == f"{float(e['common_end_decimal']):.6f}" and str(L.dec2date(common_end)) == e["common_end"]
    assert common_end <= own_end + 1e-12
    meta.update(common_end_decimal=common_end, common_end_date=str(L.dec2date(common_end)), set_of_the_common_end=cset,
                most_recent_collection_date_of_the_set_of_the_common_end=alt[cset]["most_recent_collection_date"])
    variants = []                                                 # (key, period, start date, end date text, end_used, t0, t1)
    for p in per:
        if p["period"] == "P7":
            continue
        t0 = L.decimal_year(p["start"])
        if p["end"] == "END":
            variants.append((f"{p['period']}_own_end", p["period"], p["start"], str(L.dec2date(own_end)), "own", t0, own_end))
            variants.append((f"{p['period']}_common_end", p["period"], p["start"], str(L.dec2date(common_end)), "common", t0, common_end))
        else:
            variants.append((p["period"], p["period"], p["start"], p["end"], "fixed date", t0, L.decimal_year(p["end"])))
    for v in variants:
        assert v[5] < v[6] <= own_end + 1e-12, f"period {v[0]} does not lie inside the grid of the set"
    json.dump(meta, open(f"{a.outdir}/meta_{a.analysis}_{SUFFIX}.json", "w"), indent=1)

    est = []

    def put(quantity, ps, pe, eu, median, mean, lo, hi, ess, rhat, states, unit, conv=None):
        r = dict(base)
        r.update(quantity=quantity, period_start=ps, period_end=pe, end_used=eu, median=median, mean=mean, hpd95_lower=lo, hpd95_upper=hi,
                 ess_pooled=ess, rhat_split=rhat, meets_convergence_criterion=conv if conv is not None else crit(ess, rhat),
                 states_used=states, unit=unit)
        r.update(tail_cols)
        est.append(r)

    cols_est = ["analysis", "program", "setting", "genome_set", "genomes", "most_recent_collection_date", "quantity", "period_start",
                "period_end", "end_used", "median", "mean", "hpd95_lower", "hpd95_upper", "ess_pooled", "rhat_split",
                "meets_convergence_criterion", "chains_run", "chains_retained", "states_used", "unit"]
    if not an.kept:
        reason = f"{NA}: no chain retained"
        vP = {v[0]: v for v in variants}
        for key, pname, ps, pe, eu, t0, t1 in variants:
            for q, u in ((f"growth rate, {pname}", "per day"), (f"reproduction number, {pname}", "")) + ((("doubling time, P1", "days"),) if pname == "P1" else ()):
                put(q, ps, pe, eu, reason, reason, reason, reason, reason, reason, 0, u, conv=reason)
        for q, u in (("evolutionary rate", "substitutions per site and year"), ("date of the most recent common ancestor (date)", "date"),
                     ("date of the most recent common ancestor (decimal year)", "decimal year (year of 365 days)"), ("kappa", "")) \
                + ((("precision of the Skygrid (estimated)", ""),) if estimated else ()):
            put(q, "", "", "", reason, reason, reason, reason, reason, reason, 0, u, conv=reason)
        for q, k in (("P(R of P2 < R of P1)", None), ("P(R of P1 > 1)", "P1"), ("P(R of P2 > 1)", "P2"), ("P(R of P3 < 1)", "P3_own_end")):
            put(q, vP[k][2] if k else "", vP[k][3] if k else "", vP[k][4] if k else "fixed date", reason, reason, reason, reason, reason, reason, 0,
                "probability", conv=reason)
        pd.DataFrame(est)[cols_est].to_csv(f"{a.outdir}/estimates_{a.analysis}_{SUFFIX}.csv", index=False)
        print(f"{a.analysis}: {len(sched)} chains, 0 retained; only the chain table was written")
        return 0

    # ------------------------------------------------------------------ tables of the functions of analysis.py and rt_lib.py
    tabs = A.skygrid_tables(an)
    head = A.headline(an, tabs)
    diag = A.per_chain_diagnostics(an)
    head.to_csv(f"{a.outdir}/stored_headline_{a.analysis}_{SUFFIX}.csv", index=False)
    diag.to_csv(f"{a.outdir}/stored_chain_diagnostics_{a.analysis}_{SUFFIX}.csv", index=False)
    tabs["windows"].to_csv(f"{a.outdir}/stored_windows_{a.analysis}_{SUFFIX}.csv", index=False)
    tabs["rt"].to_csv(f"{a.outdir}/stored_Rt_pairs_{a.analysis}_{SUFFIX}.csv", index=False)
    tabs["intervals"].to_csv(f"{a.outdir}/stored_intervals_{a.analysis}_{SUFFIX}.csv", index=False)
    lp = an.lnNe()
    gen = an.gen
    n_states = int(sum(len(x) for x in lp))
    kept = an.kept

    # ------------------------------------------------------------------ periods: growth rate and reproduction number
    own = []                                                       # what this script computes itself, beside the functions of analysis.py and rt_lib.py
    rs_of, Wprim = {}, tabs["windows"][tabs["windows"]["is_primary_definition"] == True].set_index("window")
    for key, pname, ps, pe, eu, t0, t1 in variants:
        c = L.window_contrast(grid, t0, t1, "interp")
        rs = [x @ c for x in lp]                                   # growth rate per year of the log, one array per retained chain
        rs_of[key] = rs
        s = A.r_summary(rs, gen, key)
        allr = np.concatenate(rs)
        if pname in STORED_NAME and eu != "common":                # the same numbers as the table of windows of analysis.py
            w = Wprim.loc[STORED_NAME[pname]]
            assert s["R_median"] == w["R_median"] and s["R_hpd_lo"] == w["R_hpd_lo"] and s["R_hpd_hi"] == w["R_hpd_hi"] \
                and s["growth_per_day_median"] == w["growth_per_day_median"], (key, s["R_median"], w["R_median"])
        put(f"growth rate, {pname}", ps, pe, eu, s["growth_per_day_median"], float(np.mean(allr / L.YEAR_DAYS)), s["growth_per_day_hpd_lo"],
            s["growth_per_day_hpd_hi"], s["ess_pooled"], s["rhat_split"], s["n_samples"], "per day")
        put(f"reproduction number, {pname}", ps, pe, eu, s["R_median"], float(np.mean(gen.r2R(allr))), s["R_hpd_lo"], s["R_hpd_hi"],
            s["ess_pooled"], s["rhat_split"], s["n_samples"], "")
        m_st = L.summarise(rs, key, 1.0 / L.YEAR_DAYS)["mean"]
        own.append(dict(analysis=a.analysis, quantity=f"growth rate, {pname}", end_used=eu, statistic="mean", value_of_this_script=float(np.mean(allr / L.YEAR_DAYS)),
                        stored_function="rt_lib.summarise (mean of the series times 1/365)", value_of_the_stored_function=m_st,
                        difference=abs(float(np.mean(allr / L.YEAR_DAYS)) - m_st)))
        own.append(dict(analysis=a.analysis, quantity=f"reproduction number, {pname}", end_used=eu, statistic="P(R > 1)", value_of_this_script=float(np.mean(gen.r2R(allr) > 1)),
                        stored_function="analysis.r_summary (share of states with a growth rate above 0)", value_of_the_stored_function=s["P_R_gt_1"],
                        difference=abs(float(np.mean(gen.r2R(allr) > 1)) - s["P_R_gt_1"])))
        m_R = L.summarise([gen.r2R(r) for r in rs], key)
        own.append(dict(analysis=a.analysis, quantity=f"reproduction number, {pname}", end_used=eu, statistic="mean", value_of_this_script=float(np.mean(gen.r2R(allr))),
                        stored_function="rt_lib.summarise (mean of the series of reproduction numbers)", value_of_the_stored_function=m_R["mean"],
                        difference=abs(float(np.mean(gen.r2R(allr))) - m_R["mean"])))
        own.append(dict(analysis=a.analysis, quantity=f"reproduction number, {pname}", end_used=eu, statistic="ESS of the series of reproduction numbers (the table gives that of the growth rate)",
                        value_of_this_script=s["ess_pooled"], stored_function="rt_lib.summarise on the reproduction numbers", value_of_the_stored_function=m_R["ess_pooled"],
                        difference=abs(s["ess_pooled"] - m_R["ess_pooled"])))
        own.append(dict(analysis=a.analysis, quantity=f"reproduction number, {pname}", end_used=eu, statistic="split R-hat of the series of reproduction numbers (the table gives that of the growth rate)",
                        value_of_this_script=s["rhat_split"], stored_function="rt_lib.summarise on the reproduction numbers", value_of_the_stored_function=m_R["rhat_split"],
                        difference=abs(s["rhat_split"] - m_R["rhat_split"])))
        if pname == "P1":
            nd = lambda v: v if v == v else f"{NA}: not defined by the stored function (the growth rate or its lower HPD bound is not above 0)"
            put("doubling time, P1", ps, pe, eu, nd(s["doubling_time_d"]),
                f"{NA}: the stored function gives no mean (the doubling time is ln 2 / median growth rate)",
                nd(s["doubling_time_d_lo"]), nd(s["doubling_time_d_hi"]), s["ess_pooled"], s["rhat_split"], s["n_samples"], "days")
        meta.setdefault("samples_with_R_set_to_0", {})[key] = int(s["n_samples_R_set_to_0"])

    # ------------------------------------------------------------------ rate, tMRCA, kappa, precision
    s = L.summarise(an.series("clock.rate"), "clock rate", 1.0)
    put("evolutionary rate", "", "", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"], s["n_samples"],
        "substitutions per site and year")
    tm = A.tmrca_summary(an)
    mean_days = float(np.mean(np.concatenate(an.root_height()) * L.YEAR_DAYS))
    fdate = lambda d: str((an.last_tip - pd.Timedelta(days=int(round(d)))).date())
    m_st = L.summarise(an.root_height(), "root height", L.YEAR_DAYS)["mean"]
    own.append(dict(analysis=a.analysis, quantity="root height in days", end_used="", statistic="mean", value_of_this_script=mean_days,
                    stored_function="rt_lib.summarise (mean of the series times 365)", value_of_the_stored_function=m_st, difference=abs(mean_days - m_st)))
    put("date of the most recent common ancestor (date)", "", "", "", tm["tmrca_median"], fdate(mean_days), tm["tmrca_hpd_early"], tm["tmrca_hpd_late"],
        tm["ess_pooled"], tm["rhat_split"], tm["n_samples"], "date (day resolution, from the root height of the log)")
    put("date of the most recent common ancestor (decimal year)", "", "", "", tm["tmrca_decimal_median"], an.t_last - mean_days / L.YEAR_DAYS,
        tm["tmrca_decimal_early"], tm["tmrca_decimal_late"], tm["ess_pooled"], tm["rhat_split"], tm["n_samples"], "decimal year (year of 365 days)")
    s = L.summarise(an.series("kappa"), "kappa")
    put("kappa", "", "", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"], s["n_samples"], "")
    assert bool(tabs["precision_fixed"]) == (not estimated), "precision fixed / estimated does not agree with the name of the analysis"
    if estimated:
        s = L.summarise(an.series("skygrid.precision"), "Skygrid precision")
        put("precision of the Skygrid (estimated)", "", "", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"],
            s["n_samples"], "")
    else:
        meta["skygrid_precision_fixed_at"] = float(tabs["precision"])

    # ------------------------------------------------------------------ posterior probabilities
    Rc = {k: [gen.r2R(r) for r in v] for k, v in rs_of.items()}
    na_p = "not applicable (probability: the value is in the column 'mean')"

    def prob(quantity, ind, ps, pe, eu):
        ind = [np.asarray(x, float) for x in ind]
        allx = np.concatenate(ind)
        with np.errstate(all="ignore"):
            ess = L.ess_multichain(ind) if all(x.std() > 0 for x in ind) else float("nan")
            rh = L.rhat_split(ind) if all(x.std() > 0 for x in ind) else float("nan")
        put(quantity, ps, pe, eu, na_p, float(allx.mean()), na_p, na_p, ess, rh, int(len(allx)), "probability",
            conv=crit(ess, rh) if ess == ess and rh == rh else f"{NA}: the indicator is constant in at least one chain")

    vP = {v[0]: v for v in variants}
    prob("P(R of P2 < R of P1)", [x2 < x1 for x1, x2 in zip(Rc["P1"], Rc["P2"])], "", "", "fixed date")
    prob("P(R of P1 > 1)", [x > 1 for x in Rc["P1"]], vP["P1"][2], vP["P1"][3], "fixed date")
    prob("P(R of P2 > 1)", [x > 1 for x in Rc["P2"]], vP["P2"][2], vP["P2"][3], "fixed date")
    prob("P(R of P3 < 1)", [x < 1 for x in Rc["P3_own_end"]], vP["P3_own_end"][2], vP["P3_own_end"][3], "own")
    E = pd.DataFrame(est)[cols_est]
    E.to_csv(f"{a.outdir}/estimates_{a.analysis}_{SUFFIX}.csv", index=False)
    pd.DataFrame(own).to_csv(f"{a.outdir}/own_code_against_stored_{a.analysis}_{SUFFIX}.csv", index=False)

    # ------------------------------------------------------------------ intervals
    iv = tabs["intervals"]
    I = pd.DataFrame(dict(analysis=a.analysis, program=PROGRAM, setting=setting, genome_set=gset, interval=iv["interval"],
                          older_edge=iv["older_edge"], newer_edge=iv["newer_edge"], mid_point=iv["mid_point"],
                          mid_point_decimal=[float(grid.mids[i - 1]) if i < grid.K else np.nan for i in iv["interval"]],
                          genomes_collected=iv["n_genomes_collected"], ln_Ne_tau_median=iv["lnNe_median"],
                          ln_Ne_tau_hpd95_lower=iv["lnNe_hpd_lo"], ln_Ne_tau_hpd95_upper=iv["lnNe_hpd_hi"],
                          ess_pooled=iv["ess_pooled"], rhat_split=iv["rhat_split"],
                          meets_convergence_criterion=[crit(x, y) for x, y in zip(iv["ess_pooled"], iv["rhat_split"])],
                          status=iv["status"], chains_retained=len(kept), states_used=n_states, unit="ln(years)"))
    I.to_csv(f"{a.outdir}/intervals_{a.analysis}_{SUFFIX}.csv", index=False)

    # ------------------------------------------------------------------ trajectories
    T = []
    tb = dict(analysis=a.analysis, program=PROGRAM, setting=setting, genome_set=gset)
    for r in iv.itertuples():
        i = int(r.interval)
        T.append(dict(tb, series="Ne tau at the mid-point of an interval", date=(str(L.dec2date(grid.mids[i - 1])) if i < grid.K else ""),
                      date_decimal=(float(grid.mids[i - 1]) if i < grid.K else np.nan), interval=i, older_date="", newer_date="",
                      median=r.Ne_tau_years_median, hpd95_lower=r.Ne_tau_years_hpd_lo, hpd95_upper=r.Ne_tau_years_hpd_hi, unit="years",
                      ess_pooled=r.ess_pooled, rhat_split=r.rhat_split,
                      note=(r.status if i < grid.K else "interval 20 is open towards the past (everything before the cut-off); it has no mid-point")))
    allp = np.concatenate(lp)
    jan1 = L.decimal_year("2026-01-01")
    first_mid = min(m for m in grid.mids[:grid.K - 1] if m > jan1)
    d = pd.Timestamp(L.dec2date(first_mid))                       # calendar date of the first mid-point after 1 January 2026
    end_date = pd.Timestamp(L.dec2date(own_end))                  # calendar date of the end of the set
    meta.update(first_date_of_the_series_of_every_7th_day=str(d.date()), last_possible_date_of_that_series=str(end_date.date()))
    while d <= end_date:
        t = L.decimal_year(d)
        w = grid.weights(t, "interp")
        x = [q @ w for q in lp]
        Ne = np.exp(np.concatenate(x))
        lo, hi = L.hpd(np.log(Ne))
        T.append(dict(tb, series="Ne tau at every 7th day", date=str(d.date()), date_decimal=t, interval=grid.interval_of(t), older_date="",
                      newer_date="", median=float(np.median(Ne)), hpd95_lower=float(np.exp(lo)), hpd95_upper=float(np.exp(hi)), unit="years",
                      ess_pooled=L.ess_multichain(x), rhat_split=L.rhat_split(x),
                      note="ln(Ne tau) read by linear interpolation between interval mid-points"))
        d = d + pd.Timedelta(days=7)
    for r in tabs["rt"].itertuples():
        T.append(dict(tb, series="reproduction number between successive mid-points", date=r.time_boundary,
                      date_decimal=float(grid.t_anchor - r.interval_pair * grid.D), interval=f"{int(r.older_interval)} to {int(r.newer_interval)}",
                      older_date=r.older_interval_mid, newer_date=r.newer_interval_mid, median=r.R_median, hpd95_lower=r.R_hpd_lo,
                      hpd95_upper=r.R_hpd_hi, unit="", ess_pooled=r.ess_pooled, rhat_split=r.rhat_split, note=r.status))
        T.append(dict(tb, series="growth rate between successive mid-points", date=r.time_boundary,
                      date_decimal=float(grid.t_anchor - r.interval_pair * grid.D), interval=f"{int(r.older_interval)} to {int(r.newer_interval)}",
                      older_date=r.older_interval_mid, newer_date=r.newer_interval_mid, median=r.growth_per_day_median,
                      hpd95_lower=r.growth_per_day_hpd_lo, hpd95_upper=r.growth_per_day_hpd_hi, unit="per day", ess_pooled=r.ess_pooled,
                      rhat_split=r.rhat_split, note=r.status))
    T = pd.DataFrame(T)
    T["meets_convergence_criterion"] = [crit(x, y) for x, y in zip(T["ess_pooled"], T["rhat_split"])]
    T["chains_retained"] = len(kept); T["states_used"] = n_states
    T.to_csv(f"{a.outdir}/trajectories_{a.analysis}_{SUFFIX}.csv", index=False)

    # ------------------------------------------------------------------ posterior states of the retained chains
    parts = []
    for j, c in enumerate(kept):
        p = c.post
        rd = L.root_dates_from_height(p[c.rh_col].values, an.last_tip)
        df = pd.DataFrame(dict(seed=c.seed, state=p["state"].values, evolutionary_rate=p["clock.rate"].values,
                               root_height_years=p[c.rh_col].values, root_date=rd.astype(str), kappa=p["kappa"].values))
        for i, col in enumerate(c.lp_cols, start=1):
            df[f"ln_Ne_tau_{i}"] = p[col].values
        df["skygrid_precision"] = p["skygrid.precision"].values
        for key in rs_of:
            df[f"growth_rate_per_day_{key}"] = rs_of[key][j] / L.YEAR_DAYS
        for key in rs_of:
            df[f"reproduction_number_{key}"] = Rc[key][j]
        parts.append(df)
    P = pd.concat(parts, ignore_index=True)
    assert len(P) == n_states
    P.to_csv(f"{a.outdir}/posterior_{a.analysis}_{SUFFIX}.csv.gz", index=False, compression=dict(method="gzip", mtime=0))
    json.dump(meta, open(f"{a.outdir}/meta_{a.analysis}_{SUFFIX}.json", "w"), indent=1)
    print(f"{a.analysis}: {len(sched)} chains run, {len(all_chains)} read, {len(kept)} retained, {n_states} states used"
          + (" [TEST ON A SUBSET OF CHAINS]" if test_only else ""))
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc()
        sys.exit(2)
