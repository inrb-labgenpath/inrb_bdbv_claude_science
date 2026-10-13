#!/usr/bin/env python3
"""derive_beast_update.py (file label 20261007): criteria, estimates, intervals, trajectories,
diagnostics and posterior files of BEAST X on the set C00 of the snapshot of 7 October 2026, settings A and B, computed with the
functions of work/pipeline_ext/ (rt_lib.py, analysis.py, beast_analysis.py, beast_first_delivery.py). No estimator is
written anew here.

usage: python derive/derive_beast_update.py --root <working root> --tree <root of the working tree> --out <dir> --analyses a,b
       --periods <periods_of_the_update_20261007.csv> --periods-by-set <ends_by_set_20261007.csv> --schedule <chain_schedule_round_20261002.csv> [--test]

Constants of the functions of work/pipeline_ext/ and the values that are passed in their place: see the list PASSED below; it is written
to <out>/passed_constants.csv (every place is listed).
--test writes every file with the prefix TEST_ and does nothing else (a test of the code on short test chains; the values are not results).
"""
import argparse, csv, glob, gzip, hashlib, json, math, os, re, sys
import numpy as np
import pandas as pd

SUFFIX = "20261007"
PROGRAM = "BEAST X 10.5.0"
SETTING = {"A": "A: Skygrid, precision fixed at 4.06335", "B": "B: Skygrid, precision estimated (Gamma(0.001, 1000) prior)"}
LABEL = {"A": "A: Skygrid (20 parameters, cutoff 0.8 yr), precision fixed at 4.06335, HKY, strict clock",
         "B": "B: Skygrid (20 parameters, cutoff 0.8 yr), precision estimated (Gamma(0.001, 1000) prior), HKY, strict clock"}
P6E = "P6_to_END (defined on 7 October 2026 after the collection dates of the new snapshot were seen and before any estimate)"
# period -> key.  W1, W2, W3, W3b, X1, X2 are keys of the function beast_first_delivery.windows_of (X2 = 1 August to the
# mid-point of the most recent interval, which is the quantity P6_to_END).  W3c, P6o, P6c and P7o are in no function of work/pipeline_ext/: their growth rates are
# computed with the functions rt_lib.window_contrast (mode 'interp') and analysis.r_summary, as windows_of computes its own periods.
PERIOD_KEYS = [("P1", "W1", "fixed date"), ("P2", "W2", "fixed date"), ("P3", "W3", "own"), ("P3", "W3c", "common"), ("P4", "W3b", "fixed date"),
               ("P5", "X1", "fixed date"), ("P6", "P6o", "own"), ("P6", "P6c", "common"), ("P7", "P7o", "own"), (P6E, "X2", "own")]


def tag_of(pname, end):
    """name of the period in the column quantity"""
    return pname if (end == "fixed date" or pname in ("P7", P6E)) else f"{pname} ({end} end)"


def col_of(pname, end):
    """name of the period in the columns of the posterior file"""
    return "P6_to_END" if pname == P6E else pname if (end == "fixed date" or pname == "P7") else f"{pname}_{end}_end"
ZERO_COLS = ["analysis", "quantity", "period_start", "period_end", "end_used", "states_with_the_reproduction_number_set_to_0", "states_used", "share_of_the_states",
             "hpd95_lower_of_the_table_of_estimates", "period", "program", "setting", "genome_set", "column_of_the_posterior_file", "note"]
EST_COLS = ["analysis", "program", "setting", "genome_set", "genomes", "most_recent_collection_date", "quantity", "period_start", "period_end", "end_used",
            "median", "mean", "hpd95_lower", "hpd95_upper", "ess_pooled", "rhat_split", "meets_convergence_criterion", "chains_run", "chains_retained",
            "states_used", "unit"]
CHAIN_COLS = ["analysis", "program", "setting", "genome_set", "seed", "command_line", "started_utc", "ended_utc", "wall_time_s", "exit_code", "last_state_logged",
              "reached_last_step", "started_again", "root_date_rule_met", "first_state_root_date_rule", "density_criterion_met", "n_sampled_trees_evaluated",
              "retained", "why", "n_states_after_burn_in", "attempts", "latest_root_date_logged", "largest_difference_of_the_density_log_units"]
INT_COLS = ["analysis", "program", "setting", "genome_set", "interval", "older_edge", "newer_edge", "mid_point", "mid_point_decimal", "genomes_collected", "ln_Ne_tau_median",
            "ln_Ne_tau_hpd95_lower", "ln_Ne_tau_hpd95_upper", "ess_pooled", "rhat_split", "meets_convergence_criterion", "status", "chains_retained", "states_used", "unit"]
TRJ_COLS = ["analysis", "program", "setting", "genome_set", "series", "date", "date_decimal", "interval", "older_date", "newer_date", "median", "hpd95_lower", "hpd95_upper", "unit",
            "ess_pooled", "rhat_split", "note", "meets_convergence_criterion", "chains_retained", "states_used"]
PASSED = []


def stem_of(analysis):
    _, s, g = analysis.split("_")
    return ("beastA_skyfix_" if s == "A" else "beastB_skyest_") + g + "_" + SUFFIX


def setup(tree, root, periods_csv, periods_by_set_csv, schedule_csv, analyses):
    """imports the modules of work/pipeline_ext/ (they address each other relative to the root of the tree) and passes the values of this analysis"""
    os.chdir(tree)
    sys.path.insert(0, "work/pipeline_ext")
    import rt_lib as L, analysis as A, beast_analysis as B, beast_first_delivery as F
    per = pd.read_csv(periods_csv, dtype=str).set_index("period")            # periods_of_the_update_20261007.csv
    pbs = pd.read_csv(periods_by_set_csv, dtype=str).set_index("set")        # ends_by_set_20261007.csv
    sched = pd.read_csv(schedule_csv, dtype=str)
    sched = sched[sched["track"] == "BEAST sets"]
    sched = sched[sched["analysis"].isin(analyses)]
    assert len(sched) == 8 * len(analyses) or os.environ.get('FILES_OF_THE_TEST_CHAIN') == 'yes', len(sched)
    # ---- periods: the dates fixed in analysis.py are those of the periods table (checked, not changed)
    assert A.W_DATES["W1"] == (per.loc["P1", "start"], per.loc["P1", "end"]) and A.W_DATES["W2"] == (per.loc["P2", "start"], per.loc["P2", "end"])
    assert A.W_DATES["W3"] == (per.loc["P3 (own end)", "start"], None) and per.loc["P3 (own end)", "end"].startswith("END") and per.loc["P3 (common end)", "start"] == "2026-06-15"
    assert A.W_DATES["W3b"] == (per.loc["P4", "start"], per.loc["P4", "end"])
    assert (per.loc["P5", "start"], per.loc["P5", "end"]) == ("2026-06-15", "2026-08-01")
    assert per.loc["P6 (own end)", "start"] == "2026-08-01" and per.loc["P6 (own end)", "end"].startswith("2026-09-01 if 1 September 2026 does not lie beyond END; otherwise END")
    assert per.loc["P6 (common end)", "start"] == "2026-08-01" and (per.loc["P7", "start"], per.loc["P7", "end"]) == ("2026-09-01", "END")
    assert (per.loc[P6E, "start"], per.loc[P6E, "end"]) == ("2026-08-01", "own END")
    PASSED.append(("analysis.py", "W_DATES (W1, W2, W3, W3b)", "periods of the earlier rounds", "unchanged: equal to P1, P2, P3 (own end) and P4 of the table of the periods of the update (checked by program)"))
    PASSED.append(("beast_first_delivery.py", "windows_of: X1 (15 June - 1 August), X2 (1 August - mid-point of the most recent interval)", "exploratory periods of the earlier rounds",
                   "unchanged: X1 is P5; X2 is the quantity P6_to_END of this round (1 August to the own END), NOT P6 (checked against the table of the periods of the update)"))
    PASSED.append(("beast_first_delivery.py", "windows_of: W3p, X2p (periods with the end date of another analysis, given only if the own end is more than a day later than that date)",
                   "rows with the end date of the Delphy primary analysis of the first round; in the round of 2 October the rows with the common end",
                   "not given by the stored function in this round (the own end of the set is the common end); the rows with the common end are computed with rt_lib.window_contrast "
                   "(mode 'interp') and analysis.r_summary, the two stored functions that windows_of itself calls"))
    PASSED.append(("no stored script", "P6 (1 August - 1 September if 1 September does not lie beyond END, otherwise END), with the own end and with the common end; P7 (1 September - END)",
                   "not in a stored script (in the round of 2 October P6 ended at END)",
                   "computed with rt_lib.window_contrast (mode 'interp') and analysis.r_summary; not evaluated under the two conditions of the stored function windows_of "
                   "(end after the last collection date; shorter than 1 day), P7 not evaluated if shorter than beast_analysis.MIN_WINDOW_DAYS"))
    # ---- ends: mid-point of the most recent interval of every set; the common end is the earliest of them (table of the ends by set)
    ends = {s_: float(L.Grid.from_cutoff(20, 0.8, pd.Timestamp(pbs.loc[s_, "most_recent_collection_date"])).mids[0]) for s_ in pbs.index}
    for s_ in pbs.index:
        assert str(L.dec2date(ends[s_])) == pbs.loc[s_, "own_end_mid_point_of_the_most_recent_interval"], (s_, ends[s_])
    t_common = min(ends.values())
    assert set(pbs["common_end_earliest_of_the_four_own_ends"]) == {str(L.dec2date(t_common))}, t_common
    old = B.T_END_PRIMARY
    B.T_END_PRIMARY = t_common
    PASSED.append(("beast_analysis.py", "T_END_PRIMARY", f"mid-point of the most recent interval for the most recent collection date 2026-09-05 ({old!r})",
                   f"common end of this round: the earliest of the own ends of the four sets of the table of the ends by set, {t_common!r} ({L.dec2date(t_common)}), computed with the "
                   f"stored function rt_lib.Grid.from_cutoff(20, 0.8, most recent collection date).mids[0]"))
    oldsfx = B.SFX
    B.SFX = SUFFIX; F.SFX = SUFFIX
    PASSED.append(("beast_analysis.py, beast_first_delivery.py", "SFX", oldsfx, SUFFIX))
    if os.environ.get('FILES_OF_THE_TEST_CHAIN') == 'yes':      # test of the code on the files of the short test chain only: the chain length is that of the schedule of the test
        B.CHAIN_LENGTH = int(sched["steps"].iloc[0])
    assert B.CHAIN_LENGTH == int(sched["steps"].iloc[0]) and set(sched["steps"]) == {str(B.CHAIN_LENGTH)} and abs(B.BURNIN - float(sched["burn_in_share"].iloc[0])) < 1e-12
    PASSED.append(("beast_analysis.py", "CHAIN_LENGTH, BURNIN", "40,000,000 states; 0.30", "unchanged: equal to the schedule of chains (checked by program)"))
    B.MODELS = {a: dict(xml=stem_of(a), short=a.split("_")[1], dataset=a.split("_")[2], kind="skygrid", label=LABEL[a.split("_")[1]], delphy="", delphy_alt="", like_for_like="")
                for a in analyses}
    PASSED.append(("beast_analysis.py", "MODELS", "the 5 analyses of the earlier round (names of the XML files, genome set, text of the model)", "the analyses of this round: " + ", ".join(
        f"{a} = {stem_of(a)}.xml on {a.split('_')[2]}" for a in analyses)))
    B.ALN = {"C00": os.path.join(root, "aln", f"C00_{SUFFIX}.fasta")}
    PASSED.append(("beast_analysis.py", "ALN", "inputs/aln/primary_20260924.fasta, inputs/aln/inrbRetained_20260924.fasta", "aln/C00_20261007.fasta (SHA-256 checked against the stated value before the XML files were written)"))
    B.xml_params = lambda m: json.load(open(os.path.join(root, "xml", f"beast_xml_parameters_{B.MODELS[m]['xml']}.json")))
    PASSED.append(("beast_analysis.py", "xml_params (reads work/beast/xml/beast_xml_parameters_<name>.json)", "parameter files of the XML files of the earlier round",
                   "the same reader on xml/beast_xml_parameters_<name>.json of this round (written by beast_make_xml.py with the XML files)"))
    F.SHORT = {a: a.split("_")[1] + " on " + a.split("_")[2] for a in analyses}
    PASSED.append(("beast_first_delivery.py", "SHORT", "short names of the 5 analyses of the earlier round", "short names of the analyses of this round"))
    g = L.GenTime(); assert (g.mean_d, g.sd_d) == (15.3, 9.3)
    PASSED.append(("rt_lib.py", "GenTime (mean 15.3 d, SD 9.3 d; year of 365.25 d for the generation time, 365 d for the dates)", "generation time of the article", "unchanged"))
    PASSED.append(("beast_analysis.py", "MIN_WINDOW_DAYS = 14", "shortest exploratory period that is evaluated", "unchanged; used for P7 (1 September to END), which is shorter and is not evaluated"))
    return L, A, B, F, per, pbs, sched, t_common


def hpd_date(L, last_tip, days):
    return str((pd.Timestamp(last_tip) - pd.Timedelta(days=int(round(days)))).date())


def crit(ess, rhat):
    if ess != ess or rhat != rhat:
        return "not assessed"
    return "yes" if (ess >= 200 and rhat <= 1.05) else "no"


def read_ledger(root):
    led = pd.read_csv(os.path.join(root, "sched", "ledger.csv"), dtype=str, keep_default_na=False)
    return led


def read_txt(p):
    return open(p).read().strip() if os.path.exists(p) else ""


def judged_by_files(part, stem, steps, tree_every):
    """(log reaches the last state, tree file complete, standard output ends with the closing lines of BEAST X) of one part of a chain"""
    lg = os.path.join(part, stem + ".log"); tr = os.path.join(part, stem + ".trees"); so = os.path.join(part, "stdout.txt")
    last = None
    if os.path.exists(lg):
        for l in open(lg, errors="replace"):
            f = l.split("\t", 1)
            if f[0].isdigit():
                last = int(f[0])
    trees, end = [], False
    if os.path.exists(tr):
        for l in open(tr, errors="replace"):
            m = re.match(r"\s*tree STATE_(\d+)\b", l)
            if m:
                trees.append(int(m.group(1)))
            elif re.match(r"\s*End;", l, re.I):
                end = True
    tail = open(so, errors="replace").read()[-6000:] if os.path.exists(so) else ""
    closing = bool(re.search(r"Operator analysis", tail)) and bool(re.search(r"^\s*[0-9.]+ (seconds|minutes|hours|days)\s*$", tail.strip().splitlines()[-1] if tail.strip() else "", re.M))
    return dict(log_last_state=last, log_ok=(last == steps), tree_ok=(bool(trees) and trees[-1] == steps and end and all(b - a == tree_every for a, b in zip(trees, trees[1:]))),
                n_trees=len(trees), closing_lines=closing)


def chain_rows(an, analysis, root, led, B):
    rows = []
    s_, g_ = analysis.split("_")[1], analysis.split("_")[2]
    stem = stem_of(analysis)
    jr = {}
    for c in an.chains:
        p = os.path.join(os.path.dirname(c.dir), "join_report.csv")
        if os.path.exists(p) and not jr:
            jr = {r["chain"]: r for r in csv.DictReader(open(p))}
    for c in an.chains:
        lab = os.path.basename(c.dir)
        rd = os.path.join(root, "runs", analysis, lab)                      # the chain directory (what the program and its wrapper wrote)
        parts = [rd]; k = 1
        while os.path.isdir(os.path.join(rd, f"continuation_{k}")):
            parts.append(os.path.join(rd, f"continuation_{k}")); k += 1
        if c.info.get("cut_at") is not None:
            parts = parts[:1]                                # summary cut at the saved states: only the first part of a cut chain is described (command, times, exit code)
        l = led[led["chain"] == lab]
        again = []
        for _, r in l.iterrows():
            if r["event"] == "started_again_from_state_0":
                prev = l[(l["utc"] < r["utc"]) & l["event"].isin(["started", "started_again_from_state_0", "continued_from_saved_state"])]
                it = l[(l["utc"] <= r["utc"]) & (l["event"] == "interrupted")]
                again.append(f"yes: the start of {prev['utc'].iloc[-1] if len(prev) else 'unknown time'} was ended from outside (last state logged "
                             f"{it['last_state_logged'].iloc[-1] if len(it) else 'not recorded'}); started again from state 0 at {r['utc']} with the same command and seed")
        if lab in jr and jr[lab].get("continued", "no") not in ("no", ""):
            again.append(jr[lab]["continued"])
        for kk, pth in enumerate(parts[1:], 1):          # every continuation: seed, the two log-likelihoods, times, where the saved state came from (rule of 8 October 2026, point 5)
            ci = read_txt(os.path.join(pth, "continuation_info.txt"))
            def g_(rx, ci=ci):
                m_ = re.search(rx, ci)
                return m_.group(1).strip() if m_ else "not recorded"
            seed_k = g_(r"seed of the continuation (\d+)"); lnl_k = g_(r"log-likelihoods of the line of BEAST X: ([^\n]+)"); dif_k = g_(r"difference of the log-likelihoods (\S+)")
            src_k = g_(r"saved state tried: ([^\n]+?) \(SHA-256")
            again.append("continuation %d: seed %s; log-likelihoods of the line of BEAST X: %s; difference %s; started %s, ended %s; saved state: %s" % (
                kk, seed_k, lnl_k, dif_k, read_txt(os.path.join(pth, "started.txt")) or "not recorded", read_txt(os.path.join(pth, "finished.txt")) or "not recorded", src_k))
        cmds = [" ".join(read_txt(os.path.join(p, "command.txt")).split()) for p in parts]
        st = read_txt(os.path.join(rd, "started.txt")); en = read_txt(os.path.join(parts[-1], "finished.txt")); ec = read_txt(os.path.join(parts[-1], "exit_code.txt"))
        wall, wall_known = 0.0, True
        for p in parts:
            a_, b_ = read_txt(os.path.join(p, "started.txt")), read_txt(os.path.join(p, "finished.txt"))
            if a_ and b_:
                wall += (pd.Timestamp(b_) - pd.Timestamp(a_)).total_seconds()
            else:
                wall_known = False
        jf = judged_by_files(parts[-1], stem, B.CHAIN_LENGTH, int(an.P["tree_every"]))
        why = ("regular: reached the last step, root-date rule not met" if (c.retained and c.info["complete"]) else c.why)
        if (not c.info["complete"] and c.info.get("cut_at") is None
                and read_txt(os.path.join(parts[-1], "exit_code.txt")) not in ("", "129", "130", "131", "137", "143")):      # not a signal from outside
            why += (f"; the program ended the chain with exit code {read_txt(os.path.join(parts[-1], 'exit_code.txt'))} before its last state: "
                    "the chain is recorded and discarded, it is not started again and not replaced (rule of the update)")
        if len(parts) > 1:
            why += f"; continued chain ({len(parts) - 1} continuation(s)); log and tree file joined from {len(parts)} parts"
        if not ec and not jf["log_ok"]:
            ec = "none: the chain has not ended regularly (its log has not reached the last state)"
        elif not ec:
            ec = "exit code not recorded"
            why += ("; judged by its files: last state of the log " + ("equals" if jf["log_ok"] else "DOES NOT equal") + " the length of the chain, tree file "
                    + ("complete" if jf["tree_ok"] else "NOT complete") + f" ({jf['n_trees']} trees), standard output of BEAST X " + ("ends" if jf["closing_lines"] else "DOES NOT end")
                    + " with its closing lines")
        rows.append(dict(analysis=analysis, program=PROGRAM, setting=SETTING[s_], genome_set=g_, seed=c.seed, command_line=" || ".join(x for x in cmds if x), started_utc=st, ended_utc=en,
                         wall_time_s=(wall if wall_known else ""), exit_code=ec, last_state_logged=c.info["last_state"],
                         reached_last_step=("yes" if c.info["complete"] else "no" if c.info.get("cut_at") is None else f"cut at the saved state {c.info['cut_at']}, which is taken as the last step of the cut chain"),
                         started_again=("; ".join(again) if again else "no"), root_date_rule_met="yes" if c.info["collapsed"] else "no",
                         first_state_root_date_rule=("" if not c.info["collapsed"] else int(c.info["first_collapsed_state"])), density_criterion_met="not applicable",
                         n_sampled_trees_evaluated="not applicable", retained="yes" if c.retained else "no", why=why, n_states_after_burn_in=len(c.post),
                         attempts=int(l["event"].isin(["started", "started_again_from_state_0", "continued_from_saved_state"]).sum()), latest_root_date_logged=c.info["latest_root_date"],
                         largest_difference_of_the_density_log_units="not applicable"))
    return rows


def derive(analysis, root, out, L, A, B, F, t_common, test=False, exclude=(), tag="", inputs=None, cut=None):
    s_, g_ = analysis.split("_")[1], analysis.split("_")[2]
    pre = "TEST_" if test else ""
    indir = os.path.join(inputs, analysis) if inputs else os.path.join(root, "runs", analysis)
    an = B.load_model(analysis, indir, False)        # the loader accepts no chain that has not reached its last step, with --test or without
    assert len(an.chains) > 0
    if cut:
        # summary cut at the saved states: the log of a continued chain ends at the state of the saved state from which
        # the chain was continued (tools/cut_continued.py). The class Chain discards a chain that has not reached state CHAIN_LENGTH; for a cut chain the saved state is taken as its
        # last step, and the chain is retained or not by the rule on the root date alone, as the class Chain decided it. The burn-in is that of the class Chain
        # (the first int(BURNIN x logged states) of the log as read). Chains that are not in the list are treated as in every other summary.
        for c in an.chains:
            lab = os.path.basename(c.dir)
            if lab in cut:
                assert int(c.info["last_state"]) == int(cut[lab]), (lab, c.info["last_state"], cut[lab])
                assert not c.info["complete"] and c.why.startswith("incomplete"), (lab, c.why)
                c.retained = not bool(c.info["collapsed"])
                c.why = (("regular" if c.retained else "entered the degenerate state (root date >= first tip - 1 d)")
                         + f"; chain CUT at the saved state {int(cut[lab])} from which it was continued (taken as the last step of the cut chain; nothing that the continuation logged is used)")
                c.info["cut_at"] = int(cut[lab])
                # burn-in: as for every chain, the first 30 % of
                # the states of a WHOLE chain (the first int(BURNIN x logged states of a whole chain) logged states: 600 of 2001, the states before 12,000,000) - not the first 30 % of
                # the logged states of the cut chain, which is what the class Chain does with a log as it reads it.
                nb_whole = int((B.CHAIN_LENGTH // int(an.P["log_every"]) + 1) * B.BURNIN)
                c.nb = nb_whole
                c.post = c.full.iloc[nb_whole:].reset_index(drop=True)
                c.info["n_post_burnin"] = len(c.post)
                if int(cut[lab]) <= nb_whole * int(an.P["log_every"]):
                    c.retained = False
                    c.why = (f"left out of the summary cut at the saved states: the chain was cut at the saved state {int(cut[lab])}, within the burn-in (at or before state "
                             f"{nb_whole * int(an.P['log_every'])}); nothing that the continuation logged is used")
            else:
                assert c.info["complete"], (lab, "a chain that is not cut has not reached its last step")
        assert sum(1 for c in an.chains if os.path.basename(c.dir) in cut) == len(cut), "a chain of the list of the cut chains was not found"
        if not any(c.info.get("cut_at") is not None and c.retained for c in an.chains):
            # no continued chain holds a state after the burn-in and before its saved state: there is no comparison
            res_none = dict(chains=pd.DataFrame(chain_rows(an, analysis, root, read_ledger(root), B))[CHAIN_COLS], an=an,
                            no_comparison="no comparison: every chain was cut within the burn-in" if all(c.info.get("cut_at") is not None for c in an.chains)
                            else "no comparison: every continued chain was cut within the burn-in")
            return res_none
    for c in an.chains:
        if os.path.basename(c.dir) in exclude and c.retained:
            c.retained = False; c.why = "left out of this summary: continued chain (second summary)"
    P = an.P
    led = read_ledger(root)
    CH = pd.DataFrame(chain_rows(an, analysis, root, led, B))[CHAIN_COLS]
    base = dict(analysis=analysis, program=PROGRAM, setting=SETTING[s_], genome_set=g_, genomes=an.n_tips, most_recent_collection_date=str(an.last_tip.date()),
                chains_run=len(an.chains), chains_retained=len(an.kept))
    res = dict(chains=CH, an=an)
    if not an.kept:
        return res
    n_used = int(sum(len(c.post) for c in an.kept))
    base["states_used"] = n_used
    g = an.grid
    S = B.series_from_analysis(an)
    tabs = A.skygrid_tables(an)
    W = {w["key"]: w for w in F.windows_of(S)}
    Q = B.quantities(S)
    rows = []

    def add(q, unit, med, mean, lo, hi, ess, rhat, ps="", pe="", end="", cr=None):
        rows.append(dict(base, quantity=q, period_start=ps, period_end=pe, end_used=end, median=med, mean=mean, hpd95_lower=lo, hpd95_upper=hi, ess_pooled=ess, rhat_split=rhat,
                         meets_convergence_criterion=(cr if cr is not None else crit(ess, rhat)), unit=unit))

    ser_r, ser_R = {}, {}
    zero = []          # states in which the conversion sets the reproduction number to 0
    t_end = float(g.mids[0]); t_sep1 = L.decimal_year("2026-09-01"); t_aug1 = L.decimal_year("2026-08-01"); t_jun15 = L.decimal_year("2026-06-15")
    FIXED_END = {"W1": "2026-05-15", "W2": "2026-06-15", "W3b": "2026-08-15", "X1": "2026-08-01"}
    res["P7_days_decimal"] = (t_end - t_sep1) * L.YEAR_DAYS
    res["P7_days_by_dates"] = int((pd.Timestamp(str(L.dec2date(t_end))) - pd.Timestamp("2026-09-01")).days)
    res["own_end"], res["common_end"] = (t_end, str(L.dec2date(t_end))), (float(B.T_END_PRIMARY), str(L.dec2date(B.T_END_PRIMARY)))

    def window(key):
        """start (text), end (text), t0, t1, the period of the function windows_of or None, reason if the period is not evaluated"""
        if key in W:                                         # W1, W2, W3, W3b, X1, X2: periods of windows_of
            w = W[key]
            t1 = t_end if key in ("W3", "X2") else L.decimal_year(FIXED_END[key])
            return w["start"], w["end"], L.decimal_year(w["start"]), t1, w, (w["note"] if w["rs"] is None else "")
        if key == "W3c":                                     # P3 with the common end
            a_, t0, t1 = "2026-06-15", t_jun15, float(B.T_END_PRIMARY)
        elif key in ("P6o", "P6c"):                          # P6: 1 August to 1 September if 1 September does not lie beyond END (own or common), otherwise to END
            e_ = t_end if key == "P6o" else float(B.T_END_PRIMARY)
            a_, t0, t1 = "2026-08-01", t_aug1, (t_sep1 if t_sep1 <= e_ else e_)
        elif key == "P7o":                                   # P7: 1 September to END
            a_, t0, t1 = "2026-09-01", t_sep1, t_end
        else:
            raise KeyError(key)
        days = (t1 - t0) * L.YEAR_DAYS
        why = ""
        if t1 > S.t_last + 1e-9:                            # the two conditions under which windows_of does not evaluate a period
            why = f"the period ends after the last collection date ({S.last_tip.date()})"
        elif days < 1.0:
            why = f"the period has no length (it would run from {a_} to {L.dec2date(t1)})"
        elif key == "P7o" and days < B.MIN_WINDOW_DAYS:
            why = (f"P7 runs from 2026-09-01 to END ({L.dec2date(t_end)}, the mid-point of the most recent interval) and is {days:.2f} days long by the convention of the stored code "
                   f"(difference of the decimal dates x 365; {res['P7_days_by_dates']} days counted by dates); a period shorter than {B.MIN_WINDOW_DAYS:.0f} days is not evaluated")
        return a_, str(L.dec2date(t1)), t0, t1, None, why

    for pname, key, end in PERIOD_KEYS:
        tagp = tag_of(pname, end)
        zb = dict(analysis=analysis, quantity=f"reproduction number {tagp}", end_used=end, period=pname, program=PROGRAM, setting=SETTING[s_], genome_set=g_,
                  column_of_the_posterior_file="R_" + col_of(pname, end))
        ps, pe, t0, t1, w, why = window(key)
        if why:
            for q, u in ((f"growth rate {tagp}", "per day"), (f"reproduction number {tagp}", "")):
                add(q, u, "", "", "", "", "", "", ps, pe, end, cr="not assessed: " + why)
            zero.append(dict(zb, period_start=ps, period_end=pe, states_with_the_reproduction_number_set_to_0="", states_used="", share_of_the_states="", hpd95_lower_of_the_table_of_estimates="",
                             column_of_the_posterior_file="", note="not assessed: " + why))
            continue
        c = L.window_contrast(g, t0, t1, "interp")
        cr_ = [x @ c for x in S.lnNe]                              # growth rate per year, one series per retained chain
        rs = A.r_summary(cr_, S.gen, key)
        if w is not None:
            for k in ("R_median", "R_hpd_lo", "R_hpd_hi", "growth_per_day_median", "ess_pooled", "rhat_split"):
                assert rs[k] == w["rs"][k] or (rs[k] != rs[k] and w["rs"][k] != w["rs"][k]), (analysis, key, k)       # the series are those of windows_of
        sr = L.summarise(cr_, key, 1.0 / L.YEAR_DAYS)
        sR = L.summarise([S.gen.r2R(x) for x in cr_], key)
        ser_r[(pname, end)] = cr_; ser_R[(pname, end)] = [S.gen.r2R(x) for x in cr_]
        assert rs["n_samples"] == n_used and rs["n_samples_R_set_to_0"] == int(sum(int((x == 0).sum()) for x in ser_R[(pname, end)])), (analysis, key)
        zero.append(dict(zb, period_start=ps, period_end=pe, states_with_the_reproduction_number_set_to_0=int(rs["n_samples_R_set_to_0"]), states_used=int(rs["n_samples"]),
                         share_of_the_states=rs["n_samples_R_set_to_0"] / rs["n_samples"], hpd95_lower_of_the_table_of_estimates=rs["R_hpd_lo"],
                         note="count of the stored function (n_samples_R_set_to_0 of analysis.r_summary; rt_lib.GenTime.n_undefined: 1 + r / B <= 0)"))
        add(f"growth rate {tagp}", "per day", rs["growth_per_day_median"], sr["mean"], rs["growth_per_day_hpd_lo"], rs["growth_per_day_hpd_hi"], rs["ess_pooled"], rs["rhat_split"], ps, pe, end)
        add(f"reproduction number {tagp}", "", rs["R_median"], sR["mean"], rs["R_hpd_lo"], rs["R_hpd_hi"], rs["ess_pooled"], rs["rhat_split"], ps, pe, end)
        if pname == "P1":
            add("doubling time P1", "days", rs["doubling_time_d"], "", rs["doubling_time_d_lo"], rs["doubling_time_d_hi"], rs["ess_pooled"], rs["rhat_split"], ps, pe, end)
    s = L.summarise(S.clock, "clock", 1.0)
    assert abs(s["median"] * 1e3 - Q["clock rate (1e-3 subst/site/yr)"]["median"]) < 1e-12
    add("evolutionary rate", "substitutions per site and year", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"])
    tm = A.tmrca_summary(an)
    srh = L.summarise(S.rh, "rh", 1.0)
    add("date of the most recent common ancestor (date)", "date; hpd95_lower is the earlier date", tm["tmrca_median"], hpd_date(L, an.last_tip, srh["mean"] * L.YEAR_DAYS), tm["tmrca_hpd_early"], tm["tmrca_hpd_late"],
        tm["ess_pooled"], tm["rhat_split"])
    add("date of the most recent common ancestor (decimal year)", "decimal year (year of 365 days); hpd95_lower is the earlier date", tm["tmrca_decimal_median"], an.t_last - srh["mean"], tm["tmrca_decimal_early"],
        tm["tmrca_decimal_late"], tm["ess_pooled"], tm["rhat_split"])
    s = L.summarise(an.series("kappa"), "kappa")
    add("kappa", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"])
    prec_estimated = not tabs["precision_fixed"]
    assert prec_estimated == (s_ == "B"), (analysis, tabs["precision"], tabs["precision_fixed"])
    if prec_estimated:
        s = L.summarise(an.series("skygrid.precision"), "precision")
        add("precision of the Skygrid", "", s["median"], s["mean"], s["hpd_lo"], s["hpd_hi"], s["ess_pooled"], s["rhat_split"])
    else:
        assert abs(tabs["precision"] - 4.06335) < 1e-12
    # posterior probabilities: share of the states used (as P_R_gt_1 of r_summary: share of the states with a growth rate above 0)
    r1, r2, r3 = np.concatenate(ser_r[("P1", "fixed date")]), np.concatenate(ser_r[("P2", "fixed date")]), np.concatenate(ser_r[("P3", "own")])
    R1, R2 = np.concatenate(ser_R[("P1", "fixed date")]), np.concatenate(ser_R[("P2", "fixed date")])
    assert float((r1 > 0).mean()) == W["W1"]["rs"]["P_R_gt_1"] and float((r2 > 0).mean()) == W["W2"]["rs"]["P_R_gt_1"]
    na = "not applicable (share of the states used; the criterion is that of the growth rates)"
    add("P(R of P2 < R of P1)", "probability (share of the states used), in the column mean", "", float((R2 < R1).mean()), "", "", "", "", cr=na)
    add("P(R of P1 > 1)", "probability (share of the states used), in the column mean", "", float((r1 > 0).mean()), "", "", "", "", W["W1"]["start"], W["W1"]["end"], "fixed date", cr=na)
    add("P(R of P2 > 1)", "probability (share of the states used), in the column mean", "", float((r2 > 0).mean()), "", "", "", "", W["W2"]["start"], W["W2"]["end"], "fixed date", cr=na)
    add("P(R of P3 < 1) (own end)", "probability (share of the states used), in the column mean", "", float((r3 < 0).mean()), "", "", "", "", W["W3"]["start"], W["W3"]["end"], "own", cr=na)
    EST = pd.DataFrame(rows)[EST_COLS]
    # ---------------- intervals
    iv = tabs["intervals"]
    INT = pd.DataFrame(dict(analysis=analysis, program=PROGRAM, setting=SETTING[s_], genome_set=g_, interval=iv["interval"], older_edge=iv["older_edge"], newer_edge=iv["newer_edge"],
                            mid_point=iv["mid_point"], mid_point_decimal=[(float(g.mids[int(i) - 1]) if int(i) < g.K else "") for i in iv["interval"]],
                            genomes_collected=iv["n_genomes_collected"], ln_Ne_tau_median=iv["lnNe_median"], ln_Ne_tau_hpd95_lower=iv["lnNe_hpd_lo"],
                            ln_Ne_tau_hpd95_upper=iv["lnNe_hpd_hi"], ess_pooled=iv["ess_pooled"], rhat_split=iv["rhat_split"],
                            meets_convergence_criterion=[crit(a, b) for a, b in zip(iv["ess_pooled"], iv["rhat_split"])], status=iv["status"], chains_retained=len(an.kept), states_used=n_used,
                            unit="ln(years)"))[INT_COLS]
    # ---------------- trajectories
    tr = []
    tb = dict(analysis=analysis, program=PROGRAM, setting=SETTING[s_], genome_set=g_, chains_retained=len(an.kept), states_used=n_used)
    mids_after = [i for i in range(g.K - 1) if g.mids[i] > 2026.0]
    first_mid = g.mids[max(mids_after)]
    for i in range(g.K - 1, 0, -1):                                # interval mid-points, oldest first (interval K is open and has no mid-point)
        r = iv[iv["interval"] == i].iloc[0]
        tr.append(dict(tb, series="Ne tau at the mid-point of an interval", date=str(L.dec2date(g.mids[i - 1])), date_decimal=float(g.mids[i - 1]), interval=i, older_date="", newer_date="",
                       median=float(np.exp(r["lnNe_median"])), hpd95_lower=float(np.exp(r["lnNe_hpd_lo"])), hpd95_upper=float(np.exp(r["lnNe_hpd_hi"])), unit="years", ess_pooled=r["ess_pooled"],
                       rhat_split=r["rhat_split"], note=r["status"], meets_convergence_criterion=crit(r["ess_pooled"], r["rhat_split"])))
    day = pd.Timestamp(L.dec2date(first_mid))
    while L.decimal_year(day) <= float(g.mids[0]) + 1e-12:
        t = L.decimal_year(day)
        wv = g.weights(t, "interp")
        sv = L.summarise([x @ wv for x in S.lnNe], "lnNe(t)")
        tr.append(dict(tb, series="Ne tau at every 7th day", date=str(day.date()), date_decimal=t, interval=g.interval_of(t), older_date="", newer_date="", median=float(np.exp(sv["median"])),
                       hpd95_lower=float(np.exp(sv["hpd_lo"])), hpd95_upper=float(np.exp(sv["hpd_hi"])), unit="years", ess_pooled=sv["ess_pooled"], rhat_split=sv["rhat_split"],
                       note="ln(Ne tau) read by linear interpolation between interval mid-points; first date: the first mid-point after 1 January 2026",
                       meets_convergence_criterion=crit(sv["ess_pooled"], sv["rhat_split"])))
        day += pd.Timedelta(days=7)
    rt = tabs["rt"]
    for ser, cols, unit in (("growth rate between successive mid-points", ("growth_per_day_median", "growth_per_day_hpd_lo", "growth_per_day_hpd_hi"), "per day"),
                            ("reproduction number between successive mid-points", ("R_median", "R_hpd_lo", "R_hpd_hi"), "")):
        for _, r in rt.sort_values("interval_pair", ascending=False).iterrows():

            nt = "date: edge between the two intervals; interval: the newer of the two" + ("; ESS and R-hat are those of the growth rate" if unit == "" else "") + (("; " + r["status"]) if r["status"] else "")
            if not r["older_interval_mid"]:
                nt += "; the older of the two intervals is the oldest interval, which is open towards the past and has no mid-point (older_date is empty)"
            tr.append(dict(tb, series=ser, date=r["time_boundary"], date_decimal=float(g.t_anchor - int(r["interval_pair"]) * g.D), interval=int(r["newer_interval"]),
                           older_date=r["older_interval_mid"], newer_date=r["newer_interval_mid"], median=r[cols[0]], hpd95_lower=r[cols[1]], hpd95_upper=r[cols[2]], unit=unit,
                           ess_pooled=r["ess_pooled"], rhat_split=r["rhat_split"], note=nt, meets_convergence_criterion=crit(r["ess_pooled"], r["rhat_split"])))
    TRJ = pd.DataFrame(tr)[TRJ_COLS]
    # ---------------- diagnostics per chain and pooled (functions of work/pipeline_ext/; periods P5, P6 and the common end added with the same functions)
    D = pd.concat([A.per_chain_diagnostics(an), B.extra_diagnostics(an)], ignore_index=True)
    extra = []
    for (pname, end), key in (((p, e), k) for p, k, e in PERIOD_KEYS):
        if key in ("W1", "W2", "W3", "W3b") or (pname, end) not in ser_r:
            continue
        ps, pe, t0, t1, w, why = window(key)
        c = L.window_contrast(g, t0, t1, "interp")
        allc = [x @ c for x in an.lnNe(an.chains)]
        name = f"growth rate {tag_of(pname, end)}"
        for j, ch in enumerate(an.chains):
            extra.append(dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, parameter=name, scope="chain", chain=ch.label, seed=ch.seed, retained=ch.retained,
                              n_samples=len(allc[j]), ess=L.ess_single(allc[j]), median=float(np.median(allc[j])), rhat_split=np.nan, rhat_rank=np.nan))
        kk = [allc[j] for j, ch in enumerate(an.chains) if ch.retained]
        sm = L.summarise(kk, name)
        extra.append(dict(analysis=an.name, engine=an.engine, model=an.model, dataset=an.dataset, parameter=name, scope="pooled retained chains", chain=f"{len(kk)} chains", seed="", retained=True,
                          n_samples=sm["n_samples"], ess=sm["ess_pooled"], ess_sum_of_chains=sm["ess_sum_of_chains"], median=sm["median"], rhat_split=sm["rhat_split"], rhat_rank=sm["rhat_rank"],
                          criterion_ess200_rhat105_met=bool(sm["ess_pooled"] >= 200 and sm["rhat_split"] <= 1.05)))
    D = pd.concat([D, pd.DataFrame(extra)], ignore_index=True)
    ren = {"growth rate W1": "growth rate P1", "growth rate W2": "growth rate P2", "growth rate W3": "growth rate P3 (own end)", "growth rate W3b": "growth rate P4", "clock.rate": "evolutionary rate",
           "root height (tMRCA)": "root height (date of the most recent common ancestor)", "skygrid.precision": "precision of the Skygrid"}
    D["parameter"] = D["parameter"].replace(ren)
    D["analysis"] = analysis
    D.insert(1, "program", PROGRAM); D.insert(2, "setting", SETTING[s_]); D.insert(3, "genome_set", g_)
    D = D.drop(columns=[c for c in ("engine", "model", "dataset") if c in D.columns])
    # ---------------- posterior file
    PS = B.samples_frame_beast(an)
    keep = ["seed", "state", "clock_rate", "root_height_years", "root_age_decimal_year", "root_date", "kappa"] + [f"lnNe_{i:02d}" for i in range(1, g.K + 1)] + ["skygrid_precision"]
    X = PS[keep].copy()
    X.insert(0, "analysis", analysis)
    for extra_col in ("joint", "prior", "likelihood", "treeLength"):
        if extra_col in PS.columns:
            X[extra_col] = PS[extra_col].values
    for pname, key, end in PERIOD_KEYS:
        if (pname, end) not in ser_r:
            continue
        nm = col_of(pname, end)
        X[f"growth_rate_per_day_{nm}"] = np.concatenate(ser_r[(pname, end)]) / L.YEAR_DAYS
        X[f"R_{nm}"] = np.concatenate(ser_R[(pname, end)])
    assert len(X) == n_used and list(X["seed"].unique()) == [c.seed for c in an.kept]
    for k in ("W1", "W2", "W3", "W3b"):                             # the growth rates of the frame of the functions of work/pipeline_ext/ are those of the columns written here
        nm = {"W1": "P1", "W2": "P2", "W3": "P3_own_end", "W3b": "P4"}[k]
        assert np.allclose(PS[f"growth_rate_per_year_{k}"].values / L.YEAR_DAYS, X[f"growth_rate_per_day_{nm}"].values, rtol=0, atol=1e-15)
    ZERO = pd.DataFrame(zero)[ZERO_COLS]
    res.update(estimates=EST, intervals=INT, trajectories=TRJ, diagnostics=D, posterior=X, tabs=tabs, S=S, series_r=ser_r, series_R=ser_R, zero=ZERO)
    return res


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True); ap.add_argument("--tree", required=True); ap.add_argument("--out", required=True)
    ap.add_argument("--analyses", required=True)
    ap.add_argument("--periods", required=True); ap.add_argument("--periods-by-set", required=True); ap.add_argument("--schedule", required=True)
    ap.add_argument("--test", action="store_true")
    ap.add_argument("--exclude-chains", default=""); ap.add_argument("--tag", default="")
    ap.add_argument("--cut-chains", default="", help="label=state,label=state: chains whose log (in --inputs, written by tools/cut_continued.py) ends at the saved state from which the chain was continued")
    ap.add_argument("--inputs", default="", help="directory written by tools/join_continued.py (logs of continued chains cut and joined); default: the chain directories")
    ap.add_argument("--constants-only", action="store_true", help="write the two tables of the passed constants for the analyses given and derive nothing (used when the tables of the track hold analyses of several calls)")
    a = ap.parse_args()
    root, tree, out = os.path.abspath(a.root), os.path.abspath(a.tree), os.path.abspath(a.out)
    # corrected on 4 October 2026: the directory of the joined files was made absolute only after the script had changed into the working tree, so that a relative
    # path was looked for there and no chain was found (the script stopped on its assertion; nothing was written)
    a.inputs = os.path.abspath(a.inputs) if a.inputs else ""
    periods, pbs, schedule = os.path.abspath(a.periods), os.path.abspath(a.periods_by_set), os.path.abspath(a.schedule)
    os.makedirs(out, exist_ok=True)
    analyses = a.analyses.split(",")
    L, A, B, F, per, pbs_, sched, t_common = setup(tree, root, periods, pbs, schedule, analyses)
    PC = pd.DataFrame(PASSED, columns=["stored_script", "constant_or_function", "value_in_the_stored_script", "value_of_this_round"])
    PC.to_csv(os.path.join(out, ("TEST_" if a.test else "") + "INTERNAL_passed_constants_with_stored_values.csv"), index=False)
    PC2 = PC[["stored_script", "constant_or_function", "value_of_this_round"]].copy()
    PC2["value_of_this_round"] = PC2["value_of_this_round"].str.replace(r"\(.*?stored function gives.*?\)", "", regex=True)
    PC2.to_csv(os.path.join(out, ("TEST_" if a.test else "") + "passed_constants.csv"), index=False)
    if a.constants_only:
        print(f"tables of the passed constants written for {len(analyses)} analyses: {', '.join(analyses)}")
        return 0
    pre = "TEST_" if a.test else ""
    excl = tuple(x for x in a.exclude_chains.split(",") if x)
    for an_ in analyses:
        cut = {x.split("=")[0]: int(x.split("=")[1]) for x in a.cut_chains.split(",") if x and x.split("=")[0].startswith(an_ + "_s")}
        assert not a.cut_chains or (cut and a.tag and a.inputs), "--cut-chains needs --tag and --inputs and at least one chain of the analysis"
        r = derive(an_, root, out, L, A, B, F, t_common, test=a.test, exclude=excl, tag=a.tag, inputs=(os.path.abspath(a.inputs) if a.inputs else None), cut=cut or None)
        d = os.path.join(out, an_ + (("_" + a.tag) if a.tag else "")); os.makedirs(d, exist_ok=True)
        r["chains"].to_csv(os.path.join(d, f"{pre}chains_part.csv"), index=False)
        if r.get("no_comparison"):
            with open(os.path.join(d, f"{pre}no_comparison.txt"), "w") as fh:
                fh.write(r["no_comparison"] + "\n")
        if "estimates" in r:
            r["estimates"].to_csv(os.path.join(d, f"{pre}estimates_part.csv"), index=False)
            r["intervals"].to_csv(os.path.join(d, f"{pre}intervals_part.csv"), index=False)
            r["trajectories"].to_csv(os.path.join(d, f"{pre}trajectories_part.csv"), index=False)
            r["diagnostics"].to_csv(os.path.join(d, f"{pre}diagnostics_part.csv"), index=False)
            r["zero"].to_csv(os.path.join(d, f"{pre}zero_part.csv"), index=False)
            r["posterior"].to_csv(os.path.join(d, f"{pre}posterior_{an_}_{SUFFIX}.csv.gz"), index=False, compression="gzip")
        an = r["an"]
        print(f"{pre}{an_}: chains {len(an.chains)}, retained {len(an.kept)}, states used {sum(len(c.post) for c in an.kept)}; "
              f"rows: estimates {len(r.get('estimates', []))}, intervals {len(r.get('intervals', []))}, trajectories {len(r.get('trajectories', []))}, "
              f"diagnostics {len(r.get('diagnostics', []))}, posterior {len(r.get('posterior', []))}", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
