#!/usr/bin/env python3
"""make_config_v4.py - make_config_v3.py with one addition: every chain of the configuration says to which attempt its start, end, wall time,
exit code and command line belong (key start_end_wall_time_exit_code_and_command_line_are_those_of); the configuration is for derive_analysis_v5.py.
make_config_v3.py was make_config_v2.py with three additions: the configuration holds a note on the analysis that says how many of its chains
were continued from a saved state; every continued chain holds the saved state from which it was continued; with --cut-at-the-saved-state TAG
the configuration is that of a COMPARISON under the name TAG: every continued chain is read up to and including its saved state and nothing of its
continuation is read (a chain that was not continued is read whole).
make_config_v2.py was make_config.py with two additions: the rule by which the dates of the series 'Ne tau at every 7th day' are set is written
into the configuration; with --without-continued-chains yes the chains that were continued from a saved state are left out and the analysis is
named <analysis>_without_continued_chains.
The configuration is for ONE analysis, read from the EXTRACTED archive of its raw output (the derivation reads the files of the archive, not the
files of the run directory), from track_config.py and from the table of the periods of this round.  Generation time: that of the article
(gamma distribution, mean 15.3 days, SD 9.3 days), passed to the function rt_lib.GenTime.
Usage: make_config_v4.py --analysis NAME --extracted DIR --alignment A.fasta --periods P.csv --pipeline work/pipeline_ext --out config.json
       [--without-continued-chains yes | --cut-at-the-saved-state TAG]"""
import argparse
import csv
import glob
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import analyses_config as C  # noqa: E402

GENERATION_TIME_MEAN_DAYS, GENERATION_TIME_SD_DAYS = 15.3, 9.3      # generation time of the article
FIRST_DATE_NOT_BEFORE = "2026-01-01"                                # series 'Ne tau at every 7th day'

ap = argparse.ArgumentParser()
for k in ("analysis", "extracted", "alignment", "periods", "pipeline", "out"):
    ap.add_argument("--" + k, required=True)
ap.add_argument("--without-continued-chains", default="no", choices=("yes", "no"))
ap.add_argument("--cut-at-the-saved-state", default="", help="name (tag) of the comparison")
a = ap.parse_args()
A = next(x for x in C.ANALYSES if x["analysis"] == a.analysis)
ex = os.path.abspath(a.extracted)
if A["kind"] == "beast":
    top = os.path.join(ex, f"{a.analysis}_raw_output_{C.SUFFIX}")
    ledger = list(csv.DictReader(open(os.path.join(top, f"ledger_{a.analysis}.csv"), newline="")))
    par = json.load(open(os.path.join(top, "xml", f"beast_xml_parameters_{A['stem']}.json")))
else:
    top = os.path.join(ex, a.analysis)
    ledger = list(csv.DictReader(open(os.path.join(top, f"ledger_rows_{a.analysis}.csv"), newline="")))
    par = {}
queue = {r["chain"]: r for r in ledger}
seeds = sorted({int(r["seed"]) for r in ledger})
chains = []
for seed in seeds:
    ch = C.chain_name(a.analysis, seed)
    att = [r for r in ledger if r["chain"] == ch]
    last = att[-1]
    e = dict(chain=ch, seed=seed, command_line=last["command"], started_utc=last["started_utc"], ended_utc=last["ended_utc"],
             wall_time_s=last["wall_time_s"], exit_code=last["exit_code"], attempts=len(att),
             started_again=("yes" if any(not r["note"].startswith("continued") for r in att[1:]) else "no"),   # from state 0; a continuation is the next column
             end_of_the_last_attempt_by_the_ledger=last["outcome"],
             continued_from_a_saved_state=("yes" if any(r["note"].startswith("continued") for r in att) else "no"))
    if A["kind"] == "beast":
        d = os.path.join(top, "runs", a.analysis, ch)
        cont = sorted(glob.glob(os.path.join(d, "continued_from_state_*")))
        if cont:
            first_cut = os.path.join(cont[-1], "first_log_cut_at_the_saved_state.log")
            e["saved_state_from_which_the_chain_was_continued"] = int(json.load(open(os.path.join(cont[-1], "continuation_record.json")))["state"])
            assert cont[-1].endswith("continued_from_state_%d" % e["saved_state_from_which_the_chain_was_continued"]) and len(cont) == 1
            e.update(log=first_cut, extra_logs=[os.path.join(cont[-1], A["stem"] + ".log")], trees=os.path.join(d, A["stem"] + ".trees"),
                     stdout=os.path.join(d, "stdout.txt"))
        else:
            e.update(log=os.path.join(d, A["stem"] + ".log"), trees=os.path.join(d, A["stem"] + ".trees"), stdout=os.path.join(d, "stdout.txt"))
        e["files_of_the_chain"] = os.path.relpath(d, ex)
    else:
        e.update(log=os.path.join(top, ch + ".log"), trees=os.path.join(top, ch + ".trees"), stdout=os.path.join(top, ch + ".stdout"))
        e["files_of_the_chain"] = os.path.relpath(top, ex) + "/" + ch + ".*"
    k_att = len(att)
    if e["continued_from_a_saved_state"] == "yes":
        e["start_end_wall_time_exit_code_and_command_line_are_those_of"] = (f"the continuation (attempt {k_att} of {k_att}), which ran only the states after the saved state; "
                      "start, end and wall time of the attempt before it are in the table of the continued chains")
    elif k_att > 1:
        e["start_end_wall_time_exit_code_and_command_line_are_those_of"] = f"the last attempt (attempt {k_att} of {k_att}), which began at state 0"
    else:
        e["start_end_wall_time_exit_code_and_command_line_are_those_of"] = "the only attempt"
    chains.append(e)
RULE = ("dates of the series: the first date is the date of the first interval mid-point after 1 January 2026, the dates follow in "
        "steps of 7 days, the last date is the date of the mid-point of the most recent interval (rule of the track; it gives "
        "the dates of the given table of the primary analysis)")
n_all = len(chains)
n_cont = sum(1 for c in chains if c["continued_from_a_saved_state"] == "yes")
NOTE = ""
if n_cont == n_all and n_all:
    NOTE = f"every chain of this analysis was continued from a saved state after the program had been ended from outside ({n_cont} of {n_all} chains)"
elif n_cont:
    NOTE = f"{n_cont} of the {n_all} chains of this analysis were continued from a saved state after the program had been ended from outside"
if a.cut_at_the_saved_state:
    assert a.without_continued_chains == "no" and n_cont > 0 and A["kind"] == "beast"
    first = {}
    for r in ledger:
        first.setdefault(r["chain"], r)
    for c in chains:
        if c["continued_from_a_saved_state"] != "yes":
            continue
        f = first[c["chain"]]
        c.update(extra_logs=[], steps_of_this_chain=c["saved_state_from_which_the_chain_was_continued"], command_line=f["command"], started_utc=f["started_utc"],
                 ended_utc="", wall_time_s="", exit_code="", end_of_the_last_attempt_by_the_ledger=f["outcome"],
                 continued_from_a_saved_state="not read: the chain is cut at its saved state")
        c["start_end_wall_time_exit_code_and_command_line_are_those_of"] = (f"the first attempt (attempt 1 of {c['attempts']} in the ledger; the column attempts counts the attempts of the ledger); "
                      "the chain is cut at its saved state and has no end, wall time or exit code of its own")
    NOTE = (f"COMPARISON, not the summary of the analysis {a.analysis}: the chains of {a.analysis}, each read up to and including the saved state from which it was continued; "
            f"nothing of a continuation is read ({n_cont} of {n_all} chains cut); command line and start are those of the first attempt; a cut chain has no end, wall time or exit code")
name = a.analysis
if a.without_continued_chains == "yes":
    chains = [c for c in chains if c["continued_from_a_saved_state"] != "yes"]
    name = a.analysis + "_without_continued_chains"
    NOTE = f"COMPARISON, not the summary of the analysis {a.analysis}: the analysis without the chains that were continued from a saved state ({n_cont}) ({n_all - n_cont} of {n_all} chains read)"
if a.cut_at_the_saved_state:
    name = a.cut_at_the_saved_state
    assert a.analysis not in name, "the tag of the comparison must not hold the name of the analysis"
cfg = dict(analysis=name, engine=A["kind"], model=("skygrid" if par.get("num_parameters") else "exponential"), program=A["program"],
           setting=A["setting"], genome_set=A["genome_set"], alignment=os.path.abspath(a.alignment), pipeline=os.path.abspath(a.pipeline),
           periods=os.path.abspath(a.periods), steps=A["steps"], num_parameters=par.get("num_parameters"), cutoff=par.get("cutoff"),
           generation_time_mean_days=GENERATION_TIME_MEAN_DAYS, generation_time_sd_days=GENERATION_TIME_SD_DAYS,
           first_date_not_before=FIRST_DATE_NOT_BEFORE, chains=chains, chains_of_the_schedule=n_all, note_on_the_analysis=NOTE,
           chains_continued_from_a_saved_state=n_cont)
if par.get("num_parameters"):
    cfg["rule_of_the_dates_of_the_series"] = RULE
json.dump(cfg, open(a.out, "w"), indent=1)
print(a.analysis, cfg["engine"], cfg["model"], "chains", len(chains), "| parameters", cfg["num_parameters"], "| cut-off", cfg["cutoff"])
