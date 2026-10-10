#!/usr/bin/env python3
"""continued_chains_table.py - one row for every chain of an analysis that was continued from a saved state, with BOTH attempts:
what the table of the chains (common form) cannot hold, because it gives the last attempt of a chain.
It reads the EXTRACTED archive of the raw output of an analysis (ledger, records of the continuation, logs: first column only)
and the evidence of the incident; it reads no estimate.
The end of a first attempt is not known to the second: the row gives the last write of the first log (time of the file in the
archive) and, apart from it, the time at which the scheduler found the process gone (the 'ended_utc' of the ledger).
Usage: continued_chains_table.py --extracted DIR --analysis NAME --out continued_chains_NAME_20261002.csv"""
import argparse
import csv
import datetime
import glob
import json
import os

ap = argparse.ArgumentParser()
for k in ("extracted", "analysis", "out"):
    ap.add_argument("--" + k, required=True)
a = ap.parse_args()
A = a.analysis
top = os.path.join(a.extracted, f"{A}_raw_output_20261002")
L = [r for r in csv.DictReader(open(os.path.join(top, f"ledger_{A}.csv"), newline="")) if r["analysis"] == A]
Q = json.load(open(os.path.join(top, f"queue_{A}.json")))["chains"]
utc = lambda t: datetime.datetime.fromtimestamp(t, datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def states(path):
    out = []
    with open(path) as fh:
        for line in fh:
            k = line.split("\t", 1)[0]
            if k.isdigit() and line.endswith("\n"):
                out.append(int(k))
    return out


rows = []
for c in Q:
    att = [r for r in L if r["chain"] == c["chain"]]
    d = os.path.join(top, c["dir"])
    cds = sorted(glob.glob(os.path.join(d, "continued_from_state_*")))
    if not cds:
        continue
    assert len(cds) == 1 and len(att) == 2, (c["chain"], cds, len(att))
    rec = json.load(open(os.path.join(cds[0], "continuation_record.json")))
    first, cut, cont = states(os.path.join(d, c["log"])), states(os.path.join(cds[0], "first_log_cut_at_the_saved_state.log")), states(os.path.join(cds[0], c["log"]))
    n = rec["state"]
    joined = cut + [s for s in cont if s > cut[-1]]
    want = list(range(0, int(c["steps"]) + 1, int(c["log_every"])))
    expected_command = c["command"].replace("-save_every", f"-load_state {c['state_file']} -force_resume -save_every").replace("../../../xml/" + os.path.basename(c["xml"]), os.path.basename(c["xml"]))
    rows.append(dict(
        analysis=A, chain=c["chain"], seed=c["seed"], attempts=len(att),
        first_attempt_started_utc=att[0]["started_utc"], first_attempt_last_write_of_the_log_utc=utc(os.path.getmtime(os.path.join(d, c["log"]))),
        first_attempt_end_by_the_ledger_utc_which_is_the_time_at_which_the_scheduler_found_the_process_gone=att[0]["ended_utc"],
        first_attempt_exit_code=att[0]["exit_code"] or "none", first_attempt_outcome=att[0]["outcome"], first_attempt_note_of_the_ledger=att[0]["note"],
        first_attempt_last_logged_state=first[-1], first_attempt_logged_states=len(first),
        saved_state_from_which_the_chain_was_continued=n, saved_state_file=rec["saved_state_file"],
        logged_states_of_the_first_attempt_kept=len(cut), logged_states_of_the_first_attempt_after_the_saved_state_discarded=len(first) - len(cut),
        states_of_the_first_attempt_after_the_saved_state_discarded=first[-1] - n,
        difference_of_the_saved_and_the_recomputed_log_likelihood=rec["difference_of_the_two_log_likelihoods"], limit_of_the_rule="1e-6", released=rec["released"],
        chain_length_of_the_input_file_of_the_continuation=rec["chain_length_of_the_input_file_of_the_continuation"],
        continuation_started_utc=att[1]["started_utc"], continuation_ended_utc=att[1]["ended_utc"], continuation_wall_time_s=att[1]["wall_time_s"],
        continuation_exit_code=att[1]["exit_code"], continuation_outcome=att[1]["outcome"], continuation_first_logged_state=cont[0], continuation_last_logged_state=cont[-1],
        continuation_logged_states=len(cont), the_continuation_logs_the_saved_state_again=("yes" if cont[0] == n else "no"),
        joined_logged_states=len(joined), joined_states_are_every_logged_state_once_from_0_to_the_last_step=("yes" if joined == want else "no"),
        command_of_the_first_attempt=att[0]["command"], command_of_the_continuation=att[1]["command"],
        command_of_the_continuation_is_the_command_of_the_schedule_with_load_state_and_force_resume=("yes" if att[1]["command"] == expected_command else "no"),
        directory_of_the_continuation=os.path.relpath(cds[0], top)))
with open(a.out, "w", newline="") as fh:
    cols = list(rows[0].keys()) if rows else ["analysis", "chain"]
    w = csv.DictWriter(fh, fieldnames=cols, lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
print(A, "| chains continued:", len(rows), "of", len(Q), "| joined states complete for:", sum(1 for r in rows if r["joined_states_are_every_logged_state_once_from_0_to_the_last_step"] == "yes"),
      "| commands as expected for:", sum(1 for r in rows if r["command_of_the_continuation_is_the_command_of_the_schedule_with_load_state_and_force_resume"] == "yes"))
