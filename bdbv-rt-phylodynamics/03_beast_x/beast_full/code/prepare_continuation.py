#!/usr/bin/env python3
"""prepare_continuation.py - prepares the continuation of ONE chain of BEAST X that was ended from outside at or beyond
half of its states.

What BEAST X 10.5.0 does after -load_state was established by a test on a small input (test_resume/):
the program runs the chain length of the input file IN ADDITION to the loaded state, logs the loaded state again as the
first row of the new log, and prints the saved and the recomputed log-likelihood.

Steps
  1. the newest saved state of the chain (copy made by the scheduler in saved_states/): state N and saved log-likelihood
  2. a directory of its own, <directory of the chain>/continued_from_state_<N>/, with
       the input file under its name, in which ONLY the chain length differs (steps - N); checked line by line
       the saved state under the name of the state file
       first_log_cut_at_the_saved_state.log and first_trees_cut_at_the_saved_state.trees: the log and the tree file of the
       attempt that was ended, without the states after N (the files of the attempt themselves are left as they are)
  3. a start of the program on a copy of the saved state (directory check_of_the_saved_state/), which is ended as soon as the
     program has printed the recomputed log-likelihood; the continuation is released only if the saved and the recomputed
     log-likelihood differ by less than 1e-6
  4. CONTINUE.json in the directory of the chain: the scheduler then starts the continuation as the next attempt
The chain is marked as continued in the table of chains (make_config.py finds the directory).
Standard library only.
Usage: prepare_continuation.py --root ROOT --control control --chain CHAIN [--check-seconds 600]
"""
import argparse
import csv
import json
import os
import re
import shutil
import signal
import subprocess
import time


def cut_log(src, dst, n):
    kept = last = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            first = line.split("\t", 1)[0]
            if first.isdigit():
                if not line.endswith("\n") or int(first) > n:
                    continue
                kept += 1
                last = int(first)
            fo.write(line)
    return kept, last


def cut_trees(src, dst, n):
    kept = 0
    with open(src) as fi, open(dst, "w") as fo:
        for line in fi:
            m = re.match(r"tree STATE_(\d+)", line)
            if m:
                if int(m.group(1)) > n or not line.rstrip().endswith(";"):
                    continue
                kept += 1
            elif line.strip().lower() == "end;":
                continue
            fo.write(line)
        fo.write("End;\n")
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True)
    ap.add_argument("--control", required=True)
    ap.add_argument("--chain", required=True)
    ap.add_argument("--check-seconds", type=int, default=600)
    a = ap.parse_args()
    root = os.path.abspath(a.root)
    c = next(x for x in json.load(open(os.path.join(a.control, "queue.json")))["chains"] if x["chain"] == a.chain)
    assert c["kind"] == "beast"
    att = [r for r in csv.DictReader(open(os.path.join(a.control, "ledger.csv"), newline="")) if r["chain"] == a.chain]
    assert att and att[-1]["outcome"] == "ended from outside" and att[-1]["what_follows"].startswith("held"), "the chain is not held for a continuation"
    d = os.path.join(root, c["dir"])
    sv = os.path.join(d, "saved_states")
    states = sorted((int(x.rsplit("_state_", 1)[1].split(".")[0]), x) for x in os.listdir(sv) if x.endswith(".state"))
    assert states, "no saved state"
    n, name = states[-1]
    with open(os.path.join(sv, name)) as fh:
        fh.readline()
        l2, l3 = fh.readline().split(), fh.readline().split()
    assert l2 == ["state", str(n)] and l3[0] == "lnL", (l2, l3)
    saved_lnl = float(l3[1])
    steps = int(c["steps"])
    assert n >= steps // 2 and n < steps, (n, steps)
    stem = os.path.basename(c["xml"])[:-4]
    cd = os.path.join(d, f"continued_from_state_{n}")
    os.makedirs(cd, exist_ok=False)
    # ---- input file: only the chain length differs
    old = open(os.path.join(root, c["xml"])).read().split("\n")
    new, changed = [], []
    for i, line in enumerate(old):
        m = re.search(r'(<mcmc id="mcmc" chainLength=")(\d+)(")', line)
        if m:
            assert int(m.group(2)) == steps
            line2 = line[:m.start()] + m.group(1) + str(steps - n) + m.group(3) + line[m.end():]
            changed.append((i + 1, line.strip(), line2.strip()))
            line = line2
        new.append(line)
    assert len(changed) == 1, changed
    open(os.path.join(cd, stem + ".xml"), "w").write("\n".join(new))
    differ = [i + 1 for i, (x, y) in enumerate(zip(old, new)) if x != y]
    assert differ == [changed[0][0]] and len(old) == len(new)
    shutil.copy2(os.path.join(sv, name), os.path.join(cd, c["state_file"]))
    kept, last = cut_log(os.path.join(d, c["log"]), os.path.join(cd, "first_log_cut_at_the_saved_state.log"), n)
    assert last == n, (last, n)
    trees = cut_trees(os.path.join(d, stem + ".trees"), os.path.join(cd, "first_trees_cut_at_the_saved_state.trees"), n)
    command = c["command"].replace("-save_every", f"-load_state {c['state_file']} -force_resume -save_every").replace("../../../xml/" + stem + ".xml", stem + ".xml")
    assert command.count("-load_state") == 1 and command.endswith(" " + stem + ".xml")
    # ---- check of the saved state: saved against recomputed log-likelihood
    ck = os.path.join(cd, "check_of_the_saved_state")
    os.makedirs(ck)
    shutil.copy2(os.path.join(cd, stem + ".xml"), ck)
    shutil.copy2(os.path.join(cd, c["state_file"]), ck)
    env = dict(os.environ, PATH=c["launcher_directory"] + ":" + os.environ.get("PATH", ""))
    env["JAVA_TOOL_OPTIONS"] = c["prelude"].split("=", 1)[1].strip("'")
    out = open(os.path.join(ck, "stdout.txt"), "w")
    p = subprocess.Popen(command.split(), cwd=ck, stdout=out, stderr=subprocess.STDOUT, env=env, start_new_session=True)
    found, t0 = None, time.time()
    while time.time() - t0 < a.check_seconds and p.poll() is None and found is None:
        time.sleep(2)
        txt = open(os.path.join(ck, "stdout.txt"), errors="replace").read()
        m = re.search(r"recomputed likelihood values \(([-0-9.eE+]+) vs\. ([-0-9.eE+]+)\)", txt)
        if m:
            found = (float(m.group(1)), float(m.group(2)))
    try:
        os.killpg(os.getpgid(p.pid), signal.SIGKILL)
    except Exception:
        pass
    p.wait()
    out.close()
    rec = dict(chain=a.chain, state=n, saved_state_file=os.path.join("saved_states", name), log_likelihood_in_the_saved_state=saved_lnl,
               log_likelihoods_printed_by_the_program=found, states_of_the_first_log_kept=kept, trees_of_the_first_tree_file_kept=trees,
               chain_length_of_the_first_input_file=steps, chain_length_of_the_input_file_of_the_continuation=steps - n,
               line_of_the_input_file_that_differs=changed[0], command=command)
    ok = found is not None and abs(found[0] - found[1]) < 1e-6 and min(abs(saved_lnl - found[0]), abs(saved_lnl - found[1])) < 1e-6
    rec["difference_of_the_two_log_likelihoods"] = (abs(found[0] - found[1]) if found else None)
    rec["released"] = "yes" if ok else "no"
    json.dump(rec, open(os.path.join(cd, "continuation_record.json"), "w"), indent=1)
    if not ok:
        print("NOT released:", json.dumps(rec, indent=1))
        return 1
    rel = os.path.relpath(cd, root)
    json.dump(dict(state=n, command=command, run_from=rel, dir=rel, log=stem + ".log",
                   outputs=[stem + ".log", stem + ".trees", stem + ".ops", stem + ".citations.txt", c["state_file"]],
                   files=dict(started="started.txt", ended="finished.txt", exit_code="exit_code.txt", stdout="stdout.txt", command="command.txt", pid="pid.txt")),
              open(os.path.join(d, "CONTINUE.json"), "w"), indent=1)
    print("released:", json.dumps(rec, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
