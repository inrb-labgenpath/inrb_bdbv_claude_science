#!/usr/bin/env python3
"""build_derive_root.py - working root for derive/derive_beast_update.py, in the form that the programs of work/pipeline_ext/ expect:
runs/<analysis>/<chain>/ (UNPACKED FROM THE SAVED ARCHIVE of raw output), sched/ledger.csv (the events of the ledger of the queue runner written in the
form that these programs read), xml/ (records of the parameters of the XML files), aln/ (the alignment), checks/ (periods, ends by set, schedule),
tree/ (the working tree pipeline_round_20261002).
usage: build_derive_root.py <working dir> <root to write> <analysis>[,<analysis>]
(Nothing is derived from the files of a chain, or from a copy of them, before all chains of its analysis have ended and its raw output is saved.)"""
import sys, os, csv, shutil, tarfile, glob
U = os.path.abspath(sys.argv[1]); W = os.path.abspath(sys.argv[2]); ANS = sys.argv[3].split(",")
S = os.path.join(U, "stored")
for d in ("sched", "xml", "aln", "checks", "runs", "tree"): os.makedirs(os.path.join(W, d), exist_ok=True)
WORD = {"started": "started", "ended": "ended", "cut": "interrupted", "started_again_from_state_0": "started_again_from_state_0", "ended_by_the_program": "ended",
        "continued_from_saved_state": "continued_from_saved_state"}
HEAD = ["utc", "event", "chain", "analysis", "genome_set", "setting", "seed", "pid", "exit_code", "last_state_logged", "command", "note"]
def round_events(ev, setting):
    """the events of the ledger of the queue runner in the form that the tools read.  The queue runner (not in the repository) writes a second event 'started' when it starts a cut
    chain again from state 0; the tools know that event as 'started_again_from_state_0'."""
    out, started = [], set()
    for e in ev:
        if e["event"] not in WORD: continue
        word = WORD[e["event"]]
        note = e["note"] + (" (event 'cut' of the ledger of the update)" if e["event"] == "cut" else " (event 'ended_by_the_program' of the ledger of the update)" if e["event"] == "ended_by_the_program" else "")
        if e["event"] == "started":
            if e["chain"] in started:
                word = "started_again_from_state_0"; note += " (second event 'started' of the chain in the ledger of the update: the attempt before was cut from outside and the chain was started again from state 0)"
            started.add(e["chain"])
        out.append([e["utc"], word, e["chain"], e["analysis"], "C00", setting[e["chain"]], e["seed"], e["pid"], e["exit_code"], e["last_state_logged"], e["command"], note])
    return out
rows = []
if os.path.exists(os.path.join(W, "sched", "ledger.csv")):
    rows = [r for r in list(csv.reader(open(os.path.join(W, "sched", "ledger.csv"))))[1:] if r[3] not in ANS]
for an in ANS:
    dst = os.path.join(W, "runs", an)
    if os.path.exists(dst): shutil.rmtree(dst)
    arc = os.path.join(U, "deliver", "%s%s_raw_output_20261007.tar.gz" % ("TEST_" if os.environ.get("FILES_OF_THE_TEST_CHAIN") == "yes" else "", an))
    tmp = os.path.join(W, ".unpacked"); 
    if os.path.exists(tmp): shutil.rmtree(tmp)
    os.makedirs(tmp)
    with tarfile.open(arc) as tar: tar.extractall(tmp)
    for d in sorted(glob.glob(os.path.join(tmp, an, an + "_s*"))): shutil.move(d, os.path.join(dst, os.path.basename(d)))
    ev = list(csv.DictReader(open(os.path.join(tmp, an, "ledger_events_rows.csv"))))
    setting = {c["label"]: c["setting"] for c in csv.DictReader(open(os.path.join(tmp, an, "queue_rows.csv")))}
    for f in glob.glob(os.path.join(tmp, an, "xml", "*.json")): shutil.copy2(f, os.path.join(W, "xml"))
    shutil.rmtree(tmp)
    assert len(glob.glob(os.path.join(dst, an + "_s*"))) == (len(setting) if os.environ.get("FILES_OF_THE_TEST_CHAIN") == "yes" else 8)      # test of the code on short test chains: any number
    rows += round_events(ev, setting)
with open(os.path.join(W, "sched", "ledger.csv"), "w", newline="") as h:
    w = csv.writer(h); w.writerow(HEAD); w.writerows(rows)
for f in glob.glob(os.path.join(U, "xml", "*.json")):
    if not os.path.exists(os.path.join(W, "xml", os.path.basename(f))): shutil.copy2(f, os.path.join(W, "xml"))
for a, b in ((os.path.join(S, "C00_20261007.fasta"), os.path.join(W, "aln", "C00_20261007.fasta")), (os.path.join(U, "tree", "pipeline_round_20261002"), os.path.join(W, "tree", "pipeline_round_20261002"))):
    if not os.path.lexists(b): os.symlink(a, b)
for f in ("periods_of_the_update_20261007.csv", "ends_by_set_20261007.csv", "chain_schedule_round_20261002.csv"): shutil.copy2(os.path.join(S, f), os.path.join(W, "checks", f))
os.makedirs(os.path.join(W, "tools"), exist_ok=True)
for f in ("join_continued_v2.py", "cut_continued.py", "chain_status.py"): shutil.copy2(os.path.join(U, "code_sets", "tools", f), os.path.join(W, "tools", f))      # tools, copied unchanged
qq = [c for c in csv.DictReader(open(os.path.join(U, "sched", "queue.csv")))]
with open(os.path.join(W, "sched", "queue.tsv"), "w") as h:                  # the queue in the form that the tools read
    h.write("order\tlabel\tanalysis\tgenome_set\tsetting\tseed\tthreads\txml_stem\n")
    for c in qq: h.write("%s\t%s\t%s\t%s\t%s\t%s\t1\t%s\n" % (c["order"], c["label"], c["analysis"], c["genome_set"], c["setting"], c["seed"], c["log"][:-4]))
print("root written:", os.path.relpath(W, U), "| analyses:", ",".join(ANS), "| from the archives of raw output", "| ledger rows:", len(rows))
