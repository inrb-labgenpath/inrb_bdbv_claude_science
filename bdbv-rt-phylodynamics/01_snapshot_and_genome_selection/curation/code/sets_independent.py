#!/usr/bin/env python3
"""sets_independent.py (round of October 2026)

Second computation of the genome sets C95, C90, C80 and C00 by code that does not call the pipeline and shares no function
with it: own reader of the alignment, own count of called positions, own date rule, own rule for records that share a sample
identifier, own search for ADAR-type clusters. The class of the quality checks is READ from the screen table (the screen is
not recomputed here).

Rules (README of the frozen pipeline):
  date        collection date given to the day (text of the form YYYY-MM-DD that is a date of the calendar)
  called      called positions (A, C, G, T among positions 1 to 18,900) >= ceil(threshold x 18,900); C00: no threshold
  identifier  for every pair of the table of pairs: the two records are compared at the positions 1 to 18,900 that are called
              in both; if they differ at one position or more, both records are left out; if they do not differ, the record
              with fewer called positions is left out (equal numbers: the later release, then the later accession)
  restricted  the genomes named as restricted use with agreement pending are left out
  class       classes E, E-, D, M and H are left out
A genome is in a set if it meets all five rules. The rules are evaluated for every genome, each on its own.

Writes (OUT = --outdir):
  independent_membership_<sfx>.csv     one row per genome: the five rules, membership, first reason in the order of the
                                       pipeline and first step in the order of the article
usage: python sets_independent.py --metadata META --alignment ALN --screen SCREEN --pairs PAIRS --restricted A,B --suffix SFX --outdir OUT
"""
import argparse
import datetime
import gzip
import math
import os
import re
import csv

L_SCORED = 18900
SETS = {"C95": 0.95, "C90": 0.90, "C80": 0.80, "C00": None}
LEFT_OUT = {"E", "E-", "D", "M", "H"}
ACGT = set("ACGT")


def fasta(path):
    out, name, buf = [], None, []
    op = gzip.open if path.endswith(".gz") else open
    with op(path, "rt") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line[:1] == ">":
                if name is not None:
                    out.append((name, "".join(buf)))
                name, buf = line[1:].split()[0], []
            elif line:
                buf.append(line)
    if name is not None:
        out.append((name, "".join(buf)))
    return out


def need(threshold):
    """smallest whole number n with n / 18,900 >= threshold, in exact arithmetic of whole numbers"""
    if threshold is None:
        return 0
    num = round(threshold * 100)                 # thresholds are given in hundredths
    return -((-num * L_SCORED) // 100)           # ceiling of num * 18900 / 100


def is_day(text):
    if not re.fullmatch(r"\d{4}-\d{2}-\d{2}", text or ""):
        return False
    try:
        datetime.date(int(text[:4]), int(text[5:7]), int(text[8:10]))
        return True
    except ValueError:
        return False


def adar_masks(seqs, called_n):
    """positions (1-based) that the pipeline sets to N: all positions of a cluster except the first.
    Cluster: 3 or more derived T>C changes (or 3 or more A>G changes) of one genome, each at most 300 nt from the next."""
    pol = [k for k, n in enumerate(called_n) if n >= need(0.95)]
    maj = []
    for p in range(L_SCORED):
        cnt = {"A": 0, "C": 0, "G": 0, "T": 0}
        for k in pol:
            b = seqs[k][p]
            if b in cnt:
                cnt[b] += 1
        best = max(cnt.values())
        maj.append(next(b for b in "ACGT" if cnt[b] == best) if best > 0 else "N")
    masks, nclu = [], []
    for s in seqs:
        m, n = [], 0
        for frm, to in (("T", "C"), ("A", "G")):
            pos = [p for p in range(L_SCORED) if s[p] == to and maj[p] == frm]
            if len(pos) < 3:
                continue
            grp = [pos[0]]
            groups = []
            for q in pos[1:]:
                if q - grp[-1] <= 300:
                    grp.append(q)
                else:
                    groups.append(grp)
                    grp = [q]
            groups.append(grp)
            for g in groups:
                if len(g) >= 3:
                    n += 1
                    m += [q + 1 for q in g[1:]]
        masks.append(m)
        nclu.append(n)
    return masks, nclu, "".join(maj)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--screen", required=True)
    ap.add_argument("--pairs", required=True)
    ap.add_argument("--restricted", required=True)
    ap.add_argument("--suffix", required=True)
    ap.add_argument("--outdir", required=True)
    a = ap.parse_args()
    recs = fasta(a.alignment)
    names = [n for n, _ in recs]
    seqs = [s.upper() for _, s in recs]
    assert len(set(names)) == len(names) and all(len(s) == 18940 for s in seqs)
    pos = {n: k for k, n in enumerate(names)}
    with open(a.metadata, newline="", encoding="utf-8") as fh:
        meta = {r["accessionVersion"]: r for r in csv.DictReader(fh)}
    with open(a.screen, newline="", encoding="utf-8") as fh:
        rd = csv.reader(fh)
        head = next(rd)
        ci = head.index("class")
        cls = {row[0]: row[ci] for row in rd}
    assert set(cls) == set(names), "screen table and alignment hold different genomes"
    called = [sum(1 for c in s[:L_SCORED] if c in ACGT) for s in seqs]
    date = [meta[n]["sampleCollectionDate"] for n in names]
    release = [meta[n]["earliestReleaseDate"] for n in names]
    r_date = [is_day(d) for d in date]
    r_class = [cls[n] not in LEFT_OUT for n in names]
    restricted = set(x.strip().rsplit(".", 1)[0] for x in a.restricted.split(",") if x.strip())
    r_use = [n.rsplit(".", 1)[0] not in restricted for n in names]
    assert sum(1 for x in r_use if not x) == len(restricted)
    # records that share a sample identifier
    r_id = [True] * len(names)
    id_why = [""] * len(names)
    partner = [""] * len(names)
    ndiff = [""] * len(names)
    with open(a.pairs, newline="", encoding="utf-8") as fh:
        prs = [(r["accessionVersion_1"], r["accessionVersion_2"]) for r in csv.DictReader(fh)]
    seen = set()
    for x, y in prs:
        assert x in pos and y in pos and x not in seen and y not in seen
        seen |= {x, y}
        i, j = pos[x], pos[y]
        d = sum(1 for p in range(L_SCORED) if seqs[i][p] in ACGT and seqs[j][p] in ACGT and seqs[i][p] != seqs[j][p])
        partner[i], partner[j] = y, x
        ndiff[i] = ndiff[j] = d
        if d > 0:
            r_id[i] = r_id[j] = False
            id_why[i] = id_why[j] = "shares a sample identifier with a discordant record"
        else:
            keep, drop = sorted([i, j], key=lambda k: (-called[k], release[k], names[k]))
            r_id[drop] = False
            id_why[drop] = "shares a sample identifier with a concordant, more complete record"
    masks, nclu, maj = adar_masks(seqs, called)
    os.makedirs(a.outdir, exist_ok=True)
    cols = ["accessionVersion", "collection_date", "release", "called_positions", "class", "partner", "differences_with_partner",
            "rule_date", "rule_identifier", "rule_restricted_use", "rule_class", "adar_clusters", "masked_positions_1based"]
    for s in SETS:
        cols += [f"rule_called_{s}", f"in_{s}", f"reason_pipeline_order_{s}", f"step_article_order_{s}"]
    with open(os.path.join(a.outdir, f"independent_membership_{a.suffix}.csv"), "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(cols)
        for k, n in enumerate(names):
            row = [n, date[k], release[k], called[k], cls[n], partner[k], ndiff[k], r_date[k], r_id[k], r_use[k], r_class[k], nclu[k],
                   ",".join(str(q) for q in masks[k])]
            for s, t in SETS.items():
                r_called = called[k] >= need(t)
                member = r_date[k] and r_called and r_id[k] and r_use[k] and r_class[k]
                if member:
                    why_p = why_a = ""
                else:
                    why_p = ("collection date not given to the day" if not r_date[k] else
                             "called fraction below threshold" if not r_called else
                             "class " + cls[n] if not r_class[k] else
                             id_why[k] if not r_id[k] else
                             "restricted use, agreement of the submitters pending")
                    why_a = ("1 no collection date to the day" if not r_date[k] else
                             "2 called fraction below the threshold" if not r_called else
                             "3 " + id_why[k] if not r_id[k] else
                             "4 restricted use, agreement of the submitters pending" if not r_use[k] else
                             "5 did not pass the quality checks (class " + cls[n] + ")")
                row += [r_called, member, why_p, why_a]
            w.writerow(row)
    with open(os.path.join(a.outdir, f"independent_majority_base_{a.suffix}.txt"), "w") as fh:
        fh.write(maj + "\n")
    print("thresholds in called positions:", {s: need(t) for s, t in SETS.items()})


if __name__ == "__main__":
    main()
