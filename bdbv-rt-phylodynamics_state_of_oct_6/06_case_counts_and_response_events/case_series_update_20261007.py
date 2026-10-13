#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Case counts, download of 7 October 2026.

Runs the stored script (build_case_series_20261002.py, checked by its SHA-256)
on the files of the repository INRB-UMIE/BDBV2026-Data that the download holds, and writes

  1. the tables of step 1 under their names, with the new date (the tables of step 2, growth
     rates and reproduction numbers, are not written: step 2 is not repeated);
  2. three lists: every report date of the files, every situation
     report of the folder of the repository, every fall of a cumulative count from one report to the next;
  3. the regression against the stored series (not for the post);
  4. the numbers that the post needs (case_counts_for_the_post_<date>.csv);
  5. a file of facts about the run.

No request to a network is made. Nothing is typed by hand: every number is computed from the files.

Example:
  python case_series_update_20261007.py --new-root inputs_20261007/bdbv_inputs_20261007_2056Z \
      --old-root inputs_20261002/bdbv_inputs_20261002_1114Z --stored stored_20261002 --outdir out_20261007 \
      --workdir work --suffix 20261007 --earlier-suffix 20261002
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")          # at most two processor cores

import argparse
import csv
import hashlib
import json
import re
import shutil
import sys

import numpy as np
import pandas as pd

STORED_SCRIPT = "build_case_series_20261002.py"
STORED_SCRIPT_SHA256 = "ea2a05e5b99055c64605e568efe3fbecc4e7d71e9a2f15c1c962b2d23b70bd75"      # the stored script as it stands in this repository
STORED_SCRIPT_VERSION = "build_case_series_20261002.py"
FZ = "data/insp_sitrep/processed/insp_sitrep__cumulative_confirmed_cases__daily.csv"
FN = "data/insp_sitrep/processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv"
FD = "data/shapefiles/DRC_Health_zones.dbf"
RAW_DIR = "data/insp_sitrep/raw/"
# the tables of step 1 and how the script of 2 October writes them
STEP1 = [("weekly_by_province", "weekly_confirmed_cases_by_province", dict(index=False), None),
         ("cumulative_by_province", "cumulative_confirmed_cases_by_province", {}, None),
         ("zone_province_map", "zone_province_map", dict(index=False), None),
         ("mean_weekly_by_zone_group", "mean_weekly_cases_by_zone_group", {}, None),
         ("mean_weekly_zone_set", "mean_weekly_cases_seven_zones", {}, None),
         ("mean_weekly_zone_set_by_treatment", "mean_weekly_cases_seven_zones_by_treatment", {}, 3),
         ("weekly_zone_set", "weekly_confirmed_cases_seven_zones", {}, 3),
         ("crosscheck", "case_series_crosscheck", {}, None),
         ("set_aside", "case_series_values_set_aside", dict(index=False), None)]
KEYS = {"weekly_confirmed_cases_by_province": ["week_ending_sunday"], "cumulative_confirmed_cases_by_province": 0, "zone_province_map": ["name_in_source"],
        "mean_weekly_cases_by_zone_group": ["four_week_period_2026"], "mean_weekly_cases_seven_zones": ["four_week_period_2026"],
        "mean_weekly_cases_seven_zones_by_treatment": ["four_week_period_2026"], "weekly_confirmed_cases_seven_zones": ["week_ending_sunday"],
        "case_series_crosscheck": ["report_date"], "case_series_values_set_aside": None}
STEP2_KEYS_OF_CHECKS = ("third_range_weeks", "third_range_given")
DAY_ZERO = "2026-05-14"         # day from which the numbers of the situation reports count, if they count days


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def git_blob(path):
    b = open(path, "rb").read()
    return hashlib.sha1(b"blob " + str(len(b)).encode() + b"\0" + b).hexdigest()


def write_step1(T, outdir, sfx):
    """Writes the tables of step 1 with the calls of the function write() of the script of 2 October.
    The control run shows that the files are those of that script, byte for byte."""
    os.makedirs(outdir, exist_ok=True)
    files = {}
    for key, name, kw, digits in STEP1:
        t = T[key] if digits is None else T[key].round(digits)
        t.to_csv(os.path.join(outdir, f"{name}_{sfx}.csv"), **kw)
        files[key] = f"{name}_{sfx}.csv"
    return files


def report_number(fname):
    m = re.match(r"(?i)SitRep[_ ]MVE[_ ]?0*(\d{1,3})[_\-]2026(.*)\.pdf$", fname)
    if not m:
        sys.exit(f"STOP: file name of a situation report not understood: {fname}")
    return int(m.group(1))


def provenance(root):
    """Files of case_counts/ against the fetch log and against the list of files of the commit."""
    cc = os.path.join(root, "case_counts")
    log = [r for r in csv.DictReader(open(os.path.join(root, "fetch_log.csv"), encoding="utf-8")) if r["part"] == "cases"]
    files_dir = [d for d in os.listdir(cc) if d.startswith("files_at_")]
    tree_f = [f for f in os.listdir(cc) if f.startswith("tree_")]
    if len(files_dir) != 1 or len(tree_f) != 1:
        sys.exit("STOP: the folder case_counts does not hold exactly one commit")
    tree = json.load(open(os.path.join(cc, tree_f[0])))
    commits = json.load(open(os.path.join(cc, "commits.json")))
    repo = json.load(open(os.path.join(cc, "repository.json")))
    ok_log = sum((sha256(os.path.join(root, r["outfile"])) == r["sha256"]) and (str(os.path.getsize(os.path.join(root, r["outfile"]))) == r["n_bytes"])
                 and r["status"] == "200" for r in log)
    blobs = {e["path"]: e["sha"] for e in tree["tree"] if e["type"] == "blob"}
    fdir = os.path.join(cc, files_dir[0])
    rel = sorted(os.path.relpath(os.path.join(d, f), fdir) for d, _, fs in os.walk(fdir) for f in fs)
    ok_git = sum(blobs.get(r) == git_blob(os.path.join(fdir, r)) for r in rel)
    raw = sorted(e["path"].split("/")[-1] for e in tree["tree"] if e["type"] == "blob" and e["path"].startswith(RAW_DIR))
    return dict(files_dir=fdir, log_rows=len(log), files_equal_log=int(ok_log), fetch_first_utc=min(r["query_start_utc"] for r in log),
                fetch_last_utc=max(r["query_end_utc"] for r in log), files_of_commit=len(rel), git_hash_equal=int(ok_git), tree_sha=tree["sha"],
                tree_truncated=bool(tree.get("truncated")), tree_is_first_commit=bool(tree["sha"] == commits[0]["sha"]),
                commit=commits[0]["sha"], commit_date_utc=commits[0]["commit"]["committer"]["date"], commit_message=commits[0]["commit"]["message"].splitlines()[0],
                commits=[dict(sha=c["sha"], date_utc=c["commit"]["committer"]["date"], message=c["commit"]["message"].splitlines()[0]) for c in commits],
                repository=repo["full_name"], url=repo["html_url"], pushed_at=repo["pushed_at"], default_branch=repo["default_branch"],
                license=(repo.get("license") or {}).get("spdx_id"), raw_files=raw)


def compare_tables(new_path, old_path, name, run):
    """Every difference between a table of 2 October and the same table of a run; one row per value."""
    key = KEYS[name]
    out = []
    if key is None:          # list without a key: compared as a whole
        a = pd.read_csv(old_path, dtype=str, keep_default_na=False)
        b = pd.read_csv(new_path, dtype=str, keep_default_na=False)
        ka = a.apply(lambda r: " | ".join(r.values), axis=1)
        kb = b.apply(lambda r: " | ".join(r.values), axis=1)
        for v in sorted(set(ka) - set(kb)):
            out.append(dict(run=run, table=name, row=v, column="(whole row)", old="present", new="absent", kind="row only in the table of the earlier round"))
        for v in sorted(set(kb) - set(ka)):
            out.append(dict(run=run, table=name, row=v, column="(whole row)", old="absent", new="present", kind="row added"))
        return out, len(a), len(b), int(len(a) * a.shape[1])
    a = pd.read_csv(old_path, dtype=str, keep_default_na=False, index_col=key)
    b = pd.read_csv(new_path, dtype=str, keep_default_na=False, index_col=key)
    n_cmp = 0
    for c in a.columns:
        if c not in b.columns:
            out.append(dict(run=run, table=name, row="(all)", column=c, old="column present", new="column absent", kind="column only in the table of the earlier round"))
    for c in b.columns:
        if c not in a.columns:
            out.append(dict(run=run, table=name, row="(all)", column=c, old="column absent", new="column present", kind="column added"))
    for i in a.index:
        if i not in b.index:
            out.append(dict(run=run, table=name, row=str(i), column="(whole row)", old="present", new="absent", kind="row only in the table of the earlier round"))
            continue
        for c in a.columns:
            if c in b.columns:
                n_cmp += 1
                if a.loc[i, c] != b.loc[i, c]:
                    out.append(dict(run=run, table=name, row=str(i), column=c, old=a.loc[i, c], new=b.loc[i, c], kind="value differs"))
    for i in b.index:
        if i not in a.index:
            for c in b.columns:
                out.append(dict(run=run, table=name, row=str(i), column=c, old="", new=b.loc[i, c], kind="row added"))
    return out, len(a), len(b), n_cmp


def compare_source(old_dir, new_dir):
    rows, summary = [], {}
    for f, col in [(FZ, "cumulative_confirmed_cases"), (FN, "national_cumulative_confirmed_cases")]:
        a = pd.read_csv(os.path.join(old_dir, f), dtype=str, keep_default_na=False)
        b = pd.read_csv(os.path.join(new_dir, f), dtype=str, keep_default_na=False)
        for x in (a, b):
            x["n"] = x.groupby(["nom", "date"]).cumcount()
        ka = a.set_index(["nom", "date", "n"])[col]
        kb = b.set_index(["nom", "date", "n"])[col]
        both = ka.index.intersection(kb.index)
        changed = [i for i in both if ka[i] != kb[i]]
        only_old = sorted(ka.index.difference(kb.index))
        only_new = sorted(kb.index.difference(ka.index), key=lambda i: (i[1], i[0]))
        last_old = max(a.date.str.slice(0, 10))
        name = os.path.basename(f)
        dates_new = sorted(set(i[1] for i in only_new))
        summary[name] = dict(rows_old=len(a), rows_new=len(b), rows_in_both=len(both), values_changed=len(changed), rows_removed=len(only_old), rows_added=len(only_new),
                             report_dates_added=dates_new, report_dates_added_inside_the_range_of_the_earlier_series=[d for d in dates_new if d <= last_old],
                             report_dates_added_after_the_earlier_series=[d for d in dates_new if d > last_old],
                             names_added=sorted(set(b.nom) - set(a.nom)), names_removed=sorted(set(a.nom) - set(b.nom)))
        for i in changed:
            rows.append(dict(run="source tables, row by row", table=name, row=i[1], column=i[0], old=ka[i], new=kb[i], kind="count changed by the source"))
        for i in only_old:
            rows.append(dict(run="source tables, row by row", table=name, row=i[1], column=i[0], old=ka[i], new="", kind="row removed by the source"))
        for i in only_new:
            rows.append(dict(run="source tables, row by row", table=name, row=i[1], column=i[0], old="", new=kb[i],
                             kind="row added by the source, report date inside the range of the earlier series" if i[1] <= last_old
                             else "row added by the source, report date after the earlier series"))
    summary["DRC_Health_zones.dbf"] = dict(identical=bool(sha256(os.path.join(old_dir, FD)) == sha256(os.path.join(new_dir, FD))))
    return rows, summary


def d_iso(x):
    return str(pd.Timestamp(x).date())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--new-root", required=True, help="folder of the new download (holds case_counts/ and fetch_log.csv)")
    ap.add_argument("--old-root", required=True, help="folder of the download of the earlier round (for the control and the comparison of the source tables)")
    ap.add_argument("--stored", required=True, help="folder with copies of the stored script, the stored list of health zones and the stored tables of the earlier round")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--workdir", default="work")
    ap.add_argument("--suffix", default="20261007")
    ap.add_argument("--earlier-suffix", default="20261002")
    ap.add_argument("--zones-file", default=None, help="default: zones_of_analysed_set_<earlier suffix>.txt of the folder --stored")
    ap.add_argument("--report-of-the-draft", default=144, type=int, help="number of the situation report that the draft of the authors names")
    ap.add_argument("--draft-sentence", default=None, help="the sentence of the draft with the four numbers to check")
    a = ap.parse_args(argv)
    sfx, esfx, O, W, S1 = a.suffix, a.earlier_suffix, a.outdir, a.workdir, a.stored
    os.makedirs(O, exist_ok=True)
    os.makedirs(W, exist_ok=True)

    # ---- the script of 2 October, checked by its SHA-256
    script = os.path.join(S1, STORED_SCRIPT)
    if sha256(script) != STORED_SCRIPT_SHA256:
        sys.exit("STOP: the script of 2 October is not the stored version " + STORED_SCRIPT_VERSION)
    sys.path.insert(0, os.path.abspath(S1))
    import build_case_series_20261002 as B
    zones_file = a.zones_file or os.path.join(S1, f"zones_of_analysed_set_{esfx}.txt")
    zones = B.read_zone_list(zones_file)

    # ---- 1. inputs
    pn, po = provenance(a.new_root), provenance(a.old_root)
    new_dir, old_dir = pn["files_dir"], po["files_dir"]
    for p_ in (pn, po):
        if p_["files_equal_log"] != p_["log_rows"] or p_["git_hash_equal"] != p_["files_of_commit"] or p_["tree_truncated"] or not p_["tree_is_first_commit"]:
            sys.exit("STOP: the supplied files do not agree with the fetch log or with the list of files of the commit")

    # ---- 2. the series from the new files (script of 2 October; all report dates)
    T = B.build(new_dir, zone_list=zones)
    files = write_step1(T, O, sfx)
    I_ = T["_internal"]
    S, P, Dz, zprov = I_["S"], I_["P"], I_["D"], I_["zprov"]
    wk = T["weekly_by_province"].copy()
    w7 = T["weekly_zone_set"].copy()
    chk = {k: v for k, v in T["checks"].items() if k not in STEP2_KEYS_OF_CHECKS}

    # ---- 2a. every report date of the files
    Zraw = pd.read_csv(os.path.join(new_dir, FZ), dtype=str, keep_default_na=False)
    Nraw = pd.read_csv(os.path.join(new_dir, FN), dtype=str, keep_default_na=False)
    Zraw["d"] = pd.to_datetime(Zraw.date.str.extract(r"(\d{4}-\d{2}-\d{2})")[0], errors="coerce")
    Zraw["numeric"] = Zraw.cumulative_confirmed_cases.str.fullmatch(r"\d+")
    raw_files = pn["raw_files"]
    numbers = {}
    for f in raw_files:
        numbers.setdefault(report_number(f), []).append(f)
    day0 = pd.Timestamp(DAY_ZERO)
    all_dates = pd.DatetimeIndex(sorted(set(Zraw.d.dropna()) | set(pd.to_datetime(Nraw.date))))
    nat = Nraw.assign(d=pd.to_datetime(Nraw.date)).set_index("d").national_cumulative_confirmed_cases
    sundays = pd.date_range(B.FIRST_SUNDAY, S.index.max(), freq="W-SUN")
    rd = []
    for d in all_dates:
        z = Zraw[Zraw.d == d]
        n = int((d - day0).days)
        rd.append(dict(report_date=d_iso(d), weekday=d.day_name(), national_count=(nat[d] if d in nat.index else ""),
                       rows_in_table_of_health_zones=len(z), rows_with_a_health_zone_and_a_number=int(((z.nom != "NA") & z.numeric).sum()),
                       rows_without_health_zone=int((z.nom == "NA").sum()), rows_with_a_value_that_is_not_a_number=int((~z.numeric).sum()),
                       date_string_repaired=int((z.date != z.d.dt.strftime("%Y-%m-%d")).sum()),
                       date_of_the_series=bool(d in S.index), sunday=bool(d.dayofweek == 6),
                       sum_of_health_zones=(float(S.sum_of_zones[d]) if d in S.index else ""), national_minus_health_zones=(float(S.not_attributed[d]) if d in S.index else ""),
                       days_after_14_May_2026=n, report_of_that_number_in_the_folder="; ".join(numbers.get(n, []))))
    RD = pd.DataFrame(rd)
    RD.to_csv(os.path.join(O, f"case_series_report_dates_{sfx}.csv"), index=False)
    files["report_dates"] = f"case_series_report_dates_{sfx}.csv"
    # ---- 2b. every situation report of the folder of the repository (names of the PDF files; the files themselves were not supplied)
    sr = []
    for n in sorted(numbers):
        d = day0 + pd.Timedelta(days=n)
        sr.append(dict(report_number=n, files_in_the_folder="; ".join(numbers[n]), date_if_the_number_counts_days_after_14_May_2026=d_iso(d),
                       national_count_on_that_date=(nat[d] if d in nat.index else ""), health_zones_with_a_count_on_that_date=int(((Zraw.d == d) & (Zraw.nom != "NA") & Zraw.numeric).sum()),
                       that_date_is_a_report_date_of_the_tables=bool(d in set(all_dates))))
    SR = pd.DataFrame(sr)
    SR.to_csv(os.path.join(O, f"case_series_situation_reports_{sfx}.csv"), index=False)
    files["situation_reports"] = f"case_series_situation_reports_{sfx}.csv"
    # ---- 2c. every fall of a cumulative count from one report to the next
    dec = []
    for z in P.columns:
        s = P[z].dropna()
        df_ = s.diff()
        for d in df_.index[df_ < 0]:
            prev = s.index[s.index.get_loc(d) - 1]
            dec.append(dict(level="health zone", place=z, province=zprov[z], previous_report=d_iso(prev), previous_count=int(s[prev]), report=d_iso(d), count=int(s[d]),
                            change=int(s[d] - s[prev]), treatment="kept as transcribed"))
    for c, lab in [("Ituri", "sum of the health zones of Ituri"), ("Nord-Kivu", "sum of the health zones of Nord-Kivu"), ("other provinces", "sum of the health zones of the other provinces"),
                   ("sum_of_zones", "sum of all health zones"), ("national", "national count")]:
        df_ = S[c].diff()
        for d in df_.index[df_ < 0]:
            prev = S.index[S.index.get_loc(d) - 1]
            dec.append(dict(level="sum", place=lab, province="", previous_report=d_iso(prev), previous_count=int(S[c][prev]), report=d_iso(d), count=int(S[c][d]),
                            change=int(S[c][d] - S[c][prev]), treatment="follows from the counts of the health zones; kept"))
    DEC = pd.DataFrame(dec, columns=["level", "place", "province", "previous_report", "previous_count", "report", "count", "change", "treatment"]).sort_values(
        ["level", "report", "place"]).reset_index(drop=True)
    DEC.to_csv(os.path.join(O, f"case_series_decreases_{sfx}.csv"), index=False)
    files["decreases"] = f"case_series_decreases_{sfx}.csv"
    neg_weeks = wk[(wk[["total_drc", "ituri", "nord_kivu", "other_provinces", "raw_ituri", "raw_nord_kivu", "raw_other"]] < 0).any(axis=1)]

    # ---- 2d. names of the source and the health-zone table
    Zt = B.read_dbf(os.path.join(new_dir, FD))
    Zt["k"] = Zt.Nom.map(B.norm)
    zmap = Zt.groupby("k").PROVINCE.agg(lambda x: sorted(set(x)))
    znames = Zt.groupby("k").Nom.agg(lambda x: sorted(set(x)))
    zm = T["zone_province_map"]
    names = []
    for n_ in zm.name_in_source:
        pv = zmap.get(B.norm(n_))
        names.append(dict(name=n_, n_provinces_in_table=0 if pv is None else len(pv), provinces="" if pv is None else ";".join(pv),
                          names_in_table="" if pv is None else ";".join(znames[B.norm(n_)])))
    NM = pd.DataFrame(names)
    zdist = zm.drop_duplicates("zone")

    # ---- 3. control and regression
    stored = lambda name: os.path.join(S1, f"{name}_{esfx}.csv")
    runs = []
    Tc = B.build(old_dir, zone_list=zones)
    write_step1(Tc, os.path.join(W, "control"), "CONTROL")
    runs.append(("control: stored script on the files of the earlier download", os.path.join(W, "control"), "CONTROL"))
    Zo = pd.read_csv(os.path.join(old_dir, FZ), dtype=str, keep_default_na=False)
    No = pd.read_csv(os.path.join(old_dir, FN), dtype=str, keep_default_na=False)
    ra = os.path.join(W, "regression_A_inputs")
    shutil.rmtree(ra, ignore_errors=True)
    for f in (FZ, FN, FD):
        os.makedirs(os.path.dirname(os.path.join(ra, f)), exist_ok=True)
    shutil.copyfile(os.path.join(new_dir, FD), os.path.join(ra, FD))
    Zraw[Zraw.date.isin(set(Zo.date))][["nom", "date", "cumulative_confirmed_cases"]].to_csv(os.path.join(ra, FZ), index=False)
    Nraw[Nraw.date.isin(set(No.date))].to_csv(os.path.join(ra, FN), index=False)
    Ta = B.build(ra, zone_list=zones)
    write_step1(Ta, os.path.join(W, "regression_A"), "REGA")
    runs.append(("regression A: new files, only the report dates that the earlier download held", os.path.join(W, "regression_A"), "REGA"))
    last_old = Tc["checks"]["last_report_date"]
    Tb = B.build(new_dir, zone_list=zones, max_report_date=last_old)
    write_step1(Tb, os.path.join(W, "regression_B"), "REGB")
    runs.append((f"regression B: new files, every report date up to {last_old}", os.path.join(W, "regression_B"), "REGB"))
    runs.append(("full series from the new files", O, sfx))
    DIFF, REG = [], []
    for run, folder, s_ in runs:
        for _, name, _, _ in STEP1:
            new_path = os.path.join(folder, f"{name}_{s_}.csv")
            dl, n_old, n_new, n_cmp = compare_tables(new_path, stored(name), name, run)
            DIFF += dl
            REG.append(dict(run=run, table=name, rows_earlier=n_old, rows_new=n_new, values_compared=n_cmp,
                            values_different=sum(d["kind"] == "value differs" for d in dl), rows_added=len(set(d["row"] for d in dl if d["kind"] == "row added")),
                            rows_missing=sum(d["kind"].startswith("row only") for d in dl), file_identical_byte_for_byte=bool(sha256(new_path) == sha256(stored(name)))))
    src_rows, src = compare_source(old_dir, new_dir)
    # checks file of the control against the stored checks file
    ck_old = json.load(open(os.path.join(S1, f"case_series_checks_{esfx}.json")))
    ck_ctrl_same = all(Tc["checks"].get(k) == v for k, v in ck_old.items() if k in Tc["checks"] and k != "input_files") and Tc["checks"]["input_files"] == ck_old["input_files"]
    REGT = pd.DataFrame(REG)
    D_all = pd.DataFrame(src_rows + DIFF, columns=["run", "table", "row", "column", "old", "new", "kind"])
    # the cause of each difference, from the report dates of the two downloads
    old_dates, new_dates = set(Tc["_internal"]["S"].index), set(S.index)
    added_dates = sorted(new_dates - old_dates)
    last_old_ts = pd.Timestamp(last_old)

    def cause(r):
        if r["run"].startswith("source tables"):
            return "the source added the report of this date; no count of an earlier report was changed"
        weekly = r["table"] in ("weekly_confirmed_cases_by_province", "weekly_confirmed_cases_seven_zones")
        if r["kind"] == "value differs" and weekly:
            w = pd.Timestamp(r["row"])
            inside = [d for d in added_dates if w - pd.Timedelta(days=6) <= d <= w]
            txt = [f"the source added {len(inside)} report dates in this week ({', '.join(d_iso(d) for d in inside)})"]
            if r["column"] != "report_dates_in_week":
                for b in (w - pd.Timedelta(days=7), w):
                    if b not in old_dates and b in new_dates:
                        lo, hi = max(d for d in old_dates if d < b), min(d for d in old_dates if d > b)
                        txt.append(f"the count of Sunday {d_iso(b)} was interpolated between the reports of {d_iso(lo)} and {d_iso(hi)} in the earlier series and is the count of a report in the new files")
            return "; ".join(txt)
        if r["kind"] == "row added" and weekly:
            return "week after the last complete week of the earlier series"
        if r["kind"] == "row added" and r["table"] == "zone_province_map":
            return "name of a health zone that the earlier files did not hold"
        if r["kind"] == "row added":
            return ("report date added by the source inside the range of the earlier series" if pd.Timestamp(r["row"]) <= last_old_ts
                    else "report date after the last report of the earlier series")
        return ""
    D_all["cause"] = [cause(r) for r in D_all.to_dict("records")] if len(D_all) else []
    D_all = D_all.rename(columns={"old": f"value_{esfx}", "new": f"value_{sfx}"})
    D_all.to_csv(os.path.join(O, f"INTERNAL_case_series_regression_{esfx}_{sfx}.csv"), index=False)
    REGT.to_csv(os.path.join(O, f"INTERNAL_case_series_regression_summary_{esfx}_{sfx}.csv"), index=False)

    # ---- 4. the numbers that the post needs
    f_wk, f_cum, f_w7, f_mb, f_rd, f_sr = (files["weekly_by_province"], files["cumulative_by_province"], files["weekly_zone_set"], files["mean_weekly_zone_set_by_treatment"],
                                           files["report_dates"], files["situation_reports"])
    POST = []

    def add(part, statement, what, value, unit, source, how):
        POST.append(dict(no=len(POST) + 1, part=part, statement_of_the_post=statement, what=what, value=value, unit=unit, source_file=source, how_obtained=how))
    g = lambda v: f"{float(v):g}"
    wki = wk.set_index("week_ending_sunday")
    w7i = w7.copy()
    w7i.index = pd.to_datetime(w7i.index)
    tag = "seven_zones"
    first_date, last_date = all_dates.min(), all_dates.max()
    last_sunday = wki.index.max()
    first_sunday = pd.Timestamp(B.FIRST_SUNDAY)
    zone_dates = pd.DatetimeIndex(sorted(set(Zraw.d[(Zraw.nom != "NA") & Zraw.numeric])))
    first_zone_used = zone_dates[zone_dates <= first_sunday].max()
    first_nat_used = S.index[S.index <= first_sunday].max()
    last_used = S.index[S.index >= last_sunday].min()
    n_of = lambda d: int((pd.Timestamp(d) - day0).days)
    in_folder = lambda d: "yes" if n_of(d) in numbers else "no"
    rule_txt = (f"number = days after {d_iso(day0)}; the files do not state which report belongs to which date; the rule fits "
                f"{int(SR.that_date_is_a_report_date_of_the_tables.sum())} of the {len(SR)} report numbers of the folder")
    A, Bp, C = "a: weekly table and range of situation reports", "b: seven health zones and Nord-Kivu", "c: check of four numbers"
    add(A, "", "first week of the weekly table (week ending)", d_iso(wki.index.min()), "date", f_wk, "first row")
    add(A, "", "last complete week of the weekly table (week ending)", d_iso(last_sunday), "date", f_wk, "last row; the last Sunday on or before the last report date")
    add(A, "", "number of weeks of the weekly table", len(wki), "weeks", f_wk, "rows")
    add(A, "", "weeks whose split by province is pooled over several weeks (weeks ending)", "; ".join(d_iso(d) for d in wki.index[wki.weeks_pooled_for_province_split > 1]), "dates", f_wk,
        "rows with weeks_pooled_for_province_split above 1")
    fewest = int(wki.report_dates_in_week.min())
    prov = wki.index[(wki.report_dates_in_week == fewest) & (wki.index == last_sunday)]
    add(A, "", "provisional week (the last week, if it rests on fewer report dates than any other week)",
        "; ".join(d_iso(d) for d in prov) if len(prov) else "none", "week ending", f_wk,
        f"the last week rests on {int(wki.report_dates_in_week.iloc[-1])} report dates; the fewest report dates of a week are {fewest} "
        f"(week ending {'; '.join(d_iso(d) for d in wki.index[wki.report_dates_in_week == fewest])}); the count on the last day of the last week is "
        + ("that of a report" if last_sunday in S.index else "interpolated between two reports"))
    add(A, "", "report dates in the last complete week", int(wki.report_dates_in_week.iloc[-1]), "report dates", f_wk, "column report_dates_in_week, last row")
    add(A, "", "Sundays of the series whose national count is interpolated between two reports", "; ".join(d_iso(d) for d in sundays if d not in S.index) or "none", "dates", f_rd,
        "Sundays from the first to the last that are not a report date of the national table")
    add(A, "", "first report date of the tables", d_iso(first_date), "date", f_rd, "first row")
    add(A, "", "last report date of the tables", d_iso(last_date), "date", f_rd, "last row")
    add(A, "", "report dates with a national count", int(len(S)), "report dates", f_rd, "rows with date_of_the_series true")
    add(A, "", "report dates with counts of health zones", int(len(zone_dates)), "report dates", f_rd, "rows with rows_with_a_health_zone_and_a_number above 0")
    add(A, "", "first report on which the weekly table rests: national count", d_iso(first_nat_used), "date", f_rd, "the report date on or before the first Sunday (" + d_iso(first_sunday) + ")")
    add(A, "", "first report on which the weekly table rests: counts of health zones", d_iso(first_zone_used), "date", f_rd,
        f"the last report date with counts of health zones on or before the first Sunday; its counts are carried forward to that Sunday; the table of the health zones itself begins on {d_iso(zone_dates.min())} "
        f"({int(((Zraw.d == zone_dates.min()) & (Zraw.nom != 'NA') & Zraw.numeric).sum())} row), a date that enters the cumulative table only")
    add(A, "", "last report on which the weekly table rests", d_iso(last_used), "date", f_rd, "the report date on or after the last Sunday (" + d_iso(last_sunday) + ") that is used for the count of that Sunday")
    add(A, "", "situation reports in the folder of the repository: PDF files", len(raw_files), "files", f_sr, "names of the files under data/insp_sitrep/raw/ in the list of files of the commit")
    add(A, "", "situation reports in the folder of the repository: distinct report numbers", len(numbers), "reports", f_sr, "rows")
    add(A, "", "lowest report number in the folder", min(numbers), "number", f_sr, "first row")
    add(A, "", "highest report number in the folder", max(numbers), "number", f_sr, "last row")
    for lab, d in [("first report date with counts of health zones on which the weekly table rests", first_zone_used),
                   ("first report date with a national count on which the weekly table rests", first_nat_used),
                   ("last report on which the weekly table rests", last_used), ("last report date of the tables", last_date)]:
        add(A, "", f"report number of the {lab}, if the numbers count days", n_of(d), "number (not stated in the files)", f_rd + "; " + f_sr, rule_txt + f"; a report of this number is in the folder: {in_folder(d)}")
    # (b) seven health zones: mean weekly cases by four-week period
    st1 = "reported cases in the same health zones rose in June and July"
    MB = T["mean_weekly_zone_set_by_treatment"]
    for per, r in MB.iterrows():
        add(Bp, st1, f"mean weekly confirmed cases of the seven health zones, four-week period {per} (treatment A)", round(float(r.treatment_A), 3), "cases per week", f_mb,
            "mean of the weekly counts of the period; weekly count = national weekly count multiplied by the share of the seven health zones among the cases attributed to health zones in the pool of weeks")
    for per, r in MB.iterrows():
        add(Bp, st1, f"mean weekly confirmed cases of the seven health zones, four-week period {per} (counts as attributed)", round(float(r.as_attributed), 3), "cases per week", f_mb,
            "mean of the weekly differences of the cumulative counts of the seven health zones")
    jj = w7i[(w7i.index >= "2026-06-01") & (w7i.index <= "2026-07-31")]

    # The words of a row say what the split did in the week of that row; each sentence is tested against the numbers of the row before it is written.
    def pool_of(d):
        return w7i.index[w7i.pool_of_weeks == w7i.pool_of_weeks[d]]

    def words_seven(d, r):
        nat, allz, att, tA = float(r.national_confirmed), float(r.all_zones_as_attributed), float(r[f"{tag}_as_attributed"]), float(r[f"{tag}_treatment_A"])
        pw = pool_of(d)
        share = float(w7i.loc[pw, f"{tag}_as_attributed"].sum() / w7i.loc[pw, "all_zones_as_attributed"].sum())
        if abs(tA - nat * share) > 1e-6 or len(pw) != int(wki.weeks_pooled_for_province_split[d]):
            sys.exit(f"STOP: treatment A of the week ending {d_iso(d)} is not the national count multiplied by the share of its pool")
        if len(pw) > 1:
            return (f"column {tag}_treatment_A: the national count of the week ({g(nat)}) multiplied by the share of the seven health zones among the cases attributed to health zones "
                    f"in the pool of {len(pw)} weeks to which the week belongs (weeks ending {d_iso(pw.min())} to {d_iso(pw.max())}; share {share:.3f}); count as attributed {g(att)}")
        if abs(nat - allz) < 1e-9:
            if abs(tA - att) > 1e-6:
                sys.exit(f"STOP: week ending {d_iso(d)}: treatment A differs from the count as attributed although nothing is unattributed")
            return (f"column {tag}_treatment_A; the week is not pooled and its national count ({g(nat)}) equals the sum of the health zones, "
                    f"so that treatment A equals the count as attributed ({g(att)})")
        if abs(nat - allz) > B.POOL_TOLERANCE * nat + 1e-9 or abs(tA - att * nat / allz) > 1e-6 or abs(tA - att) < 1e-9:
            sys.exit(f"STOP: week ending {d_iso(d)}: the numbers of the week do not fit the rule of the pools")
        return (f"column {tag}_treatment_A; the week is not pooled; its national count ({g(nat)}) differs from the sum of the health zones ({g(allz)}) by {g(round(nat - allz, 1))} cases not attributed to a health zone, "
                f"which is within 5 % of the national count; treatment A is the count as attributed ({g(att)}) multiplied by {g(nat)}/{g(allz)} and differs from it")

    def words_province(d, col, raw):
        r = wki.loc[d]
        nat, v, a_ = float(r.total_drc), float(r[col]), float(r[raw])
        allz = round(float(r.raw_ituri + r.raw_nord_kivu + r.raw_other), 1)
        pw = pool_of(d)
        if len(pw) != int(r.weeks_pooled_for_province_split):
            sys.exit(f"STOP: week ending {d_iso(d)}: the two tables do not give the same pool")
        if len(pw) > 1:
            tot = wki.loc[pw, ["raw_ituri", "raw_nord_kivu", "raw_other"]].sum().sum()
            if abs(v - round(nat * float(wki.loc[pw, raw].sum()) / float(tot), 1)) > 0.05 + 1e-9:
                sys.exit(f"STOP: week ending {d_iso(d)}: the column {col} is not the split of the pool")
            return (f"column {col}: the national count of the week ({g(nat)}) split by province in proportion to the counts attributed to health zones in the pool of {len(pw)} weeks "
                    f"to which the week belongs (weeks ending {d_iso(pw.min())} to {d_iso(pw.max())}); count as attributed {g(a_)}")
        if abs(nat - allz) < 1e-9:
            if abs(v - a_) > 1e-9:
                sys.exit(f"STOP: week ending {d_iso(d)}: the column {col} differs from the count as attributed although nothing is unattributed")
            return f"column {col}; the week is not pooled and its national count ({g(nat)}) equals the sum of the health zones, so that the column equals the count as attributed ({g(a_)})"
        if abs(nat - allz) > B.POOL_TOLERANCE * nat + 1e-9 or abs(v - round(a_ * nat / allz, 1)) > 0.05 + 1e-9:
            sys.exit(f"STOP: week ending {d_iso(d)}: the numbers of the week do not fit the rule of the pools")
        return (f"column {col}; the week is not pooled; its national count ({g(nat)}) differs from the sum of the health zones ({g(allz)}) by {g(round(nat - allz, 1))} cases not attributed to a health zone, "
                f"which is within 5 % of the national count; the column is the count as attributed ({g(a_)}) multiplied by {g(nat)}/{g(allz)}, rounded to one decimal"
                + (", and differs from it" if abs(v - a_) > 1e-9 else ", and equals it after rounding"))

    for d, r in jj.iterrows():
        add(Bp, st1, f"weekly confirmed cases of the seven health zones, week ending {d_iso(d)} (treatment A)", round(float(r[f"{tag}_treatment_A"]), 3), "cases", f_w7, words_seven(d, r))
    for d, r in jj.iterrows():
        add(Bp, st1, f"weekly confirmed cases of the seven health zones, week ending {d_iso(d)} (counts as attributed)", round(float(r[f"{tag}_as_attributed"]), 3), "cases", f_w7,
            f"column {tag}_as_attributed")
    rest = w7i[w7i.index > (pd.Timestamp(B.BLOCK_BASE) + pd.Timedelta(days=28 * len(MB)))]
    for d, r in rest.iterrows():
        add(Bp, "", f"weekly confirmed cases of the seven health zones, week ending {d_iso(d)} (after the last complete four-week period; treatment A)", round(float(r[f"{tag}_treatment_A"]), 3), "cases", f_w7,
            words_seven(d, r))
    # (b) Nord-Kivu
    st2 = "cases in Nord-Kivu rose until mid-September"
    nk = wki.nord_kivu
    top = nk[nk == nk.max()]
    if len(top) != 1:
        sys.exit("STOP: the highest weekly count of Nord-Kivu occurs in more than one week; no rule says which week to take")
    peak = top.index[0]
    add(Bp, st2, "week with the highest weekly count of Nord-Kivu (week ending)", d_iso(peak), "date", f_wk, "row with the largest value of the column nord_kivu")
    add(Bp, st2, f"weekly confirmed cases of Nord-Kivu, week ending {d_iso(peak)} (the highest)", g(nk[peak]), "cases", f_wk, words_province(peak, "nord_kivu", "raw_nord_kivu"))
    for k in range(-4, 5):
        if k == 0:
            continue
        d = peak + pd.Timedelta(days=7 * k)
        lab = f"{abs(k)} week{'s' if abs(k) > 1 else ''} {'before' if k < 0 else 'after'} the highest"
        if d in nk.index:
            add(Bp, st2, f"weekly confirmed cases of Nord-Kivu, week ending {d_iso(d)} ({lab})", g(nk[d]), "cases", f_wk, words_province(d, "nord_kivu", "raw_nord_kivu"))
        else:
            add(Bp, st2, f"weekly confirmed cases of Nord-Kivu, week ending {d_iso(d)} ({lab})", "not in the series", "cases", f_wk,
                "the series ends with the week ending " + d_iso(last_sunday) if d > last_sunday else "the series begins with the week ending " + d_iso(wki.index.min()))
    # (c) the most recent situation report of the files, the report that the article names, and the stored report
    hi = max(numbers)
    d_hi = day0 + pd.Timedelta(days=hi)
    draft_no = int(a.report_of_the_draft)
    d_draft = day0 + pd.Timedelta(days=draft_no)
    yn = lambda b: "yes" if b else "no"
    add(C, "", "most recent situation report in the files: date", d_iso(last_date), "date", f_cum,
        "last report date of the table of the health zones and of the national table of the repository (last row)")
    add(C, "", "most recent situation report in the files: number", hi, "number", f_sr,
        "highest report number of a PDF file in the folder data/insp_sitrep/raw/ (name in the list of files of the commit; the PDF itself was not supplied and was not read); the tables carry no report numbers")
    add(C, "", f"date of report {hi}, if the numbers count days", d_iso(d_hi), "date (not stated in the files)", f_sr, rule_txt)
    add(C, "", f"the last report date of the tables is the date that the count of days gives for report {hi}", yn(d_hi == last_date), "", f_sr + "; " + f_cum, "comparison of the two dates")
    named = [c for c in pn["commits"] if re.search(rf"(?i)sit\s*rep\s*0*{draft_no}\b", c["message"])]
    add(C, "", f"commits of the repository whose message names report {draft_no}", "; ".join(f"{c['sha'][:8]} ({c['date_utc']}): {c['message']}" for c in named) or "none", "", "commits.json",
        "messages of the commits of the download that name the report number")
    held = (draft_no in numbers) and (d_draft in set(S.index)) and (d_draft in set(zone_dates))
    add(C, "", f"the files hold report {draft_no} (the report that the draft names)", yn(held), "", f_sr + "; " + f_rd,
        f"a PDF file of this number is in the folder: {yn(draft_no in numbers)}; the national table holds a count for the date that the count of days gives for it ({d_iso(d_draft)}): {yn(d_draft in set(S.index))}; "
        f"the table of the health zones holds counts for that date: {yn(d_draft in set(zone_dates))}; report {draft_no} is the most recent report of the files: {yn(draft_no == hi and d_draft == last_date)}")
    add(C, "", "date of the commit of the repository", pn["commit_date_utc"], "date and time (UTC)", "commits.json", f"commit {pn['commit']}, the most recent commit of the list; message: {pn['commit_message']}")
    draft = None
    if a.draft_sentence:
        m_ = re.search(r"N°\s*(\d+)\), ([\d,]+) confirmed cases had been reported in the DRC, ([\d,]+) of them in health zones of Ituri, ([\d,]+) in health zones of Nord-Kivu, and ([\d,]+) in other provinces", a.draft_sentence)
        if not m_:
            sys.exit("STOP: the sentence of the draft is not in the expected form")
        if int(m_.group(1)) != draft_no:
            sys.exit("STOP: the sentence of the draft names another report than --report-of-the-draft")
        draft = dict(zip(["national", "Ituri", "Nord-Kivu", "other provinces"], [int(x.replace(",", "")) for x in m_.groups()[1:]]))
        d_p = S.index[S.index <= pd.Timestamp(last_old)].max()
        names_c = {"national": "the DRC", "Ituri": "health zones of Ituri", "Nord-Kivu": "health zones of Nord-Kivu", "other provinces": "other provinces"}
        for k, lab in names_c.items():
            add(C, "", f"the number of the draft for {lab} ({draft[k]}) equals the count of the most recent situation report in the files", yn(draft[k] == int(S.loc[last_date, k])), "", f_cum,
                f"draft (sentence given in the call): {draft[k]}; files, report of {d_iso(last_date)}: {int(S.loc[last_date, k])}; files, report of {d_iso(d_p)}: {int(S.loc[d_p, k])} "
                f"(equal to the number of the draft: {yn(draft[k] == int(S.loc[d_p, k]))})")
        s3 = draft["Ituri"] + draft["Nord-Kivu"] + draft["other provinces"]
        add(C, "", "sum of the three parts of the draft", s3, "cases", "sentence of the draft given in the call",
            f"{draft['Ituri']} + {draft['Nord-Kivu']} + {draft['other provinces']}; equal to the total of the draft ({draft['national']}): {yn(s3 == draft['national'])}; "
            f"equal to the national count of the report of {d_iso(d_p)} ({int(S.loc[d_p, 'national'])}): {yn(s3 == int(S.loc[d_p, 'national']))}")
    how_parts = ("sum of the cumulative counts of the health zones of the province in the report of that date; health zone to province by the health-zone table of the repository; "
                 "the count of a health zone that is absent from the report is carried forward from its last report")
    d_prev = S.index[S.index <= pd.Timestamp(last_old)].max()
    n_prev = n_of(d_prev)
    add(C, "", f"report number of the report date {d_iso(d_prev)}, if the numbers count days", n_prev, "number (not stated in the files)", f_rd + "; " + f_sr,
        rule_txt + f"; a report of this number is in the folder: {in_folder(d_prev)}")
    between_n = [n for n in range(n_prev + 1, hi) if n not in numbers]
    between_d = [d for d in pd.date_range(d_prev + pd.Timedelta(days=1), last_date - pd.Timedelta(days=1)) if d not in set(all_dates)]
    add(C, "", f"report numbers between {n_prev} and {hi} of which the folder holds no PDF file", "; ".join(str(n) for n in between_n) or "none", "numbers", f_sr,
        f"numbers {n_prev + 1} to {hi - 1} that are no row of the table; in the folder: {', '.join(str(n) for n in range(n_prev + 1, hi) if n in numbers) or 'none'}")
    add(C, "", f"days between {d_iso(d_prev)} and {d_iso(last_date)} for which the tables hold no report", "; ".join(d_iso(d) for d in between_d) or "none", "dates", f_rd,
        "days between the two report dates that are no row of the table"
        + ("; by the count of days: " + ", ".join(f"{d_iso(d)} = report {n_of(d)}" for d in between_d) + "; the files do not say whether a report of that number was issued" if between_d else ""))
    add(C, "", "last complete week of the series (week ending)", d_iso(last_sunday), "date", f_wk, "last row; the last Sunday on or before the last report date")
    reports_c = [(f"most recent situation report in the files ({d_iso(last_date)}; report {hi} by the count of days)", last_date)]
    if held and d_draft != last_date:
        reports_c.append((f"report {draft_no} by the count of days ({d_iso(d_draft)})", d_draft))
    reports_c.append((f"report of {d_iso(d_prev)} (report {n_prev} by the count of days)", d_prev))
    for lab, d in reports_c:
        r = S.loc[d]
        z = Zraw[(Zraw.d == d) & (Zraw.nom != "NA") & Zraw.numeric]
        add(C, "", f"cumulative confirmed cases of the DRC, {lab}", int(r.national), "cases", f_cum, "national table of the repository, row of that date")
        add(C, "", f"cumulative confirmed cases in health zones of Ituri, {lab}", int(r.Ituri), "cases", f_cum, how_parts)
        add(C, "", f"cumulative confirmed cases in health zones of Nord-Kivu, {lab}", int(r["Nord-Kivu"]), "cases", f_cum, how_parts)
        add(C, "", f"cumulative confirmed cases in health zones of other provinces, {lab}", int(r["other provinces"]), "cases", f_cum, how_parts)
        add(C, "", f"sum of the three parts, {lab}", int(r.Ituri + r["Nord-Kivu"] + r["other provinces"]), "cases", f_cum, "sum of the three columns")
        add(C, "", f"the three parts sum to the national count, {lab}", "yes" if int(r.Ituri + r["Nord-Kivu"] + r["other provinces"]) == int(r.national) else "no", "", f_cum, "comparison")
        zprov_map = zm.set_index("name_in_source").province
        zz = z.assign(p=z.nom.map(zprov_map), v=z.cumulative_confirmed_cases.astype(int))
        second = (int(zz[zz.p == "Ituri"].v.sum()), int(zz[zz.p == "Nord-Kivu"].v.sum()), int(zz[~zz.p.isin(["Ituri", "Nord-Kivu"])].v.sum()))
        add(C, "", f"the three parts computed a second way agree, {lab}", "yes" if second == (int(r.Ituri), int(r["Nord-Kivu"]), int(r["other provinces"])) else "no", "",
            os.path.basename(FZ), f"rows of that report date only ({len(zz)} health zones), summed by province: {second[0]}, {second[1]}, {second[2]}")
        cf = Dz[Dz.date <= d].sort_values("date").groupby("zone").tail(1)
        cf = cf[(cf.date < d)]
        add(C, "", f"health zones whose count is carried forward, {lab}", "; ".join(f"{r_.zone} ({int(r_.v)} cases, report of {d_iso(r_.date)})" for r_ in cf.itertuples()) or "none", "", f_cum,
            "health zones with a count in an earlier report and none in the report of that date")
    PT = pd.DataFrame(POST)
    PT.to_csv(os.path.join(O, f"case_counts_for_the_post_{sfx}.csv"), index=False)

    # ---- 5. checks and facts
    chk["files_written"] = dict(files)
    chk["source_note"] = f"INRB-UMIE/BDBV2026-Data, commit {pn['commit']}"
    chk["decreases_of_a_health_zone"] = int((DEC.level == "health zone").sum())
    chk["weeks_with_a_negative_count"] = [d_iso(d) for d in neg_weeks.week_ending_sunday]
    chk["sundays_with_interpolated_national_count"] = [d_iso(d) for d in sundays if d not in S.index]
    chk["report_dates_in_the_tables"] = int(len(all_dates))
    json.dump(chk, open(os.path.join(O, f"case_series_checks_{sfx}.json"), "w"), indent=1, ensure_ascii=False)
    X = T["crosscheck"]
    nzd = X[X.national_minus_health_zones != 0]
    lic = [l for l in open(os.path.join(new_dir, "LICENSE.md"), encoding="utf-8").read().split("\n")[:4] if l.strip()]
    meta = open(os.path.join(new_dir, "data/insp_sitrep/metadata.yaml"), encoding="utf-8").read()
    rdm = open(os.path.join(new_dir, "README.md"), encoding="utf-8").read()
    facts = dict(suffix=sfx, earlier_suffix=esfx, new=pn, old={k: v for k, v in po.items() if k != "raw_files"}, zones=zones, zones_file=os.path.basename(zones_file),
                 stored_script=dict(file=STORED_SCRIPT, version=STORED_SCRIPT_VERSION, sha256=STORED_SCRIPT_SHA256),
                 license_first_lines=lic, metadata_license=re.search(r'license: "(.*?)"', meta).group(1), metadata_period=re.search(r'period: "(.*?)"', meta).group(1),
                 readme_license_sentence=re.search(r"The repository code is licensed under the terms in LICENSE\.[^\n]*", rdm).group(0),
                 readme_notice_preliminary=bool(re.search(r"should be considered preliminary", rdm)),
                 dbf_records=int(len(Zt)), names_in_source=int(len(NM)), names_with_unique_province_in_table=int((NM.n_provinces_in_table == 1).sum()),
                 names_fixed=NM[NM.n_provinces_in_table != 1].to_dict("records"), merged=[(r.name_in_source, r.zone) for r in zm.itertuples() if r.name_in_source != r.zone],
                 zones_after_merging=int(len(zdist)), zones_by_group={"Ituri": int((zdist.province == "Ituri").sum()), "Nord-Kivu": int((zdist.province == "Nord-Kivu").sum()),
                                                                      "other provinces": int((~zdist.province.isin(["Ituri", "Nord-Kivu"])).sum())},
                 other_provinces_list=sorted(set(zdist.province) - {"Ituri", "Nord-Kivu"}), zones_without_any_numeric_value=chk["health_zones_without_any_numeric_value"],
                 names_new_since_earlier=src[os.path.basename(FZ)]["names_added"],
                 pools=[(d_iso(g_.index.min()), d_iso(g_.index.max()), int(len(g_))) for _, g_ in w7i.groupby("pool_of_weeks") if len(g_) > 1],
                 report_dates_in_last_weeks={d_iso(w): [d_iso(d) for d in S.index[(S.index > w - pd.Timedelta(days=7)) & (S.index <= w)]] for w in wki.index[-3:]},
                 report_dates_after_last_complete_week=[d_iso(d) for d in S.index[S.index > last_sunday]],
                 days_without_report=[d_iso(d) for d in pd.date_range(first_date, last_date) if d not in set(all_dates)],
                 dates_without_zone_values=[d_iso(d) for d in S.index if d not in set(zone_dates)],
                 crosscheck=dict(report_dates=int(len(X)), equal=int((X.national_minus_health_zones == 0).sum()), different=int(len(nzd)), last_difference=d_iso(nzd.index.max()),
                                 equal_since=int((X.index > nzd.index.max()).sum()), first_equal_since=d_iso(X.index[X.index > nzd.index.max()].min()),
                                 largest_positive=float(nzd.national_minus_health_zones.max()), largest_positive_date=d_iso(nzd.national_minus_health_zones.idxmax()),
                                 largest_negative=float(nzd.national_minus_health_zones.min()), largest_negative_date=d_iso(nzd.national_minus_health_zones.idxmin())),
                 day_rule=dict(day_zero=DAY_ZERO, report_numbers=int(len(SR)), numbers_with_a_report_date=int(SR.that_date_is_a_report_date_of_the_tables.sum()),
                               numbers_without_a_report_date=[int(n) for n in SR.report_number[~SR.that_date_is_a_report_date_of_the_tables]],
                               report_dates=int(len(RD)), report_dates_with_a_file=int((RD.report_of_that_number_in_the_folder != "").sum()),
                               report_dates_without_a_file=[(r.report_date, int(r.days_after_14_May_2026)) for r in RD[RD.report_of_that_number_in_the_folder == ""].itertuples()],
                               last_number_without_date=max([int(n) for n in SR.report_number[~SR.that_date_is_a_report_date_of_the_tables]] or [0]),
                               last_date_without_file=max([r.report_date for r in RD[RD.report_of_that_number_in_the_folder == ""].itertuples()] or [""])),
                 raw=dict(files=len(raw_files), numbers=len(numbers), lowest=min(numbers), highest=max(numbers), numbers_with_two_files=sorted(n for n in numbers if len(numbers[n]) > 1),
                          absent_numbers=[n for n in range(min(numbers), max(numbers) + 1) if n not in numbers],
                          files_added_since_earlier=sorted(set(raw_files) - set(po["raw_files"])), files_removed_since_earlier=sorted(set(po["raw_files"]) - set(raw_files))),
                 regression=dict(control_checks_file_equal=bool(ck_ctrl_same), last_report_date_of_earlier_series=last_old,
                                 by_run={run: dict(tables=int((REGT.run == run).sum()), values_compared=int(REGT[REGT.run == run].values_compared.sum()),
                                                   values_different=int(REGT[REGT.run == run].values_different.sum()), rows_added=int(REGT[REGT.run == run].rows_added.sum()),
                                                   rows_missing=int(REGT[REGT.run == run].rows_missing.sum()),
                                                   files_identical=int(REGT[REGT.run == run].file_identical_byte_for_byte.sum())) for run, _, _ in runs},
                                 source=src),
                 part_c=dict(report_of_the_draft=draft_no, numbers_of_the_draft=draft, date_of_the_report_of_the_draft_by_count_of_days=d_iso(d_draft), held=bool(held), highest_number=int(hi),
                             last_report_date=d_iso(last_date), commits_naming_the_report=named, earlier_report_date=d_iso(d_prev), earlier_report_number_by_count_of_days=int(n_prev),
                             numbers_between_without_file=[int(n) for n in between_n], numbers_between_with_file=[int(n) for n in range(n_prev + 1, hi) if n in numbers],
                             days_between_without_report=[d_iso(d) for d in between_d],
                             cumulative={d_iso(d): {k: int(S.loc[d, k]) for k in ("national", "Ituri", "Nord-Kivu", "other provinces", "sum_of_zones", "not_attributed")} for d in (last_date, d_prev)},
                             health_zones_in_report={d_iso(d): int(((Zraw.d == d) & (Zraw.nom != "NA") & Zraw.numeric).sum()) for d in (last_date, d_prev)}),
                 first_zone_date_used=d_iso(first_zone_used), first_national_date_used=d_iso(first_nat_used), last_date_used=d_iso(last_used),
                 first_sunday=d_iso(first_sunday), last_sunday=d_iso(last_sunday),
                 post_rows=int(len(PT)), runs=[r[0] for r in runs])
    json.dump(facts, open(os.path.join(O, f"INTERNAL_case_series_facts_{sfx}.json"), "w"), indent=1, ensure_ascii=False, default=str)
    print(json.dumps(dict(commit=pn["commit"][:8], last_report_date=chk["last_report_date"], last_complete_week=chk["last_complete_week_ending"], weeks=len(wk),
                          regression={k: (v["values_compared"], v["values_different"], v["rows_added"]) for k, v in facts["regression"]["by_run"].items()},
                          post_rows=len(PT), files=sorted(os.listdir(O))), indent=1))


if __name__ == "__main__":
    main()
