#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Weekly confirmed cases by province and by group of health zones, and reproduction numbers from
the growth of the weekly counts (2026 outbreak of Bundibugyo virus, Democratic Republic of the Congo).

Source: situation reports of the Institut National de Sante Publique (INSP), as transcribed in the public
repository INRB-UMIE/BDBV2026-Data. The script reads three files of one commit of that repository and
nothing else; it makes no request to a network:

    <indir>/data/insp_sitrep/processed/insp_sitrep__cumulative_confirmed_cases__daily.csv
    <indir>/data/insp_sitrep/processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv
    <indir>/data/shapefiles/DRC_Health_zones.dbf

Rules of processing: those of the round of 24 September 2026, taken from the code that built
weekly_confirmed_cases_by_province_20260924.csv, mean_weekly_cases_*_20260924.csv and
case_growth_comparators_20260924.csv and from compute_case_growth_seven_zones_20260925.py
(treatments of the weeks with delayed attribution).
What is new here: no date, no list of health zones and no file name is fixed in the code; they are
arguments (see `--help`). The list of health zones is read from a text file (`--zones-file`), so that
it can be replaced without a change of the code.

Example:
    python build_case_series_20261002.py --indir files_at_8eb57154 --outdir out --suffix 20261002 \
        --zones-file zones_of_analysed_set_20261002.txt --source-note "commit 8eb57154"
"""
import os

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "2")          # at most two processor cores

import argparse
import hashlib
import json
import re
import struct
import sys
import unicodedata

import numpy as np
import pandas as pd
import statsmodels.api as sm

# ----------------------------------------------------------------------------------------------
# Fixed rules of the round of 24 September 2026
# ----------------------------------------------------------------------------------------------
GT_MEAN, GT_SD = 15.3, 9.3                    # generation time of the article: gamma, mean and SD in days
FIRST_SUNDAY = "2026-05-17"                   # the first weekly bin ends one week later (24 May 2026)
BLOCK_BASE = "2026-05-24"                     # four-week periods are counted from 25 May 2026
POOL_TOLERANCE = 0.05                         # a pool of weeks is closed when |not attributed| <= 5 % of national
DEFAULT_ZONES = ["Bunia", "Rwampara", "Nizi", "Mongbwalu", "Lita", "Mangala", "Bambu"]
# spellings of one health zone that are merged, and provinces that the health-zone table does not give uniquely
PROVINCE_FIX = {"Gethy": ("Gety", "Ituri"), "Lubunga": ("Lubunga", "Tshopo"),
                "Lubunga (Tshopo)": ("Lubunga", "Tshopo"), "Rumba": ("Rimba", "Ituri")}
PROVINCE_SET = {"Gety": "Ituri", "Rimba": "Ituri"}
SPELLING_FIX = {"Nia Nia": "Nia-Nia", "Miti Murhesa": "Miti-Murhesa", "Makiso Kisangani": "Makiso-Kisangani"}
# treatment M: the two weeks that the table of 24 September replaced by their mean in the series of the groups
# of health zones (kept so that the rows of that table can be reproduced; the article uses treatment A)
MEAN_OF_WEEKS = ["2026-08-09", "2026-08-16"]
# ranges of weeks: (label, last Sunday before the range, last Sunday of the range)
RANGES_PROVINCES = [("W2-like: weeks ending 24 May-14 Jun", "2026-05-17", "2026-06-14"),
                    ("W3-like: weeks ending 21 Jun-30 Aug", "2026-06-14", "2026-08-30"),
                    ("weeks ending 21 Jun-2 Aug", "2026-06-14", "2026-08-02"),
                    ("weeks ending 9 Aug-13 Sep", "2026-08-02", "2026-09-13")]
RANGES_ZONE_GROUPS = [("weeks ending 21 Jun-2 Aug", "2026-06-14", "2026-08-02"),
                      ("weeks ending 9 Aug-13 Sep", "2026-08-02", "2026-09-13"),
                      ("weeks ending 21 Jun-30 Aug", "2026-06-14", "2026-08-30")]
# third range, defined on 2 October 2026 before the new counts were seen
THIRD_RANGE_AFTER = "2026-09-13"
THIRD_RANGE_MIN_WEEKS = 3
METHOD_24SEP = ("log-linear regression of weekly confirmed cases on time; R with generation time "
                "gamma(mean 15.3 d, SD 9.3 d); 95 % confidence interval of the slope")
METHOD_25SEP = ("log-linear regression of weekly confirmed cases on time (ordinary least squares); R with generation "
                "time gamma(mean 15.3 d, SD 9.3 d); 95 % confidence interval of the slope")
F_ZONES = "data/insp_sitrep/processed/insp_sitrep__cumulative_confirmed_cases__daily.csv"
F_NATIONAL = "data/insp_sitrep/processed/insp_sitrep__national_cumulative_confirmed_cases__daily.csv"
F_DBF = "data/shapefiles/DRC_Health_zones.dbf"


# ----------------------------------------------------------------------------------------------
# Helpers (unchanged from the earlier round)
# ----------------------------------------------------------------------------------------------
def R_from_r(r, mu=GT_MEAN, sd=GT_SD):
    """R = (1 + r sd^2/mu)^(mu^2/sd^2) for a gamma-distributed generation time."""
    k = (mu / sd) ** 2
    th = sd ** 2 / mu
    x = 1 + np.asarray(r, float) * th
    return np.where(x > 0, np.abs(x) ** k, 0.0)


def read_dbf(path, enc="utf-8"):
    b = open(path, "rb").read()
    nrec, hlen, rlen = struct.unpack("<xxxxIHH", b[:12])
    fields = []
    pos = 32
    while b[pos] != 0x0D:
        name = b[pos:pos + 11].split(b"\x00")[0].decode("ascii", "ignore")
        ftype = chr(b[pos + 11])
        flen = b[pos + 16]
        fields.append((name, ftype, flen))
        pos += 32
    rows = []
    off = hlen
    for _ in range(nrec):
        rec = b[off:off + rlen]
        off += rlen
        if rec[:1] == b"*":
            continue
        p = 1
        row = {}
        for name, ftype, flen in fields:
            raw = rec[p:p + flen]
            p += flen
            try:
                row[name] = raw.decode(enc).strip()
            except UnicodeDecodeError:
                row[name] = raw.decode("latin-1").strip()
        rows.append(row)
    return pd.DataFrame(rows)


def norm(s):
    s = unicodedata.normalize("NFKD", str(s)).encode("ascii", "ignore").decode().lower()
    s = re.sub(r"\(.*?\)", "", s)
    return re.sub(r"[^a-z0-9]", "", s)


def growth(series, a, b, digits=True):
    """Log-linear regression of the weekly counts of the weeks ending after `a` and up to `b`.
    Weeks with a count of 0 or less are left out (rule of the earlier round)."""
    s = series[(series.index > pd.Timestamp(a)) & (series.index <= pd.Timestamp(b))]
    n_in_range = len(s)
    s = s[s > 0]
    x = (s.index - s.index[0]).days.values.astype(float)
    y = np.log(s.values.astype(float))
    m = sm.OLS(y, sm.add_constant(x)).fit()
    r = float(m.params[1])
    ci = [float(v) for v in m.conf_int()[1]]
    out = dict(n_weeks=len(s), first_week_ending=str(s.index[0].date()), last_week_ending=str(s.index[-1].date()),
               r_per_day=r, r_lo=ci[0], r_hi=ci[1], R=float(R_from_r(r)), R_lo=float(R_from_r(ci[0])),
               R_hi=float(R_from_r(ci[1])), mean_weekly=float(s.mean()))
    if digits:                                 # rounding of the table of 24 September
        for k_ in ("r_per_day", "r_lo", "r_hi"):
            out[k_] = round(out[k_], 4)
        for k_ in ("R", "R_lo", "R_hi"):
            out[k_] = round(out[k_], 2)
        out["mean_weekly"] = round(out["mean_weekly"], 0)
    out["_weeks_in_range"] = n_in_range
    return out


def sha256_of(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def read_zone_list(path):
    if path is None:
        return list(DEFAULT_ZONES)
    zones = []
    for line in open(path, encoding="utf-8"):
        line = line.split("#")[0].strip()
        if line:
            zones.append(line)
    if not zones:
        sys.exit(f"STOP: the list of health zones {path} is empty")
    return zones


# ----------------------------------------------------------------------------------------------
# The series
# ----------------------------------------------------------------------------------------------
def build(indir, zone_list=None, max_report_date=None, zone_set_label="seven zones", zone_set_tag="seven_zones"):
    """Returns a dictionary of tables. Raises SystemExit with a message if a rule cannot be applied."""
    T = {}
    chk = {}
    fz, fn, fd = (os.path.join(indir, f) for f in (F_ZONES, F_NATIONAL, F_DBF))
    for f in (fz, fn, fd):
        if not os.path.exists(f):
            sys.exit(f"STOP: input file is missing: {f}")
    chk["input_files"] = {os.path.basename(f): dict(bytes=os.path.getsize(f), sha256=sha256_of(f)) for f in (fz, fn, fd)}

    # -- 1. read; dates and values as in the earlier round
    C = pd.read_csv(fz, dtype=str)
    Craw = pd.read_csv(fz, dtype=str, keep_default_na=False)        # only to report the literal strings
    N = pd.read_csv(fn, dtype=str)
    C["date_raw"] = C.date
    C["date"] = pd.to_datetime(C.date.str.extract(r"(\d{4}-\d{2}-\d{2})")[0], errors="coerce")
    C["v"] = pd.to_numeric(C.cumulative_confirmed_cases, errors="coerce")
    C["nom_literal"] = Craw.nom.values
    N["date_ts"] = pd.to_datetime(N.date, errors="coerce")
    N["v"] = pd.to_numeric(N.national_cumulative_confirmed_cases, errors="coerce")
    if N.date_ts.isna().any() or N.v.isna().any() or N.date_ts.duplicated().any():
        sys.exit("STOP: the national table holds a date or a value that cannot be read, or a date twice; "
                 "the rules of the earlier round do not cover this case")
    if C.date.isna().any():
        sys.exit("STOP: the table of the health zones holds a row without a readable date; "
                 "the rules of the earlier round do not cover this case")
    chk["rows_zone_table_all_dates"] = int(len(C))
    chk["rows_national_table_all_dates"] = int(len(N))
    chk["last_report_date_in_files"] = str(max(C.date.max(), N.date_ts.max()).date())
    if max_report_date is not None:
        C = C[C.date <= pd.Timestamp(max_report_date)].copy()
        N = N[N.date_ts <= pd.Timestamp(max_report_date)].copy()
    chk["restricted_to_report_dates_up_to"] = None if max_report_date is None else str(pd.Timestamp(max_report_date).date())
    chk["rows_zone_table_used_dates"] = int(len(C))
    chk["rows_national_table_used_dates"] = int(len(N))
    chk["report_dates_zone_table"] = int(C.date.nunique())
    chk["report_dates_national_table"] = int(N.date_ts.nunique())
    chk["first_report_date"] = str(min(C.date.min(), N.date_ts.min()).date())
    chk["last_report_date"] = str(max(C.date.max(), N.date_ts.max()).date())

    # -- 2. what is set aside or repaired (listed, not silently dropped)
    aside = []
    for _, r in C[C.date_raw != C.date.dt.strftime("%Y-%m-%d")].iterrows():
        aside.append(dict(kind="date string repaired, row used", name_in_source=r.nom_literal, date_in_source=r.date_raw,
                          date_used=str(r.date.date()), value_in_source=r.cumulative_confirmed_cases))
    for _, r in C[C.v.isna()].iterrows():
        aside.append(dict(kind="value not numeric, row set aside", name_in_source=r.nom_literal, date_in_source=r.date_raw,
                          date_used=str(r.date.date()), value_in_source=r.cumulative_confirmed_cases))
    for _, r in C[C.nom.isna() & C.v.notna()].iterrows():
        aside.append(dict(kind="no health zone named ('NA'), row not attributed to a province", name_in_source=r.nom_literal,
                          date_in_source=r.date_raw, date_used=str(r.date.date()), value_in_source=r.cumulative_confirmed_cases))
    T["set_aside"] = pd.DataFrame(aside, columns=["kind", "name_in_source", "date_in_source", "date_used", "value_in_source"])
    chk["date_strings_repaired"] = int((C.date_raw != C.date.dt.strftime("%Y-%m-%d")).sum())
    chk["values_not_numeric"] = {str(k): int(v) for k, v in C.cumulative_confirmed_cases[C.v.isna()].value_counts(dropna=False).items()}
    chk["rows_without_zone_name"] = int(C.nom.isna().sum())
    chk["duplicate_zone_date_rows"] = int(C[C.nom.notna()].duplicated(["nom", "date"]).sum())

    # -- 3. health zones to provinces
    Z = read_dbf(fd)
    Z["k"] = Z.Nom.map(norm)
    zmap = Z.groupby("k").PROVINCE.agg(lambda x: sorted(set(x)))
    noms = sorted(C.nom.dropna().unique())
    rowsm = []
    for n_ in noms:
        k = norm(n_)
        pv = zmap.get(k)
        rowsm.append(dict(nom=n_, key=k, province=(pv[0] if pv is not None and len(pv) == 1 else None),
                          candidates=(";".join(pv) if pv is not None else "")))
    MAP = pd.DataFrame(rowsm)
    prov = dict(zip(MAP.nom, MAP.province))
    canon = {n_: n_ for n_ in noms}
    for k_, (cn, pv) in PROVINCE_FIX.items():
        prov[k_] = pv
        canon[k_] = cn
    prov.update(PROVINCE_SET)
    canon.update(SPELLING_FIX)
    D = C[C.nom.notna() & C.v.notna()].copy()
    D["zone"] = D.nom.map(canon)
    D["province"] = D.nom.map(prov)
    if D.province.isna().any():
        sys.exit("STOP: no province for the health zone(s) " + ", ".join(sorted(D[D.province.isna()].nom.unique()))
                 + "; the health-zone table gives none or more than one, and the list of the earlier round does not hold them")
    names_used = sorted(D.nom.unique())
    T["zone_province_map"] = pd.DataFrame(dict(name_in_source=list(canon.keys()), zone=[canon[k] for k in canon],
                                               province=[prov.get(k) for k in canon]))
    chk["names_in_source"] = int(len(noms))
    chk["health_zones_after_merging"] = int(len(set(canon[n_] for n_ in noms)))
    chk["health_zones_with_a_count"] = int(D.zone.nunique())
    chk["health_zones_without_any_numeric_value"] = sorted(set(canon[n_] for n_ in noms) - set(D.zone))
    chk["spellings_merged"] = sorted(f"{k} = {v}" for k, v in canon.items() if k != v and k in noms)

    # -- 4. one value per zone and report date (largest over spellings), carried forward over the report dates
    P = D.pivot_table(index="date", columns="zone", values="v", aggfunc="max").sort_index()
    dates = pd.DatetimeIndex(sorted(set(P.index) | set(N.date_ts)))
    P = P.reindex(dates)
    viol = {z: int((P[z].dropna().diff() < 0).sum()) for z in P.columns if (P[z].dropna().diff() < 0).any()}
    chk["zones_with_a_decrease_of_the_cumulative_count"] = viol
    Pf = P.ffill().fillna(0)
    zprov = D.sort_values("date").drop_duplicates("zone", keep="last").set_index("zone").province
    grp = zprov.map(lambda p: "Ituri" if p == "Ituri" else "Nord-Kivu" if p == "Nord-Kivu" else "other provinces")
    S = pd.DataFrame({g: Pf[[z for z in Pf.columns if grp[z] == g]].sum(1) for g in ["Ituri", "Nord-Kivu", "other provinces"]})
    Nn = N.set_index("date_ts").v.astype("int64")
    Nn.index.name = "date"
    S["sum_of_zones"] = S.sum(1)
    S["national"] = Nn.reindex(S.index)
    S["not_attributed"] = S.national - S.sum_of_zones
    S = S[S.national.notna()]
    T["cumulative_by_province"] = S.copy()
    chk["report_dates_of_the_series"] = int(len(S))
    chk["report_dates_of_zone_table_without_national_value"] = sorted(str(d.date()) for d in set(P.dropna(how="all").index) - set(S.index))
    chk["national_cumulative_decreases_at"] = [str(d.date()) for d in S.index[S.national.diff() < 0]]

    # -- 5. daily grid by linear interpolation between report dates; differences at Sundays
    g = pd.date_range(S.index.min(), S.index.max(), freq="D")
    cols = ["Ituri", "Nord-Kivu", "other provinces", "not_attributed", "national"]
    I = pd.DataFrame({c: np.interp(g.values.astype("i8"), S.index.values.astype("i8"), S[c].values) for c in cols}, index=g)
    sund = pd.date_range(FIRST_SUNDAY, I.index.max(), freq="W-SUN")
    Wc = I.reindex(sund).diff().dropna()
    Wc.index.name = "week_ending_sunday"
    Wc["report_dates_in_week"] = [int(((S.index > d - pd.Timedelta(days=7)) & (S.index <= d)).sum()) for d in Wc.index]
    Wc = Wc.round(1)
    last_complete = Wc.index.max()
    chk["last_complete_week_ending"] = str(last_complete.date())
    chk["days_after_last_complete_week"] = int((S.index.max() - last_complete).days)

    # -- 6. pools of weeks and split of the national weekly count by province
    Wd = Wc.copy()
    grp_id = []
    g_ = 0
    acc_na = 0
    acc_nat = 0
    for d, r in Wd.iterrows():
        acc_na += r.not_attributed
        acc_nat += r.national
        grp_id.append(g_)
        if abs(acc_na) <= POOL_TOLERANCE * acc_nat:
            g_ += 1
            acc_na = 0
            acc_nat = 0
    Wd["pool"] = grp_id
    chk["last_pool_closed"] = bool(acc_nat == 0)
    out = []
    for gid, gdf in Wd.groupby("pool"):
        att = gdf[["Ituri", "Nord-Kivu", "other provinces"]].sum()
        sh = att / att.sum()
        for d, r in gdf.iterrows():
            out.append(dict(week_ending_sunday=d, total_drc=r.national, ituri=r.national * sh["Ituri"],
                            nord_kivu=r.national * sh["Nord-Kivu"], other_provinces=r.national * sh["other provinces"],
                            weeks_pooled_for_province_split=len(gdf), report_dates_in_week=int(r.report_dates_in_week),
                            raw_ituri=r.Ituri, raw_nord_kivu=r["Nord-Kivu"], raw_other=r["other provinces"],
                            raw_not_attributed=r.not_attributed))
    WKP = pd.DataFrame(out).round(1)
    WKP["nord_kivu_share"] = (WKP.nord_kivu / WKP.total_drc).round(3)
    WKP["complete_week"] = WKP.week_ending_sunday <= last_complete
    T["weekly_by_province"] = WKP.copy()
    pool_of_week = pd.Series([p + 1 for p in grp_id], index=Wc.index)

    # -- 7. groups of health zones
    zl = read_zone_list(None) if zone_list is None else list(zone_list)
    key2zone = {}
    for z in Pf.columns:
        key2zone.setdefault(norm(z), z)
    for k_, v_ in canon.items():
        key2zone.setdefault(norm(k_), v_)
    core, unknown = [], []
    for z in zl:
        zz = key2zone.get(norm(z))
        (core if zz is not None else unknown).append(zz if zz is not None else z)
    if unknown:
        sys.exit("STOP: health zone(s) of the list not found in the table of the health zones: " + ", ".join(unknown))
    if len(set(core)) != len(core):
        sys.exit("STOP: the list of health zones names a health zone twice")
    chk["zone_set"] = core
    chk["zone_set_label"] = zone_set_label
    chk["zone_set_provinces"] = {z: zprov[z] for z in core}
    chk["zone_set_all_in_ituri"] = bool(all(zprov[z] == "Ituri" for z in core))
    lab_other_it = "other Ituri zones"
    G2 = pd.DataFrame({zone_set_label: Pf[core].sum(1),
                       lab_other_it: Pf[[z for z in Pf.columns if grp[z] == "Ituri" and z not in core]].sum(1),
                       "Nord-Kivu": Pf[[z for z in Pf.columns if grp[z] == "Nord-Kivu" and z not in core]].sum(1),
                       "other provinces": Pf[[z for z in Pf.columns if grp[z] == "other provinces" and z not in core]].sum(1)})
    G2 = G2.loc[S.index]
    G2["national"] = S.national
    G2["not_attributed"] = G2.national - G2[[zone_set_label, lab_other_it, "Nord-Kivu", "other provinces"]].sum(1)
    gI = pd.DataFrame({c: np.interp(g.values.astype("i8"), G2.index.values.astype("i8"), G2[c].values) for c in G2.columns}, index=g)
    W2 = gI.reindex(sund).diff().dropna()
    W2.index.name = "week_ending_sunday"
    # four-week periods counted from 25 May; complete periods only
    ends = []
    e = pd.Timestamp(BLOCK_BASE) + pd.Timedelta(days=28)
    while e <= last_complete:
        ends.append(e)
        e = e + pd.Timedelta(days=28)
    ends = pd.DatetimeIndex(ends)
    block_labels = [f"{(e - pd.Timedelta(days=27)).strftime('%d %b')}-{e.strftime('%d %b')}" for e in ends]
    lev = gI.reindex(pd.to_datetime([BLOCK_BASE]).append(ends))
    B = lev.diff().dropna()
    B.index = block_labels
    Bout = (B / 4).round(1)
    Bout.insert(0, "weeks", 4)
    Bout.index.name = "four_week_period_2026"
    T["mean_weekly_by_zone_group"] = Bout
    lz = pd.DataFrame({z: np.interp(g.values.astype("i8"), S.index.values.astype("i8"), Pf.loc[S.index, z].values) for z in core},
                      index=g).reindex(pd.to_datetime([BLOCK_BASE]).append(ends)).diff().dropna()
    lz.index = block_labels
    T["mean_weekly_zone_set"] = (lz / 4).round(1).rename_axis("four_week_period_2026")
    chk["weeks_after_last_complete_four_week_period"] = int(round((last_complete - (ends[-1] if len(ends) else pd.Timestamp(BLOCK_BASE))).days / 7))

    # -- 8. weekly series of the set of health zones, with the treatments of the weeks with delayed attribution
    Sz = pd.DataFrame({"set": Pf[core].sum(1), "all_zones": Pf.sum(1)})
    Sz["national"] = Nn.reindex(Sz.index)
    Sz = Sz[Sz.national.notna()]
    Iz = pd.DataFrame({c: np.interp(g.values.astype("i8"), Sz.index.values.astype("i8"), Sz[c].values) for c in Sz.columns}, index=g)
    Wz = Iz.reindex(sund).diff().dropna()
    assert (Wz.index == Wc.index).all() and float((Wz.national - WKP.set_index("week_ending_sunday").total_drc).abs().max()) < 0.051
    _wk = WKP.set_index("week_ending_sunday")
    # Consistency of the unrounded sum of all health zones with the three province columns of the weekly table.
    # Each of the three columns is rounded to one decimal, so that their sum may differ from the unrounded sum
    # by up to 0.15. The script of 25 September allowed 0.051, which held as long as all weekly counts were
    # multiples of 0.5; it does not hold for a Sunday that lies between two reports three days apart.
    chk["max_difference_all_zones_unrounded_vs_rounded_province_columns"] = round(
        float((Wz.all_zones - (_wk.raw_ituri + _wk.raw_nord_kivu + _wk.raw_other)).abs().max()), 4)
    assert chk["max_difference_all_zones_unrounded_vs_rounded_province_columns"] < 0.151
    tag = zone_set_tag
    CASES = pd.DataFrame({"national_confirmed": Wz.national, f"{tag}_as_attributed": Wz["set"],
                          "all_zones_as_attributed": Wz.all_zones, "pool_of_weeks": pool_of_week.values}, index=Wz.index)
    _s = CASES.groupby("pool_of_weeks")[[f"{tag}_as_attributed", "all_zones_as_attributed"]].transform("sum")
    CASES[f"share_of_{tag}_in_pool"] = _s[f"{tag}_as_attributed"] / _s.all_zones_as_attributed
    CASES[f"{tag}_treatment_A"] = CASES.national_confirmed * CASES[f"share_of_{tag}_in_pool"]
    CASES[f"{tag}_treatment_B"] = CASES.groupby("pool_of_weeks")[f"{tag}_as_attributed"].transform("mean")
    CASES[f"{tag}_treatment_M"] = CASES[f"{tag}_as_attributed"]
    adw = [pd.Timestamp(d) for d in MEAN_OF_WEEKS if pd.Timestamp(d) in CASES.index]
    if len(adw) == len(MEAN_OF_WEEKS):
        CASES.loc[adw, f"{tag}_treatment_M"] = CASES.loc[adw, f"{tag}_as_attributed"].mean()
    # treatment A, applied to the attributed cases of Ituri, gives the series of Ituri of the weekly table
    _ci = pd.DataFrame({"raw": _wk.raw_ituri, "all": _wk.raw_ituri + _wk.raw_nord_kivu + _wk.raw_other, "p": pool_of_week.values})
    _cs = _ci.groupby("p")[["raw", "all"]].transform("sum")
    assert float((_wk.total_drc * _cs.raw / _cs["all"] - _wk.ituri).abs().max()) < 0.06
    CASES.index.name = "week_ending_sunday"
    T["weekly_zone_set"] = CASES.copy()
    # the same for the other health zones of Ituri (treatment A and the treatment of 24 September)
    OTH = pd.DataFrame({"as_attributed": W2[lab_other_it], "all": Wz.all_zones.values, "national": Wz.national.values,
                        "p": pool_of_week.values}, index=W2.index)
    _so = OTH.groupby("p")[["as_attributed", "all"]].transform("sum")
    OTH["treatment_A"] = OTH.national * _so.as_attributed / _so["all"]
    OTH["treatment_M"] = OTH.as_attributed
    if len(adw) == len(MEAN_OF_WEEKS):
        OTH.loc[adw, "treatment_M"] = OTH.loc[adw, "as_attributed"].mean()
    # mean weekly cases of the set by four-week period, for each treatment
    mb = []
    for lab_, e in zip(block_labels, ends):
        w = CASES[(CASES.index > e - pd.Timedelta(days=28)) & (CASES.index <= e)]
        mb.append({"four_week_period_2026": lab_, "weeks": len(w), "as_attributed": w[f"{tag}_as_attributed"].mean(),
                   "treatment_A": w[f"{tag}_treatment_A"].mean(), "treatment_B": w[f"{tag}_treatment_B"].mean(),
                   "national": w.national_confirmed.mean()})
    T["mean_weekly_zone_set_by_treatment"] = pd.DataFrame(mb).set_index("four_week_period_2026")

    # -- 9. growth of the weekly counts and reproduction numbers
    Wk = WKP.set_index("week_ending_sunday")
    third_weeks = Wk.index[Wk.index > pd.Timestamp(THIRD_RANGE_AFTER)]
    n3 = len(third_weeks)
    chk["third_range_weeks"] = [str(d.date()) for d in third_weeks]
    chk["third_range_given"] = bool(n3 >= THIRD_RANGE_MIN_WEEKS)
    third = None
    if n3 >= 1:
        f3, l3 = third_weeks[0], third_weeks[-1]
        third = (f"weeks ending {f3.day} {f3.strftime('%b')}-{l3.day} {l3.strftime('%b')}", THIRD_RANGE_AFTER, str(l3.date()))
    area_label = zone_set_label

    def not_estimated(period, area, treatment, n_):
        return dict(period=period, area=area, n_weeks=n_, first_week_ending=str(third_weeks[0].date()) if n_ else "",
                    last_week_ending=str(third_weeks[-1].date()) if n_ else "", r_per_day=np.nan, r_lo=np.nan, r_hi=np.nan,
                    R=np.nan, R_lo=np.nan, R_hi=np.nan, mean_weekly=np.nan, treatment=treatment,
                    status=f"not estimated: {n_} complete week(s) after {pd.Timestamp(THIRD_RANGE_AFTER).strftime('%d %b %Y')}, "
                           f"fewer than {THIRD_RANGE_MIN_WEEKS}")

    def est(series, lab, a, b, area, treatment):
        o = growth(series, a, b)
        n_range = o.pop("_weeks_in_range")
        st = "estimated" if o["n_weeks"] == n_range else f"estimated; {n_range - o['n_weeks']} week(s) with a count of 0 or less left out"
        return dict(period=lab, area=area, **o, treatment=treatment, status=st)

    TR_PROV = "split by province within pools of weeks (weekly table of the provinces)"
    TR_NAT = "national weekly count (no split)"
    TR_A = "A: national weekly count multiplied by the share of the group among the attributed cases of the pool of weeks"
    TR_DEL = "M: weeks ending 9 and 16 August replaced by their mean; all other weeks as attributed"
    rows = []
    for lab, a, b in RANGES_PROVINCES + ([third] if third else []):
        is_third = third is not None and (lab, a, b) == third
        for col, nm in [("total_drc", "DRC"), ("ituri", "Ituri"), ("nord_kivu", "Nord-Kivu")]:
            tr = TR_NAT if nm == "DRC" else TR_PROV
            if is_third and n3 < THIRD_RANGE_MIN_WEEKS:
                rows.append(not_estimated(lab, nm, tr, n3))
            else:
                rows.append(est(Wk[col], lab, a, b, nm, tr))
    for lab, a, b in RANGES_ZONE_GROUPS + ([third] if third else []):
        is_third = third is not None and (lab, a, b) == third
        for nm, sA, sD in [(area_label, CASES[f"{tag}_treatment_A"], CASES[f"{tag}_treatment_M"]),
                           (lab_other_it, OTH.treatment_A, OTH.treatment_M)]:
            for tr, s in [(TR_A, sA), (TR_DEL, sD)]:
                if is_third and n3 < THIRD_RANGE_MIN_WEEKS:
                    rows.append(not_estimated(lab, nm, tr, n3))
                else:
                    rows.append(est(s, lab, a, b, nm, tr))
    CG = pd.DataFrame(rows)
    CG["method"] = METHOD_24SEP
    CG = CG[["period", "area", "n_weeks", "first_week_ending", "last_week_ending", "r_per_day", "r_lo", "r_hi", "R", "R_lo", "R_hi",
             "mean_weekly", "method", "treatment", "status"]]
    T["growth"] = CG

    TRT = [("A", f"{tag}_treatment_A", "national weekly count multiplied by the share of the seven health zones among the attributed cases of the pool "
            "of weeks (pools as in the weekly table of the provinces); the treatment of the series of the provinces"),
           ("B", f"{tag}_treatment_B", "every week of a pool replaced by the mean of the attributed cases of the seven health zones in the pool"),
           ("M", f"{tag}_treatment_M", "weeks ending 9 and 16 August replaced by their mean; all other weeks as attributed"),
           ("none", f"{tag}_as_attributed", "weekly differences of the attributed cumulative counts, no treatment")]
    if zone_set_label != "seven zones":
        TRT = [(a_, b_, c_.replace("the seven health zones", "the health zones of the set")) for a_, b_, c_ in TRT]
    area_long = (f"seven health zones ({', '.join(core)})" if zone_set_label == "seven zones"
                 else f"{zone_set_label} ({', '.join(core)})")
    rows7 = []
    for nm, col, desc in TRT:
        for lab, a, b in RANGES_ZONE_GROUPS + ([third] if third else []):
            is_third = third is not None and (lab, a, b) == third
            if is_third and n3 < THIRD_RANGE_MIN_WEEKS:
                r_ = not_estimated(lab, area_long, nm, n3)
                r_.pop("treatment")
                rows7.append(dict(treatment=nm, treatment_description=desc, **r_))
            else:
                o = growth(CASES[col], a, b, digits=False)
                n_range = o.pop("_weeks_in_range")
                st = "estimated" if o["n_weeks"] == n_range else f"estimated; {n_range - o['n_weeks']} week(s) with a count of 0 or less left out"
                rows7.append(dict(treatment=nm, treatment_description=desc, period=lab, area=area_long, **o, status=st))
    G7 = pd.DataFrame(rows7)
    G7["method"] = METHOD_25SEP
    G7 = G7[["treatment", "treatment_description", "period", "area", "n_weeks", "first_week_ending", "last_week_ending", "r_per_day",
             "r_lo", "r_hi", "R", "R_lo", "R_hi", "mean_weekly", "method", "status"]]
    T["growth_zone_set_treatments"] = G7

    # -- 10. cross-check of the national totals within the repository
    X = S[["national", "sum_of_zones", "not_attributed"]].copy()
    X.columns = ["national_table", "sum_of_health_zones", "national_minus_health_zones"]
    X["difference_percent_of_national"] = (100 * X.national_minus_health_zones / X.national_table).round(2)
    nz = D.groupby("date").zone.nunique()
    X["health_zones_with_a_value_in_the_report"] = nz.reindex(X.index).fillna(0).astype(int)
    X["health_zones_carried_forward"] = [int((Pf.loc[d] > 0).sum() - (P.loc[d].notna() & (Pf.loc[d] > 0)).sum()) for d in X.index]
    X.index.name = "report_date"
    T["crosscheck"] = X
    wk_sum = float(WKP.total_drc.sum())
    cum_diff = float(I.national.loc[last_complete] - I.national.loc[pd.Timestamp(FIRST_SUNDAY)])
    chk["sum_of_weekly_national_counts"] = round(wk_sum, 1)
    chk["national_cumulative_last_sunday_minus_first_sunday"] = round(cum_diff, 1)
    chk["weekly_province_columns_sum_to_national_max_abs_difference"] = round(
        float((WKP.ituri + WKP.nord_kivu + WKP.other_provinces - WKP.total_drc).abs().max()), 2)
    chk["report_dates_with_national_equal_to_sum_of_zones"] = int((X.national_minus_health_zones == 0).sum())
    chk["report_dates_with_national_different_from_sum_of_zones"] = int((X.national_minus_health_zones != 0).sum())
    chk["last_report_date_with_a_difference"] = (str(X.index[X.national_minus_health_zones != 0].max().date())
                                                 if (X.national_minus_health_zones != 0).any() else None)
    chk["weekly_counts_nonnegative"] = bool((WKP[["total_drc", "ituri", "nord_kivu", "other_provinces"]] >= 0).all().all())
    last = S.index.max()
    chk["cumulative_at_last_report_date"] = dict(date=str(last.date()), national=int(S.national.loc[last]), ituri=int(S.Ituri.loc[last]),
                                                 nord_kivu=int(S["Nord-Kivu"].loc[last]), other_provinces=int(S["other provinces"].loc[last]),
                                                 not_attributed=int(S.not_attributed.loc[last]),
                                                 health_zones_with_a_value_in_the_last_report=int(nz.reindex([last]).fillna(0).iloc[0]))
    T["checks"] = chk
    T["_internal"] = dict(Pf=Pf, P=P, S=S, I=I, W2=W2, grp=grp, zprov=zprov, D=D, OTH=OTH)
    return T


def write(T, outdir, suffix, zone_set_tag="seven_zones", source_note=""):
    os.makedirs(outdir, exist_ok=True)
    p = lambda name: os.path.join(outdir, name)
    files = {}

    def w(key, name, **kw):
        T[key].to_csv(p(name), **kw)
        files[key] = name
    w("weekly_by_province", f"weekly_confirmed_cases_by_province_{suffix}.csv", index=False)
    w("cumulative_by_province", f"cumulative_confirmed_cases_by_province_{suffix}.csv")
    w("zone_province_map", f"zone_province_map_{suffix}.csv", index=False)
    w("mean_weekly_by_zone_group", f"mean_weekly_cases_by_zone_group_{suffix}.csv")
    w("mean_weekly_zone_set", f"mean_weekly_cases_{zone_set_tag}_{suffix}.csv")
    # three decimals, so that a value printed without decimals is not rounded twice
    T["mean_weekly_zone_set_by_treatment"].round(3).to_csv(p(f"mean_weekly_cases_{zone_set_tag}_by_treatment_{suffix}.csv"))
    files["mean_weekly_zone_set_by_treatment"] = f"mean_weekly_cases_{zone_set_tag}_by_treatment_{suffix}.csv"
    T["weekly_zone_set"].round(3).to_csv(p(f"weekly_confirmed_cases_{zone_set_tag}_{suffix}.csv"))
    files["weekly_zone_set"] = f"weekly_confirmed_cases_{zone_set_tag}_{suffix}.csv"
    G = T["growth"].copy()
    G7 = T["growth_zone_set_treatments"].copy()
    if source_note:
        G7["source"] = source_note
    G.to_csv(p(f"case_growth_comparators_{suffix}.csv"), index=False)
    files["growth"] = f"case_growth_comparators_{suffix}.csv"
    G7.to_csv(p(f"case_growth_{zone_set_tag}_treatments_{suffix}.csv"), index=False)
    files["growth_zone_set_treatments"] = f"case_growth_{zone_set_tag}_treatments_{suffix}.csv"
    w("crosscheck", f"case_series_crosscheck_{suffix}.csv")
    w("set_aside", f"case_series_values_set_aside_{suffix}.csv", index=False)
    chk = dict(T["checks"])
    chk["source_note"] = source_note
    chk["files_written"] = files
    json.dump(chk, open(p(f"case_series_checks_{suffix}.json"), "w"), indent=1, ensure_ascii=False)
    files["checks"] = f"case_series_checks_{suffix}.json"
    return files


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--indir", required=True, help="folder with the files of one commit of INRB-UMIE/BDBV2026-Data")
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--suffix", default="20261002", help="date in the names of the tables")
    ap.add_argument("--zones-file", default=None, help="text file, one health zone per line; default: the seven zones of the earlier round")
    ap.add_argument("--zone-set-label", default="seven zones", help="name of the set of health zones in the tables")
    ap.add_argument("--zone-set-tag", default="seven_zones", help="name of the set of health zones in file and column names")
    ap.add_argument("--max-report-date", default=None, help="use report dates up to this date only (regression against an earlier round)")
    ap.add_argument("--source-note", default="", help="text of the column 'source', for example the commit")
    a = ap.parse_args(argv)
    T = build(a.indir, zone_list=read_zone_list(a.zones_file), max_report_date=a.max_report_date,
              zone_set_label=a.zone_set_label, zone_set_tag=a.zone_set_tag)
    files = write(T, a.outdir, a.suffix, zone_set_tag=a.zone_set_tag, source_note=a.source_note)
    c = T["checks"]
    print(json.dumps(dict(last_report_date=c["last_report_date"], last_complete_week_ending=c["last_complete_week_ending"],
                          weeks=len(T["weekly_by_province"]), third_range_weeks=c["third_range_weeks"],
                          third_range_given=c["third_range_given"], files=sorted(files.values())), indent=1))


if __name__ == "__main__":
    main()
