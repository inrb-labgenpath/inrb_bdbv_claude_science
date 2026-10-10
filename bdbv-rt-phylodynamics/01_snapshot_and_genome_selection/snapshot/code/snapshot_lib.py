"""Derivations for the Pathoplexus BDBV snapshots.

The definitions follow the code that produced the snapshot of 24 September 2026
(lineage of pp_snapshot_20260924_metadata.csv and snapshot_validation_20260924.csv).
Differences against that code:
  * health zones are assigned by lookup in the map of 24 September (zone_harmonisation_map_20260924.csv),
    key = (geoLocCountry, geoLocAdmin1, geoLocAdmin2) as recorded; a key that the map does not hold
    is flagged and receives no zone;
  * a sample identifier that occurs more than twice gives all pairs (the earlier code asserted two);
  * two columns are added: n_acgt_1_18900 and n_acgt_1_18940;
  * the baseline of the comparison columns is passed in (None = columns left empty);
  * the label of a genome without collection month names the reason (range / missing / other);
  * for a record whose spelling the map does not hold, the four flags of the map are left empty (not assessed);
  * rule on blanks: the character U+0020 at the beginning or end of
    geoLocCountry, geoLocAdmin1 and geoLocAdmin2 is removed before the map is applied; a record on which the rule acted
    carries the rule in zone_rule and the flag flag_zone_interpreted; the database fields themselves stay as recorded.
No network access. No random numbers.
"""
import gzip
import hashlib
import itertools
import re

import numpy as np
import pandas as pd

L_CALL = 18900
L_ALN = 18940

DERIVED16 = ["snapshotDatetimeUtc", "lapisDataVersion", "outbreak2026", "latest_nonrevoked", "called_frac",
             "collection_date_iso8601", "collection_date_precision", "zone", "zone_rule", "zone_ref_province",
             "flag_province_conflict", "flag_zone_not_in_reference", "flag_zone_interpreted", "flag_zone_missing",
             "duplicate_sample_id", "flags"]
NEW_COLS = ["n_acgt_1_18900", "n_acgt_1_18940"]

RULE_NOT_IN_MAP = "raw value not in the zone map (flagged, no zone assigned)"
LABEL_NOT_IN_MAP = "(spelling not in the zone map)"
RULE_BLANK = "blank at the beginning or end removed, then "
BLANK = " "            # U+0020; the rule removes this character and no other
PLACE_FIELDS = ["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "geoLocCity", "geoLocSite", "geoLocLatitude", "geoLocLongitude",
                "hostOriginCountry", "ncbiSubmitterCountry"]


def outer_blanks(df, fields=None):
    """Every record and place field with the character U+0020 at the beginning or end, and every record and place field
    with another white-space character there (which the rule does not remove). One row per record and field."""
    rows = []
    for c in (fields or [f for f in PLACE_FIELDS if f in df.columns]):
        for av, v in zip(df["accessionVersion"], df[c]):
            if v != v.strip():
                lead, trail = v[:len(v) - len(v.lstrip())], v[len(v.rstrip()):]
                rows.append(dict(accessionVersion=av, field=c, recorded=repr(v), after_removal=repr(v.strip(BLANK)),
                                 characters_at_beginning=";".join(f"U+{ord(ch):04X}" for ch in lead),
                                 characters_at_end=";".join(f"U+{ord(ch):04X}" for ch in trail),
                                 rule_acts=v != v.strip(BLANK),
                                 other_white_space_left=v.strip(BLANK) != v.strip()))
    return pd.DataFrame(rows, columns=["accessionVersion", "field", "recorded", "after_removal", "characters_at_beginning",
                                       "characters_at_end", "rule_acts", "other_white_space_left"])


# ---------------------------------------------------------------- files
def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_str(s):
    return hashlib.sha256(s.encode()).hexdigest()


def read_fasta(path):
    """Return (records, info). records = list of (header, sequence) in file order, nothing dropped."""
    opener = gzip.open if str(path).endswith(".gz") else open
    recs, name, buf = [], None, []
    info = dict(n_lines=0, n_crlf=0, n_blank=0, n_lines_before_first_header=0)
    with opener(path, "rb") as fh:
        for raw in fh:
            info["n_lines"] += 1
            if raw.endswith(b"\r\n"):
                info["n_crlf"] += 1
            line = raw.rstrip(b"\r\n")
            if line.startswith(b">"):
                if name is not None:
                    recs.append((name, b"".join(buf).decode("ascii")))
                name, buf = line[1:].decode("ascii"), []
            elif line == b"":
                info["n_blank"] += 1
            else:
                if name is None:
                    info["n_lines_before_first_header"] += 1
                else:
                    buf.append(line.strip())
        if name is not None:
            recs.append((name, b"".join(buf).decode("ascii")))
    return recs, info


def write_fasta_gz(path, seqs):
    """Header = accessionVersion only, sorted, one line per sequence; gzip without time stamp or file name."""
    with open(path, "wb") as raw_fh:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw_fh, mtime=0) as gz:
            for av in sorted(seqs):
                gz.write(f">{av}\n{seqs[av]}\n".encode("ascii"))


def read_csv_str(path):
    return pd.read_csv(path, dtype=str, keep_default_na=False)


# ---------------------------------------------------------------- dates
def norm_date(s):
    if pd.isna(s) or s == "":
        return pd.NaT
    s = str(s)[:10]
    if len(s) == 4:
        s += "-07-01"
    elif len(s) == 7:
        s += "-15"
    return pd.to_datetime(s, errors="coerce")


def precision(s):
    if s == "":
        return "missing"
    if "/" in s:
        return "range"
    return {4: "year", 7: "month", 10: "day"}.get(len(s), "other")


def iso(series):
    return pd.to_datetime(series, format="ISO8601", errors="coerce")


# ---------------------------------------------------------------- outbreak definition (section 2 of the report of 24 Sep)
def outbreak_definition(d, parser="first10"):
    """Rows of d (all versions) that belong to the 2026-outbreak set, in the order of d.
    1. versionStatus == LATEST_VERSION
    2. collection date >= 2026-01-01, OR collection date not interpretable as a single date and earliestReleaseDate >= 2026-05-01
    3. not a revocation entry
    parser='first10': collection date from the first 10 characters (year -> 1 July, year-month -> 15th, range -> lower bound)
    parser='iso8601': pandas.to_datetime(format='ISO8601', errors='coerce')
    Returns (rows, arm) with arm = 'collection' or 'release' per row."""
    if parser == "first10":
        coll = d["sampleCollectionDate"].map(norm_date)
    elif parser == "iso8601":
        coll = iso(d["sampleCollectionDate"])
    else:
        raise ValueError(parser)
    rel = d["earliestReleaseDate"].map(norm_date)
    latest = d["versionStatus"] == "LATEST_VERSION"
    arm_coll = coll >= pd.Timestamp("2026-01-01")
    arm_rel = coll.isna() & (rel >= pd.Timestamp("2026-05-01"))
    not_rev = d["isRevocation"].str.lower() != "true"
    m = latest & (arm_coll | arm_rel) & not_rev
    arm = pd.Series(np.where(arm_coll, "collection", np.where(arm_rel, "release", "")), index=d.index)
    return d[m], arm[m]


# ---------------------------------------------------------------- sequences
def count_acgt(s, upto):
    seg = s[:upto].upper()
    return seg.count("A") + seg.count("C") + seg.count("G") + seg.count("T")


def called_frac(s):
    return count_acgt(s, L_CALL) / L_CALL


def seq_counts(s):
    seg = s[:L_CALL].upper()
    n = count_acgt(s, L_CALL)
    n_n, n_gap = seg.count("N"), seg.count("-")
    return dict(n_called=n, n_N=n_n, n_gap=n_gap, n_ambiguous=L_CALL - n - n_n - n_gap)


ACGT = np.array([b"A", b"C", b"G", b"T"])


def pair_dist(sa, sb, upto=L_CALL):
    x = np.frombuffer(sa[:upto].upper().encode(), dtype="S1")
    y = np.frombuffer(sb[:upto].upper().encode(), dtype="S1")
    ok = np.isin(x, ACGT) & np.isin(y, ACGT)
    return int((x[ok] != y[ok]).sum()), int(ok.sum())


# ---------------------------------------------------------------- health zones by the zone map
ZONE_FLAGS = ["flag_province_conflict", "flag_zone_not_in_reference", "flag_zone_missing", "flag_zone_interpreted"]


def load_zone_map(path):
    zm = read_csv_str(path)
    key = list(zip(zm.geoLocCountry, zm.geoLocAdmin1, zm.geoLocAdmin2))
    assert len(set(key)) == len(key), "stored zone map: key not unique"
    for c in ZONE_FLAGS:
        assert set(zm[c]) <= {"True", "False"}, c
    look = {}
    for k, (_, r) in zip(key, zm.iterrows()):
        look[k] = dict(zone=r.zone, zone_rule=r.zone_rule, zone_ref_province=r.zone_ref_province,
                       **{c: r[c] == "True" for c in ZONE_FLAGS})
    return zm, look


def harmonise_by_map(o, look):
    rows = []
    for c, a1, a2 in zip(o.geoLocCountry, o.geoLocAdmin1, o.geoLocAdmin2):
        key = (c.strip(BLANK), a1.strip(BLANK), a2.strip(BLANK))      # rule on blanks, applied before the map
        acted = key != (c, a1, a2)
        hit = look.get(key)
        if hit is not None and acted:
            hit = dict(hit, zone_rule=RULE_BLANK + hit["zone_rule"], flag_zone_interpreted=True)
        if hit is None:
            # the map has no entry: no zone, and the four flags of the map are NOT ASSESSED (left empty), not False
            rows.append(dict(zone="", zone_rule=RULE_NOT_IN_MAP, zone_ref_province="",
                             flag_province_conflict=None, flag_zone_not_in_reference=None,
                             flag_zone_missing=None, flag_zone_interpreted=None, flag_zone_not_in_map=True))
        else:
            rows.append(dict(hit, flag_zone_not_in_map=False))
        rows[-1]["rule_on_blanks_acted"] = acted
    hz = pd.DataFrame(rows, index=o.index)
    for c in ZONE_FLAGS:
        hz[c] = hz[c].astype(object).where(hz[c].notna(), None)
    return hz


def is_true(v):
    """True only for a real True; a flag that was not assessed (None/NaN) is not counted as set."""
    return (v is True) or (isinstance(v, (bool, np.bool_)) and bool(v))


# ---------------------------------------------------------------- rule of the code of 24 September (information only)
REF_ITURI = ["Adi", "Adja", "Angumu", "Ariwara", "Aru", "Aungba", "Bambu", "Biringi", "Boga", "Bunia", "Damas", "Drodro", "Fataki",
             "Gety", "Jiba", "Kambala", "Kilo", "Komanda", "Laybo", "Linga", "Lita", "Logo", "Lolwa", "Mahagi", "Mambasa", "Mandima",
             "Mangala", "Mongbwalu", "Nia-Nia", "Nizi", "Nyankunde", "Nyarambe", "Rethy", "Rimba", "Rwampara", "Tchomia"]
REF_NK = ["Alimbongo", "Bambo", "Beni", "Biena", "Binza", "Birambizo", "Butembo", "Goma", "Itebero", "Kalunguta", "Kamango",
          "Karisimbi", "Katoyi", "Katwa", "Kayna", "Kibirizi", "Kibua", "Kirotshe", "Kyondo", "Lubero", "Mabalako", "Manguredjipa",
          "Masereka", "Masisi", "Musienene", "Mutwanga", "Mweso", "Nyiragongo", "Oicha", "Pinga", "Rutshuru", "Rwanguba", "Vuhovi",
          "Walikale"]
REF_SK = ["Bagira", "Bunyakiri", "Fizi", "Hauts-Plateaux", "Ibanda", "Idjwi", "Itombwe", "Kabare", "Kadutu", "Kalehe", "Kalole",
          "Kalonge", "Kamituga", "Kaniola", "Katana", "Kaziba", "Kimbi-Lulenge", "Kitutu", "Lemera", "Lulingu", "Minembwe", "Minova",
          "Miti-Murhesa", "Mubumbano", "Mulungu", "Mwana", "Mwenga", "Nundu", "Nyangezi", "Nyantende", "Ruzizi", "Shabunda", "Uvira",
          "Walungu"]


def _key(s):
    return re.sub(r"[^A-Z]", "", s.upper())


_REF = {}
for _prov, _lst in (("Ituri", REF_ITURI), ("Nord-Kivu", REF_NK), ("Sud-Kivu", REF_SK)):
    for _z in _lst:
        _REF[_key(_z)] = (_z, _prov)
_VARIANTS = {"MONGWALU": ("Mongbwalu", "spelling variant"), "MUNGWALU": ("Mongbwalu", "spelling variant"),
             "NYAKUNDE": ("Nyankunde", "spelling variant"), "GETHY": ("Gety", "spelling variant"),
             "NIA": ("Nia-Nia", "truncated value (interpreted)")}
_NOT_PLACE = {"NOTPROVIDED", "UNKNOWN", "NA", "MISSING", ""}
_OTHER_COUNTRY_REF = {"DRC", "BUNYACONGO"}


def rule_of_20260924_code(country, admin1, admin2):
    """The function harmonise() of the code that produced zone_harmonisation_map_20260924.csv, transcribed from the
    code of the earlier round. Used here ONLY to validate the transcription against the map of 24 September and to report,
    for a spelling that the map does not hold, what that rule would give. Its result is NOT written to the column zone."""
    k = _key(admin2)
    res = dict(zone="", zone_rule="", zone_ref_province="", flag_province_conflict=False, flag_zone_not_in_reference=False,
               flag_zone_missing=False, flag_zone_interpreted=False)
    if k in _NOT_PLACE:
        res.update(zone_rule="missing" if admin2 == "" else "placeholder text -> missing", flag_zone_missing=True)
        return res
    if country == "Democratic Republic of the Congo":
        if k in _REF:
            z, prov = _REF[k]
            rule = "exact" if admin2 == z else "case/separator variant"
        elif k in _VARIANTS:
            z, rule = _VARIANTS[k]
            prov = _REF[_key(z)][1]
            res["flag_zone_interpreted"] = True
        else:
            res.update(zone=admin2.strip().title(), zone_rule="not in reference list (raw value kept, title-cased)",
                       flag_zone_not_in_reference=True)
            return res
        res.update(zone=z, zone_rule=rule, zone_ref_province=prov)
        if admin1 != "" and _key(admin1) != _key(prov):
            res["flag_province_conflict"] = True
        return res
    if country == "Uganda":
        if k in _OTHER_COUNTRY_REF:
            res.update(zone_rule="admin2 names another country (raw value kept in geoLocAdmin2) -> missing", flag_zone_missing=True)
            return res
        res.update(zone=admin2.strip().title(), zone_rule="Ugandan district (raw value, case-normalised)")
        return res
    res.update(zone_rule="outside DRC/Uganda -> missing", flag_zone_missing=True)
    return res


# ---------------------------------------------------------------- pairs of records that share a sample identifier
def norm_sid(s):
    return re.sub(r"[^A-Z0-9]", "", s.upper())


DUP_COLS = ["specimenCollectorSampleId", "n_accessions", "accessionVersion_1", "accessionVersion_2",
            "earliestReleaseDate_1", "earliestReleaseDate_2", "sampleCollectionDate_1", "sampleCollectionDate_2",
            "same_collection_date", "zone_1", "zone_2", "same_zone", "sequencingProtocol_1", "sequencingProtocol_2",
            "sequencingDate_1", "sequencingDate_2", "called_frac_1", "called_frac_2", "n_sites_both_called",
            "n_differences", "identical_aligned_sequence", "groupId", "same_group", "dataUseTerms_1", "dataUseTerms_2"]


def duplicate_pairs(o, seqs):
    rows = []
    sub = o[(o.specimenCollectorSampleId != "") & o.specimenCollectorSampleId.duplicated(keep=False)]
    for sidv, g in sub.groupby("specimenCollectorSampleId"):
        g = g.sort_values(["release_date", "accession"])
        for i, j in itertools.combinations(range(len(g)), 2):
            a, b = g.index[i], g.index[j]
            if a in seqs and b in seqs:
                nd_, nc_ = pair_dist(seqs[a], seqs[b])
                ident = seqs[a] == seqs[b]
            else:
                nd_, nc_, ident = "", "", ""
            rows.append(dict(
                specimenCollectorSampleId=sidv, n_accessions=len(g),
                accessionVersion_1=a, accessionVersion_2=b,
                earliestReleaseDate_1=g.earliestReleaseDate.iloc[i], earliestReleaseDate_2=g.earliestReleaseDate.iloc[j],
                sampleCollectionDate_1=g.sampleCollectionDate.iloc[i], sampleCollectionDate_2=g.sampleCollectionDate.iloc[j],
                same_collection_date=g.sampleCollectionDate.iloc[i] == g.sampleCollectionDate.iloc[j],
                zone_1=g.zone.iloc[i], zone_2=g.zone.iloc[j],
                same_zone=("" if (g.flag_zone_not_in_map.iloc[i] or g.flag_zone_not_in_map.iloc[j])
                           else g.zone.iloc[i] == g.zone.iloc[j]),
                sequencingProtocol_1=g.sequencingProtocol.iloc[i], sequencingProtocol_2=g.sequencingProtocol.iloc[j],
                sequencingDate_1=g.sequencingDate.iloc[i], sequencingDate_2=g.sequencingDate.iloc[j],
                called_frac_1=round(g.called_frac.iloc[i], 6), called_frac_2=round(g.called_frac.iloc[j], 6),
                n_sites_both_called=nc_, n_differences=nd_, identical_aligned_sequence=ident,
                groupId=g.groupId.iloc[i], same_group=g.groupId.nunique() == 1,
                dataUseTerms_1=g.dataUseTerms.iloc[i], dataUseTerms_2=g.dataUseTerms.iloc[j]))
    dups = pd.DataFrame(rows, columns=DUP_COLS)
    if len(dups):
        dups = dups.sort_values("specimenCollectorSampleId", kind="stable").reset_index(drop=True)
    return dups


# ---------------------------------------------------------------- comparison of two metadata tables
def norm_val(x):
    x = "" if x is None else str(x)
    if x.lower() in ("nan", "none", "null"):
        return ""
    if x.lower() in ("true", "false"):
        return x.lower()
    try:
        f = float(x)
        return repr(round(f, 9))
    except Exception:
        return x


FLAGCOLS = ["flag_province_conflict", "flag_zone_not_in_reference", "flag_zone_interpreted", "flag_zone_missing",
            "flag_admin1_repeats_country", "flag_collection_date_not_day_precision", "flag_received_before_collection",
            "flag_sequenced_before_collection", "flag_sample_id_missing", "duplicate_sample_id", "flag_lab_passaged",
            "flag_sample_id_on_revoked_accession"]
SHORT = {"flag_province_conflict": "province_conflict", "flag_zone_not_in_reference": "zone_not_in_reference",
         "flag_zone_interpreted": "zone_spelling_interpreted", "flag_zone_missing": "zone_missing",
         "flag_admin1_repeats_country": "admin1_repeats_country",
         "flag_collection_date_not_day_precision": "date_not_day_precision",
         "flag_received_before_collection": "received_before_collection",
         "flag_sequenced_before_collection": "sequenced_before_collection",
         "flag_sample_id_missing": "sample_id_missing", "duplicate_sample_id": "duplicate_sample_id",
         "flag_lab_passaged": "lab_passaged", "flag_sample_id_on_revoked_accession": "sample_id_on_revoked_accession",
         "flag_zone_not_in_map": "zone_spelling_not_in_map"}

VAL_COLS = ["accession", "version", "status_vs_baseline", "baseline_accessionVersion", "baseline_seq_sha256", "seq_sha256",
            "metadata_changed_vs_baseline", "called_frac", "called_frac_ge_0.95", "n_called", "n_N", "n_ambiguous", "n_gap",
            "completeness", "length", "sampleCollectionDate", "collection_date_iso8601", "collection_date_precision",
            "collection_date_exact", "sampleCollectionDateRangeLower", "sampleCollectionDateRangeUpper", "collection_month",
            "earliestReleaseDate", "releasedDate", "submittedDate", "sampleReceivedDate", "sequencingDate",
            "geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "geoLocCity", "geoLocSite", "zone", "zone_rule",
            "zone_ref_province", "sequencingProtocol", "sequencingInstrument", "sequencingAssayType", "purposeOfSequencing",
            "purposeOfSampling", "ampliconPcrPrimerScheme", "ampliconSize", "specimenCollectorSampleId",
            "duplicate_sample_id", "duplicate_partner_accessionVersion", "identical_aligned_seq_partners", "groupId",
            "groupName", "dataUseTerms", "dataUseTermsRestrictedUntil", "dataUseTermsUrl", "isLabHost", "passageNumber",
            "cellLine", "hostNameScientific", "outbreak", "flag_province_conflict", "flag_zone_not_in_reference",
            "flag_zone_interpreted", "flag_zone_missing", "flag_admin1_repeats_country",
            "flag_collection_date_not_day_precision", "flag_received_before_collection",
            "flag_sequenced_before_collection", "flag_sample_id_missing", "flag_lab_passaged",
            "flag_sample_id_on_revoked_accession", "flags"]


def month_label(month, prec):
    if month != "":
        return month
    return {"range": "(no month: date range)", "missing": "(no month: date missing)"}.get(prec, "(no month: date not readable)")


def build(full, seqs, zone_lookup, snapshot_query_utc, data_version, baseline_full=None, baseline_seqs=None,
          add_new_columns=True):
    """full: all records of all versions, database fields only (strings). seqs: dict accessionVersion -> aligned sequence.
    Returns dict with o (outbreak set, derived), M (metadata table), V (validation table), dups, zmap, byzone, zm."""
    raw_cols = list(full.columns)
    ob_full, arm = outbreak_definition(full, "first10")
    ob_iso, arm_iso = outbreak_definition(full, "iso8601")

    o = ob_full.set_index("accessionVersion").copy()
    o["definition_arm_first10"] = arm.values
    o["definition_arm_iso8601"] = (pd.Series(arm_iso.values, index=ob_iso.accessionVersion.values)
                                   .reindex(o.index).fillna("(not in set)"))
    o["called_frac"] = pd.Series({k: called_frac(seqs[k]) for k in o.index if k in seqs}, dtype=float).reindex(o.index)
    rawd = o.sampleCollectionDate
    o["collection_date"] = iso(rawd)
    o["collection_date_precision"] = rawd.map(precision)
    o["release_date"] = iso(o.earliestReleaseDate)
    o["submitted_date"] = iso(o.submittedDate)
    o["collection_month"] = o.collection_date.dt.strftime("%Y-%m").fillna("")
    o["called_frac_ge_0.95"] = o.called_frac >= 0.95
    o["sid_norm"] = o.specimenCollectorSampleId.map(norm_sid)

    hz = harmonise_by_map(o, zone_lookup)
    for c in hz.columns:
        o[c] = hz[c]

    dups = duplicate_pairs(o, seqs)
    dpart = {}
    for _, r in dups.iterrows():
        dpart.setdefault(r.accessionVersion_1, []).append(r.accessionVersion_2)
        dpart.setdefault(r.accessionVersion_2, []).append(r.accessionVersion_1)
    o["duplicate_sample_id"] = o.index.isin(dpart.keys())
    o["duplicate_partner_accessionVersion"] = o.index.map(lambda k: ";".join(sorted(dpart.get(k, []))))

    o["seq_sha256"] = o.index.map(lambda k: sha256_str(seqs[k]) if k in seqs else "")
    # comparison with a baseline
    if baseline_full is not None:
        b_out, _ = outbreak_definition(baseline_full, "first10")
        b_acc = {r.accession: r.accessionVersion for r in b_out.itertuples()}
        st, bav_, bsha_ = {}, {}, {}
        for k, acc in zip(o.index, o.accession):
            if acc not in b_acc:
                st[k], bav_[k], bsha_[k] = "added-since-baseline", "", ""
                continue
            bav = b_acc[acc]
            bsha = sha256_str(baseline_seqs[bav]) if bav in baseline_seqs else ""
            bav_[k], bsha_[k] = bav, bsha
            same_seq = bsha != "" and o.at[k, "seq_sha256"] == bsha
            if bsha == "" or o.at[k, "seq_sha256"] == "":
                st[k] = "sequence-not-available"
            elif k == bav and same_seq:
                st[k] = "identical"
            elif k != bav and same_seq:
                st[k] = "re-versioned"
            else:
                st[k] = "sequence-changed"
        o["status_vs_baseline"] = pd.Series(st)
        o["baseline_accessionVersion"] = pd.Series(bav_)
        o["baseline_seq_sha256"] = pd.Series(bsha_)
        shared = [c for c in baseline_full.columns if c in full.columns and c != "accessionVersion"]
        B = baseline_full.set_index("accessionVersion")[shared]
        N = full.set_index("accessionVersion")[shared]
        # compared with the version that the baseline holds of the SAME ACCESSION (also when the version number changed);
        # empty only for a genome whose accession the baseline set does not hold
        changed, fields_changed = {}, {}
        for k in o.index:
            bav = bav_.get(k, "")
            if bav == "" or bav not in B.index:
                changed[k], fields_changed[k] = "", ""
            else:
                d = [c for c in shared if norm_val(B.at[bav, c]) != norm_val(N.at[k, c])]
                changed[k], fields_changed[k] = bool(d), ";".join(d)
        o["metadata_changed_vs_baseline"] = pd.Series(changed, dtype=object)
        o["metadata_fields_changed_vs_baseline"] = pd.Series(fields_changed, dtype=object)
    else:
        for c in ("status_vs_baseline", "baseline_accessionVersion", "baseline_seq_sha256", "metadata_changed_vs_baseline",
                  "metadata_fields_changed_vs_baseline"):
            o[c] = ""

    o["collection_date_iso8601"] = o.collection_date.dt.strftime("%Y-%m-%d").fillna("")
    o["collection_date_exact"] = np.where(o.collection_date_precision == "day", o.collection_date_iso8601, "")

    recv = iso(o.sampleReceivedDate)
    sqd = iso(o.sequencingDate)
    o["flag_collection_date_not_day_precision"] = o.collection_date_precision != "day"
    o["flag_received_before_collection"] = (recv < o.collection_date).fillna(False)
    o["flag_sequenced_before_collection"] = (sqd < o.collection_date).fillna(False)
    o["flag_sample_id_missing"] = o.specimenCollectorSampleId == ""
    o["flag_lab_passaged"] = (o.isLabHost.str.lower() == "true") | (o.cellLine != "") | (o.passageNumber != "")
    o["flag_admin1_repeats_country"] = (o.geoLocAdmin1 != "") & (o.geoLocAdmin1 == o.geoLocCountry)

    sc = pd.DataFrame({k: seq_counts(seqs[k]) for k in o.index if k in seqs}).T.reindex(o.index)
    for c in ["n_called", "n_N", "n_gap", "n_ambiguous"]:
        o[c] = sc[c] if c in sc.columns else np.nan
    o["n_acgt_1_18900"] = pd.Series({k: count_acgt(seqs[k], L_CALL) for k in o.index if k in seqs}).reindex(o.index)
    o["n_acgt_1_18940"] = pd.Series({k: count_acgt(seqs[k], L_ALN) for k in o.index if k in seqs}).reindex(o.index)

    partners = {}
    with_seq = o[o.seq_sha256 != ""]
    for h, idx in with_seq.groupby("seq_sha256").groups.items():
        if len(idx) > 1:
            for k in idx:
                partners[k] = ";".join(sorted(set(idx) - {k}))
    o["identical_aligned_seq_partners"] = o.index.map(lambda k: partners.get(k, ""))

    nonlatest = full[~full.accessionVersion.isin(o.index)]
    nl = nonlatest[nonlatest.specimenCollectorSampleId != ""]
    hit = o[o.specimenCollectorSampleId.isin(set(nl.specimenCollectorSampleId)) & (o.specimenCollectorSampleId != "")]
    x = hit.reset_index()[["accessionVersion", "accession", "specimenCollectorSampleId"]].merge(
        nl[["accession", "accessionVersion", "versionStatus", "isRevocation", "specimenCollectorSampleId"]],
        on="specimenCollectorSampleId", suffixes=("", "_other"))
    x_other = x[x.accession != x.accession_other]
    o["flag_sample_id_on_revoked_accession"] = o.index.isin(set(x_other.accessionVersion))
    detail = {}
    for r in x_other.itertuples():
        tag = f"{r.accessionVersion_other}[{r.versionStatus}{',revocation entry' if r.isRevocation.lower() == 'true' else ''}]"
        detail.setdefault(r.accessionVersion, []).append(tag)
    o["sample_id_on_other_accession_detail"] = o.index.map(lambda k: ";".join(sorted(detail.get(k, []))))

    flagcols = FLAGCOLS + ["flag_zone_not_in_map"]
    o["flags"] = o[flagcols].apply(lambda r: ";".join(SHORT[c] for c in flagcols if is_true(r[c])), axis=1)

    # ---------- metadata table: all records of all versions
    M = full[raw_cols].copy()
    M.insert(0, "snapshotDatetimeUtc", snapshot_query_utc)
    M.insert(1, "lapisDataVersion", data_version)
    M["outbreak2026"] = M.accessionVersion.isin(o.index)
    M["latest_nonrevoked"] = (M.versionStatus == "LATEST_VERSION") & (M.isRevocation.str.lower() == "false")
    cf_all = pd.Series({k: called_frac(v) for k, v in seqs.items()}, dtype=float)
    M["called_frac"] = M.accessionVersion.map(cf_all)
    cd_all = iso(M.sampleCollectionDate)
    M["collection_date_iso8601"] = cd_all.dt.strftime("%Y-%m-%d").fillna("")
    M["collection_date_precision"] = M.sampleCollectionDate.map(precision)
    for c in ["zone", "zone_rule", "zone_ref_province"]:
        M[c] = M.accessionVersion.map(o[c]).fillna("")
    for c in ["flag_province_conflict", "flag_zone_not_in_reference", "flag_zone_interpreted", "flag_zone_missing",
              "duplicate_sample_id"]:
        M[c] = M.accessionVersion.map(o[c])
    M["flags"] = M.accessionVersion.map(o["flags"]).fillna("")
    if add_new_columns:
        n1 = {k: str(count_acgt(v, L_CALL)) for k, v in seqs.items()}
        n2 = {k: str(count_acgt(v, L_ALN)) for k, v in seqs.items()}
        M["n_acgt_1_18900"] = M.accessionVersion.map(n1).fillna("")
        M["n_acgt_1_18940"] = M.accessionVersion.map(n2).fillna("")

    # ---------- validation table: one row per genome of the outbreak set
    i_gap = VAL_COLS.index("n_gap") + 1
    vcols = VAL_COLS[:i_gap] + (NEW_COLS if add_new_columns else []) + VAL_COLS[i_gap:]
    V = o[vcols].copy()
    if add_new_columns:
        V.insert(len(V.columns) - 1, "metadata_fields_changed_vs_baseline", o["metadata_fields_changed_vs_baseline"])
        V.insert(len(V.columns) - 1, "flag_zone_not_in_map", o["flag_zone_not_in_map"])
        V.insert(len(V.columns) - 1, "rule_on_blanks_acted", o["rule_on_blanks_acted"])
        V.insert(len(V.columns) - 1, "definition_arm", o["definition_arm_first10"])
    for c in ["n_called", "n_N", "n_gap", "n_ambiguous"] + (NEW_COLS if add_new_columns else []):
        V[c] = V[c].map(lambda v: "" if pd.isna(v) else str(int(v)))
    V = V.rename(columns={"completeness": "completeness_pathoplexus", "length": "length_unaligned",
                          "outbreak": "outbreak_field_pathoplexus"})
    V.index.name = "accessionVersion"
    V = V.sort_index()
    V.insert(0, "lapisDataVersion", data_version)
    V.insert(1, "snapshotQueryUtc", snapshot_query_utc)
    assert V.index.is_unique and V.columns.is_unique

    # ---------- zone map and counts
    zcols = ["geoLocCountry", "geoLocAdmin1", "geoLocAdmin2", "zone", "zone_rule", "zone_ref_province",
             "flag_province_conflict", "flag_zone_not_in_reference", "flag_zone_missing", "flag_zone_interpreted"]
    zsrc = o[zcols + ["flag_zone_not_in_map"]].copy()
    for c in ZONE_FLAGS:
        zsrc[c] = zsrc[c].map(lambda v: "" if v is None or (isinstance(v, float) and np.isnan(v)) else str(bool(v)))
    zsrc["flag_zone_not_in_map"] = zsrc["flag_zone_not_in_map"].map(lambda v: str(bool(v)))
    zmap = (zsrc.groupby(zcols + (["flag_zone_not_in_map"] if add_new_columns else [])).size().rename("n_genomes").reset_index()
             .sort_values(["geoLocCountry", "zone", "geoLocAdmin2"]))

    o["zone_label"] = np.where(o.flag_zone_not_in_map, LABEL_NOT_IN_MAP, np.where(o.zone == "", "(not recorded)", o.zone))
    for c in ZONE_FLAGS:
        o["_cnt_" + c] = o[c].map(is_true)
    o["admin1_label"] = np.where(o.geoLocAdmin1 == "", "(not recorded)", o.geoLocAdmin1)
    o["_dayd"] = o.collection_date.where(o.collection_date_precision == "day")
    byzone = (o.groupby(["geoLocCountry", "admin1_label", "zone_label"])
               .agg(n_genomes=("accession", "size"), n_called_frac_ge_095=("called_frac_ge_0.95", "sum"),
                    first_collection_day_precision=("_dayd", "min"), last_collection_day_precision=("_dayd", "max"),
                    n_not_day_precision=("flag_collection_date_not_day_precision", "sum"),
                    n_raw_spellings=("geoLocAdmin2", "nunique"), n_province_conflict=("_cnt_flag_province_conflict", "sum"),
                    n_not_in_reference=("_cnt_flag_zone_not_in_reference", "sum"),
                    n_spelling_interpreted=("_cnt_flag_zone_interpreted", "sum"))
               .reset_index().rename(columns={"geoLocCountry": "country", "admin1_label": "admin1_raw", "zone_label": "zone"}))
    for c in ("first_collection_day_precision", "last_collection_day_precision"):
        byzone[c] = byzone[c].dt.strftime("%Y-%m-%d").fillna("")
    byzone = byzone.sort_values(["country", "admin1_raw", "n_genomes", "zone"], ascending=[True, True, False, True]).reset_index(drop=True)
    assert byzone.n_genomes.sum() == len(o)
    if (byzone.zone == LABEL_NOT_IN_MAP).any():
        for c in ("n_province_conflict", "n_not_in_reference", "n_spelling_interpreted"):
            byzone[c] = byzone[c].astype(object)
            byzone.loc[byzone.zone == LABEL_NOT_IN_MAP, c] = ""          # not assessed

    o["month_label"] = [month_label(m, p) for m, p in zip(o.collection_month, o.collection_date_precision)]
    zm = pd.crosstab([o.geoLocCountry, o.admin1_label, o.zone_label], o.month_label, margins=True, margins_name="total")
    zm.index.names = ["country", "admin1_raw", "zone"]

    return dict(o=o, M=M, V=V, dups=dups, zmap=zmap, byzone=byzone, zm=zm, ob_iso=ob_iso)
