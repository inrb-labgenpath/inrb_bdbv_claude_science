"""Step 3 - terms of use: table by submitting group (for the authors), the list of the restricted-use records of other
groups (rule_list below), and the three pages of the terms against those of the download of 2 October 2026 (control).
  * the table by group is written by a function, applied (a) to the metadata of 2 October as a control against the
    stored table of that download and (b) to the metadata of 7 October;
  * the rule list, with one row per listed record in the table by group;
  * the pages are compared with those of the download of 2 October, as text and as raw HTML.
Input:  out/pp_snapshot_20261007_metadata.csv, raw/pathoplexus_terms/*.html, out/snapshot_fetch_log_20261007.csv,
        metadata and table by group of 2 October, control_on_20261002/raw/pathoplexus_terms/*.html and
        control_on_20261002/out/snapshot_fetch_log_20261002.csv (both from the zip of 2 October).
Output: out/INTERNAL_data_use_by_submitting_group_20261007.csv, out/restricted_use_records_of_other_groups_20261007.txt,
        out/INTERNAL_terms_pages_differences_20261002_20261007.txt, out/terms_text_20261007/*.txt,
        out/terms_text_20261002/*.txt, out/terms_facts.json
Run from snapshot_20261007/ :  python code/05_terms_of_use.py
"""
import datetime as dt
import difflib
import json
import pathlib
import re
import sys
from html.parser import HTMLParser

import pandas as pd

sys.path.insert(0, "code")
import snapshot_lib as sl  # noqa: E402

CFG = json.load(open("code/inputs.json"))
P = CFG["paths"]
OUT = pathlib.Path("out")
A, TAG = "20261002", "20261007"
SNAP_DAY = {"20261002": "2026-10-02", "20261007": "2026-10-07"}
FROM_SETS_TABLE = ["n_in_analysed_set_of_round_20260924", "n_restricted_in_analysed_set_of_round_20260924",
                   "n_left_out_in_round_20260924_agreement_pending", "accessions_left_out_in_round_20260924"]


def dates_with_counts(s):
    vc = s[s != ""].value_counts().sort_index()
    return "; ".join(f"{d}: {n}" for d, n in vc.items())


def rule_list(ob):
    """The rule by which the list is made:
       rule = sorted(o[(o.dataUseTerms == "RESTRICTED") & (o.groupName != main)].accessionVersion),
       main = the group with the most records of the outbreak set."""
    counts = ob.groupName.value_counts()
    main = counts.idxmax()
    assert int((counts == counts.max()).sum()) == 1, "two groups share the largest number of records: the rule does not name one main group"
    rule = sorted(ob[(ob.dataUseTerms == "RESTRICTED") & (ob.groupName != main)].accessionVersion)
    return main, rule


def group_table(M, snap_day):
    ob = M[M.outbreak2026 == "True"].copy()
    ob["cf"] = ob.called_frac.astype(float)
    assert set(ob.dataUseTerms) <= {"OPEN", "RESTRICTED"}
    assert ob.groupby("groupId").groupName.nunique().max() == 1 and ob.groupby("groupName").groupId.nunique().max() == 1
    main, rule = rule_list(ob)
    rows = []
    for gid, g in sorted(ob.groupby("groupId"), key=lambda kv: -len(kv[1])):
        r_ = g[g.dataUseTerms == "RESTRICTED"]
        o_ = g[g.dataUseTerms == "OPEN"]
        listed = sorted(set(g.accessionVersion) & set(rule))
        rows.append(dict(
            groupId=gid, groupName=g.groupName.iloc[0], n_genomes_outbreak_set=len(g), n_open=len(o_), n_restricted=len(r_),
            restricted_until_first=r_.dataUseTermsRestrictedUntil.min() if len(r_) else "",
            restricted_until_last=r_.dataUseTermsRestrictedUntil.max() if len(r_) else "",
            restricted_until_dates_with_number_of_genomes=dates_with_counts(r_.dataUseTermsRestrictedUntil),
            n_restricted_without_date=int((r_.dataUseTermsRestrictedUntil == "").sum()),
            n_restricted_with_restriction_in_force_on_snapshot_day=int((r_.dataUseTermsRestrictedUntil > snap_day).sum()),
            earliest_release_of_restricted_first=r_.earliestReleaseDate.min() if len(r_) else "",
            earliest_release_of_restricted_last=r_.earliestReleaseDate.max() if len(r_) else "",
            open_since_first=o_.dataBecameOpenAt.min() if len(o_) else "", open_since_last=o_.dataBecameOpenAt.max() if len(o_) else "",
            n_open_with_insdc_accession=int((o_.insdcAccessionFull != "").sum()),
            countries="; ".join(sorted(set(v for v in g.geoLocCountry if v))),
            n_submitter_accounts=int(g.submitter.nunique()),
            n_restricted_collection_date_to_the_day=int((r_.collection_date_precision == "day").sum()),
            n_restricted_called_fraction_ge_0_95=int((r_.cf >= 0.95).sum()), n_restricted_called_fraction_ge_0_90=int((r_.cf >= 0.90).sum()),
            n_restricted_called_fraction_ge_0_80=int((r_.cf >= 0.80).sum()),
            dataUseTermsUrl="; ".join(sorted(set(g.dataUseTermsUrl))),
            agreement_of_the_group_as_stated_by_the_authors="",
            group_with_the_most_records_of_the_outbreak_set=str(g.groupName.iloc[0] == main),
            n_records_listed_by_the_rule=len(listed), records_listed_by_the_rule=";".join(listed)))
    T = pd.DataFrame(rows)
    tot = {c: "" for c in T.columns}
    tot.update(groupId="(all groups)", groupName="(all groups)")
    for c in T.columns:
        if c.startswith("n_"):
            tot[c] = int(T[c].sum())
    r_all = ob[ob.dataUseTerms == "RESTRICTED"]
    tot.update(restricted_until_first=r_all.dataUseTermsRestrictedUntil.min(), restricted_until_last=r_all.dataUseTermsRestrictedUntil.max(),
               restricted_until_dates_with_number_of_genomes=dates_with_counts(r_all.dataUseTermsRestrictedUntil),
               n_submitter_accounts=int(ob.submitter.nunique()), countries="; ".join(sorted(set(v for v in ob.geoLocCountry if v))),
               records_listed_by_the_rule=";".join(rule))
    T = pd.concat([T, pd.DataFrame([tot])], ignore_index=True)
    assert int(T.n_genomes_outbreak_set.iloc[:-1].sum()) == len(ob) == int(T.n_genomes_outbreak_set.iloc[-1])
    assert int(T.n_records_listed_by_the_rule.iloc[-1]) == len(rule)
    return ob, T, main, rule


F = {}
# ---------------------------------------------------------------- control: the function on the metadata of 2 October
Mo = sl.read_csv_str(P["metadata_20261002"])
ob_o, T_o, main_o, rule_o = group_table(Mo, SNAP_DAY[A])
S = sl.read_csv_str(P["terms_table_20261002"])
shared = [c for c in S.columns if c in T_o.columns]
ctl = dict(rows_stored=len(S), rows_today=len(T_o), columns_stored=len(S.columns), columns_compared=len(shared),
           columns_of_stored_table_not_computed=[c for c in S.columns if c not in T_o.columns],
           columns_new_today=[c for c in T_o.columns if c not in S.columns])
To_s = T_o.astype(str)
ctl["cells_compared"] = int(len(S) * len(shared))
ctl["cells_equal"] = int(sum((S[c].values == To_s[c].values).sum() for c in shared)) if len(S) == len(T_o) else 0
ctl["columns_different"] = [c for c in shared if not (S[c].values == To_s[c].values).all()] if len(S) == len(T_o) else shared
assert ctl["columns_of_stored_table_not_computed"] == FROM_SETS_TABLE, ctl["columns_of_stored_table_not_computed"]
F["control_table_of_20261002"] = ctl
F["rule_on_20261002"] = dict(n_groups=int(ob_o.groupId.nunique()), main_groupId=ob_o.groupId[ob_o.groupName == main_o].iloc[0],
                             n_records_of_main_group=int((ob_o.groupName == main_o).sum()), rule=rule_o,
                             equals_list_named_in_the_brief=bool(rule_o == ["PP_006XXY5.1", "PP_007L9XG.1"]))

# ---------------------------------------------------------------- the table and the rule list of 7 October
M = sl.read_csv_str(OUT / f"pp_snapshot_{TAG}_metadata.csv")
ob, T, main, rule = group_table(M, SNAP_DAY[TAG])
T["row_type"] = ["submitting group"] * (len(T) - 1) + ["all groups"]
rec_cols = ["record_accessionVersion", "record_groupId", "record_geoLocCountry", "record_sampleCollectionDate",
            "record_collection_date_precision", "record_dataUseTerms", "record_dataUseTermsRestrictedUntil"]
for c in rec_cols:
    T[c] = ""
R = ob[ob.accessionVersion.isin(rule)].sort_values("accessionVersion")
rec_rows = []
for r in R.itertuples():
    row = {c: "" for c in T.columns}
    row.update(groupId=r.groupId, groupName=r.groupName, row_type="record listed by the rule (restricted use, other group)",
               record_accessionVersion=r.accessionVersion, record_groupId=r.groupId, record_geoLocCountry=r.geoLocCountry,
               record_sampleCollectionDate=r.sampleCollectionDate, record_collection_date_precision=r.collection_date_precision,
               record_dataUseTerms=r.dataUseTerms, record_dataUseTermsRestrictedUntil=r.dataUseTermsRestrictedUntil)
    rec_rows.append(row)
T_out = pd.concat([T, pd.DataFrame(rec_rows, columns=T.columns)], ignore_index=True) if rec_rows else T
T_out.to_csv(OUT / f"INTERNAL_data_use_by_submitting_group_{TAG}.csv", index=False)
with open(OUT / f"restricted_use_records_of_other_groups_{TAG}.txt", "w", encoding="utf-8", newline="\n") as fh:
    fh.write("".join(av + "\n" for av in rule))
back = [x for x in open(OUT / f"restricted_use_records_of_other_groups_{TAG}.txt", encoding="utf-8").read().split("\n") if x]
assert back == rule and list(R.accessionVersion) == rule
assert all(re.fullmatch(r"PP_[0-9A-Z]+\.[0-9]+", x) for x in back)

F.update(n_outbreak=len(ob), n_groups=int(ob.groupId.nunique()), terms_counts=ob.dataUseTerms.value_counts().to_dict(),
         main_groupId=ob.groupId[ob.groupName == main].iloc[0], n_records_of_main_group=int((ob.groupName == main).sum()),
         n_records_of_second_largest_group=int(sorted(ob.groupName.value_counts().tolist())[-2]) if ob.groupName.nunique() > 1 else 0,
         rule=rule, rule_n=len(rule), rule_same_as_20261002=bool(rule == rule_o),
         rule_by_groupId_gives_same_list=bool(rule == sorted(ob[(ob.dataUseTerms == "RESTRICTED") & (ob.groupId != ob.groupId[ob.groupName == main].iloc[0])].accessionVersion)),
         rule_records=[dict(accessionVersion=r.accessionVersion, groupId=r.groupId, geoLocCountry=r.geoLocCountry, geoLocAdmin1=r.geoLocAdmin1,
                            sampleCollectionDate=r.sampleCollectionDate, collection_date_precision=r.collection_date_precision,
                            dataUseTerms=r.dataUseTerms, dataUseTermsRestrictedUntil=r.dataUseTermsRestrictedUntil,
                            earliestReleaseDate=r.earliestReleaseDate, called_frac=r.called_frac) for r in R.itertuples()],
         restricted_records_of_main_group=int(((ob.dataUseTerms == "RESTRICTED") & (ob.groupName == main)).sum()),
         open_records_of_main_group=int(((ob.dataUseTerms == "OPEN") & (ob.groupName == main)).sum()),
         open_records_of_other_groups=int(((ob.dataUseTerms == "OPEN") & (ob.groupName != main)).sum()),
         table_rows=len(T_out), table_columns=len(T_out.columns), table_rows_by_type=T_out.row_type.value_counts().to_dict())
F["by_group_without_names"] = T.drop(columns=["groupName"] + rec_cols).to_dict("records")
F["by_group_20261002_without_names"] = T_o.drop(columns=["groupName"]).to_dict("records")
r_all = ob[ob.dataUseTerms == "RESTRICTED"]
d = (pd.to_datetime(r_all.dataUseTermsRestrictedUntil) - pd.to_datetime(r_all.earliestReleaseDate)).dt.days
F["restricted_until_minus_earliest_release_days"] = {int(k): int(v) for k, v in d.value_counts().items()}
d2 = (pd.to_datetime(r_all.dataUseTermsRestrictedUntil) - pd.to_datetime(r_all.submittedDate)).dt.days
F["restricted_until_minus_submitted_date_of_current_version_days"] = {int(k): int(v) for k, v in d2.value_counts().items()}
F["restricted_until_not_365_days_after_submitted_date"] = sorted(r_all.accessionVersion[d2.values != 365])
F["restricted_without_date"] = sorted(r_all.accessionVersion[r_all.dataUseTermsRestrictedUntil == ""])
F["restricted_with_date_not_after_snapshot_day"] = sorted(r_all.accessionVersion[(r_all.dataUseTermsRestrictedUntil != "") & (r_all.dataUseTermsRestrictedUntil <= SNAP_DAY[TAG])])
F["open_with_restricted_until_date"] = sorted(ob.accessionVersion[(ob.dataUseTerms == "OPEN") & (ob.dataUseTermsRestrictedUntil != "")])
F["open_with_insdc_accession"] = int(((ob.dataUseTerms == "OPEN") & (ob.insdcAccessionFull != "")).sum())
F["restricted_with_insdc_accession"] = int(((ob.dataUseTerms == "RESTRICTED") & (ob.insdcAccessionFull != "")).sum())
F["restricted_until_distinct_dates"] = int(r_all.dataUseTermsRestrictedUntil.nunique())
lnr_other = M[(M.latest_nonrevoked == "True") & (M.outbreak2026 != "True")]
F["other_latest_nonrevoked_terms"] = lnr_other.dataUseTerms.value_counts().to_dict()
# terms of the records held on 2 October: any change of terms, date or group?
both = sorted(set(ob.accessionVersion) & set(ob_o.accessionVersion))
a_ = ob_o.set_index("accessionVersion").loc[both]
b_ = ob.set_index("accessionVersion").loc[both]
F["terms_fields_of_genomes_held_on_20261002"] = {c: int((a_[c].values != b_[c].values).sum()) for c in
                                                 ["dataUseTerms", "dataUseTermsRestrictedUntil", "dataUseTermsUrl", "dataBecameOpenAt", "groupId", "groupName", "submitter"]}
F["n_genomes_held_on_20261002"] = len(both)
new_ob = ob[~ob.accession.isin(set(Mo.accession))]
F["new_genomes_by_group_and_terms"] = {f"groupId={k[0]}|{k[1]}|restricted until {k[2]}": int(v) for k, v in
                                       new_ob.groupby(["groupId", "dataUseTerms", "dataUseTermsRestrictedUntil"]).size().items()}
F["new_genomes_n"] = len(new_ob)
F["new_genomes_of_main_group"] = int((new_ob.groupName == main).sum())


# ---------------------------------------------------------------- text of the supplied pages
class TextExtractor(HTMLParser):
    BLOCK = {"p", "div", "li", "ul", "ol", "h1", "h2", "h3", "h4", "h5", "h6", "br", "tr", "section", "article", "header", "footer",
             "nav", "main", "table", "blockquote", "hr", "title"}
    SKIP = {"script", "style", "noscript", "svg", "template"}

    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.out, self.skip = [], 0

    def handle_starttag(self, tag, attrs):
        if tag in self.SKIP:
            self.skip += 1
        if tag in self.BLOCK:
            self.out.append("\n")
        if tag in ("h1", "h2", "h3", "h4", "h5", "h6"):
            self.out.append("#" * int(tag[1]) + " ")
        if tag == "li":
            self.out.append("- ")

    def handle_endtag(self, tag):
        if tag in self.SKIP and self.skip > 0:
            self.skip -= 1
        if tag in self.BLOCK:
            self.out.append("\n")

    def handle_data(self, data):
        if self.skip == 0:
            self.out.append(data)


def page_lines(path):
    raw = open(path, "rb").read().decode("utf-8", errors="replace")
    p = TextExtractor()
    p.feed(raw)
    t = re.sub(r"[ \t\r\f\v]+", " ", "".join(p.out))
    return raw, [x.strip() for x in t.split("\n") if x.strip()]


PAGES = ["data-use-terms.html", "restricted-data.html", "open-data.html"]
DIRS = {A: pathlib.Path("control_on_20261002/raw/pathoplexus_terms"), TAG: pathlib.Path("raw/pathoplexus_terms")}
raws, texts = {}, {}
for tag in (A, TAG):
    TT = OUT / f"terms_text_{tag}"
    TT.mkdir(exist_ok=True)
    for f in PAGES:
        raw, lines = page_lines(DIRS[tag] / f)
        raws[(tag, f)], texts[(tag, f)] = raw, lines
        (TT / f.replace(".html", ".txt")).write_text("\n".join(f"{i:03d}  {x}" for i, x in enumerate(lines)) + "\n", encoding="utf-8")
pages = {f: texts[(TAG, f)] for f in PAGES}


def find(page, start, n=1):
    """the line(s) of the page text of 7 October that begin with `start` (verbatim quote, with its line number)"""
    hits = [(i, x) for i, x in enumerate(pages[page]) if x.startswith(start)]
    assert len(hits) == n, (page, start, len(hits))
    return [dict(page=page, line=i, text=x) for i, x in hits]


Q = {}
Q["version_of_terms"] = find("data-use-terms.html", "These Terms were updated on")
Q["duration_full"] = find("data-use-terms.html", "Data can remain")
Q["duration_summary"] = find("restricted-data.html", "Restricted Data on Pathoplexus is subject to")
Q["publication_permission_full"] = find("data-use-terms.html", "In scientific publications and preprints, Restricted-Use data")
Q["publication_permission_summary"] = find("restricted-data.html", "- Use of restricted data often requires explicit permission")
Q["seqset_full"] = find("data-use-terms.html", "Users utilizing data from Pathoplexus in publications and preprints must create a SeqSet")
Q["focal_definition"] = find("data-use-terms.html", "- The “Focal Set” is defined as")
Q["background_definition"] = find("data-use-terms.html", "- The “Background Set” is defined as")
i0 = find("data-use-terms.html", "#### Requirements for use of Restricted-Use Data in publications and preprints:")[0]["line"]
i1 = find("data-use-terms.html", "#### 4.2.3 Deciding if Restricted-Use Data is part of a Focal or Background Set")[0]["line"]
Q["requirements_block"] = [dict(page="data-use-terms.html", line=i, text=pages["data-use-terms.html"][i]) for i in range(i0, i1)]
Q["if_unsure"] = find("data-use-terms.html", "If unsure whether data should be in the Focal or Background Set")
Q["focal_geographic_area"] = find("data-use-terms.html", "- Are these sequences from a specific geographic area")
Q["focal_accessions_mentioned"] = find("data-use-terms.html", "- Do I mention these accession IDs")
h0 = find("data-use-terms.html", "Hint: In outbreak settings")[0]["line"]
h1 = find("data-use-terms.html", "In some cases there are a large number of data submitters")[0]["line"]
Q["outbreak_hint"] = [dict(page="data-use-terms.html", line=h0, text=" ".join(pages["data-use-terms.html"][h0:h1]))]
Q["few_groups"] = find("data-use-terms.html", "In some cases there are a large number of data submitters")
Q["third_party"] = find("data-use-terms.html", "Restricted-Use Data from Pathoplexus can be shared onward")
Q["unpublished"] = find("data-use-terms.html", "- For Restricted-Use Data, it must be acknowledged")
Q["summary_note"] = find("restricted-data.html", "This summary is provided as a general guideline")
s0 = find("restricted-data.html", "- Usage requires")[0]["line"]
Q["summary_focal"] = [dict(page="restricted-data.html", line=i, text=pages["restricted-data.html"][i]) for i in range(s0, s0 + 3)]
assert Q["summary_focal"][2]["text"].startswith("- Obtaining explicit written permission")
F["quotes"] = Q

# ---------------------------------------------------------------- the pages against those of the download of 2 October
olog = sl.read_csv_str("control_on_20261002/out/snapshot_fetch_log_20261002.csv")
nlog = sl.read_csv_str(OUT / f"snapshot_fetch_log_{TAG}.csv")


def tokens(raw):
    """raw HTML cut into tags and text between tags; long tags are cut further into pieces between quotation marks and commas.
    Returns the pieces and, for each piece, whether it is text between tags ('text') or part of a tag ('tag')."""
    out, kinds = [], []
    for t in re.findall(r"<[^>]*>|[^<]+", raw):
        kind = "tag" if t.startswith("<") else "text"
        if len(t) > 200:
            for x in re.split(r"(&quot;|,|\"|\s+)", t):
                if x != "":
                    out.append(x)
                    kinds.append(kind)
        else:
            out.append(t)
            kinds.append(kind)
    return out, kinds


def classify(old, new, ctx, kinds):
    """class of a passage of the raw HTML that differs; a passage that fits no class is 'OTHER' and is listed"""
    both_ = old + " " + new
    if kinds and all(k == "text" for k in kinds):
        return "text of the page"
    if re.fullmatch(r"\d{13}\]\}", old) and re.fullmatch(r"\d{13}\]\}", new) and "&quot;:[0," in ctx:
        return "time of the server when the page was made"
    if re.search(r"Astro v\d", both_):
        return "version of the program that makes the pages"
    if re.search(r"/_astro/[A-Za-z]+\.[A-Za-z0-9_-]+\.(css|js)", both_):
        return "name of a style or script file"
    if ctx.endswith("uid=\"") or re.search(r"uid=", both_):
        return "identifier of a page component"
    if re.search(r"hikungunya", both_):
        return "menu of organisms (Chikungunya Virus added)"
    if re.search(r"class=\"group", both_):
        return "style class of the menu of organisms"
    return "OTHER"


cmp_pages, passages, raw_rows = [], [], []
for f in PAGES:
    o = olog[olog.outfile.str.endswith("/" + f)]
    n = nlog[nlog.outfile.str.endswith("/" + f)]
    assert len(o) == 1 and len(n) == 1
    la, lb = texts[(A, f)], texts[(TAG, f)]
    sm = difflib.SequenceMatcher(None, la, lb, autojunk=False)
    n_pass = 0
    for op, i1, i2, j1, j2 in sm.get_opcodes():
        if op == "equal":
            continue
        n_pass += 1
        old_t, new_t = la[i1:i2], lb[j1:j2]
        menu = all(x.startswith("Andes Virus") for x in old_t + new_t)
        passages.append(dict(page=f, lines_20261002=f"{i1}-{i2 - 1}" if i2 > i1 else "", lines_20261007=f"{j1}-{j2 - 1}" if j2 > j1 else "",
                             text_20261002=old_t, text_20261007=new_t,
                             kind="menu of organisms at the head of the page (not part of the terms)" if menu else "wording of the terms",
                             context_before=la[i1 - 1] if i1 > 0 else "", context_after=la[i2] if i2 < len(la) else ""))
    (ta, ka), (tb, kb) = tokens(raws[(A, f)]), tokens(raws[(TAG, f)])
    assert "".join(ta) == raws[(A, f)] and "".join(tb) == raws[(TAG, f)]
    sm2 = difflib.SequenceMatcher(None, ta, tb, autojunk=False)
    n_raw = 0
    for op, i1, i2, j1, j2 in sm2.get_opcodes():
        if op == "equal":
            continue
        n_raw += 1
        old_r, new_r = "".join(ta[i1:i2]), "".join(tb[j1:j2])
        ctx = "".join(ta[max(0, i1 - 3):i1])[-80:]
        kind = classify(old_r, new_r, ctx, ka[i1:i2] + kb[j1:j2])
        raw_rows.append(dict(page=f, kind=kind, context_before=ctx, raw_20261002=old_r, raw_20261007=new_r))
    st = re.search(r"serverTime&quot;:\[0,(\d{13})\]", raws[(TAG, f)])
    st_o = re.search(r"serverTime&quot;:\[0,(\d{13})\]", raws[(A, f)])

    def as_time(m):
        return dt.datetime.fromtimestamp(int(m.group(1)) / 1000, dt.timezone.utc).isoformat(timespec="milliseconds") if m else ""

    cmp_pages.append(dict(page=f, n_bytes_20261002=o.n_bytes.iloc[0], n_bytes_20261007=n.n_bytes.iloc[0],
                          difference_bytes=int(n.n_bytes.iloc[0]) - int(o.n_bytes.iloc[0]),
                          sha256_20261002=o.sha256.iloc[0], sha256_20261007=n.sha256.iloc[0],
                          sha256_equal=bool(o.sha256.iloc[0] == n.sha256.iloc[0]),
                          url_20261002=o.url.iloc[0], url_20261007=n.url.iloc[0], url_equal=bool(o.url.iloc[0] == n.url.iloc[0]),
                          last_modified_20261002=o.last_modified.iloc[0], last_modified_20261007=n.last_modified.iloc[0],
                          text_lines_20261002=len(la), text_lines_20261007=len(lb), text_equal=bool(la == lb),
                          n_text_passages_different=n_pass, n_raw_passages_different=n_raw,
                          server_time_in_page_20261002=as_time(st_o), server_time_in_page_20261007=as_time(st),
                          request_start_20261007=n.query_start_utc.iloc[0], request_end_20261007=n.query_end_utc.iloc[0],
                          server_time_within_request_20261007=bool(st and n.query_start_utc.iloc[0][:23] <= as_time(st)[:23] <= n.query_end_utc.iloc[0][:23]),
                          request_start_20261002=o.query_start_utc.iloc[0], request_end_20261002=o.query_end_utc.iloc[0],
                          server_time_within_request_20261002=bool(st_o and o.query_start_utc.iloc[0][:23] <= as_time(st_o)[:23] <= o.query_end_utc.iloc[0][:23])))
F["pages_vs_download_of_20261002"] = cmp_pages
F["text_passages_different"] = passages
F["n_text_passages_wording_of_the_terms"] = int(sum(p["kind"] == "wording of the terms" for p in passages))
RR = pd.DataFrame(raw_rows)
F["raw_passages_by_page_and_kind"] = {f"{k[0]}|{k[1]}": int(v) for k, v in RR.groupby(["page", "kind"]).size().items()}
F["raw_passages_of_kind_OTHER"] = RR[RR.kind == "OTHER"].to_dict("records")
F["n_raw_passages_different"] = len(RR)
with open(OUT / f"INTERNAL_terms_pages_differences_{A}_{TAG}.txt", "w", encoding="utf-8", newline="\n") as fh:
    fh.write(f"Pages of the terms of use: download of 7 October 2026 against download of 2 October 2026 (INTERNAL check)\n"
             f"Part 1: every passage of the TEXT of the pages that differs. Part 2: every passage of the raw HTML that differs.\n\n")
    fh.write("PART 1 - TEXT (page text without tags, scripts and styles; line numbers as in terms_text_<date>/)\n")
    for i, p in enumerate(passages, 1):
        fh.write(f"\n[{i}] page {p['page']} | {p['kind']} | lines of 2 October: {p['lines_20261002']} | lines of 7 October: {p['lines_20261007']}\n")
        for x in p["text_20261002"]:
            fh.write(f"    2 October: {x}\n")
        for x in p["text_20261007"]:
            fh.write(f"    7 October: {x}\n")
    fh.write("\nPART 2 - RAW HTML (cut into tags and, inside long tags, into pieces between quotation marks and commas)\n")
    for i, r in enumerate(raw_rows, 1):
        fh.write(f"\n[{i}] page {r['page']} | {r['kind']}\n    after:     {r['context_before']!r}\n    2 October: {r['raw_20261002']!r}\n    7 October: {r['raw_20261007']!r}\n")
F["terms_urls_in_metadata"] = ob.dataUseTermsUrl.value_counts().to_dict()
F["sha256"] = {p.name: sl.sha256_file(p) for p in [OUT / f"INTERNAL_data_use_by_submitting_group_{TAG}.csv",
                                                   OUT / f"restricted_use_records_of_other_groups_{TAG}.txt",
                                                   OUT / f"INTERNAL_terms_pages_differences_{A}_{TAG}.txt"]}
json.dump(F, open(OUT / "terms_facts.json", "w"), indent=1, ensure_ascii=False)
print(T.drop(columns=["groupName", "restricted_until_dates_with_number_of_genomes", "dataUseTermsUrl"] + rec_cols).T.to_string())
print(json.dumps({k: F[k] for k in ["control_table_of_20261002", "rule_on_20261002", "rule", "rule_records", "rule_by_groupId_gives_same_list",
                                    "restricted_until_minus_earliest_release_days", "restricted_until_not_365_days_after_submitted_date",
                                    "terms_fields_of_genomes_held_on_20261002", "new_genomes_by_group_and_terms", "other_latest_nonrevoked_terms",
                                    "raw_passages_by_page_and_kind", "n_text_passages_wording_of_the_terms", "open_with_restricted_until_date",
                                    "restricted_without_date", "restricted_with_date_not_after_snapshot_day"]}, indent=1)[:6000])
print(json.dumps(cmp_pages, indent=1)[:3000])
print("OTHER raw passages:", json.dumps(F["raw_passages_of_kind_OTHER"], indent=1)[:3000])
