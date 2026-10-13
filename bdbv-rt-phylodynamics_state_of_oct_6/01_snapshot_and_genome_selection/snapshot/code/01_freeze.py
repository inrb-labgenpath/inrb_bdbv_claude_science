"""Step 1 - freeze a supplied download (part pathoplexus/ and pathoplexus_terms/ only).
The script of 2 October (01_freeze.py) with the names of the download as parameters;
nothing else is changed. Without arguments it freezes the download of 7 October 2026.
Reads:  the zip of the authors (path in code/inputs.json), nothing else.
Writes: <workdir>/raw/ (unchanged copies), <workdir>/out/snapshot_fetch_log_<tag>.csv, <workdir>/out/freeze_facts.json,
        <workdir>/out/INTERNAL_pp_snapshot_<tag>_raw_files.tar.gz
Run from the directory snapshot_20261007/ :
    python code/01_freeze.py                                                        (download of 7 October)
    python code/01_freeze.py 20261002 zip_20261002 bdbv_inputs_20261002_1114Z/ control_on_20261002   (control)
"""
import csv
import datetime as dt
import gzip
import io
import json
import os
import pathlib
import shutil
import sys
import tarfile
import zipfile
import zlib
from collections import Counter

import pandas as pd

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent) if "__file__" in globals() else "code")
from snapshot_lib import read_fasta, sha256_file  # noqa: E402

TAG = sys.argv[1] if len(sys.argv) > 1 else "20261007"
ZIPKEY = sys.argv[2] if len(sys.argv) > 2 else "zip_20261007"
ROOT = sys.argv[3] if len(sys.argv) > 3 else "bdbv_inputs_20261007_2056Z/"
HERE = pathlib.Path(sys.argv[4]) if len(sys.argv) > 4 else pathlib.Path(".")
CFG = json.load(open(pathlib.Path("code") / "inputs.json"))
ZIP = CFG["paths"][ZIPKEY]
RAW = HERE / "raw"
OUT = HERE / "out"
OUT.mkdir(parents=True, exist_ok=True)
MY_PARTS = ("pathoplexus/", "pathoplexus_terms/")

facts = {"zip_sha256": sha256_file(ZIP), "zip_bytes": os.path.getsize(ZIP), "zip_version_id": CFG["version_ids"][ZIPKEY]}

# ---- 1. copy my part out of the zip, byte for byte
zf = zipfile.ZipFile(ZIP)
members = [n for n in zf.namelist() if not n.endswith("/")]
facts["zip_members_total"] = len(members)
facts["zip_members_by_part"] = dict(Counter(n[len(ROOT):].split("/")[0] for n in members))
mine = [n for n in members if n[len(ROOT):].startswith(MY_PARTS) or n == ROOT + "fetch_log.csv"]
# integrity inside the zip: CRC-32 and size of the members of pathoplexus/ and pathoplexus_terms/ only, computed while they are copied.
# The other members of the zip are not read and therefore NOT tested here.
crc_bad = []
for n in mine:
    zi = zf.getinfo(n)
    data = zf.read(n)
    if (zlib.crc32(data) & 0xFFFFFFFF) != zi.CRC or len(data) != zi.file_size:
        crc_bad.append(n[len(ROOT):])
    out = RAW / n[len(ROOT):]
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "wb") as dst:
        dst.write(data)
facts["zip_crc32_and_size_checked_members"] = len(mine)
facts["zip_crc32_or_size_mismatch"] = crc_bad
facts["zip_members_not_read_and_not_tested"] = len(members) - len(mine)
facts["members_copied"] = sorted(n[len(ROOT):] for n in mine)
facts["members_not_opened"] = sorted(set(n[len(ROOT):].split("/")[0] for n in members if n not in mine))

# ---- 2. verify against the log
log = pd.read_csv(RAW / "fetch_log.csv", dtype=str, keep_default_na=False)
facts["log_rows_total"] = len(log)
facts["log_rows_by_part"] = log["part"].value_counts().to_dict()
facts["fetch_log_sha256"] = sha256_file(RAW / "fetch_log.csv")
facts["fetch_log_bytes"] = os.path.getsize(RAW / "fetch_log.csv")
my = log[log["part"].isin(["pathoplexus", "terms"])].copy().reset_index(drop=True)
assert len(my) == 11, len(my)
ver = []
for r in my.itertuples():
    p = RAW / r.outfile
    exists = p.is_file()
    ver.append(dict(file_found=exists,
                    verified_n_bytes=os.path.getsize(p) if exists else "",
                    verified_sha256=sha256_file(p) if exists else ""))
ver = pd.DataFrame(ver)
my["file_found"] = ver.file_found
my["verified_n_bytes"] = ver.verified_n_bytes.astype(str)
my["verified_sha256"] = ver.verified_sha256
my["size_matches_log"] = my.verified_n_bytes == my.n_bytes
my["sha256_matches_log"] = my.verified_sha256 == my.sha256
my["url_equals_final_url"] = my.url == my.final_url
in_log = set(my.outfile)
on_disk = set(str(p.relative_to(RAW)) for p in RAW.rglob("*") if p.is_file()) - {"fetch_log.csv"}
facts["files_in_my_part_without_log_row"] = sorted(on_disk - in_log)
facts["log_rows_of_my_part_without_file"] = sorted(in_log - on_disk)

# ---- 3. what the answers say about themselves
body_version, body_lapis, body_silo, body_reqinfo, n_rec = [], [], [], [], []
for r in my.itertuples():
    p = RAW / r.outfile
    bv = bl = bs = bi = ""
    n = ""
    if r.outfile.endswith(".json"):
        j = json.load(open(p))
        inf = j.get("info", j)
        bv, bl, bs, bi = inf.get("dataVersion", ""), inf.get("lapisVersion", ""), inf.get("siloVersion", ""), inf.get("requestInfo", "")
        if "data" in j:
            assert len(j["data"]) == 1 and list(j["data"][0]) == ["count"]
            n = j["data"][0]["count"]
    elif r.outfile.endswith(".csv"):
        with open(p, newline="", encoding="utf-8") as fh:
            rd = csv.reader(fh)
            header = next(rd)
            rows = list(rd)
        assert all(len(x) == len(header) for x in rows), "ragged csv"
        n = len(rows)
        facts[f"n_columns::{r.outfile}"] = len(header)
    elif r.outfile.endswith(".fasta"):
        recs, info = read_fasta(p)
        n = len(recs)
        facts[f"fasta_info::{r.outfile}"] = dict(info, n_records=len(recs), n_unique_headers=len(set(h for h, _ in recs)),
                                                 n_empty_sequences=sum(len(s) == 0 for _, s in recs),
                                                 lengths=dict(sorted(Counter(len(s) for _, s in recs).items())))
    body_version.append(bv); body_lapis.append(bl); body_silo.append(bs); body_reqinfo.append(bi); n_rec.append(str(n))
my["data_version_in_body"] = body_version
my["lapis_version_in_body"] = body_lapis
my["silo_version_in_body"] = body_silo
my["request_info_in_body"] = body_reqinfo
my["n_records_in_answer"] = n_rec
my.to_csv(OUT / f"snapshot_fetch_log_{TAG}.csv", index=False)

lap = my[my.part == "pathoplexus"]
versions_header = sorted(set(v for v in lap.lapis_data_version if v != ""))
versions_body = sorted(set(v for v in lap.data_version_in_body if v != ""))
facts.update(
    n_files_my_part=len(my), all_found=bool(my.file_found.all()), all_sizes_match=bool(my.size_matches_log.all()),
    all_sha256_match=bool(my.sha256_matches_log.all()), all_status_200=bool((my.status == "200").all()),
    all_url_equal_final_url=bool(my.url_equals_final_url.all()),
    lapis_requests=len(lap),
    lapis_first_request_utc=lap.query_start_utc.min(), lapis_last_answer_utc=lap.query_end_utc.max(),
    all_first_request_utc=my.query_start_utc.min(), all_last_answer_utc=my.query_end_utc.max(),
    details_allversions_query_start_utc=lap.loc[lap.outfile == "pathoplexus/details_allversions.csv", "query_start_utc"].iloc[0],
    data_versions_in_log_header=versions_header,
    lapis_answers_without_version_in_log=lap.loc[lap.lapis_data_version == "", "outfile"].tolist(),
    data_versions_in_json_bodies=versions_body,
    lapis_answers_with_version_in_body=lap.loc[lap.data_version_in_body != "", "outfile"].tolist(),
    lapis_versions=sorted(set(v for v in lap.lapis_version_in_body if v)), silo_versions=sorted(set(v for v in lap.silo_version_in_body if v)),
    counts={r.outfile: r.n_records_in_answer for r in lap.itertuples()},
)
dv = versions_body[0] if len(versions_body) == 1 else None
facts["data_version"] = dv
facts["data_version_as_unix_time_utc"] = dt.datetime.fromtimestamp(int(dv), dt.timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC") if dv else None
t0 = dt.datetime.fromisoformat(facts["lapis_first_request_utc"]); t1 = dt.datetime.fromisoformat(facts["lapis_last_answer_utc"])
facts["lapis_window_seconds"] = round((t1 - t0).total_seconds(), 3)
facts["age_of_data_version_at_first_request_hours"] = round((t0 - dt.datetime.fromtimestamp(int(dv), dt.timezone.utc)).total_seconds() / 3600, 2) if dv else None

# ---- 4. archive the raw files unchanged (deterministic archive), then read the archive back
arch = OUT / f"INTERNAL_pp_snapshot_{TAG}_raw_files.tar.gz"
with open(arch, "wb") as fh:
    with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=6) as gz:
        with tarfile.open(fileobj=gz, mode="w", format=tarfile.PAX_FORMAT) as tar:
            for rel in ["fetch_log.csv"] + sorted(in_log):
                ti = tar.gettarinfo(str(RAW / rel), arcname=f"pp_snapshot_{TAG}_raw/" + rel)
                ti.mtime = 0; ti.uid = ti.gid = 0; ti.uname = ti.gname = ""; ti.mode = 0o644
                with open(RAW / rel, "rb") as src:
                    tar.addfile(ti, src)
back = {}
with tarfile.open(arch, "r:gz") as tar:
    for ti in tar.getmembers():
        import hashlib
        back[ti.name.split("/", 1)[1]] = hashlib.sha256(tar.extractfile(ti).read()).hexdigest()
facts["raw_archive"] = arch.name
facts["raw_archive_sha256"] = sha256_file(arch)
facts["raw_archive_members"] = len(back)
facts["raw_archive_readback_matches_log"] = bool(all(back[r.outfile] == r.sha256 for r in my.itertuples())
                                                 and back["fetch_log.csv"] == facts["fetch_log_sha256"])
json.dump(facts, open(OUT / "freeze_facts.json", "w"), indent=1)
print(json.dumps({k: v for k, v in facts.items() if not k.startswith("fasta_info::pathoplexus/unaligned")}, indent=1)[:6000])
