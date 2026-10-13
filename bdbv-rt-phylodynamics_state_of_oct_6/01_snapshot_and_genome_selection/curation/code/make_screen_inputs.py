#!/usr/bin/env python3
"""make_screen_inputs.py (snapshot of 7 October 2026) - writes the table of metadata and the alignment that the
program of the quality checks reads: those of the snapshot WITHOUT the restricted-use records of other groups.

Rule: the records of the list are removed from the
table of metadata and from the alignment; nothing else is changed. The list is checked against the rule as coded on 2 October
(records of the outbreak set whose terms are RESTRICTED and whose group is not the group with the most records).

The metadata is filtered line by line (one record per physical line; checked), so that every byte of the records that stay is
that of the snapshot. The alignment keeps the order and the text of the records that stay; it is written as gzip with the
modification time 0, no file name in the header, and compression level 9.

usage: make_screen_inputs.py METADATA ALIGNMENT.fasta.gz LIST.txt SUFFIX OUTDIR
writes OUTDIR/INTERNAL_screen_input_metadata_<SUFFIX>.csv, OUTDIR/screen_input_outbreak_aligned_<SUFFIX>.fasta.gz and
       OUTDIR/screen_inputs_<SUFFIX>.csv (one row per file: SHA-256, records, records removed)
"""
import csv
import gzip
import hashlib
import io
import os
import sys

import pandas as pd

meta_path, aln_path, list_path, sfx, outdir = sys.argv[1:6]
os.makedirs(outdir, exist_ok=True)


def sha(b):
    return hashlib.sha256(b).hexdigest()


remove = [l.strip() for l in open(list_path) if l.strip()]
assert len(remove) == len(set(remove))
M = pd.read_csv(meta_path, low_memory=False, dtype=str, keep_default_na=False)
o = M[M.outbreak2026 == "True"]
main = o.groupName.value_counts().index[0]
rule = sorted(o[(o.dataUseTerms == "RESTRICTED") & (o.groupName != main)].accessionVersion)
assert rule == sorted(remove), ("the list is not what the rule as coded gives", rule, remove)
nov = lambda a: a.rsplit(".", 1)[0]
assert not M[M.accessionVersion.map(nov).isin([nov(a) for a in remove]) & ~M.accessionVersion.isin(remove)].shape[0], "another version of a removed accession is in the table"

# ---- metadata, line by line ----
raw = open(meta_path, "rb").read()
lines = raw.split(b"\n")
assert lines[-1] == b"" and len(lines) - 2 == len(M), "a record of the metadata is not one physical line"
col = list(M.columns).index("accessionVersion")
keep, gone = [lines[0]], []
for ln, av in zip(lines[1:-1], M.accessionVersion):
    first = next(csv.reader(io.StringIO(ln.decode("utf-8"))))
    assert first[col] == av
    (gone if av in remove else keep).append(ln if av not in remove else av)
assert sorted(gone) == sorted(remove)
out_meta = b"\n".join(keep) + b"\n"
f_meta = os.path.join(outdir, f"INTERNAL_screen_input_metadata_{sfx}.csv")
open(f_meta, "wb").write(out_meta)
M2 = pd.read_csv(f_meta, low_memory=False, dtype=str, keep_default_na=False)
assert len(M2) == len(M) - len(remove) and M2.equals(M[~M.accessionVersion.isin(remove)].reset_index(drop=True))

# ---- alignment ----
txt = gzip.open(aln_path, "rb").read()
recs = txt.split(b">")
assert recs[0] == b""
kept, n_all, removed_aln = [], 0, []
for r in recs[1:]:
    name = r.split(b"\n", 1)[0].decode().strip()
    n_all += 1
    if name in remove:
        removed_aln.append(name)
    else:
        kept.append(b">" + r)
assert sorted(removed_aln) == sorted(remove)
out_txt = b"".join(kept)
f_aln = os.path.join(outdir, f"screen_input_outbreak_aligned_{sfx}.fasta.gz")
with open(f_aln, "wb") as fh:
    with gzip.GzipFile(filename="", mode="wb", fileobj=fh, mtime=0, compresslevel=9) as gz:
        gz.write(out_txt)
assert gzip.open(f_aln, "rb").read() == out_txt
n_out = int((M2.outbreak2026 == "True").sum())
rows = [dict(file=os.path.basename(f_meta), sha256=sha(open(f_meta, "rb").read()), bytes=os.path.getsize(f_meta), records=len(M2), records_of_the_outbreak_set=n_out,
             sha256_of_the_uncompressed_text="", records_in_the_file_of_the_snapshot=len(M), records_removed=len(remove), removed="; ".join(sorted(remove)),
             file_of_the_snapshot=os.path.basename(meta_path), sha256_of_the_file_of_the_snapshot=sha(raw)),
        dict(file=os.path.basename(f_aln), sha256=sha(open(f_aln, "rb").read()), bytes=os.path.getsize(f_aln), records=len(kept), records_of_the_outbreak_set=len(kept),
             sha256_of_the_uncompressed_text=sha(out_txt), records_in_the_file_of_the_snapshot=n_all, records_removed=len(remove), removed="; ".join(sorted(remove)),
             file_of_the_snapshot=os.path.basename(aln_path), sha256_of_the_file_of_the_snapshot=sha(open(aln_path, "rb").read()))]
pd.DataFrame(rows).to_csv(os.path.join(outdir, f"screen_inputs_{sfx}.csv"), index=False)
print(pd.DataFrame(rows)[["file", "records", "records_of_the_outbreak_set", "records_removed", "removed"]].to_string())
