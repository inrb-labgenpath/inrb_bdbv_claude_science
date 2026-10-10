import io
import os
import re
import tarfile

import pandas as pd

# Load genome sets by called fraction
GSv = 'inputs/genome_sets_by_called_fraction_20261002.csv'
GS_ = pd.read_csv(GSv, dtype=str, keep_default_na=False)
acc_col = "accession_version"
assert len(GS_) == 810 and GS_[acc_col].is_unique

MEMB = {}
SRC = {}

for s_ in ("C95", "C90", "C80", "C00"):
    MEMB[s_] = set(GS_[acc_col][GS_[f"in_{s_}"] == "True"])
    SRC[s_] = ("Curation", f"genome_sets_by_called_fraction_20261002.csv")

def heads(fh):
    return [l[1:].strip() for l in io.TextIOWrapper(fh, encoding="utf-8", errors="replace") if l.startswith(">")]

for trk, vid, pref in (
    ("Sensitivity", 'inputs/alignments_sensitivity_20261002.tar.gz', "sens"),
    ("Spatial", 'inputs/alignments_spatial_20261002.tar.gz', "spat"),
):
    with tarfile.open(vid) as t:
        for m in t.getmembers():
            if not (m.isfile() and re.search(r"\.(fasta|fa)$", m.name)):
                continue
            h = heads(t.extractfile(m))
            ids = [x.split("|")[0] for x in h]
            nm_ = f"{pref}:{os.path.basename(m.name).replace('_20261002.fasta', '').replace('.fasta', '')}"
            assert len(ids) == len(set(ids)), (m.name, "a genome twice")
            MEMB[nm_] = set(ids)
            SRC[nm_] = (trk, f"{m.name}")

allacc = set(GS_[acc_col])

MT = pd.DataFrame({"accession_version": sorted(allacc)})
for k in MEMB:
    MT[k] = MT.accession_version.isin(MEMB[k])

setcols = [c for c in MT.columns if c != "accession_version"]
MT["in_any_analysed_set"] = MT[setcols].any(axis=1)
MT["in_any_set_beyond_C00"] = MT[[c for c in setcols if c != "C00"]].any(axis=1) & ~MT["C00"]

MT.to_csv("membership_of_every_genome_in_every_analysed_set_20261003.csv", index=False)