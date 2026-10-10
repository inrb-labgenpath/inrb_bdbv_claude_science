#!/usr/bin/env python3
"""write_restricted_list.py (round of October 2026) - the restricted-use genomes that are left out.
Until the authors state for which groups they have agreement, the 2 genomes that were left out in the earlier round
are left out again, and no other genome. The list is read from the column restricted_use_agreement_pending of the
table of genome sets of 24 Sep and must equal what the rule of step 0 of run_curation_20260924.sh gives on the new metadata.
usage: python write_restricted_list.py EARLIER_TABLE NEW_METADATA OUT.txt"""
import sys

import pandas as pd

stored, meta, out = sys.argv[1:4]
t = pd.read_csv(stored, dtype=str, keep_default_na=False)
d7 = sorted(t.loc[t.restricted_use_agreement_pending == "True", "accessionVersion"])
m = pd.read_csv(meta, low_memory=False)
o = m[m.outbreak2026 == True]  # noqa: E712
main = o.groupName.value_counts().index[0]
rule = sorted(o[(o.dataUseTerms == "RESTRICTED") & (o.groupName != main)].accessionVersion)
assert len(d7) == 2 and d7 == rule, "the genomes of the stored table are not those that the rule of step 0 gives on the new metadata"
assert set(d7) <= set(o.accessionVersion), "a genome of the list is not a latest version of the new snapshot"
open(out, "w").write(",".join(d7) + "\n")
print(",".join(d7))
