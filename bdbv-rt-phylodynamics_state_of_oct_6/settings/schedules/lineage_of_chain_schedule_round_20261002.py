import pandas as pd
import numpy as np
import re

MD = pd.read_csv('inputs/run_manifest_delphy_20260924.csv', dtype=str, keep_default_na=False, low_memory=False)
MB_ = pd.read_csv('inputs/run_manifest_beast_20260924.csv', dtype=str, keep_default_na=False, low_memory=False)

kD = MD[MD.analysis.isin(["ladder_cells8000", "primary_cells8000_more_seeds"]) & (MD.retained == "yes")]
te_D = sorted(set(kD.tree_every.astype(float)))

kB = {a: MB_[MB_.analysis == a] for a in ("A_skyfix_primary", "B_skyest_primary")}
setB = {a: dict(tree_every=sorted(set(v.tree_every)), threads=sorted(set(v.threads))) for a, v in kB.items()}

rows = []
def add(aid, program, setting, gset, n, first_seed, track, note=""):
    for j in range(n):
        rows.append(dict(analysis=aid, program=program, setting=setting, genome_set=gset, chain=j+1, seed=first_seed+j, track=track, note=note))

for s_, base_ in (("C95", 31000), ("C90", 32000), ("C80", 33000), ("C00", 34000)):
    add(f"delphy_fixed_{s_}", "Delphy 1.4.1", "Skygrid, smoothing fixed (program default)", s_, 24, base_+1, "Delphy sets")
for s_, base_ in (("C95", 41000), ("C90", 42000), ("C80", 43000), ("C00", 44000)):
    add(f"delphy_estimated_{s_}", "Delphy 1.4.1", "Skygrid, smoothing estimated (hyperprior at the program defaults)", s_, 16, base_+1, "Delphy sets")
for d_, base_ in ((1, 51000), (2, 52000)):
    add(f"delphy_fixed_C95masked_draw{d_}", "Delphy 1.4.1", "Skygrid, smoothing fixed (program default)", f"C95 masked, draw {d_} (seed of the draw {20261001+d_})", 12, base_+1, "Masking")
for s_, base_ in (("C90", 62000), ("C80", 63000), ("C00", 64000)):
    add(f"beast_A_{s_}", "BEAST X 10.5.0", "A: Skygrid, precision fixed at 4.06335", s_, 8, base_+1, "BEAST sets")
for s_, base_ in (("C90", 72000), ("C80", 73000), ("C00", 74000)):
    add(f"beast_B_{s_}", "BEAST X 10.5.0", "B: Skygrid, precision estimated (Gamma(0.001, 1000) prior)", s_, 8, base_+1, "BEAST sets")
add("beast_A_C95", "BEAST X 10.5.0", "A: Skygrid, precision fixed at 4.06335", "C95", 8, 61001, "Masking", "unmasked side of the masking experiment (section 4 of the plan)")
add("beast_A_C95masked_draw1", "BEAST X 10.5.0", "A: Skygrid, precision fixed at 4.06335", "C95 masked, draw 1 (seed of the draw 20261002)", 8, 65001, "Masking")

SCH = pd.DataFrame(rows)
isD = SCH.program.str.startswith("Delphy")
SCH["steps"] = np.where(isD, 1500000000, 40000000)
SCH["log_every"] = np.where(isD, 300000, 20000)
SCH["tree_every"] = np.where(isD, 30000000, np.where(SCH.setting.str.startswith("A"), int(setB["A_skyfix_primary"]["tree_every"][0]), int(setB["B_skyest_primary"]["tree_every"][0])))
SCH["threads"] = np.where(isD, 2, int(setB["A_skyfix_primary"]["threads"][0]))
SCH["burn_in_share"] = 0.30

SCH.to_csv("chain_schedule_round_20261002.csv", index=False)