import pandas as pd
import hashlib

# Load the original schedule
orig = pd.read_csv("inputs/phase5_schedule_20261003.csv", dtype=str, keep_default_na=False)
s5 = orig.copy()

h_old = hashlib.sha256(s5.to_csv(index=False).encode()).hexdigest()

# Change 4: Add 8 rows for delphy_spatial_threshold090
base = s5[s5.analysis == "delphy_spatial"].copy()
add = base.copy()
add["analysis"] = "delphy_spatial_threshold090"
add["seed"] = (add.seed.astype(int) + 100000).astype(str)
add["what_is_varied"] = "genomes: the regional set by the rule of the earlier round (threshold 0.90 kept, with its exception)"
add["description"] = "comparison for the reconstruction of spread between areas: 8 chains of 3 x 10^8 steps with 400 cells, as delphy_spatial; added on 3 October 2026 (specification, section 9, change 4)"

i_last = s5.index[s5.analysis == "delphy_spatial"].max()
s5 = pd.concat([s5.iloc[:i_last+1], add, s5.iloc[i_last+1:]], ignore_index=True)

# Fix description: change 3 -> change 4 in delphy_spatial_threshold090 rows
m = s5.analysis == "delphy_spatial_threshold090"
s5.loc[m, "description"] = s5.loc[m, "description"].str.replace("section 9, change 3)", "section 9, change 4)", regex=False)

# Change 5: Add 4 rows for delphy_area_Mangala
lita = s5[s5.analysis == "delphy_area_Lita"].copy()
mg = lita.copy()
mg["analysis"] = "delphy_area_Mangala"
mg["seed"] = [str(110001 + i) for i in range(4)]
mg["chain"] = [str(i + 1) for i in range(4)]
mg["description"] = mg.description.str.replace("Lita", "Mangala", regex=False)
mg["counterpart_in_the_round_of_24_September"] = "none: the stored rule selects the area on the regional set of this round and not on that of the earlier round (specification, section 9, changes 1 and 5)"

i_last = s5.index[s5.analysis == "delphy_area_Bunia_Rwampara"].max()
s5 = pd.concat([s5.iloc[:i_last+1], mg, s5.iloc[i_last+1:]], ignore_index=True)

s5.to_csv("phase5_schedule_20261003.csv", index=False)