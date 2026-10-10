#!/usr/bin/env bash
# build_sets_round.sh - builds genome sets from ONE snapshot with the frozen pipeline of 24 September:
# the commands of round_20261002/regress_curation.sh (steps 0, 1, 2 and 3 of run_curation_20260924.sh) with the names of the
# input files, the suffix and the set specification as arguments.  The same script is run on the snapshot of 24 September
# (reproduction test, rule 8) and on the snapshot of 2 October.  One process at a time, one thread.
# usage: build_sets_round.sh WORKDIR SUFFIX METADATA ALIGNMENT VALIDATION PAIRS SPEC
#   WORKDIR holds inputs/ (artic-pan-ebola_v2_amplicons.csv, inrb_1046_exclusions.csv, public_quality_screen_20260912.csv)
#   and pipeline_new/ (a copy of pipeline_20260924/ of the working tree); METADATA .. SPEC are paths of files, which are copied to inputs/
set -euo pipefail
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
W="$1"; SFX="$2"; META="$3"; ALN="$4"; VAL="$5"; PAIRS="$6"; SPEC="$7"
cp "$META" "$W/inputs/metadata_$SFX.csv"; cp "$ALN" "$W/inputs/aligned_$SFX.fasta.gz"; cp "$VAL" "$W/inputs/validation_$SFX.csv"
cp "$PAIRS" "$W/inputs/pairs_$SFX.csv"; cp "$SPEC" "$W/inputs/set_spec_$SFX.json"
cd "$W"
A="--metadata ../inputs/metadata_$SFX.csv --alignment ../inputs/aligned_$SFX.fasta.gz"
AMP="--amplicons ../inputs/artic-pan-ebola_v2_amplicons.csv"
PREV="--previous-screen ../inputs/public_quality_screen_20260912.csv"
V1SW="--undated-comparators exclude --threshold-precision exact --batch-rule holm-exact --loo-scan --extended-columns"
mkdir -p run_final run_v2 run_align_v2 logs
# ---- 0. restricted-use genomes of groups other than the main data provider (rule of run_curation_20260924.sh, step 0) ----
python - "inputs/metadata_$SFX.csv" <<'PY'
import sys
import pandas as pd
m = pd.read_csv(sys.argv[1], low_memory=False)
o = m[m.outbreak2026 == True]
main = o.groupName.value_counts().index[0]
d7 = o[(o.dataUseTerms == "RESTRICTED") & (o.groupName != main)]
assert len(d7) == 2 and d7.groupName.nunique() == 2, "expected two groups with one restricted-use genome each"
open("d7_accessions.txt", "w").write(",".join(d7.accessionVersion) + "\n")
PY
D7=$(tr -d '\n' < d7_accessions.txt)
# ---- 1. first-stage screen (rules v1) ----
( cd run_final && python ../pipeline_new/screen.py $A $AMP --outdir . --out-suffix $SFX --expected-n 810 $PREV $V1SW \
    --public-table --validation ../inputs/validation_$SFX.csv --withhold "$D7" --withhold-name minusD7 > ../logs/run_final.log 2>&1 )
# ---- 2. final screen (rules v2) ----
( cd run_v2 && python ../pipeline_new/screen.py $A $AMP --outdir . --out-suffix $SFX --expected-n 810 $PREV --rules v2 \
    --attach-class rules_v1=../run_final/public_quality_screen_$SFX.csv \
    --public-table --validation ../inputs/validation_$SFX.csv --withhold "$D7" --withhold-name minusD7 > ../logs/run_v2.log 2>&1 )
# ---- 3. sets and alignments ----
cp inputs/set_spec_$SFX.json run_align_v2/alignment_set_spec_$SFX.json
( cd run_align_v2 && python ../pipeline_new/build_alignments.py $A --screen ../run_v2/public_quality_screen_$SFX.csv \
    --out-suffix $SFX --outdir . --set-spec alignment_set_spec_$SFX.json --completeness-precision exact \
    --duplicate-pairs ../inputs/pairs_$SFX.csv --data-use-exclude "$D7" \
    --batch-table ../run_v2/screen_batches_$SFX.csv --validation ../inputs/validation_$SFX.csv \
    --exclusion-list ../inputs/inrb_1046_exclusions.csv --tail-days 42 > ../logs/run_align_v2.log 2>&1 )
echo "build_sets_round.sh: finished ($SFX)"
