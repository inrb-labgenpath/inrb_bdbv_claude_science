#!/usr/bin/env bash
# run_chain.sh - starts ONE Delphy chain of the schedule of the round of October 2026.
#
# usage: run_chain.sh <analysis> <genome_set> <seed> <fixed|estimated> <attempt>
# Run from the run root (the directory that holds programs/, aln/ and runs/).
#
# The command line is the same for every chain; only the alignment, the seed
# and the output names change, and --v0-skygrid-infer-prior-smoothness is added for the analyses with the smoothing
# estimated.  Every value of the command line is fixed here; none is read from a result.
#
# Environment of the process:
#   ulimit -v 8388608          limit of the address space, 8 GiB (set in this shell, which starts the program)
#   KMP_AFFINITY=disabled      as work/pipeline_ext/run_delphy.sh sets it
#   ulimit -c 0                no core file (does not act on the computation)
#
# Files written for the chain (prefix runs/<analysis>/<analysis>_s<seed>):
#   .cmd     command line        .start / .end   seconds since the epoch (UTC)
#   .stdout  standard output and standard error of Delphy (the density criterion reads it)
#   .log     log of Delphy       .trees  sampled trees
#   .exit    exit code of Delphy (written when the program has ended; absent if this script was ended from outside)
#   .env     limits and variables that were in force      .attempt  number of the attempt
#   .lock    held by this script and by Delphy while the chain runs
# No 'set -e': the exit code of Delphy has to be recorded whatever it is.
set -u

ANALYSIS="$1"; SET="$2"; SEED="$3"; SMOOTHING="$4"; ATTEMPT="$5"
STEPS=1500000000; THREADS=2; NUMPAR=20; CUTOFF=0.8; LOG_EVERY=300000; TREE_EVERY=30000000; CELLS=8000
SUFFIX=20261002
DELPHY=programs/delphy_1.4.1/bin/delphy
FASTA="aln/${SET}_${SUFFIX}.fasta"
OUT="runs/${ANALYSIS}/${ANALYSIS}_s${SEED}"

[[ -x "$DELPHY" ]] || { echo "delphy not executable: $DELPHY" >&2; exit 90; }
[[ -s "$FASTA" ]] || { echo "alignment not found: $FASTA" >&2; exit 91; }
[[ "$SMOOTHING" == "fixed" || "$SMOOTHING" == "estimated" ]] || { echo "fixed|estimated" >&2; exit 92; }
mkdir -p "runs/${ANALYSIS}" || exit 93

# one process per chain: the lock is inherited by Delphy and is free again when both have ended
exec 9> "$OUT.lock" || exit 94
flock -n 9 || { echo "chain is running already: $OUT" >&2; exit 99; }
# an attempt never writes over the files of an earlier attempt: they are moved aside before this script starts
for f in "$OUT.log" "$OUT.trees" "$OUT.stdout" "$OUT.exit"; do
  [[ ! -e "$f" ]] || { echo "output exists already: $f" >&2; exit 98; }
done

ulimit -c 0
ulimit -v 8388608 || exit 95
export KMP_AFFINITY=disabled

cmd=("$DELPHY" --v0-in-fasta "$FASTA" --v0-steps "$STEPS" --v0-seed "$SEED" --v0-threads "$THREADS"
     --v0-pop-model skygrid --v0-skygrid-num-parameters "$NUMPAR" --v0-skygrid-cutoff "$CUTOFF"
     --v0-log-every "$LOG_EVERY" --v0-tree-every "$TREE_EVERY"
     --v0-out-log-file "$OUT.log" --v0-out-trees-file "$OUT.trees"
     --v0-target-coal-prior-cells "$CELLS")
[[ "$SMOOTHING" == "fixed" ]] || cmd+=(--v0-skygrid-infer-prior-smoothness)

echo "${cmd[*]}" > "$OUT.cmd"
echo "$ATTEMPT" > "$OUT.attempt"
{ echo "address_space_limit_kib=$(ulimit -v)"; echo "KMP_AFFINITY=${KMP_AFFINITY}"; echo "core_limit=$(ulimit -c)"; } > "$OUT.env"
date -u +%s > "$OUT.start"
"${cmd[@]}" > "$OUT.stdout" 2>&1
rc=$?
date -u +%s > "$OUT.end"
echo "$rc" > "$OUT.exit.tmp" && mv "$OUT.exit.tmp" "$OUT.exit"
exit "$rc"
