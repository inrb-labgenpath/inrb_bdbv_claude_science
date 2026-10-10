#!/usr/bin/env bash
# run_chain_v2.sh - runs ONE part of ONE chain of BEAST X (version 2 of 2 October 2026). Started by a scheduler that is not part
# of this repository. THIS FILE IS NOT EDITED ONCE IT IS IN USE (bash reads a running script from its file); a change gets a new
# file name.
# usage: run_chain_v2.sh <root> <label> <analysis> <genome_set> <setting> <seed> <threads> <xml-stem> <mode> <n>
#   mode fresh     start from state 0 in runs/<analysis>/<label>/ (n = number of the start, 1 = first):
#                  beast -threads T -beagle_CPU -beagle_SSE -beagle_threads T -seed S -save_every N -save_state <label>.state -overwrite <stem>.xml
#   mode continue  continuation number n in a directory of its own, runs/<analysis>/<label>/continuation_<n>/ :
#                  beast ... -seed S -load_state loaded_state_<N>.state -force_resume -save_every N -save_state <label>.state -overwrite <stem>.xml
#     - loaded is a COPY of a saved state of the part before; no file of the part before is changed;
#     - the XML of the continuation is the XML of the chain with ONLY the chain length changed to the states that remain
#       (tools/make_continuation_xml.py; BEAST X counts the chain length from the loaded state), so that the chain ends at its last state;
#     - saved states are tried in this order: the state file of the part before, then the copies of earlier saved states kept for that part,
#       newest first. A saved state is not used if it is incomplete, if BEAST X ends before it has logged a state, or if the difference
#       between the saved and the recomputed log-likelihood that BEAST X reports is not below 1e-6. The files of an attempt that is not
#       used are kept in rejected_continuation_<n>_attempt_<a>/ . If no saved state can be used, all_saved_states_rejected.txt is written
#       and the scheduler starts the chain again from state 0.
set -u
ROOT=$1; LABEL=$2; ANALYSIS=$3; GSET=$4; SETTING=$5; SEED=$6; THREADS=$7; STEM=$8; MODE=${9:-fresh}; N=${10:-1}
STEPS=${STEPS:-40000000}; LIMIT=${LIKELIHOOD_LIMIT:-1e-6}
. "$ROOT/sched/config.sh"
CHAINDIR="$ROOT/runs/$ANALYSIS/$LABEL"
mkdir -p "$CHAINDIR"
exec 9>>"$CHAINDIR/chain.lock"
flock -n 9 || exit 97
export PATH="$BEAST_BIN:$PATH"
export JAVA_TOOL_OPTIONS="$JAVA_TOOL_OPTIONS_VALUE"
unset OMP_NUM_THREADS
now() { date -u +%Y-%m-%dT%H:%M:%SZ; }
CMD=""
ledger() {  # event exit_code last_state note
  local line
  line="$(now),$1,$LABEL,$ANALYSIS,$GSET,$SETTING,$SEED,$$,$2,$3,\"$CMD\",\"$4\""
  ( flock 6; echo "$line" >> "$ROOT/sched/ledger.csv" ) 6>>"$ROOT/sched/ledger.lock" 9>&-
}
last_state() { local f="$1"; [ -s "$f" ] && tail -1 "$f" | cut -f1 | grep -E '^[0-9]+$' || echo ""; }
state_number() { awk -F'\t' 'NR<=3 && $1=="state"{print $2; exit}' "$1" 2>/dev/null; }
complete_state() {  # file: first lines rng, state, lnL; last line names node 2 n - 2
  local f=$1 nt; nt=$(grep -c '<taxon id=' "$CHAINDIR/$STEM.xml")
  [ -s "$f" ] && [ "$(sed -n 1p "$f" | cut -f1)" = rng ] && [ "$(sed -n 3p "$f" | cut -f1)" = lnL ] && [[ "$(state_number "$f")" =~ ^[0-9]+$ ]] && [ "$(tail -1 "$f" | cut -f1)" = "$((2*nt-2))" ] && [ -z "$(tail -c1 "$f")" ]
}
stop_child() { pkill -TERM -P "$1" 2>/dev/null; kill -TERM "$1" 2>/dev/null; sleep 2; pkill -KILL -P "$1" 2>/dev/null; kill -KILL "$1" 2>/dev/null; wait "$1" 2>/dev/null; }
gone() { [ ! -e "$CHAINDIR/run_pid_v2.txt" ] || [ ! -e "$ROOT/sched/ledger.csv" ]; }
[ -s "$CHAINDIR/$STEM.xml" ] || cp "$ROOT/xml/$STEM.xml" "$CHAINDIR/$STEM.xml"
echo $$ > "$CHAINDIR/run_pid_v2.txt.tmp"

if [ "$MODE" = "fresh" ]; then
  DIR="$CHAINDIR"; cd "$DIR" || exit 99
  CMD="beast -threads $THREADS -beagle_CPU -beagle_SSE -beagle_threads $THREADS -seed $SEED -save_every $SAVE_EVERY -save_state $LABEL.state -overwrite $STEM.xml"
  if [ "$N" -gt 1 ]; then EVENT="started_again_from_state_0"; else EVENT="started"; fi
  NOTE="start $N"
  echo "$CMD" > command.txt
  echo "JAVA_TOOL_OPTIONS=$JAVA_TOOL_OPTIONS; PATH begins with $BEAST_BIN; directory ${DIR#$ROOT/}" > environment.txt
  now > started.txt
  mv -f "$CHAINDIR/run_pid_v2.txt.tmp" "$CHAINDIR/run_pid_v2.txt"
  ledger "$EVENT" "" "" "$NOTE"
  $CMD > stdout.txt 2>&1 9>&- &
  CHILD=$!
  while kill -0 "$CHILD" 2>/dev/null; do
    sleep 15
    if gone; then pkill -TERM -P "$CHILD" 2>/dev/null; kill -TERM "$CHILD" 2>/dev/null; exit 98; fi
  done
  wait "$CHILD"; RC=$?
  echo "$RC" > exit_code.txt.tmp && mv -f exit_code.txt.tmp exit_code.txt
  now > finished.txt
  LS=$(last_state "$STEM.log"); echo "$LS" > last_state_logged.txt
  ledger "ended" "$RC" "$LS" "$NOTE"
  rm -f "$CHAINDIR/run_pid_v2.txt"
  exit 0
fi

# ------------------------------------------------------------------------------------------------ mode continue
if [ "$N" -gt 1 ]; then PREV="$CHAINDIR/continuation_$((N-1))"; PART="continuation_$((N-1))"; else PREV="$CHAINDIR"; PART="first_part"; fi
ENDED_AT=$(last_state "$PREV/$STEM.log")
CANDS=("$PREV/$LABEL.state")
if [ -d "$ROOT/state_copies/$ANALYSIS/$LABEL" ]; then
  while read -r c; do [ -n "$c" ] && CANDS+=("$ROOT/state_copies/$ANALYSIS/$LABEL/$c"); done < <(ls -1 "$ROOT/state_copies/$ANALYSIS/$LABEL" | grep -E "^$LABEL\.$PART\.state\.[0-9]{9}$" | sort -r)
fi
now > "$CHAINDIR/started_continuation_$N.txt"
mv -f "$CHAINDIR/run_pid_v2.txt.tmp" "$CHAINDIR/run_pid_v2.txt"
ATT=0; ACCEPTED=0; TRIED=""
for CAND in "${CANDS[@]}"; do
  [ -e "$CAND" ] || continue
  SUM=$(sha256sum "$CAND" | cut -c1-64); case " $TRIED " in *" $SUM "*) continue;; esac; TRIED="$TRIED $SUM"
  ATT=$((ATT+1)); DIR="$CHAINDIR/continuation_$N"; REJ="$CHAINDIR/rejected_continuation_${N}_attempt_$ATT"
  mkdir -p "$DIR"; cd "$DIR" || exit 99
  CK=$(state_number "$CAND"); WHY=""
  echo "saved state tried: ${CAND#$ROOT/} (SHA-256 $SUM); state ${CK:-not readable}; the part before was ended at state $ENDED_AT" > continuation_info.txt
  if ! complete_state "$CAND"; then WHY="saved state incomplete"
  elif [ -n "$ENDED_AT" ] && [ "$CK" -gt "$ENDED_AT" ]; then WHY="saved state $CK lies after the last state logged $ENDED_AT"
  elif [ "$CK" -ge "$STEPS" ]; then WHY="saved state $CK is not below the length of the chain"
  fi
  if [ -z "$WHY" ]; then
    cp -p "$CAND" "loaded_state_$CK.state"
    REMAIN=$(python3 "$ROOT/tools/make_continuation_xml.py" "$CHAINDIR/$STEM.xml" "$DIR/$STEM.xml" "$STEPS" "$CK" 2> xml_error.txt) || WHY="XML of the continuation was not made: $(cat xml_error.txt)"
    [ -s xml_error.txt ] || rm -f xml_error.txt
  fi
  if [ -z "$WHY" ]; then
    CMD="beast -threads $THREADS -beagle_CPU -beagle_SSE -beagle_threads $THREADS -seed $SEED -load_state loaded_state_$CK.state -force_resume -save_every $SAVE_EVERY -save_state $LABEL.state -overwrite $STEM.xml"
    echo "$CMD" > command.txt
    echo "JAVA_TOOL_OPTIONS=$JAVA_TOOL_OPTIONS; PATH begins with $BEAST_BIN; directory ${DIR#$ROOT/}" > environment.txt
    echo "state of the checkpoint $CK; chain length of the XML of the continuation $REMAIN; length of the chain $STEPS" >> continuation_info.txt
    now > started.txt
    ledger "continued_from_saved_state" "" "" "state $CK; continuation $N; attempt $ATT; first part or continuation before ended at state $ENDED_AT; chain length of the XML $REMAIN"
    $CMD > stdout.txt 2>&1 9>&- &
    CHILD=$!
    CHECKED=0
    while kill -0 "$CHILD" 2>/dev/null; do
      sleep 15
      if gone; then pkill -TERM -P "$CHILD" 2>/dev/null; kill -TERM "$CHILD" 2>/dev/null; exit 98; fi
      if [ "$CHECKED" = 0 ] && [ -n "$(last_state "$STEM.log")" ]; then
        CHECKED=1
        LINE=$(grep -m1 -E "Forcing analysis to resume regardless of recomputed likelihood values|COMPARING LIKELIHOODS" stdout.txt || true)
        if [ -n "$LINE" ]; then
          DIFF=$(echo "$LINE" | grep -oE -- '-?[0-9]+\.[0-9]+(E-?[0-9]+)?' | head -2 | awk 'NR==1{a=$1} NR==2{b=$1} END{d=a-b; if(d<0)d=-d; printf "%.6e", d}')
          echo "$LINE" > likelihood_check.txt; echo "absolute difference $DIFF" >> likelihood_check.txt
          if awk -v d="$DIFF" -v l="$LIMIT" 'BEGIN{exit !(d+0 >= l+0)}'; then WHY="difference of the log-likelihoods $DIFF not below $LIMIT"; fi
        else
          DIFF=0; echo "BEAST X reported no difference between the saved and the recomputed log-likelihood" > likelihood_check.txt; echo "absolute difference 0" >> likelihood_check.txt
        fi
        echo "difference of the log-likelihoods $DIFF" >> continuation_info.txt
        [ -n "$WHY" ] && { stop_child "$CHILD"; break; }
      fi
    done
    if [ -z "$WHY" ]; then
      wait "$CHILD"; RC=$?
      if [ "$CHECKED" = 0 ]; then WHY="BEAST X ended with code $RC before it had logged a state"; fi
    fi
  fi
  if [ -n "$WHY" ]; then
    echo "$(now) attempt $ATT of continuation $N not used: $WHY" > continuation_rejected.txt
    ledger "continuation_rejected" "" "$(last_state "$STEM.log")" "continuation $N, attempt $ATT, saved state ${CK:-not readable}: $WHY"
    cd "$CHAINDIR"; mv "$DIR" "$REJ"
    continue
  fi
  ACCEPTED=1
  echo "$RC" > exit_code.txt.tmp && mv -f exit_code.txt.tmp exit_code.txt
  now > finished.txt
  LS=$(last_state "$STEM.log"); echo "$LS" > last_state_logged.txt
  ledger "ended" "$RC" "$LS" "state $CK; continuation $N; attempt $ATT"
  break
done
if [ "$ACCEPTED" = 0 ]; then
  echo "$(now) continuation $N: no saved state could be used ($ATT tried); the chain is to be started again from state 0" > "$CHAINDIR/all_saved_states_rejected.txt"
  CMD=""; ledger "no_saved_state_usable" "" "$ENDED_AT" "continuation $N: $ATT saved states tried, none used"
fi
rm -f "$CHAINDIR/run_pid_v2.txt"
exit 0
