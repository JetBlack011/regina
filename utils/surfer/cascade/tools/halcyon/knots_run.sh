#!/bin/bash
# On halcyon: the knots in <list.csv> (pool_A.csv layout) the atlas has not verified,
# each proved constructively to its literature g4 by
# the cascade (v8 build), while the Pool B baseline is paused:
#   1. pause Pool B's loop between targets (SIGSTOP its subshell by PID; the
#      target in flight finishes, and its driver.log keeps its true wall);
#   2. pass 1: the cascade alone, master withheld, 7,200 CPU-s each, hops
#      50k surfaces, a revisit up to 1M;
#   3. pass 2, for any knot still unmet: the same with the atlas's witnesses
#      as free edges (--master-witnesses, read only) -- constructive facts,
#      what verifyslicegenus itself chains through;
#   4. every certificate replayed by cascade_check.py;
#   5. Pool B resumed.
# usage: knots_run.sh <list.csv> <tag>. Log: <tag>.log.
B=$HOME/cascade-bench
D=$HOME/Projects/cobordism-atlas/data
MASTER=$HOME/Projects/cobordism-atlas/results/cobordisms.csv
V=$B/v3
LIST=$1 TAG=$2
L=$B/$TAG.log
log() { echo "$*" >> "$L"; }

# pool() as the Pool B baseline runs it.
eval "$(sed -n '/^pool() {/,/^}/p' "$B/cascade_v8.sh")"

# Whether a process is still running: not gone, and not a zombie (the
# paused loop cannot reap its finished target, and kill -0 succeeds on a
# zombie).
running() { local s; s=$(ps -o stat= -p "$1" 2>/dev/null); [ -n "$s" ] && [ "${s:0:1}" != Z ]; }

# 1. Pause Pool B: its loop subshell is the parent of the cascadesearch
# working under results/poolB_base (bracketed pattern: never matches awk).
# PAUSED=<pid> adopts a loop already stopped.
paused=${PAUSED:-}
[ -n "$paused" ] && log "Pool B loop $paused already paused; adopted $(date -Is)"
for i in $(seq 1 60); do
  [ -n "$paused" ] && break
  read -r cs sub < <(ps -eo pid,ppid,cmd | awk '/[c]ascade\/cascadesearch/ && /poolB_base/ {print $1, $2; exit}')
  if [ -n "$sub" ]; then
    kill -STOP "$sub" && paused=$sub
    log "Pool B paused: loop $sub stopped, target in flight $cs $(date -Is)"
    while running "$cs"; do sleep 2; done
    log "Pool B's target in flight finished $(date -Is)"
    break
  fi
  grep -q 'BATCH DONE' "$B/results/poolB_base.txt" 2>/dev/null && { log "Pool B already done"; break; }
  sleep 2
done

# 2. Pass 1: the cascade alone.
log "pass 1 (master withheld) starts $(date -Is)"
pool "$B/$LIST" ${TAG}_pure --hop-surfaces 50000 --max-hop-surfaces 1000000 \
  --max-expansions 100000 --cpu-budget 7200
log "pass 1 done $(date -Is)"
sed 's/ :: .*//' "$B/results/${TAG}_pure.txt" | while read -r line; do log "  $line"; done
grep -h 'GOAL MET' "$B/results/${TAG}_pure.txt" | cut -d' ' -f1 | sed 's/^/  met: /' >> "$L"

# 3. Pass 2: the unmet, with the atlas's witnesses.
head -1 "$B/$LIST" > "$B/${TAG}_unmet.csv"
tail -n +2 "$B/$LIST" | while IFS=, read -r name rest; do
  grep -q "^$name .*GOAL MET" "$B/results/${TAG}_pure.txt" || echo "$name,$rest" >> "$B/${TAG}_unmet.csv"
done
if [ "$(wc -l < "$B/${TAG}_unmet.csv")" -gt 1 ] && [ -f "$MASTER" ]; then
  log "pass 2 (with master witnesses) for $(tail -n +2 "$B/${TAG}_unmet.csv" | cut -d, -f1 | tr '\n' ' ') starts $(date -Is)"
  pool "$B/${TAG}_unmet.csv" ${TAG}_master --hop-surfaces 50000 --max-hop-surfaces 1000000 \
    --max-expansions 100000 --cpu-budget 7200 --master-witnesses "$MASTER"
  log "pass 2 done $(date -Is)"
  sed 's/ :: .*//' "$B/results/${TAG}_master.txt" | while read -r line; do log "  $line"; done
fi

# 4. Replay every certificate.
for c in "$B"/results/${TAG}_*/*/certificate.json; do
  [ -f "$c" ] || continue
  out=${c%certificate.json}check.txt
  CASCADE_ATLAS=$HOME/Projects/cobordism-atlas "$HOME/.venvs/atlas/bin/python" \
    "$V/utils/surfer/cascade/tools/cascade_check.py" \
    --farsidediagram "$V/build/utils/surfer/farsidediagram" "$c" > "$out" 2>&1
  log "check $(basename "$(dirname "$(dirname "$c")")")/$(basename "$(dirname "$c")"): $(tail -1 "$out")"
done

# 5. Resume Pool B.
if [ -n "$paused" ]; then kill -CONT "$paused"; log "Pool B resumed $(date -Is)"; fi
log "$TAG done $(date -Is)"
