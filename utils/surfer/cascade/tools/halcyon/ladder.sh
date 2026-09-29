#!/bin/bash
# On halcyon: a list of knots (pool_A.csv layout, goal = literature g4) up a
# ladder of cascade settings, each rung taking only what the rungs below left
# unmet, then every certificate replayed and every CERTIFIED proof recorded:
#   1  the cascade alone (master withheld), 50k-surface hops, B1 CPU-s
#   2  with the atlas's witnesses as free edges, B2 CPU-s
#   3  with the atlas's witnesses and 1M-surface hops (up to 4M), B3 CPU-s
# usage: ladder.sh <list.csv> <tag>   (B1, B2, B3 default 1800, 1800, 7200)
# Log: <tag>.log; results in results/<tag>_r1 .. _r3; recorded as runs
# 2026-09-28_<tag>_r1 .. _r3.
B=$HOME/cascade-bench
D=$HOME/Projects/cobordism-atlas/data
MASTER=$HOME/Projects/cobordism-atlas/results/cobordisms.csv
V=$B/v3
LIST=$1 TAG=$2
B1=${B1:-1800} B2=${B2:-1800} B3=${B3:-7200}
L=$B/$TAG.log
log() { echo "$*" >> "$L"; }
eval "$(sed -n '/^pool() {/,/^}/p' "$B/cascade_v8.sh")"

unmet_of() { # <list> <results txt> <out list>
  head -1 "$1" > "$3"
  tail -n +2 "$1" | while IFS=, read -r name rest; do
    grep -q "^$name .*GOAL MET" "$2" 2>/dev/null || echo "$name,$rest" >> "$3"
  done
}
count() { echo $(( $(wc -l < "$1") - 1 )); }

log "ladder $TAG: $(count "$B/$LIST") targets, budgets $B1/$B2/$B3 CPU-s, starts $(date -Is)"
cur=$B/$LIST
for rung in 1 2 3; do
  [ "$(count "$cur")" -eq 0 ] && break
  case $rung in
    1) flags=(--hop-surfaces 50000 --max-hop-surfaces 1000000 --cpu-budget "$B1") ;;
    2) flags=(--hop-surfaces 50000 --max-hop-surfaces 1000000 --cpu-budget "$B2"
              --master-witnesses "$MASTER") ;;
    3) flags=(--hop-surfaces 1000000 --max-hop-surfaces 4000000 --cpu-budget "$B3"
              --master-witnesses "$MASTER") ;;
  esac
  log "rung $rung: $(count "$cur") targets starts $(date -Is)"
  pool "$cur" "${TAG}_r$rung" --max-expansions 100000 "${flags[@]}"
  met=$(grep -c 'GOAL MET' "$B/results/${TAG}_r$rung.txt")
  log "rung $rung done $(date -Is): $met met"
  unmet_of "$cur" "$B/results/${TAG}_r$rung.txt" "$B/${TAG}_unmet$rung.csv"
  cur=$B/${TAG}_unmet$rung.csv
done
log "unmet after the ladder: $(tail -n +2 "$cur" | cut -d, -f1 | tr '\n' ' ')"

for c in "$B"/results/${TAG}_r*/*/certificate.json; do
  [ -f "$c" ] || continue
  out=${c%certificate.json}check.txt
  CASCADE_ATLAS=$HOME/Projects/cobordism-atlas "$HOME/.venvs/atlas/bin/python" \
    "$V/utils/surfer/cascade/tools/cascade_check.py" \
    --farsidediagram "$V/build/utils/surfer/farsidediagram" "$c" > "$out" 2>&1
  tail -1 "$out" | grep -q 'VERDICT: CERTIFIED' || log "CHECK FAILED: $c: $(tail -1 "$out")"
done
log "checks: $(grep -l 'VERDICT: CERTIFIED' "$B"/results/${TAG}_r*/*/check.txt 2>/dev/null | wc -l) CERTIFIED"

args=()
for rung in 1 2 3; do
  [ -d "$B/results/${TAG}_r$rung" ] && args+=("2026-09-28_${TAG}_r$rung=${TAG}_r$rung")
done
bash "$B/record_all.sh" "${args[@]}"
log "recorded $(date -Is)"
log "ladder $TAG done $(date -Is)"
