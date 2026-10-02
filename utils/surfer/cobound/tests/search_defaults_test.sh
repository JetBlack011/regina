#!/bin/sh
# search_defaults_test.sh <verifyslicegenus> <cascadesearch> <test data dir>
#
# One set of defaults for every search (plan divergence 3):
#   - thicken and collar layers 2/2, no cone, and the boundary condition
#     `proper`, in verifyslicegenus as in the cascade: a row searched with no
#     shape options at all is the row searched with them spelled out;
#   - resolve_unlinked has NO default: it changes which surfaces count toward
#     a surface target and is part of every frontier's fingerprint, so a run
#     that builds a search and does not say which is refused, by both
#     drivers. A solve builds no search and needs neither.
set -eu

V=$1
C=$2
DATA=$3
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
head -2 "$DATA/knots_to_6.csv" > "$T/rows.csv"
PD=$(sed -n 2p "$DATA/knots_to_6.csv" | cut -d, -f2)

run() { # dir args...
  d=$T/$1; shift
  mkdir -p "$d"
  echo "$HDR" > "$d/cobordisms.csv"
  rc=0
  "$V" --input "$T/rows.csv" --output "$d/out.csv" --cobordisms "$d/cobordisms.csv" \
       --census-db "$d/none.sqlite" --knot-table "$DATA/knots_to_6.csv" \
       --link-table "$DATA/links_to_6.csv" --no-census-updates --no-retriangulate-on-miss \
       --max-faces 3 --threads 2 "$@" > "$d/log" 2> "$d/err" || rc=$?
  return $rc
}

# 1. verifyslicegenus refuses a search without resolve_unlinked.
rc=0; run none || rc=$?
if [ "$rc" -eq 0 ] || ! grep -q 'resolve_unlinked has no default' "$T/none/err"; then
  echo "FAIL: verifyslicegenus searched without --resolve-unlinked/--no-resolve-unlinked (exit $rc)"
  exit 1
fi
# ... but a solve needs neither.
mkdir -p "$T/solve"; echo "$HDR" > "$T/solve/cobordisms.csv"
if ! "$V" --input "$T/rows.csv" --output "$T/solve/out.csv" --cobordisms "$T/solve/cobordisms.csv" \
     --knot-table "$DATA/knots_to_6.csv" --link-table "$DATA/links_to_6.csv" --solve-only \
     > "$T/solve/log" 2>&1; then
  echo "FAIL: --solve-only was refused without resolve_unlinked"; exit 1
fi

# 2. The default shape is the spelled-out 2/2, no cone, proper.
run defaults --no-resolve-unlinked
run spelled --no-resolve-unlinked --thicken-layers 2 --collar-layers 2 --no-cone \
    --boundary-condition proper
acct() { grep ': accounting:' "$1" | sed -E 's/recorded [0-9]+, duplicate [0-9]+/recorded <n>, duplicate <n>/'; }
acct "$T/defaults/log" > "$T/defaults.acct"
acct "$T/spelled/log" > "$T/spelled.acct"
cat "$T/defaults.acct"
if [ ! -s "$T/defaults.acct" ] || ! cmp -s "$T/defaults.acct" "$T/spelled.acct"; then
  echo "FAIL: the default search shape is not 2/2, no cone, proper"
  cat "$T/spelled.acct"; exit 1
fi
# 3_1 at cap 3 under that shape is the canary's (canaries.expected: accepted 1752).
if ! grep -q 'accepted 1752, described 1752' "$T/defaults.acct"; then
  echo "FAIL: 3_1 at cap 3 did not search the canary shape"; exit 1
fi

# 3. cascadesearch refuses a run without resolve_unlinked.
rc=0
"$C" --target-pd "$PD" --target-name 3_1 --work "$T/goal" --knot-table "$DATA/knots_to_6.csv" \
     --link-table "$DATA/links_to_6.csv" --census-db "$T/none.sqlite" --goal-genus 1 \
     --max-expansions 0 --threads 2 > "$T/goal.log" 2>&1 || rc=$?
if [ "$rc" -eq 0 ] || ! grep -q 'resolve_unlinked has no default' "$T/goal.log"; then
  echo "FAIL: cascadesearch ran without --resolve-unlinked/--no-resolve-unlinked (exit $rc)"
  exit 1
fi
echo "PASS: one set of search defaults, and resolve_unlinked required by both drivers"
