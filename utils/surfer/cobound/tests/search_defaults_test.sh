#!/bin/sh
# search_defaults_test.sh <cobound> <test data dir>
#
# One set of defaults for every search (plan divergence 3; the config, phase
# 5):
#   - layers 2, collared through both, and the boundary condition `proper`:
#     a run whose config states no shape is the run with it spelled out;
#   - resolve_unlinked has NO default: it changes which surfaces count toward
#     a surface target and is part of every frontier's fingerprint, so a run
#     that does not say which is refused, with a goal and without. Nor has
#     `work` (a run's pending files go there). A solve builds no search and
#     needs neither;
#   - the retired options are refused: the config replaces them;
#   - every run writes the configuration it ran with to <work>/cobound.conf,
#     and that file, read back, runs the same search.
set -eu

C=$1
DATA=$2
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
head -2 "$DATA/knots_to_6.csv" > "$T/rows.csv"
PD=$(sed -n 2p "$DATA/knots_to_6.csv" | cut -d, -f2)
cat > "$T/base.conf" <<CONF
targets = $T/rows.csv
knot_table = $DATA/knots_to_6.csv
link_table = $DATA/links_to_6.csv
census_updates = 0
retriangulate_on_miss = 0
max_faces = 3
threads = 2
CONF

run() { # dir --set... ; the run's own files under dir
  d=$T/$1; shift
  mkdir -p "$d"
  echo "$HDR" > "$d/cobordisms.csv"
  rc=0
  "$C" run --config "$T/base.conf" --set "verdicts=$d/out.csv" \
       --set "cobordisms=$d/cobordisms.csv" --set "census=$d/none.sqlite" "$@" \
       > "$d/log" 2> "$d/err" || rc=$?
  return $rc
}

# 1. A run without a goal refuses to search without resolve_unlinked, or
#    without work.
rc=0; run none --set "work=$T/none/work" || rc=$?
if [ "$rc" -eq 0 ] || ! grep -q 'resolve_unlinked is required' "$T/none/err"; then
  echo "FAIL: a run searched without resolve_unlinked (exit $rc)"; cat "$T/none/err"
  exit 1
fi
rc=0; run nowork --set resolve_unlinked=0 || rc=$?
if [ "$rc" -eq 0 ] || ! grep -q 'work is required' "$T/nowork/err"; then
  echo "FAIL: a run searched without work (exit $rc)"; exit 1
fi
# ... but a solve needs neither.
mkdir -p "$T/solve"; echo "$HDR" > "$T/solve/cobordisms.csv"
if ! "$C" solve --set "targets=$T/rows.csv" --set "verdicts=$T/solve/out.csv" \
     --set "cobordisms=$T/solve/cobordisms.csv" --set "knot_table=$DATA/knots_to_6.csv" \
     --set "link_table=$DATA/links_to_6.csv" > "$T/solve/log" 2>&1; then
  echo "FAIL: solve was refused without resolve_unlinked"; cat "$T/solve/log"; exit 1
fi
# The retired options are refused.
rc=0; "$C" run --max-faces 3 > "$T/retired.log" 2>&1 || rc=$?
if [ "$rc" -ne 2 ] || ! grep -q 'unknown option --max-faces' "$T/retired.log"; then
  echo "FAIL: a retired option was accepted (exit $rc)"; exit 1
fi

# 2. The default shape is the spelled-out 2 layers, proper.
run defaults --set resolve_unlinked=0 --set "work=$T/defaults/work"
run spelled --set resolve_unlinked=0 --set "work=$T/spelled/work" --set layers=2 \
    --set boundary_condition=proper
acct() { grep ': accounting:' "$1" | sed -E 's/recorded [0-9]+, duplicate [0-9]+/recorded <n>, duplicate <n>/'; }
acct "$T/defaults/log" > "$T/defaults.acct"
acct "$T/spelled/log" > "$T/spelled.acct"
cat "$T/defaults.acct"
if [ ! -s "$T/defaults.acct" ] || ! cmp -s "$T/defaults.acct" "$T/spelled.acct"; then
  echo "FAIL: the default search shape is not 2 layers, proper"
  cat "$T/spelled.acct"; exit 1
fi
# 3_1 at cap 3 under that shape is the canary's (canaries.expected: accepted 1752).
if ! grep -q 'accepted 1752, described 1752' "$T/defaults.acct"; then
  echo "FAIL: 3_1 at cap 3 did not search the canary shape"; exit 1
fi

# 3. The configuration it ran with, read back, runs the same search.
conf=$T/defaults/work/cobound.conf
if [ ! -s "$conf" ] || ! grep -q '^resolve_unlinked = 0$' "$conf" || ! grep -q '^layers = 2$' "$conf"; then
  echo "FAIL: the run did not write its configuration to $conf"; exit 1
fi
d=$T/again; mkdir -p "$d"; echo "$HDR" > "$d/cobordisms.csv"
"$C" run --config "$conf" --set "verdicts=$d/out.csv" --set "cobordisms=$d/cobordisms.csv" \
     --set "work=$d/work" > "$d/log" 2> "$d/err"
acct "$d/log" > "$T/again.acct"
if ! cmp -s "$T/defaults.acct" "$T/again.acct"; then
  echo "FAIL: the written configuration ran a different search"; exit 1
fi

# 4. A run with a goal refuses without resolve_unlinked too.
rc=0
"$C" run --set "target_pd=$PD" --set target_name=3_1 --set "work=$T/goal" \
     --set "knot_table=$DATA/knots_to_6.csv" --set "link_table=$DATA/links_to_6.csv" \
     --set "census=$T/none.sqlite" --set goal_genus=1 --set max_searches=0 --set threads=2 \
     > "$T/goal.log" 2>&1 || rc=$?
if [ "$rc" -eq 0 ] || ! grep -q 'resolve_unlinked is required' "$T/goal.log"; then
  echo "FAIL: a goal run ran without resolve_unlinked (exit $rc)"
  exit 1
fi
echo "PASS: one set of search defaults, resolve_unlinked and work required, the retired options refused, and the written configuration runs the same search"
