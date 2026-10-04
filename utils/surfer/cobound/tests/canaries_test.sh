#!/bin/sh
# canaries_test.sh <cobound> <test data dir>
#
# The atlas's canaries, as a cobound test pinning phase 0's numbers:
#   - canaries.sh's nine rows (data/canaries.csv: 3_1, 6_1, 8_8, 8_20 and five
#     links, their table lines verbatim), searched EXHAUSTIVELY at cap 3 by one
#     run without a goal, census-free; each row's accounting line, the
#     recorded/duplicate split summed (it follows names, which are not
#     canonical), must be data/canaries.expected's -- the atlas's
#     tools/orchestrate/canaries.expected, phase 0's numbers, byte for byte;
#   - cascade_canaries.sh's two goal runs (3_1 goal 1 and L2a1{0} goal 0, one
#     exhaustive search each): each search's accounting, the target's best
#     genus and the run's outcome must be data/goal_canaries.expected's (the
#     atlas's cascade_canaries.expected, byte for byte).
# Exhaustive and name-free, so the numbers are a function of the code alone:
# they hold with the small tables here as with the atlas's (naming moves only
# the recorded/duplicate split).
set -eu

C=$1
D=$2
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
echo "$HDR" > "$T/cobordisms.csv"
cat > "$T/depth0.conf" <<CONF
targets = $D/canaries.csv
verdicts = $T/out.csv
cobordisms = $T/cobordisms.csv
work = $T/work
census = $T/none.sqlite
knot_table = $D/knots_to_6.csv
link_table = $D/links_to_6.csv
census_updates = 0
retriangulate_on_miss = 0
layers = 2
max_faces = 3
boundary_condition = proper
resolve_unlinked = 0
threads = ${CANARY_THREADS:-4}
CONF
"$C" run --config "$T/depth0.conf" > "$T/log" 2> "$T/err" || {
  echo "FAIL: cobound run exited $?"; tail -5 "$T/err"; exit 1; }
grep ': accounting:' "$T/log" | sed -E \
  's/^\[\+\] (.*): accounting: accepted ([0-9]+), described ([0-9]+), recorded ([0-9]+), duplicate ([0-9]+), (.*)/\1 accepted=\2 described=\3 examined=\4+\5 \6/' |
  awk '{ split($4, rd, "="); split(rd[2], p, "+"); $4 = "examined=" (p[1] + p[2]); print }' |
  LC_ALL=C sort > "$T/got"
if ! diff -u "$D/canaries.expected" "$T/got"; then
  echo "FAIL: the canaries' accounting is not phase 0's"; exit 1
fi

goal() { # name goal pd
  cat > "$T/goal.conf" <<CONF
target_pd = $3
target_name = $1
work = $T/goal-$2
knot_table = $D/knots_to_6.csv
link_table = $D/links_to_6.csv
census = $T/none.sqlite
goal_genus = $2
literature = 0
threads = ${CANARY_THREADS:-4}
max_searches = 1
surface_target = 1000000000
max_surface_target = 1000000000
max_faces = 3
iddfs_iterations = 0
iddfs_start = 0
iddfs_step = 0
root_budget_start = 0
resolve_unlinked = 1
CONF
  "$C" run --config "$T/goal.conf" > "$T/goal-$2.log" 2>&1 || true
  grep -E '^\[\+\] hop [0-9]+ .*: accounting:' "$T/goal-$2.log" | sed -E \
    's/^\[\+\] hop ([0-9]+) (.*): accounting: accepted ([0-9]+), described ([0-9]+), recorded ([0-9]+), duplicate ([0-9]+), (.*)/\2 hop\1 accepted=\3 described=\4 examined=\5+\6 \7/' |
    awk '{ split($5, rd, "="); split(rd[2], p, "+"); $5 = "examined=" (p[1] + p[2]); print }'
  echo "$1 $(grep -oE 'Target best: [0-9a-z]+' "$T/goal-$2.log")" \
       "$(grep -oE 'outcome [-a-z]+' "$T/goal-$2.log" | tail -1)"
}
pd31=$(grep '^3_1,' "$D/canaries.csv" | cut -d, -f2)
pdhopf=$(grep -F 'L2a1{0},' "$D/canaries.csv" | cut -d, -f2)
{ goal 3_1 1 "$pd31"; goal 'L2a1{0}' 0 "$pdhopf"; } | LC_ALL=C sort > "$T/cascade.got"
if ! diff -u "$D/goal_canaries.expected" "$T/cascade.got"; then
  echo "FAIL: the goal runs' accounting, best genus or outcome is not phase 0's"; exit 1
fi
echo "PASS: $(wc -l < "$T/got") canary rows and $(wc -l < "$T/cascade.got") goal-run lines are phase 0's"
