#!/bin/sh
# contradiction_halt_test.sh <cobound> <test data dir>
#
# The cobordism graph's contradiction gates run in every run (plan divergence
# 6), and a contradiction halts it, each with its documented code: 2 without
# a goal (depth 0), 3 with one. The crafted case, as cobordismgraph_test
# builds its gate cases: a literature lower bound made falsely high, and a
# search that finds a surface below it. The Hopf link L2a1{0} bounds an
# annulus (a genus-0 surface for the link alone, found at once); its table
# entry is rewritten to claim 4-genus 1. The surfaces are real, so a run
# writes what it found before it halts (depth 0: the database).
set -eu

C=$1
DATA=$2
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
# The crafted tables: L2a1{0}'s 4-genus 0 becomes 1.
sed 's/^\(L2a1{0},.*\),0$/\1,1/' "$DATA/links_to_6.csv" > "$T/links.csv"
grep -q '^L2a1{0},.*,1$' "$T/links.csv" || { echo "FAIL: could not craft the table"; exit 1; }
{ head -1 "$T/links.csv"; grep -F 'L2a1{0},' "$T/links.csv"; } > "$T/rows.csv"
echo "$HDR" > "$T/cobordisms.csv"

# 1. Depth 0: exit 2, the FATAL banner, and the target's cobordisms on disk.
rc=0
cat > "$T/depth0.conf" <<CONF
targets = $T/rows.csv
verdicts = $T/out.csv
cobordisms = $T/cobordisms.csv
work = $T/work
census = $T/none.sqlite
knot_table = $DATA/knots_to_6.csv
link_table = $T/links.csv
census_updates = 0
retriangulate_on_miss = 0
resolve_unlinked = 0
max_faces = 3
threads = 2
CONF
"$C" run --config "$T/depth0.conf" > "$T/log" 2> "$T/err" || rc=$?
grep -E 'CONTRADICTION|contradict|BELOW|FATAL' "$T/log" "$T/err" | head -5 || true
if [ "$rc" -ne 2 ] || ! grep -q 'FATAL' "$T/err"; then
  echo "FAIL: depth 0 exited $rc (want 2, halted)"; exit 1; fi
if [ "$(wc -l < "$T/cobordisms.csv")" -lt 2 ]; then
  echo "FAIL: the halted row's witnesses were not written first"; exit 1; fi

# 2. With a goal: exit 3, even though the same surface meets the goal.
pd=$(grep -F 'L2a1{0},' "$T/links.csv" | cut -d, -f2)
rc=0
cat > "$T/goal.conf" <<CONF
target_pd = $pd
target_name = L2a1{0}
work = $T/goal
knot_table = $DATA/knots_to_6.csv
link_table = $T/links.csv
census = $T/none.sqlite
goal_genus = 0
literature = 0
threads = 2
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
"$C" run --config "$T/goal.conf" > "$T/goal.log" 2>&1 || rc=$?
grep -E 'CONTRADICTION|outcome' "$T/goal.log" | head -5 || true
if [ "$rc" -ne 3 ] || ! grep -q 'CONTRADICTION' "$T/goal.log"; then
  echo "FAIL: the goal run exited $rc (want 3, a contradiction)"; exit 1; fi
if [ ! -s "$T/goal/hop_0_n0/log.txt" ]; then
  echo "FAIL: the goal run's hop was not recorded"; exit 1; fi
echo "PASS: a contradiction halts every run: 2 at depth 0, 3 with a goal"
