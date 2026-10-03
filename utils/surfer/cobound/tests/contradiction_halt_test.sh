#!/bin/sh
# contradiction_halt_test.sh <verifyslicegenus> <cascadesearch> <test data dir>
#
# The cobordism graph's contradiction gates run in every run (plan divergence
# 6), and a contradiction halts it, each with its documented code: 2 without
# a goal (depth 0), 3 with one. The crafted case, as cobordismgraph_test
# builds its gate cases: a literature lower bound made falsely high, and a
# search that finds a surface below it. The Hopf link L2a1{0} bounds an
# annulus (a genus-0 surface for the link alone, found at once); its table
# entry is rewritten to claim 4-genus 1. The surfaces are real, so a run
# writes what it found before it halts (depth 0: the witness file).
set -eu

V=$1
C=$2
DATA=$3
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
# The crafted tables: L2a1{0}'s 4-genus 0 becomes 1.
sed 's/^\(L2a1{0},.*\),0$/\1,1/' "$DATA/links_to_6.csv" > "$T/links.csv"
grep -q '^L2a1{0},.*,1$' "$T/links.csv" || { echo "FAIL: could not craft the table"; exit 1; }
{ head -1 "$T/links.csv"; grep -F 'L2a1{0},' "$T/links.csv"; } > "$T/rows.csv"
echo "$HDR" > "$T/cobordisms.csv"

# 1. Depth 0: exit 2, the FATAL banner, and the row's witnesses on disk.
rc=0
"$V" --input "$T/rows.csv" --output "$T/out.csv" --cobordisms "$T/cobordisms.csv" \
     --census-db "$T/none.sqlite" --knot-table "$DATA/knots_to_6.csv" --link-table "$T/links.csv" \
     --no-census-updates --no-retriangulate-on-miss --no-resolve-unlinked \
     --max-faces 3 --threads 2 > "$T/log" 2> "$T/err" || rc=$?
grep -E 'CONTRADICTION|contradict|BELOW|FATAL' "$T/log" "$T/err" | head -5 || true
if [ "$rc" -ne 2 ] || ! grep -q 'FATAL' "$T/err"; then
  echo "FAIL: depth 0 exited $rc (want 2, halted)"; exit 1; fi
if [ "$(wc -l < "$T/cobordisms.csv")" -lt 2 ]; then
  echo "FAIL: the halted row's witnesses were not written first"; exit 1; fi

# 2. With a goal: exit 3, even though the same surface meets the goal.
pd=$(grep -F 'L2a1{0},' "$T/links.csv" | cut -d, -f2)
rc=0
"$C" --target-pd "$pd" --target-name 'L2a1{0}' --work "$T/goal" \
     --knot-table "$DATA/knots_to_6.csv" --link-table "$T/links.csv" --census-db "$T/none.sqlite" \
     --goal-genus 0 --constructive --threads 2 --max-expansions 1 \
     --hop-surfaces 1000000000 --max-hop-surfaces 1000000000 --hop-max-faces 3 \
     --hop-iddfs-iterations 0 --hop-iddfs-start 0 --hop-iddfs-step 0 --hop-root-budget 0 \
     --resolve-unlinked > "$T/goal.log" 2>&1 || rc=$?
grep -E 'CONTRADICTION|outcome' "$T/goal.log" | head -5 || true
if [ "$rc" -ne 3 ] || ! grep -q 'CONTRADICTION' "$T/goal.log"; then
  echo "FAIL: the goal run exited $rc (want 3, a contradiction)"; exit 1; fi
if [ ! -s "$T/goal/hop_0_n0/log.txt" ]; then
  echo "FAIL: the goal run's hop was not recorded"; exit 1; fi
echo "PASS: a contradiction halts every run: 2 at depth 0, 3 with a goal"
