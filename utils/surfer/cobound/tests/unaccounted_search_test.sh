#!/bin/sh
# unaccounted_search_test.sh <verifyslicegenus> <cascadesearch> <test data dir>
#
# An accounting imbalance ends only that search (plan divergence 2): its
# outcome is `unaccounted`, with no frontier and no exhaustion claim, and the
# run goes on. The exit code reports what the run was for:
#   - without a goal (depth 0) completeness is the product, so the run exits 2,
#     after searching every other row;
#   - with a goal the cascade's codes stay (0 met, 1 not met): the imbalance is
#     completeness only, and the hop is marked suspect on its `[!!]` line.
# SURFER_TEST_UNACCOUNTED=<search name> is the fixture: that search drops its
# first described surface from every bucket, as a surface lost between the
# drain and the record would be.
# An impossible state is different: it halts the process (exit 2), after
# writing what was found, at depth 0 and with a goal, even a met one.
# SURFER_TEST_IMPOSSIBLE=<search name> counts one surface in an impossible
# bucket.
set -eu

V=$1
C=$2
DATA=$3
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
head -3 "$DATA/knots_to_6.csv" > "$T/rows.csv"   # 3_1, 4_1
echo "$HDR" > "$T/cobordisms.csv"

# 1. Depth 0: 3_1 imbalanced, 4_1 searched as usual, exit 2.
rc=0
SURFER_TEST_UNACCOUNTED=3_1 "$V" --input "$T/rows.csv" --output "$T/out.csv" \
     --cobordisms "$T/cobordisms.csv" --census-db "$T/none.sqlite" \
     --knot-table "$DATA/knots_to_6.csv" --link-table "$DATA/links_to_6.csv" \
     --no-census-updates --no-retriangulate-on-miss --no-resolve-unlinked \
     --max-faces 2 --frontier-dir "$T/frontiers" --threads 2 > "$T/log" 2> "$T/err" || rc=$?
grep -E '^\[\+\] (3_1|4_1): [0-9]+ new witnesses|: accounting:|^\[!!\]|EXHAUSTIVE' "$T/log" || true
outcome() { awk -F, -v r="$1" '
  NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "search_outcome") c = i; if ($i == "exhausted_depth") x = i }; next }
  $1 == r { print $c " " $x }' "$T/out.csv"; }
if [ "$rc" -ne 2 ]; then echo "FAIL: depth 0 exited $rc, not 2"; exit 1; fi
if [ "$(outcome 3_1)" != "unaccounted -1" ]; then
  echo "FAIL: 3_1 recorded '$(outcome 3_1)', not 'unaccounted -1'"; exit 1; fi
if grep -q '3_1: EXHAUSTIVE' "$T/log" || [ -e "$T/frontiers/3_1.frontier" ]; then
  echo "FAIL: the imbalanced search claimed exhaustion or kept a frontier"; exit 1; fi
if [ "$(outcome 4_1)" != "exhausted 2" ] || [ ! -e "$T/frontiers/4_1.frontier" ]; then
  echo "FAIL: the run did not go on to 4_1 as usual ('$(outcome 4_1)')"; exit 1; fi
if ! grep -q '^\[!!\] 3_1: surface accounting failed' "$T/log"; then
  echo "FAIL: no [!!] line for the imbalanced search"; exit 1; fi

# 2. With a goal: the codes stay. 3_1, goal 1, not met at cap 3 (cascade
#    canaries): exit 1; L2a1{0}, goal 0, met by its annulus: exit 0.
goal() { # name pd goal table -> exit code
  rc=0
  SURFER_TEST_UNACCOUNTED=$1 "$C" --target-pd "$2" --target-name "$1" --work "$T/g$3" \
     --knot-table "$DATA/knots_to_6.csv" --link-table "$DATA/links_to_6.csv" \
     --census-db "$T/none.sqlite" --goal-genus "$3" --constructive --threads 1 \
     --max-expansions 1 --hop-surfaces 1000000000 --max-hop-surfaces 1000000000 \
     --hop-max-faces 3 --hop-iddfs-iterations 0 --hop-iddfs-start 0 --hop-iddfs-step 0 \
     --hop-root-budget 0 --resolve-unlinked > "$T/g$3.log" 2>&1 || rc=$?
  grep -E '^\[!!\] hop|outcome' "$T/g$3.log" >&2 || true
  if ! grep -q '^\[!!\] hop 0: surface accounting failed' "$T/g$3.log"; then
    echo "FAIL: $1's hop was not marked suspect" >&2; echo 99; return; fi
  if ! grep -q 'outcome unaccounted' "$T/g$3/hop_0_n0/log.txt"; then
    echo "FAIL: $1's hop outcome is not unaccounted" >&2; echo 99; return; fi
  echo "$rc"
}
pd31=$(sed -n 2p "$DATA/knots_to_6.csv" | cut -d, -f2)
pdhopf=$(grep -F 'L2a1{0},' "$DATA/links_to_6.csv" | cut -d, -f2)
r1=$(goal 3_1 "$pd31" 1)
r0=$(goal 'L2a1{0}' "$pdhopf" 0)
if [ "$r1" != 1 ] || [ "$r0" != 0 ]; then
  echo "FAIL: goal runs exited $r1 (want 1, not met) and $r0 (want 0, met)"; exit 1; fi
# 3. An impossible state halts, after writing what was found: at depth 0
#    before 4_1 is searched, and with a goal even when the goal is met.
rc=0
rm -rf "$T/frontiers"; echo "$HDR" > "$T/cobordisms.csv"; rm -f "$T/out.csv"
SURFER_TEST_IMPOSSIBLE=3_1 "$V" --input "$T/rows.csv" --output "$T/out.csv" \
     --cobordisms "$T/cobordisms.csv" --census-db "$T/none.sqlite" \
     --knot-table "$DATA/knots_to_6.csv" --link-table "$DATA/links_to_6.csv" \
     --no-census-updates --no-retriangulate-on-miss --no-resolve-unlinked \
     --max-faces 2 --frontier-dir "$T/frontiers" --threads 2 > "$T/ilog" 2> "$T/ierr" || rc=$?
if [ "$rc" -ne 2 ]; then echo "FAIL: impossible state at depth 0 exited $rc, not 2"; exit 1; fi
if ! grep -q 'FATAL' "$T/ierr" || grep -q '^\[+\] 4_1:' "$T/ilog"; then
  echo "FAIL: the run did not halt at 3_1"; exit 1; fi
if [ "$(wc -l < "$T/cobordisms.csv")" -lt 2 ]; then
  echo "FAIL: 3_1's witnesses were not written before the halt"; exit 1; fi
rc=0
SURFER_TEST_IMPOSSIBLE='L2a1{0}' "$C" --target-pd "$pdhopf" --target-name 'L2a1{0}' --work "$T/gi" \
   --knot-table "$DATA/knots_to_6.csv" --link-table "$DATA/links_to_6.csv" \
   --census-db "$T/none.sqlite" --goal-genus 0 --constructive --threads 1 \
   --max-expansions 1 --hop-surfaces 1000000000 --max-hop-surfaces 1000000000 \
   --hop-max-faces 3 --hop-iddfs-iterations 0 --hop-iddfs-start 0 --hop-iddfs-step 0 \
   --hop-root-budget 0 --resolve-unlinked > "$T/gi.log" 2>&1 || rc=$?
grep -E '^\[!!\] HALT|outcome' "$T/gi.log" || true
if [ "$rc" -ne 2 ] || ! grep -q '^\[!!\] HALT: hop 0' "$T/gi.log"; then
  echo "FAIL: an impossible state in a goal run exited $rc, not 2 with a HALT line"; exit 1; fi

echo "PASS: an imbalanced search ends only itself (exit 2 at depth 0, the goal's code with a goal); an impossible state halts (exit 2)"
