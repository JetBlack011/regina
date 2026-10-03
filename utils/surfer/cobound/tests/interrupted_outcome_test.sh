#!/bin/sh
# interrupted_outcome_test.sh <verifyslicegenus> <cascadesearch> <test data dir>
#
# One process-level policy for SIGINT and SIGTERM (plan divergence 8): the
# first signal ends the running search cleanly -- its drain finishes, its
# pending cobordisms are fsynced, and it records `interrupted` -- and the run
# starts no further search; the exit code reports what the run was for: 0
# without a goal (the search is thinner, not failed), 1 with one (not met).
#
# search_outcome says why a row's search stopped, and "exhausted" means it ran
# out of candidates on its own. A search stopped by a signal did not, so its
# row must not be written as exhausted: that value is what licenses reading a
# negative (paper prop:exhausted), and every merge keeps the maximum.
#
# Each case starts a search whose cap-5 run takes minutes, signals it after a
# few seconds, and reads what it recorded: at depth 0 two rows (the second
# must never start), and with a goal one hop.
set -eu

V=$1
C=$2
DATA=$3
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
# 6_1 then 6_2: rows run in crossing order (stable), so 6_2 would be next.
{ head -1 "$DATA/knots_to_6.csv"; grep '^6_1,' "$DATA/knots_to_6.csv"; grep '^6_2,' "$DATA/knots_to_6.csv"; } > "$T/rows.csv"
# The goal case: 6_2 (4-genus 1) with goal 0, which no hop can meet.
pd62=$(grep '^6_2,' "$DATA/knots_to_6.csv" | cut -d, -f2)

# Signals $1 after 10 s (stdout is block-buffered into the log, so wait a
# fixed time rather than for a line: 6_1's thickening builds in well under a
# second, and its cap-5 search takes minutes). Sets rc to its exit code, or
# to 77 when it had already ended.
signal_after() { # pid signal
  sleep 10
  if ! kill "-$2" "$1" 2>/dev/null; then rc=77; return; fi
  rc=0; wait "$1" || rc=$?
}

depth0() { # signal
  d=$T/d0-$1; mkdir -p "$d"
  echo "$HDR" > "$d/cobordisms.csv"
  "$V" --input "$T/rows.csv" --output "$d/out.csv" --cobordisms "$d/cobordisms.csv" \
       --census-db "$d/none.sqlite" --knot-table "$DATA/knots_to_6.csv" \
       --link-table "$DATA/links_to_6.csv" --no-census-updates --no-retriangulate-on-miss \
       --max-faces 5 --no-resolve-unlinked --threads 2 > "$d/log" 2> "$d/err" &
  signal_after $! "$1"
  [ "$rc" = 77 ] && { echo "SKIP: the run ended before it could be signalled"; exit 77; }
  if grep -q '6_1: EXHAUSTIVE to' "$d/log"; then
    echo "SKIP: the search finished before the signal arrived"; exit 77; fi
  outcome=$(awk -F, '
    NR == 1 { for (i = 1; i <= NF; i++) if ($i == "search_outcome") c = i; next }
    $1 == "6_1" { print $c }' "$d/out.csv")
  echo "depth 0, $1: exit $rc, 6_1 search_outcome '$outcome'"
  [ "$rc" = 0 ] || { echo "FAIL: depth 0 exited $rc after $1, not 0"; exit 1; }
  [ "$outcome" = interrupted ] || { echo "FAIL: 6_1 recorded '$outcome', not interrupted"; exit 1; }
  if grep -q 'Searching 6_2' "$d/log" || ! grep -q "SIG$1: the run searches no further row" "$d/log"; then
    echo "FAIL: the run started another search after $1, or did not say it would not"; exit 1; fi
}

goal() { # signal
  d=$T/g-$1
  "$C" --target-pd "$pd62" --target-name 6_2 --work "$d" --knot-table "$DATA/knots_to_6.csv" \
       --link-table "$DATA/links_to_6.csv" --census-db "$T/none.sqlite" --goal-genus 0 \
       --constructive --threads 2 --max-expansions 4 --hop-surfaces 1000000000 \
       --max-hop-surfaces 1000000000 --hop-max-faces 5 --hop-iddfs-iterations 0 \
       --hop-iddfs-start 0 --hop-iddfs-step 0 --hop-root-budget 0 --resolve-unlinked \
       > "$T/g-$1.log" 2>&1 &
  signal_after $! "$1"
  [ "$rc" = 77 ] && { echo "SKIP: the goal run ended before it could be signalled"; exit 77; }
  hop=$(grep -h 'outcome' "$d"/hop_0_n0/log.txt 2>/dev/null || true)
  last=$(grep -oE 'outcome [a-z-]+$' "$T/g-$1.log" | tail -1)
  echo "goal, $1: exit $rc, hop: '$hop', run: '$last'"
  if grep -q 'GOAL MET' "$T/g-$1.log"; then
    echo "SKIP: the goal was met before the signal arrived"; exit 77; fi
  [ "$rc" = 1 ] || { echo "FAIL: the goal run exited $rc after $1, not 1"; exit 1; }
  case $hop in *"outcome interrupted"*) ;; *) echo "FAIL: the hop did not record interrupted"; exit 1 ;; esac
  [ "$last" = "outcome interrupted" ] || { echo "FAIL: the run's outcome is '$last'"; exit 1; }
  if grep -qE '^\[\+\] hop 1 ' "$T/g-$1.log"; then
    echo "FAIL: the goal run started another search after $1"; exit 1; fi
}

depth0 INT
depth0 TERM
goal INT
goal TERM
echo "PASS: SIGINT and SIGTERM end the running search cleanly, at depth 0 and with a goal"
