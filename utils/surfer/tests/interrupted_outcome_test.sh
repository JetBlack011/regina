#!/bin/sh
# interrupted_outcome_test.sh <verifyslicegenus>
#
# search_outcome says why a row's search stopped, and "exhausted" means it ran
# out of candidates on its own. A search stopped by SIGINT did not, so its
# row must not be written as exhausted: that value is what licenses reading
# a negative (paper prop:exhausted), every merge keeps the maximum, and
# without --research-settled a row recorded exhausted at a face cap is never
# searched again at that cap.
#
# One SIGINT inside the search stops only that row's search (runSearch_'s
# SigintScope); the driver then writes the row and goes on. So: start a row
# whose cap-5 search takes minutes, interrupt it after a few seconds, and read
# its search_outcome.
set -eu

V=$1
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

cat > "$T/rows.csv" <<'EOF'
Name,PD Notation,Genus-4D
6_1,[[1;7;2;6];[3;10;4;11];[5;3;6;2];[7;1;8;12];[9;4;10;5];[11;9;12;8]],0
EOF
printf 'kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices\n' \
  > "$T/cobordisms.csv"

"$V" --input "$T/rows.csv" --output "$T/out.csv" \
     --cobordisms "$T/cobordisms.csv" --census-db "$T/none.sqlite" \
     --knot-table "$T/rows.csv" \
     --no-census-updates --no-retriangulate-on-miss \
     --thicken-layers 2 --collar-layers 2 --max-faces 5 \
     --no-cone --harvest --boundary-condition proper --research-settled \
     --threads 2 > "$T/log" 2> "$T/err" &
P=$!

# stdout is block-buffered into the log, so wait a fixed time rather than
# for a line: building 6_1's thickening takes well under a second, and its
# cap-5 search takes minutes.
sleep 10
if ! kill -INT "$P" 2>/dev/null; then
  echo "SKIP: the run ended before it could be interrupted"; exit 77
fi
wait "$P" || true

if grep -q 'EXHAUSTIVE to' "$T/log"; then
  echo "SKIP: the search finished before the signal arrived"; exit 77
fi

outcome=$(awk -F, '
  NR == 1 { for (i = 1; i <= NF; i++) if ($i == "search_outcome") c = i; next }
  $1 == "6_1" { print $c }' "$T/out.csv")
echo "search_outcome after SIGINT: '$outcome'"
grep -E '^\[\+\] 6_1: [0-9]+ new witnesses' "$T/log" || true

if [ -z "$outcome" ]; then
  echo "FAIL: no output row for 6_1"; exit 1
fi
if [ "$outcome" = "exhausted" ]; then
  echo "FAIL: an interrupted search was recorded as exhausted"; exit 1
fi
echo "PASS: an interrupted search is not recorded as exhausted"
