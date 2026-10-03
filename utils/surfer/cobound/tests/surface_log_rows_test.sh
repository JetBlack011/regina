#!/bin/sh
# surface_log_rows_test.sh <verifyslicegenus>
#
# --surface-log makes one CsvWriter per row. CsvWriter caches each thread's
# shard in a function-local thread_local pointer, shared by every CsvWriter
# in the process, so a thread that writes in one row and again in the next
# uses a shard of the previous, destroyed writer. Two rows on one thread
# must both be logged, and the run must finish.
set -eu

V=$1
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

cat > "$T/rows.csv" <<'EOF'
Name,PD Notation,Genus-4D
3_1,[[1;5;2;4];[3;1;4;6];[5;3;6;2]],1
4_1,[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]],1
EOF
# The tables the search's cobordism graph names by (no link is needed here).
printf 'Name,PD Notation,Genus-4D\n' > "$T/links.csv"
printf 'kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices\n' \
  > "$T/cobordisms.csv"

rc=0
"$V" --input "$T/rows.csv" --output "$T/out.csv" \
     --cobordisms "$T/cobordisms.csv" --census-db "$T/none.sqlite" \
     --knot-table "$T/rows.csv" --link-table "$T/links.csv" \
     --no-census-updates --no-retriangulate-on-miss \
     --thicken-layers 2 --collar-layers 2 --max-faces 2 \
     --no-cone --harvest --boundary-condition proper --research-settled \
     --no-resolve-unlinked \
     --surface-log "$T/surfaces.csv" --threads 1 \
     > "$T/log" 2> "$T/err" || rc=$?

if [ "$rc" -ne 0 ]; then
  echo "FAIL: verifyslicegenus exited $rc"
  grep -a 'terminate\|what()' "$T/err" || tail -3 "$T/err"
  exit 1
fi
rows=$(grep -c ': accounting:' "$T/log" || true)
if [ "$rows" -ne 2 ]; then
  echo "FAIL: expected two rows searched, got $rows"; exit 1
fi
echo "PASS: two rows logged on one thread"
