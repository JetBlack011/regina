#!/bin/bash
# On halcyon: one cascade target, constructive, master withheld, then its
# certificate replayed by cascade_check.py. Leaves Pool B alone.
# usage: one_target.sh <tag> <name> <goal genus> '<PD>' [cascadesearch flags...]
# Result line in results/<tag>.txt, the run in results/<tag>/<name>/.
B=$HOME/cascade-bench
D=$HOME/Projects/cobordism-atlas/data
V=$B/v3
S=$V/build/utils/surfer
tag=$1 name=$2 g4=$3 pd=$4
shift 4
w=$B/results/$tag/$name
rm -rf "$w"; mkdir -p "$w"
start=$(date +%s)
"$S/cascade/cascadesearch" --target-pd "$pd" --target-name "$name" --work "$w" \
  --knot-table "$D/4d_smooth_slice_genus_13_crossings_pd_codes.csv" \
  --link-table "$D/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv" \
  --knot-symmetry "$D/knot_symmetry.csv" \
  --census-db "$B/census.sqlite" --goal-genus "$g4" --constructive --threads 14 \
  --hop-surfaces 50000 --max-hop-surfaces 1000000 --max-expansions 100000 --cpu-budget 7200 \
  "$@" > "$w/driver.log" 2>&1 < /dev/null
res=$(grep -E 'GOAL MET|done:|CONTRADICTION|refused|cascadesearch:|\[!!\]' "$w/driver.log" | tr '\n' ' ')
echo "$name g4=$g4 wall=$(( $(date +%s) - start )) :: $res" >> "$B/results/$tag.txt"
if [ -f "$w/certificate.json" ]; then
  CASCADE_ATLAS=$HOME/Projects/cobordism-atlas "$HOME/.venvs/atlas/bin/python" \
    "$V/utils/surfer/cascade/tools/cascade_check.py" \
    --farsidediagram "$S/farsidediagram" "$w/certificate.json" > "$w/check.txt" 2>&1
  echo "check $name: $(tail -1 "$w/check.txt")" >> "$B/results/$tag.txt"
fi
echo "DONE $name" >> "$B/results/$tag.txt"
