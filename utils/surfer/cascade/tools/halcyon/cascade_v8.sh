#!/bin/bash
# On halcyon, detached, after cascade_v7.sh:
#   v8 = v7 + patchG_hopshape (--hop-max-faces, --hop-iddfs-start,
#        --hop-iddfs-iterations, --hop-root-budget; the defaults are the
#        campaign's shape, so Pool A should not change)
# Then the Pool B baseline: the 66 rows with a known (point) literature value
# the atlas has no constructive bound for, each run alone at a fixed budget
# of 1,200 CPU-s, constructive, master withheld. Each line carries the row's
# fixed_search tag: only 10 were ever searched by a binary without the
# row-drop bugs (c4), so the other 56 negatives are weak.
# Log: v8.log.
B=$HOME/cascade-bench
D=$HOME/Projects/cobordism-atlas/data
L=$B/v8.log
log() { echo "$*" >> "$L"; }
log "waiting for v7 $(date -Is)"
until grep -q 'v7 done' "$B/v7.log" 2>/dev/null; do sleep 20; done
if grep -q 'FAILED' "$B/v7.log"; then log "v7 failed; not starting"; exit 1; fi

eval "$(sed -n '/^poolA() {/,/^}/p; /^check() {/,/^}/p' "$B/cascade_v2v3.sh")"
V=$B/v3

start=$(date +%s)
if ( cd "$V" && git apply "$B/patchG_hopshape.patch" \
     && make -C build -j16 cascadesearch cascade_profile_test cascade_proofgraph_test \
          cascade_diagramiso_test cascade_nodes_test cascade_hopedges_test \
          cascade_hoprunner_test cascade_leaves_test > build.v8.log 2>&1 \
     && ( cd build/utils/surfer/cascade && ctest --output-on-failure ) > tests.v8.log 2>&1 ); then
  log "v8 built and tested in $(( $(date +%s) - start )) s $(date -Is): $(grep -h 'tests passed' "$V/tests.v8.log")"
else
  log "v8 FAILED $(date -Is)"
  exit 1
fi
poolA "$V" v8
log "v8 Pool A done $(date -Is)"
check "$V" v8 10_27 'L10a136{1;0}' 11a_355

# pool <csv> <tag> <cascadesearch flags...>: every row of a Pool B-style
# list, its result line tagged with status and fixed_search.
pool() {
  local csvf=$1 tag=$2; shift 2
  local S=$V/build/utils/surfer R=$B/results/$tag
  mkdir -p "$R"
  : > "$R.txt"
  tail -n +2 "$csvf" | while IFS=, read -r name g4 rows chain sweepcpu known pd status lit derived fixed; do
    local w=$R/$name
    rm -rf "$w"; mkdir -p "$w"
    local start=$(date +%s)
    "$S/cascade/cascadesearch" --target-pd "$pd" --target-name "$name" --work "$w" \
      --knot-table "$D/4d_smooth_slice_genus_13_crossings_pd_codes.csv" \
      --link-table "$D/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv" \
      --knot-symmetry "$D/knot_symmetry.csv" \
      --census-db "$B/census.sqlite" --goal-genus "$g4" --constructive \
      --threads 14 "$@" > "$w/driver.log" 2>&1 < /dev/null
    local hops cpu res
    hops=$(wc -l < "$w/cascade.jsonl" 2>/dev/null || echo 0)
    cpu=$(python3 -c "import json;print(round(sum(json.loads(l).get('cpu',0) for l in open('$w/cascade.jsonl'))))" 2>/dev/null)
    res=$(grep -E 'GOAL MET|done:|CONTRADICTION|refused|cascadesearch:|\[!!\]' "$w/driver.log" | tr '\n' ' ')
    echo "$name g4=$g4 status=$status fixed=$fixed hops=$hops cpu=$cpu wall=$(( $(date +%s) - start )) :: $res" >> "$R.txt"
  done
  echo "BATCH DONE" >> "$R.txt"
}

log "Pool B baseline starts $(date -Is)"
pool "$B/pool_B.csv" poolB_base --hop-surfaces 50000 --max-hop-surfaces 200000 \
  --max-expansions 1000 --cpu-budget 1200
log "Pool B baseline done $(date -Is)"
log "v8 done $(date -Is)"
