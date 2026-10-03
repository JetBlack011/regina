#!/usr/bin/env bash
# compare_surface_sets.sh <binaryA> <binaryB> [row...]
#
# Do two builds accept, and describe, exactly the same surfaces? Each binary
# is cobound (by the name it resolves to: `cobound run` with a config) or a
# retired verifyslicegenus (its options), so a build can be compared with the
# reference. Each row is searched EXHAUSTIVELY at a small face cap by both,
# with a surface log: one line per described surface -- orientable, genus,
# tubed genus, punctures, triangle count and pair signature. A pair
# signature is canonical for the (ambient, surface) pair, so the sorted logs
# are identical exactly when the two builds found the same surfaces and
# described each the same way, whatever order they found them in. The
# accounting lines (less the recorded/duplicate split, which follows
# non-canonical far-side names; see canaries.sh) must match too.
#
# This is the equivalence check for a change that reorders the search
# (the traversal rework) as much as for one that should change nothing.
#
# Rows default to canaries.sh's. Environment:
#   CAP     face cap (default 3, as canaries.sh)
#   THREADS default 8
#   BUDGET  --root-budget-start (default 0: unbudgeted)
#   PENDING_CAP  --pending-surface-cap (default: the binary's own). Set it
#           small to make the search pause and drain its queue many times
#           in a row that would otherwise never reach the cap.
#   RESOLVE_FLAG  how a verifyslicegenus is told resolve_unlinked, which has
#           no default (default --no-resolve-unlinked; set it empty to compare
#           builds from before that, which take no such option); cobound is
#           told resolve_unlinked = 1 exactly when it is --resolve-unlinked
#   ATLAS   cobordism-atlas checkout (default: as bench_search.sh)
#   WORK    scratch directory (default: a mktemp directory, removed on
#           success and kept on failure)
set -euo pipefail

A=${1:?usage: compare_surface_sets.sh <binaryA> <binaryB> [row...]}
B=${2:?usage: compare_surface_sets.sh <binaryA> <binaryB> [row...]}
shift 2
ROWS=("$@")
[ ${#ROWS[@]} -gt 0 ] || ROWS=(3_1 6_1 8_8 8_20 'L2a1{0}' 'L2a1{1}' 'L6a3{0}' 'L6a3{1}' 'L6a4{0;0}')
CAP=${CAP:-3}
THREADS=${THREADS:-8}
BUDGET=${BUDGET:-0}
CAPFLAG=()
[ -n "${PENDING_CAP:-}" ] && CAPFLAG=(--pending-surface-cap "$PENDING_CAP")
RESOLVE_FLAG=${RESOLVE_FLAG---no-resolve-unlinked}
[ -n "$RESOLVE_FLAG" ] && CAPFLAG+=("$RESOLVE_FLAG")
if [ -z "${ATLAS:-}" ]; then
  for d in "$HOME/Projects/cobordism-atlas" "$HOME/Projects/triangles/cobordism-atlas"; do
    [ -d "$d/data" ] && { ATLAS=$d; break; }
  done
fi
DATA=${ATLAS:?set ATLAS}/data
KEEP=${WORK:+1}
WORK=${WORK:-$(mktemp -d)}
mkdir -p "$WORK"

# Whether $1 is cobound (a symbolic link to it counts).
is_cobound() { [ "$(basename "$(readlink -f "$1")")" = cobound ]; }

run() {
  local bin=$1 tag=$2 out=$WORK/$2
  mkdir -p "$out"
  {
    echo "Name,PD Notation,Genus-4D"
    for r in "${ROWS[@]}"; do
      awk -F, -v r="$r" '$1 == r { print; found = 1; exit } END { exit !found }' \
          "$DATA/4d_smooth_slice_genus_13_crossings_pd_codes.csv" ||
      awk -F, -v r="$r" '$1 == r { print; found = 1; exit } END { exit !found }' \
          "$DATA/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv" ||
      { echo "compare_surface_sets.sh: no row $r" >&2; return 2; }
    done
  } > "$out/rows.csv"
  printf 'kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices\n' \
      > "$out/cobordisms.csv"
  # One process per row, so each row's surface log is its own file.
  local i=0
  for r in "${ROWS[@]}"; do
    i=$((i + 1))
    { head -1 "$out/rows.csv"; sed -n "$((i + 1))p" "$out/rows.csv"; } > "$out/row$i.csv"
    if is_cobound "$bin"; then
      cat > "$out/run$i.conf" <<CONF
targets = $out/row$i.csv
verdicts = $out/out$i.csv
cobordisms = $out/cobordisms.csv
work = $out/work
census = $out/none.sqlite
census_updates = 0
retriangulate_on_miss = 0
knot_table = $DATA/4d_smooth_slice_genus_13_crossings_pd_codes.csv
link_table = $DATA/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv
layers = 2
max_faces = $CAP
root_budget_start = $BUDGET
root_budget_growth = 2
boundary_condition = proper
surface_log = $out/surfaces$i.csv
threads = $THREADS
resolve_unlinked = $([ "$RESOLVE_FLAG" = --resolve-unlinked ] && echo 1 || echo 0)
CONF
      [ -n "${PENDING_CAP:-}" ] && echo "pending_surface_cap = $PENDING_CAP" >> "$out/run$i.conf"
      "$bin" run --config "$out/run$i.conf" > "$out/log$i" 2> "$out/err$i" || {
        echo "compare_surface_sets.sh: $tag failed on $r; see $out/err$i" >&2; return 1; }
    else
    "$bin" --input "$out/row$i.csv" --output "$out/out$i.csv" \
        --cobordisms "$out/cobordisms.csv" \
        --census-db "$out/none.sqlite" --no-census-updates --no-retriangulate-on-miss \
        --knot-table "$DATA/4d_smooth_slice_genus_13_crossings_pd_codes.csv" \
        --link-table "$DATA/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv" \
        --thicken-layers 2 --collar-layers 2 --max-faces "$CAP" \
        --root-budget-start "$BUDGET" --root-budget-growth 2 \
        --no-cone --harvest --boundary-condition proper --research-settled \
        --surface-log "$out/surfaces$i.csv" "${CAPFLAG[@]}" \
        --threads "$THREADS" > "$out/log$i" 2> "$out/err$i" || {
      echo "compare_surface_sets.sh: $tag failed on $r; see $out/err$i" >&2; return 1; }
    fi
    tail -n +2 "$out/surfaces$i.csv" | sort > "$out/surfaces$i.sorted"
    grep ': accounting:' "$out/log$i" |
      sed -E 's/recorded ([0-9]+), duplicate ([0-9]+)/examined \1+\2/' |
      awk '{ for (k = 1; k <= NF; k++) if ($k == "examined") { split($(k+1), p, "+"); $(k+1) = (p[1] + p[2]) "," } print }' \
      > "$out/acct$i"
  done
}

run "$A" A
run "$B" B

fail=0
i=0
for r in "${ROWS[@]}"; do
  i=$((i + 1))
  n=$(wc -l < "$WORK/A/surfaces$i.sorted")
  if ! cmp -s "$WORK/A/surfaces$i.sorted" "$WORK/B/surfaces$i.sorted"; then
    echo "DIFFER  $r: surface sets differ ($n vs $(wc -l < "$WORK/B/surfaces$i.sorted") lines)"
    diff "$WORK/A/surfaces$i.sorted" "$WORK/B/surfaces$i.sorted" | head -5
    fail=1
  elif ! cmp -s "$WORK/A/acct$i" "$WORK/B/acct$i"; then
    echo "DIFFER  $r: accounting differs"
    diff "$WORK/A/acct$i" "$WORK/B/acct$i"
    fail=1
  else
    echo "same    $r: $n surfaces, identical sets and accounting"
  fi
done

if [ "$fail" -eq 0 ]; then
  echo "compare_surface_sets.sh: IDENTICAL on ${#ROWS[@]} rows (cap $CAP, budget $BUDGET${PENDING_CAP:+, pending cap $PENDING_CAP})"
  [ -n "$KEEP" ] || rm -rf "$WORK"
else
  echo "compare_surface_sets.sh: DIFFERENT; kept $WORK"
  exit 1
fi
