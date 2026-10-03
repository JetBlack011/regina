#!/bin/sh
# goal_layout_test.sh <cobound>
#
# A goal run must explore the same links whatever the heap layout. Until
# phase 3.0 it did not: boundary edges were collected in unordered sets of
# pointers, so which edge started each curve, the order of a link's
# components and the collar seed's order followed the addresses the
# allocator happened to hand out. Those decide how each outgoing link is
# drawn, and so which links a goal run interns. The same run explored 9 or
# 10 links on 6_2 depending only on the length of its --work path.
#
# So: run 6_2's goal run (G-regression's shape: exhaustive at cap 3, four
# searches) from work directories of five path lengths, which shifts every
# later allocation, and require one normalised node_bounds.jsonl (tabulated
# links by name, untabulated ones by their bounds, as G-regression compares
# them). The run names far sides from its diagrams, which needs the atlas's
# full tables: SURFER_TEST_ATLAS_DATA, or the usual checkouts; without them
# the test is skipped.
set -eu

C=$1
DATA=${SURFER_TEST_ATLAS_DATA:-}
for d in "$HOME/Projects/cobordism-atlas/data" "$HOME/Projects/triangles/cobordism-atlas/data"; do
  [ -n "$DATA" ] || { [ -f "$d/knot_symmetry.csv" ] && DATA=$d; } || true
done
KNOTS=$DATA/4d_smooth_slice_genus_13_crossings_pd_codes.csv
LINKS=$DATA/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv
if [ -z "$DATA" ] || [ ! -f "$KNOTS" ] || [ ! -f "$LINKS" ]; then
  echo "SKIP: no atlas tables (set SURFER_TEST_ATLAS_DATA)"; exit 77
fi

T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
PD='[[1;8;2;9];[3;11;4;10];[5;1;6;12];[7;2;8;3];[9;7;10;6];[11;5;12;4]]'
X=xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

n=0
for pad in "" x "$X" "$X$X" "$X$X$X$X"; do
  n=$((n + 1))
  dir=$T/w$pad
  mkdir -p "$dir"
  cat > "$dir/run.conf" <<CONF
target_pd = $PD
target_name = 6_2
work = $dir/work
knot_table = $KNOTS
link_table = $LINKS
knot_symmetry = $DATA/knot_symmetry.csv
census = $dir/none.sqlite
goal_genus = 1
literature = 0
threads = 4
max_searches = 4
surface_target = 1000000000
max_surface_target = 1000000000
max_faces = 3
iddfs_iterations = 0
iddfs_start = 0
iddfs_step = 0
root_budget_start = 0
resolve_unlinked = 1
CONF
  "$C" run --config "$dir/run.conf" > "$dir/stdout.txt" 2>&1 || true
  if [ ! -s "$dir/work/node_bounds.jsonl" ]; then
    echo "FAIL: layout $n wrote no node_bounds.jsonl"; tail -5 "$dir/stdout.txt"; exit 1
  fi
  python3 - "$dir/work/node_bounds.jsonl" > "$T/bounds.$n" <<'EOF'
import json, sys
out = []
for line in open(sys.argv[1]):
    n = json.loads(line)
    ident = n.get('table') or 'untabulated'
    up = sorted((u['partition'], u['genus'], u.get('constructive')) for u in n.get('upper', []))
    lo = sorted(json.dumps({k: v for k, v in l.items() if k not in ('record', 'id')}, sort_keys=True)
                for l in n.get('lower', []))
    out.append(f"{ident}\t{n['components']}\tupper {up}\tlower {lo}")
print('\n'.join(sorted(out)))
EOF
  echo "layout $n: work path of $(printf %s "$dir/work" | wc -c) characters," \
       "$(wc -l < "$T/bounds.$n") links, $(grep -c '^\[+\] hop [0-9]* .*: accounting:' "$dir/stdout.txt") searches"
done

if [ "$(grep -c '^\[+\] hop [0-9]* .*: accounting:' "$T/w/stdout.txt")" -lt 2 ]; then
  echo "FAIL: the goal run made fewer than two searches, so it tests no exploration"; exit 1
fi
i=2
while [ $i -le $n ]; do
  if ! cmp -s "$T/bounds.1" "$T/bounds.$i"; then
    echo "FAIL: layout $i explored different links from layout 1"
    diff "$T/bounds.1" "$T/bounds.$i" | head -10
    exit 1
  fi
  i=$((i + 1))
done
echo "PASS: one set of links and bounds under $n heap layouts"
