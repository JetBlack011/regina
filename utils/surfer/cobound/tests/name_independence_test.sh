#!/bin/sh
# name_independence_test.sh <verifyslicegenus>
#
# No gate of the search may depend on an identified name. Names are not
# canonical -- a census hit's "#N" varies between identifications of the
# same manifold -- and a gate that compared them once silently discarded
# whole rows (the D1 bug, 2026-09-26).
#
# So: search the same rows twice, exhaustively and deterministically, once
# with SURFER_TEST_PERTURB_NAMES set (every identified name gets a fresh
# suffix, so no two identifications ever agree). Which surfaces are accepted,
# and why each is or is not recorded, must not change. Only the split between
# "recorded" and "duplicate" may move, since dedup is BY name for far sides
# that bear a bound; their sum may not.
set -eu

V=$1
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

cat > "$T/rows.csv" <<'EOF'
Name,PD Notation,Genus-4D
6_1,[[1;7;2;6];[3;10;4;11];[5;3;6;2];[7;1;8;12];[9;4;10;5];[11;9;12;8]],0
L4a1{1},PD[X[6; 2; 7; 1]; X[8; 4; 5; 3]; X[2; 8; 3; 7]; X[4; 6; 1; 5]],1
EOF
# The same rows as the tables, so that far sides are named from their
# diagrams (farsidenaming.h) -- the production path -- and the perturbation
# reaches those names too.
head -2 "$T/rows.csv" > "$T/knots.csv"
{ head -1 "$T/rows.csv"; tail -1 "$T/rows.csv"; } > "$T/links.csv"

run() {
  dir=$1
  mkdir -p "$dir"
  printf 'kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices\n' \
    > "$dir/cobordisms.csv"
  "$V" --input "$T/rows.csv" --output "$dir/out.csv" \
       --cobordisms "$dir/cobordisms.csv" --census-db "$dir/none.sqlite" \
       --knot-table "$T/knots.csv" --link-table "$T/links.csv" \
       --no-census-updates --no-retriangulate-on-miss \
       --thicken-layers 2 --collar-layers 2 --max-faces 3 \
       --no-cone --harvest --boundary-condition proper --research-settled \
     --no-resolve-unlinked \
       --threads 2 > "$dir/log" 2> "$dir/err"
  # accepted, described, other-orientation, search-side-elsewhere,
  # impossible, drain, verdict -- and recorded+duplicate as one number.
  grep ': accounting:' "$dir/log" | sed -E \
    's/^\[\+\] (.*): accounting: accepted ([0-9]+), described ([0-9]+), recorded ([0-9]+), duplicate ([0-9]+), (.*)/\1 accepted=\2 described=\3 examined=\4+\5 \6/' |
    awk '{ split($4, rd, "="); split(rd[2], p, "+"); $4 = "examined=" (p[1] + p[2]); print }'
}

SURFER_TEST_PERTURB_NAMES=0 run "$T/plain" > "$T/plain.acct"
SURFER_TEST_PERTURB_NAMES=1 run "$T/perturbed" > "$T/perturbed.acct"

echo "--- plain"; cat "$T/plain.acct"
echo "--- perturbed"; cat "$T/perturbed.acct"

if [ ! -s "$T/plain.acct" ] || [ "$(wc -l < "$T/plain.acct")" -ne 2 ]; then
  echo "FAIL: expected one accounting line per row"; exit 1
fi
if ! grep -q 'examined=[1-9]' "$T/plain.acct"; then
  echo "FAIL: the plain run examined nothing, so the comparison is vacuous"
  exit 1
fi
if ! grep -qE ': diagram naming: [1-9][0-9]* far sides drawn' "$T/plain/log"; then
  echo "FAIL: no far side was named from its diagram, so the namer went untested"
  exit 1
fi
if ! cmp -s "$T/plain.acct" "$T/perturbed.acct"; then
  echo "FAIL: perturbing names changed which surfaces were accepted or why"
  exit 1
fi
if ! grep -q 'SURFER_TEST_PERTURB_NAMES' "$T/perturbed/err"; then
  echo "FAIL: the perturbed run did not take the test hook"; exit 1
fi
echo "PASS: acceptance is independent of every identified name"
