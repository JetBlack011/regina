#!/bin/sh
# frontier_pending_test.sh <verifyslicegenus> <cascadesearch> <test data dir>
#
# One frontier rule (plan divergence 1). A search's frontier records its
# pending file and that file's fsynced byte length; `sign` records how far it
# has signed each pending file (<pending>.signed). A later search that would
# skip the frontier's prefix first requires that file signed at least that
# far -- else the prefix's cobordisms would live only in a pending file nobody
# signs -- unless the file is in its own run directory, which its own sign
# step signs. Refused, it says why and searches from the start.
#
#   1. run A (work WA) searches 3_1 to a surface target and signs at its end;
#      its frontier names WA's pending file and length;
#   2. A's signed record is removed: as if A had been killed before signing;
#   3. run B (work WB) refuses A's frontier, naming the reason;
#   4. `cascadesearch --sign-only --work WA` signs WA: the record is back;
#   5. run C (work WC) resumes A's frontier;
#   6. with the record removed again, run D in A's own work directory WA
#      resumes it (its own sign step signs WA's files).
set -eu

V=$1
C=$2
DATA=$3
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
head -2 "$DATA/knots_to_6.csv" > "$T/rows.csv"   # 3_1
echo "$HDR" > "$T/cobordisms.csv"
run() { # tag work target [resume]
  r=""; [ -n "${4:-}" ] && r="--resume-frontier-dir $T/fr"
  "$V" --input "$T/rows.csv" --output "$T/out.csv" --cobordisms "$T/cobordisms.csv" \
       --census-db "$T/none.sqlite" --knot-table "$DATA/knots_to_6.csv" \
       --link-table "$DATA/links_to_6.csv" --no-census-updates --no-retriangulate-on-miss \
       --no-resolve-unlinked --max-faces 4 --threads 1 --surface-target "$3" \
       --work "$2" --frontier-dir "$T/fr$1" $r > "$T/$1.log" 2> "$T/$1.err"
  grep -oE 'breadth: .*' "$T/$1.log" | sed -E 's/fingerprint [0-9a-f]+/fingerprint <f>/'
}

run A "$T/WA" 300
mkdir -p "$T/fr"
cp "$T/frA/3_1.frontier" "$T/fr/3_1.frontier" || { echo "FAIL: run A wrote no frontier"; exit 1; }
pending=$(sed -n 's/^pending \([0-9]*\) \(.*\)$/\1 \2/p' "$T/fr/3_1.frontier")
[ -n "$pending" ] || { echo "FAIL: the frontier records no pending file"; exit 1; }
bytes=${pending%% *}; file=${pending#* }
echo "frontier: pending $bytes bytes of $file"
[ "$bytes" -gt 0 ] || { echo "FAIL: nothing pending was recorded"; exit 1; }
[ "$(cat "$file.signed")" -ge "$bytes" ] || { echo "FAIL: run A did not sign its pending file"; exit 1; }

rm "$file.signed"
run B "$T/WB" 600 resume
grep -q "resumed no: its pending file $file is signed through 0 of the $bytes bytes" "$T/B.log" ||
  { echo "FAIL: run B resumed a frontier whose pending file is not signed"; exit 1; }

"$C" --sign-only --work "$T/WA" --witness-store "$T/cobordisms.csv" \
     --knot-table "$DATA/knots_to_6.csv" --link-table "$DATA/links_to_6.csv" > "$T/sign.log" 2>&1
[ "$(cat "$file.signed")" -ge "$bytes" ] || { echo "FAIL: sign recorded no signed length"; exit 1; }

run C "$T/WC" 600 resume
grep -q 'resumed yes' "$T/C.log" || { echo "FAIL: run C refused a signed frontier"; exit 1; }

rm "$file.signed"
run D "$T/WA" 600 resume
grep -q 'resumed yes' "$T/D.log" || { echo "FAIL: run D refused its own run's pending file"; exit 1; }
echo "PASS: a frontier is resumed only once its pending file is signed (or is the run's own)"
