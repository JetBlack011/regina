#!/bin/sh
# frontier_pending_test.sh <cobound> <test data dir>
#
# One frontier rule (plan divergence 1). A search's frontier records its
# pending file and that file's fsynced byte length; `sign` records how far it
# has signed each pending file (<pending>.signed). A later search that would
# skip the frontier's prefix first requires that file signed at least that
# far -- else the prefix's cobordisms would live only in a pending file nobody
# signs -- unless the file is in its own run directory, which its own sign
# step signs. Refused, it says why and searches from the start.
#
# The frontier records the pending file RELATIVE to its own directory
# (phase 5), so a work tree that is packed, synced or copied keeps its resume.
#
#   1. run A (work WA) searches 3_1 to a surface target and signs at its end;
#      its frontier names WA's pending file and length;
#   2. A's signed record is removed: as if A had been killed before signing;
#   3. run B (work WB) refuses A's frontier, naming the reason;
#   4. `cobound sign` with work WA signs WA: the record is back;
#   5. run C (work WC) resumes A's frontier;
#   6. with the record removed again, run D in A's own work directory WA
#      resumes it (its own sign step signs WA's files);
#   7. the whole tree is moved: run E (work WE, in the moved tree) resumes
#      A's frontier, whose pending file moved with it and is signed;
#   8. with the moved record removed, run F in the moved WA resumes it.
set -eu

C=$1
DATA=$2
R=$(mktemp -d)
trap 'rm -rf "$R"' EXIT
T=$R/tree
mkdir -p "$T"
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
head -2 "$DATA/knots_to_6.csv" > "$T/rows.csv"   # 3_1
echo "$HDR" > "$T/cobordisms.csv"
run() { # tag work target [resume]
  r=""; [ -n "${4:-}" ] && r="resume_frontier_dir=$T/fr"
  cat > "$T/$1.conf" <<CONF
targets = $T/rows.csv
verdicts = $T/out.csv
cobordisms = $T/cobordisms.csv
census = $T/none.sqlite
knot_table = $DATA/knots_to_6.csv
link_table = $DATA/links_to_6.csv
census_updates = 0
retriangulate_on_miss = 0
resolve_unlinked = 0
max_faces = 4
threads = 1
surface_target = $3
work = $2
frontier_dir = $T/fr$1
CONF
  [ -n "$r" ] && echo "$r" | sed 's/=/ = /' >> "$T/$1.conf"
  "$C" run --config "$T/$1.conf" > "$T/$1.log" 2> "$T/$1.err"
  grep -oE 'breadth: .*' "$T/$1.log" | sed -E 's/fingerprint [0-9a-f]+/fingerprint <f>/'
}
# The pending file the frontier in $T/fr names, absolute (it is recorded
# relative to the frontier's directory).
pending_file() {
  rel=$(sed -n 's/^pending [0-9]* \(.*\)$/\1/p' "$T/fr/3_1.frontier")
  case $rel in /*) echo "$rel" ;; *) echo "$(cd "$T/fr/$(dirname "$rel")" && pwd)/$(basename "$rel")" ;; esac
}

run A "$T/WA" 300
mkdir -p "$T/fr"
cp "$T/frA/3_1.frontier" "$T/fr/3_1.frontier" || { echo "FAIL: run A wrote no frontier"; exit 1; }
bytes=$(sed -n 's/^pending \([0-9]*\) .*$/\1/p' "$T/fr/3_1.frontier")
rel=$(sed -n 's/^pending [0-9]* \(.*\)$/\1/p' "$T/fr/3_1.frontier")
[ -n "$bytes" ] || { echo "FAIL: the frontier records no pending file"; exit 1; }
case $rel in /*) echo "FAIL: the pending file is recorded absolute ($rel)"; exit 1 ;; esac
file=$(pending_file)
echo "frontier: pending $bytes bytes of $rel ($file)"
[ "$bytes" -gt 0 ] || { echo "FAIL: nothing pending was recorded"; exit 1; }
[ "$(cat "$file.signed")" -ge "$bytes" ] || { echo "FAIL: run A did not sign its pending file"; exit 1; }

rm "$file.signed"
run B "$T/WB" 600 resume
grep -q "resumed no: its pending file $file is signed through 0 of the $bytes bytes" "$T/B.log" ||
  { echo "FAIL: run B resumed a frontier whose pending file is not signed"; exit 1; }

"$C" sign --set "work=$T/WA" --set "cobordisms=$T/cobordisms.csv" \
     --set "knot_table=$DATA/knots_to_6.csv" --set "link_table=$DATA/links_to_6.csv" > "$T/sign.log" 2>&1
[ "$(cat "$file.signed")" -ge "$bytes" ] || { echo "FAIL: sign recorded no signed length"; exit 1; }

run C "$T/WC" 600 resume
grep -q 'resumed yes' "$T/C.log" || { echo "FAIL: run C refused a signed frontier"; exit 1; }

rm "$file.signed"
run D "$T/WA" 600 resume
grep -q 'resumed yes' "$T/D.log" || { echo "FAIL: run D refused its own run's pending file"; exit 1; }
[ "$(cat "$file.signed")" -ge "$bytes" ] || { echo "FAIL: run D did not sign its own pending files"; exit 1; }

# The tree moves (packed and unpacked elsewhere, say): the frontier's record
# follows its pending file.
mv "$T" "$R/moved"
T=$R/moved
file=$(pending_file)
[ -f "$file" ] || { echo "FAIL: the moved frontier names no existing file ($file)"; exit 1; }
run E "$T/WE" 900 resume
grep -q 'resumed yes' "$T/E.log" ||
  { echo "FAIL: run E refused the moved tree's signed frontier"; grep -o 'resumed .*' "$T/E.log"; exit 1; }
rm "$file.signed"
run F "$T/WA" 900 resume
grep -q 'resumed yes' "$T/F.log" ||
  { echo "FAIL: run F refused its own (moved) run's pending file"; grep -o 'resumed .*' "$T/F.log"; exit 1; }
echo "PASS: a frontier is resumed only once its pending file is signed (or is the run's own), and a moved tree keeps its resume"
