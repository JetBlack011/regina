#!/bin/sh
# given_diagram_test.sh <cobound> <test data dir>
#
# Each target is searched on the diagram it was given, certified in-binary
# (plan, phase 7.2): the incoming link drawn back from the triangulation T
# built from the PD must be that PD's diagram (search::certifyIncoming()).
#   1. A goal run searches a table target on its PD as written, not on the
#      diagram its graph simplified it to: so its T is every other run's.
#   2. A table row whose given diagram does not certify is refused as an
#      error and never searched on another diagram (that would silently
#      change T): at depth 0 recorded as a build failure, with a goal the run
#      ends with exit 2. The fixture is 3_1 with a nugatory kink (a fourth
#      crossing [5;7;6;6]): the drawer (todiagram.h) cannot draw a nugatory
#      crossing, so its T never redraws as its diagram.
#   3. Only an untabulated target under a goal may fall back to its
#      simplified diagram, and the run says so.
#   4. What the run records a target's searches under -- log.txt's first
#      line, kept.csv's row_pd (column 14), the database's .rows.csv and
#      certificate.json's row_pd -- is the given PD in the frozen formats'
#      spelling, [[a;b;c;d];...], labels and crossing order as given: a link
#      table row given as LinkInfo writes it (PD[X[4; 1; 3; 2]; ...]), and a
#      knot's PD given with spaces (as the knot table writes it from 11
#      crossings), are respelt, never recorded as written.
set -eu

C=$1
DATA=$2
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
fail() { echo "FAIL: $*"; exit 1; }
pd31=$(grep '^3_1,' "$DATA/knots_to_6.csv" | cut -d, -f2)
kinked='[[1;5;2;4];[3;1;4;8];[5;7;6;6];[7;3;8;2]]'

goal() { # dir name pd [config line...] -> runs; rc in $rc (later lines win)
  mkdir -p "$1"
  cat > "$1.conf" <<CONF
target_pd = $3
target_name = $2
work = $1
knot_table = $DATA/knots_to_6.csv
link_table = $DATA/links_to_6.csv
census = $1.none.sqlite
goal_genus = 1
literature = 0
threads = 2
max_searches = 1
surface_target = 1000000000
max_surface_target = 1000000000
max_faces = 2
iddfs_iterations = 0
iddfs_start = 0
iddfs_step = 0
root_budget_start = 0
resolve_unlinked = 1
CONF
  d=$1; shift 3
  for line in "$@"; do echo "$line" >> "$d.conf"; done
  rc=0; "$C" run --config "$d.conf" > "$d.log" 2>&1 || rc=$?
}

# 1. A table target, searched on its PD as written.
goal "$T/table" 3_1 "$pd31"
[ "$rc" -eq 1 ] || fail "3_1 goal 1 at cap 2 exited $rc, not 1 (not met)"
[ "$(head -1 "$T/table/hop_0_n0/log.txt")" = "[+] 3_1 $pd31" ] ||
  fail "3_1 was not searched on its given PD: '$(head -1 "$T/table/hop_0_n0/log.txt")'"

# 2a. A table row that does not certify, with a goal: an error, exit 2, nothing searched.
goal "$T/refused" 3_1 "$kinked"
grep -E 'certif|refused' "$T/refused.log" || true
[ "$rc" -eq 2 ] || fail "a table target that does not certify exited $rc, not 2"
grep -q 'is a table row whose given diagram does not certify' "$T/refused.log" ||
  fail "no refusal message for the table target"
[ ! -e "$T/refused/hop_0_n0" ] || fail "the refused table target was searched"

# 2b. The same at depth 0: refused and recorded as a build failure.
mkdir -p "$T/d0"; echo "$HDR" > "$T/d0/cobordisms.csv"
printf 'Name,PD Notation,Genus-4D\n3_1,%s,1\n' "$kinked" > "$T/d0/rows.csv"
cat > "$T/d0/run.conf" <<CONF
targets = $T/d0/rows.csv
verdicts = $T/d0/out.csv
cobordisms = $T/d0/cobordisms.csv
work = $T/d0/work
census = $T/d0/none.sqlite
knot_table = $DATA/knots_to_6.csv
link_table = $DATA/links_to_6.csv
census_updates = 0
retriangulate_on_miss = 0
resolve_unlinked = 0
max_faces = 2
threads = 2
CONF
rc=0; "$C" run --config "$T/d0/run.conf" > "$T/d0/log" 2> "$T/d0/err" || rc=$?
grep -E 'failed to build' "$T/d0/err" || true
outcome=$(awk -F, 'NR == 1 { for (i = 1; i <= NF; i++) if ($i == "search_outcome") c = i; next }
                   $1 == "3_1" { print $c }' "$T/d0/out.csv")
[ "$outcome" = build-failed ] || fail "depth 0 recorded '$outcome' for a row that does not certify"
grep -q '3_1: failed to build' "$T/d0/err" || fail "no build-failure line at depth 0"
! grep -q ': accounting:' "$T/d0/log" || fail "the uncertified row was searched at depth 0"

# 3. An untabulated target that does not certify: its simplified diagram, said so.
goal "$T/fallback" kinked_trefoil "$kinked"
grep -E 'certif' "$T/fallback.log" || true
[ "$rc" -eq 1 ] || fail "the untabulated fallback exited $rc, not 1 (searched, not met)"
grep -q '^\[!\] target kinked_trefoil: its given diagram does not certify' "$T/fallback.log" ||
  fail "the fallback was not logged"
grep -q '^\[+\] hop 0 .*: accounting:' "$T/fallback.log" || fail "the fallback searched nothing"

# 4a. A link table row given as LinkInfo writes it: recorded in the bracket spelling.
pdL=$(grep '^L2a1{0},' "$DATA/links_to_6.csv" | cut -d, -f2)
spelt=$(printf '%s' "$pdL" | sed 's/^PD\[/[/; s/X\[/[/g; s/ //g')
case $pdL in PD\[X*) ;; *) fail "the fixture L2a1{0} is not in LinkInfo's spelling: $pdL" ;; esac
goal "$T/link" 'L2a1{0}' "$pdL" "max_faces = 3" "cobordisms = $T/link.store.csv" \
  "run_name = given_diagram_test/L2a1"
[ "$rc" -eq 0 ] || fail "L2a1{0} goal 1 at cap 3 exited $rc, not 0 (met)"
[ "$(head -1 "$T/link/hop_0_n0/log.txt")" = "[+] L2a1{0} $spelt" ] ||
  fail "L2a1{0}'s log.txt does not record its given PD respelt: '$(head -1 "$T/link/hop_0_n0/log.txt")'"
k=$(wc -l < "$T/link/hop_0_n0/kept.csv")
[ "$k" -gt 0 ] || fail "L2a1{0} kept nothing"
[ "$(grep -cF ",$spelt," "$T/link/hop_0_n0/kept.csv")" -eq "$k" ] ||
  fail "a kept.csv line of L2a1{0} does not record row_pd $spelt"
grep -qF "\"row_pd\":\"$spelt\"" "$T/link/certificate.json" || fail "certificate.json's row_pd is not $spelt"
grep -qF ",$spelt" "$T/link.store.csv.rows.csv" || fail "the .rows.csv sidecar does not record $spelt"
for f in "$T/link/hop_0_n0/log.txt" "$T/link/hop_0_n0/kept.csv" "$T/link.store.csv.rows.csv"; do
  ! grep -q 'PD\[X' "$f" || fail "$f records LinkInfo's spelling"
done
! grep -q '"row_pd":"PD\[X' "$T/link/certificate.json" || fail "certificate.json records LinkInfo's spelling"

# 4b. A knot's PD given with spaces: recorded without them.
spaced=$(printf '%s' "$pd31" | sed 's/;/; /g')
goal "$T/spaced" 3_1 "$spaced"
[ "$rc" -eq 1 ] || fail "3_1 (spaced PD) goal 1 at cap 2 exited $rc, not 1 (not met)"
[ "$(head -1 "$T/spaced/hop_0_n0/log.txt")" = "[+] 3_1 $pd31" ] ||
  fail "3_1's spaced PD was not respelt: '$(head -1 "$T/spaced/hop_0_n0/log.txt")'"

echo "PASS: targets are searched on their given diagrams, certified; a table row that does not certify is refused (exit 2 with a goal, build-failed at depth 0); an untabulated target falls back, logged; a given PD is recorded in the bracket spelling, labels and order as given"
