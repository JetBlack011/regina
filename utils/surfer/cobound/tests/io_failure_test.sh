#!/bin/sh
# io_failure_test.sh <cobound> <test data dir>
#
# A failed output write ends the search, never the process (plan, phase 7.3).
# Several writers run on the drain's and the judge's threads (the surface
# log's rows, rejection samples, the pending file's checkpoints), where a
# throw is std::terminate: phase 6 lost a whole run, SIGABRT and a core dump,
# when the disk filled under a surface log. Now the write failure is the
# search's: its outcome is `io-error`, it records no frontier and claims no
# exhaustion, what it kept reaches its pending file and the database, and the
# run reports it -- without a goal the run goes on to its other targets and
# exits 2; with a goal it halts (exit 2). Each case below must end that way,
# never in a signal:
#   1. a frontier that cannot be written (3_1's frontier path is a directory):
#      3_1 is an I/O error, 4_1 is searched as usual;
#   2. surface_stats and self_intersection_census at /dev/full (every write
#      fails, ENOSPC);
#   3. surface_log at /dev/full (its shard files beside it cannot be opened);
#   4. a disk that fills under the surface log: a file-size limit (ulimit -f,
#      SIGXFSZ ignored, so writes past it fail with EFBIG) far below the log
#      and far above every other file, so the shards' writes fail on the
#      drain threads -- phase 6's crash -- and nothing else does;
#   5. rejection_sample_log at /dev/full;
#   6. a goal run whose search directory refuses writes (its pending file),
#      and one whose frontier.txt cannot be written.
# And the writers after (or beside) the search, on the main thread (phase 7's
# review), each a halt or an io-error outcome, never a silent loss:
#   7. the judge thread's forced checkpoint at the first constructive find
#      (depth 0, L2a1{0}): the run's umask makes the search directory it creates
#      refuse files, so the pending file's first write -- that checkpoint, in a
#      run of seconds -- fails, right after the find is announced;
#   8. a goal met whose certificate cannot be written: exit 2, no GOAL MET
#      anywhere (stdout, stderr, the run directory's records);
#   9. a run record that cannot be written, at the run's end (node_bounds.jsonl:
#      outcome io-error) and mid-run (cascade.jsonl: a halt), the goal met all
#      the same: exit 2, no GOAL MET, no certificate;
#  10. the sign step at the run's end (a database that refuses appends), at
#      depth 0 (the search reported again with outcome io-error, its exhaustion
#      claim withdrawn) and with a goal (halt, no GOAL MET); and a verdicts
#      shard that cannot be written (depth 0: outcome io-error, no claim).
set -eu

C=$1
DATA=$2
T=$(mktemp -d)
trap 'chmod -R u+w "$T" 2>/dev/null; rm -rf "$T"' EXIT
HDR='kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices'
fail() { echo "FAIL: $*"; exit 1; }

conf() { # dir rows max_faces key=value...
  d=$1; rows=$2; faces=$3; shift 3
  mkdir -p "$d"; echo "$HDR" > "$d/cobordisms.csv"
  {
    echo "targets = $rows"; echo "verdicts = $d/out.csv"; echo "cobordisms = $d/cobordisms.csv"
    echo "work = $d/work"; echo "census = $d/none.sqlite"
    echo "knot_table = $DATA/knots_to_6.csv"; echo "link_table = $DATA/links_to_6.csv"
    echo "census_updates = 0"; echo "retriangulate_on_miss = 0"; echo "resolve_unlinked = 0"
    echo "max_faces = $faces"; echo "threads = 4"
    for kv in "$@"; do echo "${kv%%=*} = ${kv#*=}"; done
  } > "$d/run.conf"
}
outcome() { # dir row -> "search_outcome exhausted_depth"
  awk -F, -v r="$2" '
    NR == 1 { for (i = 1; i <= NF; i++) { if ($i == "search_outcome") c = i; if ($i == "exhausted_depth") x = i }; next }
    $1 == r { print $c " " $x }' "$1/out.csv"; }
# The run's exit code, which must be 2 and not a signal's (128 + n).
expect2() { # rc what
  [ "$1" -lt 128 ] || fail "$2: the process died on signal $(($1 - 128))"
  [ "$1" -eq 2 ] || fail "$2: exited $1, not 2"; }
signed() { [ "$(wc -l < "$1/cobordisms.csv")" -ge 2 ] || fail "$2: no cobordism reached the database"; }

{ head -1 "$DATA/knots_to_6.csv"; grep -E '^(3_1|4_1),' "$DATA/knots_to_6.csv"; } > "$T/knots.csv"
{ head -1 "$DATA/knots_to_6.csv"; grep '^6_1,' "$DATA/knots_to_6.csv"; } > "$T/6_1.csv"
{ head -1 "$DATA/links_to_6.csv"; grep -F 'L2a1{0},' "$DATA/links_to_6.csv"; } > "$T/hopf.csv"

# 1. A frontier that cannot be written.
d=$T/frontier; conf "$d" "$T/knots.csv" 2 "frontier_dir=$d/frontiers"
mkdir -p "$d/frontiers/3_1.frontier/x"
rc=0; "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err" || rc=$?
grep -E 'outcome|EXHAUSTIVE|^\[!!?\]' "$d/log" || true
expect2 $rc "an unwritable frontier"
[ "$(outcome "$d" 3_1)" = "io-error -1" ] || fail "3_1 recorded '$(outcome "$d" 3_1)', not 'io-error -1'"
! grep -q '3_1: EXHAUSTIVE' "$d/log" || fail "3_1 claimed exhaustion without its frontier"
grep -q '^\[!!\] 3_1: an output write failed -- frontier' "$d/log" || fail "no [!!] line for 3_1"
[ "$(outcome "$d" 4_1)" = "exhausted 2" ] && [ -f "$d/frontiers/4_1.frontier" ] ||
  fail "the run did not go on to 4_1 as usual ('$(outcome "$d" 4_1)')"
signed "$d" "an unwritable frontier"

# 2. Surface stats and the self-intersection census at /dev/full.
d=$T/stats; conf "$d" "$T/knots.csv" 2 surface_stats=/dev/full self_intersection_census=/dev/full
rc=0; "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err" || rc=$?
grep -E 'outcome|^\[!!?\]' "$d/log" || true
expect2 $rc "surface stats at /dev/full"
for r in 3_1 4_1; do
  [ "$(outcome "$d" $r)" = "io-error -1" ] || fail "$r recorded '$(outcome "$d" $r)' with stats at /dev/full"
done
grep -q '^\[!!\] 3_1: an output write failed -- surface_stats' "$d/log" || fail "no surface_stats failure line"
signed "$d" "surface stats at /dev/full"

# 3. The surface log at /dev/full.
d=$T/devlog; conf "$d" "$T/knots.csv" 2 surface_log=/dev/full
rc=0; "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err" || rc=$?
grep -E 'outcome|^\[!!?\]' "$d/log" || true
expect2 $rc "the surface log at /dev/full"
[ "$(outcome "$d" 3_1)" = "io-error -1" ] || fail "3_1 recorded '$(outcome "$d" 3_1)' with the log at /dev/full"
grep -q '^\[!!\] 3_1: an output write failed -- surface_log' "$d/log" || fail "no surface_log failure line"
signed "$d" "the surface log at /dev/full"

# 4. The disk fills under the surface log (6_1, cap 3: tens of MB of log; every
#    other file of the run is under 200 KB). ulimit -f counts 512-byte blocks in
#    dash and 1024-byte ones in bash: 4096 is 2-4 MB either way.
d=$T/full; conf "$d" "$T/6_1.csv" 3 "surface_log=$d/logs/surfaces.csv"
mkdir -p "$d/logs"
rc=0; (trap '' XFSZ; ulimit -f 4096; exec "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err") || rc=$?
grep -E 'outcome|^\[!!?\]' "$d/log" || true
expect2 $rc "a disk filling under the surface log"
[ "$(outcome "$d" 6_1)" = "io-error -1" ] || fail "6_1 recorded '$(outcome "$d" 6_1)' with a full disk"
grep -q '^\[!!\] 6_1: an output write failed -- surface_log' "$d/log" || fail "no surface_log failure line"
[ -z "$(ls "$d/logs")" ] || fail "the failed surface log left files behind: $(ls "$d/logs")"
signed "$d" "a disk filling under the surface log"

# 5. Rejection samples at /dev/full (L2a1{0} turns surfaces away on orientation).
d=$T/rejections; conf "$d" "$T/hopf.csv" 3 rejection_sample_log=/dev/full
rc=0; "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err" || rc=$?
grep -E 'outcome|^\[!!?\]' "$d/log" || true
expect2 $rc "rejection samples at /dev/full"
[ "$(outcome "$d" 'L2a1{0}')" = "io-error -1" ] || fail "L2a1{0} recorded '$(outcome "$d" 'L2a1{0}')'"
grep -q 'an output write failed -- rejection_sample_log' "$d/log" || fail "no rejection_sample_log failure line"

# 6. Goal runs: 3_1, goal 1, one exhaustive cap-3 search (the goal canaries'),
#    signing into a database.
goalconf() { # work
  mkdir -p "$1"; echo "$HDR" > "$1.db.csv"
  cat > "$1.conf" <<CONF
target_pd = $(grep '^3_1,' "$DATA/knots_to_6.csv" | cut -d, -f2)
target_name = 3_1
work = $1
cobordisms = $1.db.csv
run_name = io_failure/3_1
knot_table = $DATA/knots_to_6.csv
link_table = $DATA/links_to_6.csv
census = $1.none.sqlite
goal_genus = 1
literature = 0
threads = 4
max_searches = 1
surface_target = 1000000000
max_surface_target = 1000000000
max_faces = 3
iddfs_iterations = 0
iddfs_start = 0
iddfs_step = 0
root_budget_start = 0
resolve_unlinked = 1
CONF
}
#    a. The search directory refuses writes: no pending file.
g=$T/goal-pending; goalconf "$g"; mkdir -p "$g/hop_0_n0"; chmod a-w "$g/hop_0_n0"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
grep -E '^\[!!?\]|outcome' "$g.log" || true
expect2 $rc "a goal run's unwritable search directory"
grep -q '^\[!!\] HALT: hop 0: an output write failed -- pending file' "$g.log" || fail "no HALT line for the pending file"
#    b. frontier.txt cannot be written (a directory of that name).
g=$T/goal-frontier; goalconf "$g"; mkdir -p "$g/hop_0_n0/frontier.txt/x"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
grep -E '^\[!!?\]|outcome' "$g.log" || true
expect2 $rc "a goal run's unwritable frontier"
grep -q '^\[!!\] HALT: hop 0: an output write failed -- frontier' "$g.log" || fail "no HALT line for the frontier"
[ "$(wc -l < "$g.db.csv")" -ge 2 ] || fail "the halted goal run signed nothing"

# The last outcome line a run prints for a row (dispatch.py reads the last).
lastoutcome() { grep -E "^\[\+\] $2: [0-9]+ new witnesses, outcome " "$1" | tail -1 | sed 's/.*outcome \([-a-z]*\).*/\1/'; }
# No goal is reported met: not in the run's output, not in its directory.
notmet() { # log workdir what
  ! grep -q 'GOAL MET' "$1" || fail "$3: a GOAL MET line was printed"
  ! grep -rqsE 'GOAL MET|goal_met' "$2" || fail "$3: a record in $2 says the goal was met"
  [ ! -f "$2/certificate.json" ] || fail "$3: a certificate was written"; }

# 7. The judge thread's forced checkpoint fails (depth 0).
d=$T/forced; conf "$d" "$T/hopf.csv" 2; mkdir -p "$d/work"
rc=0; (umask 0222; exec "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err") || rc=$?
grep -E 'CONSTRUCTIVE|outcome|^\[!!?\]' "$d/log" || true
expect2 $rc "the forced checkpoint's failure"
grep -A1 'L2a1{0}: CONSTRUCTIVE witness found.*Checkpointing now' "$d/log" | tail -1 |
  grep -q '^\[!\] L2a1{0}: an output write failed (pending file: ' ||
  fail "the pending file's failure does not follow the first constructive find's checkpoint"
[ "$(outcome "$d" 'L2a1{0}')" = "io-error -1" ] || fail "L2a1{0} recorded '$(outcome "$d" 'L2a1{0}')'"
[ "$(lastoutcome "$d/log" 'L2a1\{0\}')" = io-error ] || fail "L2a1{0}'s outcome line is not io-error"

# 8-10 (goal runs): L2a1{0}, goal 0, one exhaustive cap-3 search, which meets it.
metconf() { # work
  mkdir -p "$1"; echo "$HDR" > "$1.db.csv"
  cat > "$1.conf" <<CONF
target_pd = $(grep -F 'L2a1{0},' "$DATA/links_to_6.csv" | cut -d, -f2)
target_name = L2a1{0}
work = $1
cobordisms = $1.db.csv
run_name = io_failure/L2a1
knot_table = $DATA/knots_to_6.csv
link_table = $DATA/links_to_6.csv
census = $1.none.sqlite
goal_genus = 0
literature = 0
threads = 4
max_searches = 1
surface_target = 1000000000
max_surface_target = 1000000000
max_faces = 3
iddfs_iterations = 0
iddfs_start = 0
iddfs_step = 0
root_budget_start = 0
resolve_unlinked = 1
CONF
}
#    The control: the goal is met, with its certificate.
g=$T/met; metconf "$g"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
[ "$rc" -eq 0 ] && grep -q '^\[+\] GOAL MET' "$g.log" && [ -s "$g/certificate.json" ] ||
  fail "the control goal run did not meet its goal with a certificate (exit $rc)"

# 8. A certificate that cannot be written (a directory of that name).
g=$T/cert; metconf "$g"; mkdir -p "$g/certificate.json/x"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
grep -E '^\[!!?\]|outcome' "$g.log" || true
expect2 $rc "an unwritable certificate"
grep -q '^\[!!\] HALT: an output write failed -- certificate: ' "$g.log" || fail "no HALT line for the certificate"
[ "$(lastoutcome "$g.log" 'L2a1\{0\}')" = io-error ] || fail "the run's outcome is not io-error"
! grep -q 'GOAL MET' "$g.log" || fail "an unwritable certificate: a GOAL MET line was printed"
! grep -rqsE 'GOAL MET|goal_met' "$g" || fail "an unwritable certificate: a record says the goal was met"
[ "$(wc -l < "$g.db.csv")" -ge 2 ] || fail "the goal run with no certificate signed nothing"

# 9a. A run record that cannot be written at the run's end (node_bounds.jsonl).
g=$T/bounds; metconf "$g"; mkdir -p "$g/node_bounds.jsonl/x"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
grep -E '^\[!!?\]|outcome' "$g.log" || true
expect2 $rc "an unwritable node_bounds.jsonl"
grep -q '^\[!!\] HALT: an output write failed -- .*node_bounds.jsonl' "$g.log" || fail "no HALT line for node_bounds.jsonl"
[ "$(lastoutcome "$g.log" 'L2a1\{0\}')" = io-error ] || fail "the run's outcome is not io-error"
notmet "$g.log" "$g" "an unwritable node_bounds.jsonl"

# 9b. A run record that cannot be written mid-run (cascade.jsonl, at the search's record).
g=$T/cascade; metconf "$g"; mkdir -p "$g/cascade.jsonl/x"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
grep -E '^\[!!?\]|outcome' "$g.log" || true
expect2 $rc "an unwritable cascade.jsonl"
grep -q '^\[!!\] HALT: an output write failed -- .*cascade.jsonl' "$g.log" || fail "no HALT line for cascade.jsonl"
[ "$(lastoutcome "$g.log" 'L2a1\{0\}')" = halted ] || fail "the run's outcome is not halted"
notmet "$g.log" "$g" "an unwritable cascade.jsonl"

# 10a. The sign step at depth 0: a database that refuses appends.
d=$T/sign0; conf "$d" "$T/knots.csv" 2; chmod a-w "$d/cobordisms.csv"
rc=0; "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err" || rc=$?
grep -E 'outcome|^\[!!?\]|not signed' "$d/log" || true
expect2 $rc "a database that refuses the sign step (depth 0)"
grep -q '^\[!!\] an output write failed at the run.s end -- database: ' "$d/log" || fail "no [!!] line for the sign step"
for r in 3_1 4_1; do
  [ "$(lastoutcome "$d/log" $r)" = io-error ] || fail "$r's last outcome line is not io-error"
  [ "$(outcome "$d" $r)" = "io-error -1" ] || fail "$r recorded '$(outcome "$d" $r)', not 'io-error -1'"
done
[ -s "$d/work/hop_0_n0/kept.csv" ] || fail "the unsigned cobordisms are not in their pending file"

# 10b. The sign step with a goal: halts, no GOAL MET.
g=$T/sign; metconf "$g"; chmod a-w "$g.db.csv"
rc=0; "$C" run --config "$g.conf" > "$g.log" 2>&1 || rc=$?
grep -E '^\[!!?\]|outcome|not signed' "$g.log" || true
expect2 $rc "a database that refuses the sign step (goal)"
grep -q '^\[!!\] HALT: an output write failed -- database: ' "$g.log" || fail "no HALT line for the sign step"
[ "$(lastoutcome "$g.log" 'L2a1\{0\}')" = io-error ] || fail "the run's outcome is not io-error"
notmet "$g.log" "$g" "a database that refuses the sign step (goal)"

# 10c. A verdicts shard that cannot be written (depth 0): its directory refuses files.
d=$T/verdicts; conf "$d" "$T/knots.csv" 2 "verdicts=$d/v/out.csv"; mkdir -p "$d/v"; chmod a-w "$d/v"
rc=0; "$C" run --config "$d/run.conf" > "$d/log" 2> "$d/err" || rc=$?
grep -E 'outcome|EXHAUSTIVE|^\[!!?\]' "$d/log" || true
expect2 $rc "an unwritable verdicts shard"
grep -q '^\[!!\] 3_1: an output write failed -- verdicts: ' "$d/log" || fail "no [!!] line for the verdicts"
! grep -q 'EXHAUSTIVE' "$d/log" || fail "a search claimed exhaustion with no verdicts written"
for r in 3_1 4_1; do
  [ "$(lastoutcome "$d/log" $r)" = io-error ] || fail "$r's last outcome line is not io-error"
done

echo "PASS: a failed output write ends its search as io-error (no frontier, no exhaustion claim, what it kept signed), never the process; exit 2; a run record, certificate, sign step or verdicts shard that cannot be written halts or reports io-error, and no goal is reported met without its certificate"
