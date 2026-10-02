#!/usr/bin/env bash
# Rough, repeatable search benchmarks for verifyslicegenus, one row per run.
# Results land in $BENCH_DIR/results.tsv; the "Search performance history"
# table in utils/surfer/README.md is written from them.
#
#   bench_search.sh run <binary> <b1|b2> <row>      one run, one result line
#   bench_search.sh ab  <binaryA> <binaryB> [reps]   B1 on every B1 row, A and B
#                                                    alternated (ABAB...), then
#                                                    a per-row comparison
#   bench_search.sh b2  <binary>                     B2 on its row
#   bench_search.sh summary [tag]                    re-print a comparison
#
# B1 (exhaust-4): one IDDFS round at cap 4, no surface target, an empty
#   witness file. Fixed work: every correct version accepts exactly the same
#   surfaces, so `accepted` must match across binaries (a free correctness
#   check) and wall time compares like for like. About a minute per row on
#   halcyon at f183dffbe. Its budget passes make each root's traversal a
#   function of the root alone, so the `search profile:` walk counters
#   (nodes, attempts, evaluated, replayed) are deterministic too: `summary`
#   flags any difference, which an output-identical change must not have.
# B2 (production): the full [campaign] shape (surface target, caps 4 then 5),
#   with the witness store in $STORE copied in if set -- round 2, the drain
#   tail and peak RSS as a campaign sees them.
#
# Every run starts from the same frozen census snapshot, taken once per
# $BENCH_DIR, so no run's census writes help a later one.
#
# Environment:
#   BENCH_DIR  runs and results.tsv (default ~/bench-runs)
#   ATLAS      cobordism-atlas checkout: data/ and tools/orchestrate/hosts.conf
#              (default: the first of ~/Projects/cobordism-atlas and
#              ~/Projects/triangles/cobordism-atlas that exists)
#   HOST       hosts.conf section for threads and cache limits
#              (default: this machine's short hostname)
#   THREADS    override the host's thread count
#   CENSUS     census db to snapshot (default: the host section's census_db)
#   STORE      witness store copied in for B2 (default: none, i.e. empty)
#   B1_ROWS    default: 10_141 L10a14{0} L10a127{1;1}
#   B2_ROW     default: 10_141
#   ROOT_BUDGET_START  override [campaign]'s (0 = unbudgeted; used to check
#              that a budgeted run's charged - replayed equals an unbudgeted
#              run's attempts)
#   TAG        tag for `run` results (default: none)
#   EXTRA_ARGS further verifyslicegenus flags, e.g. --audit-linking
set -euo pipefail

BENCH_DIR=${BENCH_DIR:-$HOME/bench-runs}
if [ -z "${ATLAS:-}" ]; then
  for d in "$HOME/Projects/cobordism-atlas" "$HOME/Projects/triangles/cobordism-atlas"; do
    [ -d "$d/data" ] && { ATLAS=$d; break; }
  done
fi
: "${ATLAS:?bench_search.sh: set ATLAS to a cobordism-atlas checkout}"
HOST=${HOST:-$(hostname -s)}
CONF=$ATLAS/tools/orchestrate/hosts.conf
DATA=$ATLAS/data
B1_ROWS=${B1_ROWS:-"10_141 L10a14{0} L10a127{1;1}"}
B2_ROW=${B2_ROW:-10_141}
RESULTS=$BENCH_DIR/results.tsv
mkdir -p "$BENCH_DIR"

# The search shape from [campaign], and the host's own limits, exactly as
# dispatch.py hands them to remote_run.sh.
eval "$(python3 - "$CONF" "$HOST" <<'EOF'
import configparser, shlex, sys
c = configparser.ConfigParser(inline_comment_prefixes=('#',), interpolation=None)
c.read(sys.argv[1])
host = sys.argv[2]
for sec, keys in (
    ('campaign', ['surface_target', 'max_faces', 'iddfs_iterations',
                  'iddfs_start', 'iddfs_step', 'root_budget_start',
                  'root_budget_growth', 'thicken_layers', 'collar_layers',
                  'max_crossings', 'resolve_unlinked', 'exact_far_side_names']),
    (host, ['threads', 'secs', 'census_db', 'pending_surface_cap',
            'petal_cache_limit', 'recognition_cache_limit',
            'boundary_signature_cache_limit', 'retriangulate_time_budget'])):
    for k in keys:
        v = c.get(sec, k, fallback='') if c.has_section(sec) else ''
        print(f"CFG_{k.upper()}={shlex.quote(v)}")
EOF
)"
[ -n "$CFG_THREADS" ] || { echo "bench_search.sh: no [$HOST] section in $CONF" >&2; exit 2; }
THREADS=${THREADS:-$CFG_THREADS}
CENSUS=${CENSUS:-$CFG_CENSUS_DB}

FROZEN=$BENCH_DIR/census.frozen.sqlite
if [ ! -s "$FROZEN" ]; then
  # The backup API, never a bare copy: worker censuses are WAL-mode.
  python3 - "$CENSUS" "$FROZEN" <<'EOF'
import sqlite3, sys
src = sqlite3.connect(f"file:{sys.argv[1]}?mode=ro", uri=True)
dst = sqlite3.connect(sys.argv[2])
src.backup(dst)
dst.close(); src.close()
EOF
fi

run_one() {
  local bin=$1 mode=$2 row=$3 tag=${4:-}
  local name safe run
  name=$(basename "$bin")
  safe=$(printf '%s' "$row" | tr -c 'A-Za-z0-9_.-' '_')
  run=$BENCH_DIR/runs/$name/$mode/$safe${ROOT_BUDGET_START:+-budget$ROOT_BUDGET_START}
  rm -rf "$run"; mkdir -p "$run"

  {
    echo "Name,PD Notation,Genus-4D"
    awk -F, -v r="$row" '$1 == r { print; found = 1; exit } END { exit !found }' \
        "$DATA/4d_smooth_slice_genus_13_crossings_pd_codes.csv" ||
    awk -F, -v r="$row" '$1 == r { print; found = 1; exit } END { exit !found }' \
        "$DATA/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv"
  } > "$run/rows.csv" || { echo "bench_search.sh: no row $row in the tables" >&2; return 2; }
  if [ "$mode" = b2 ] && [ -n "${STORE:-}" ]; then
    cp --reflink=auto "$STORE" "$run/cobordisms.csv"
  else
    printf 'kind,subject,subject_components,other,other_candidates,other_components,genus,tubed,pairsig,source_row,thicken_layers,max_faces,resolved_vertices\n' \
        > "$run/cobordisms.csv"
  fi
  cp "$FROZEN" "$run/census.sqlite"

  local cmd=(
    "$bin" --input "$run/rows.csv" --output "$run/out.csv"
    --cobordisms "$run/cobordisms.csv" --surface-stats "$run/surface_stats.csv"
    --census-db "$run/census.sqlite"
    --knot-table "$DATA/4d_smooth_slice_genus_13_crossings_pd_codes.csv"
    --link-table "$DATA/links_4d_smooth_slice_genus_11_crossings_pd_codes.csv"
    --name-aliases "$DATA/name_aliases.csv"
    --max-crossings "$CFG_MAX_CROSSINGS"
    --thicken-layers "$CFG_THICKEN_LAYERS" --collar-layers "$CFG_COLLAR_LAYERS"
    --root-budget-start "${ROOT_BUDGET_START:-$CFG_ROOT_BUDGET_START}"
    --root-budget-growth "$CFG_ROOT_BUDGET_GROWTH"
    --no-cone --harvest --boundary-condition proper --research-settled
    --threads "$THREADS"
    --retriangulate-time-budget "$CFG_RETRIANGULATE_TIME_BUDGET"
    --pending-surface-cap "$CFG_PENDING_SURFACE_CAP"
    --petal-cache-limit "$CFG_PETAL_CACHE_LIMIT"
    --recognition-cache-limit "$CFG_RECOGNITION_CACHE_LIMIT"
    --boundary-signature-cache-limit "$CFG_BOUNDARY_SIGNATURE_CACHE_LIMIT"
  )
  case $mode in
    b1) cmd+=(--max-faces 4 --iddfs-iterations 1 --iddfs-start 4 --iddfs-step 1
              --per-knot-time-limit 7200) ;;
    b2) cmd+=(--max-faces "$CFG_MAX_FACES" --iddfs-iterations "$CFG_IDDFS_ITERATIONS"
              --iddfs-start "$CFG_IDDFS_START" --iddfs-step "$CFG_IDDFS_STEP"
              --per-knot-time-limit "$CFG_SECS")
        [ -n "$CFG_SURFACE_TARGET" ] && cmd+=(--surface-target "$CFG_SURFACE_TARGET") ;;
    *) echo "bench_search.sh: mode must be b1 or b2" >&2; return 2 ;;
  esac
  if [ "$CFG_RESOLVE_UNLINKED" = 1 ]; then cmd+=(--resolve-unlinked); else cmd+=(--no-resolve-unlinked); fi
  [ "$CFG_EXACT_FAR_SIDE_NAMES" = 1 ] && cmd+=(--exact-far-side-names)
  # Word-split on purpose: extra flags, e.g. EXTRA_ARGS=--audit-linking.
  # shellcheck disable=SC2206
  [ -n "${EXTRA_ARGS:-}" ] && cmd+=($EXTRA_ARGS)

  local rc=0
  if [ -x /usr/bin/time ]; then
    /usr/bin/time -f '%e %M' -o "$run/time.txt" "${cmd[@]}" > "$run/log" 2> "$run/err" || rc=$?
  else
    # No GNU time: wall from the clock, peak RSS from the kernel's own
    # high-water mark (VmHWM), polled until the process exits.
    local t0 hwm=0 pid now
    t0=$(date +%s.%N)
    "${cmd[@]}" > "$run/log" 2> "$run/err" &
    pid=$!
    while kill -0 "$pid" 2>/dev/null; do
      now=$(awk '/^VmHWM:/ { print $2 }' "/proc/$pid/status" 2>/dev/null || true)
      [ -n "$now" ] && [ "$now" -gt "$hwm" ] && hwm=$now
      sleep 0.5
    done
    wait "$pid" || rc=$?
    echo "$(echo "$(date +%s.%N) - $t0" | bc) $hwm" > "$run/time.txt"
  fi
  [ "$rc" -eq 0 ] || {
    echo "bench_search.sh: $name exited $rc on $row ($mode); see $run" >&2; return 1; }

  python3 - "$RESULTS" "$run" "$name" "$mode" "$row" "$THREADS" "$HOST" "$tag" <<'EOF'
import os, re, sys, time
results, run, name, mode, row, threads, host, tag = sys.argv[1:]
wall, rss_kb = open(f"{run}/time.txt").read().split()[-2:]
log = re.sub(r'\x1b\[[0-9;]*[A-Za-z]', '', open(f"{run}/log").read())
def grab(pattern, default=''):
    m = re.search(pattern, log)
    return m.group(1) if m else default
f = {
    'time': time.strftime('%Y-%m-%dT%H:%M:%S'), 'host': host, 'tag': tag,
    'binary': name, 'mode': mode, 'row': row, 'threads': threads,
    'wall_s': wall, 'rss_mb': str(int(rss_kb) // 1024),
    'accepted': grab(r'accounting: accepted (\d+)'),
    'outcome': grab(r'new witnesses, outcome ([a-z-]+)'),
    'exhausted': grab(r'EXHAUSTIVE to (\d+)'),
    'prototype_s': grab(r'search profile: prototype ([\d.]+)s'),
    'proto_unknot': grab(r'prototype [\d.]+s \(unknot misses (\d+ in [\d.]+s)'),
    'rounds_s': grab(r'; rounds ([\d.s ]+);').strip(),
    'drain_tail': grab(r'drain tail (\d+ surfaces in [\d.]+s)'),
    'nodes': grab(r'; nodes (\d+)'),
    'attempts': grab(r', attempts (\d+)'),
    'evaluated': grab(r', evaluated (\d+)'),
    'charged': grab(r', charged (\d+)'),
    'replayed': grab(r', replayed (\d+)'),
    'unknot_misses': grab(r'petal misses: unknot (\d+ in [\d.]+s)'),
    'linking_misses': grab(r'petal misses: unknot \d+ in [\d.]+s, linking (\d+ in [\d.]+s)'),
}
new = not os.path.exists(results)
if not new:
    # Never append a row under a header with other columns: every later
    # reader would shift the fields.
    header = open(results).readline().rstrip('\n').split('\t')
    if header != list(f):
        sys.exit(f"bench_search.sh: {results} has columns {header}, this "
                 f"version writes {list(f)}; move it aside first")
with open(results, 'a') as out:
    if new:
        out.write('\t'.join(f) + '\n')
    out.write('\t'.join(v.replace('\t', ' ') for v in f.values()) + '\n')
print('  '.join(f"{k}={v}" for k, v in f.items() if k not in ('time', 'host')))
EOF
}

summary() {
  python3 - "$RESULTS" "${1:-}" <<'EOF'
import csv, statistics, sys
results, tag = sys.argv[1], sys.argv[2]
rows = list(csv.DictReader(open(results), delimiter='\t'))
if not tag:
    tags = [r['tag'] for r in rows if r['tag']]
    tag = tags[-1] if tags else ''
rows = [r for r in rows if r['tag'] == tag]
bins = list(dict.fromkeys(r['binary'] for r in rows))
print(f"tag {tag}: " + ' vs '.join(bins))
total = {b: 0.0 for b in bins}
for row in dict.fromkeys(r['row'] for r in rows):
    cells, accepted = [], set()
    for b in bins:
        runs = [r for r in rows if r['row'] == row and r['binary'] == b]
        walls = [float(r['wall_s']) for r in runs]
        accepted |= {r['accepted'] for r in runs}
        med = statistics.median(walls) if walls else float('nan')
        total[b] += med
        cells.append(f"{b} {med:7.1f}s ({', '.join(f'{w:.1f}' for w in walls)})")
    ratio = ''
    if len(bins) == 2:
        a = [float(r['wall_s']) for r in rows if r['row'] == row and r['binary'] == bins[0]]
        z = [float(r['wall_s']) for r in rows if r['row'] == row and r['binary'] == bins[1]]
        if a and z:
            ratio = f"  -> {statistics.median(z) / statistics.median(a):.3f}x"
    flag = '' if len(accepted) == 1 else f"  ACCEPTED COUNTS DIFFER: {sorted(accepted)}"
    # In a budgeted exhaustive round each root's traversal is a function of
    # the root alone, so these are deterministic: equal for an output-
    # identical change, and expected to differ only for a traversal change.
    walk = {tuple(r.get(k, '') for k in ('nodes', 'attempts', 'evaluated',
                                         'charged', 'replayed'))
            for r in rows if r['row'] == row and r['nodes']}
    if len(walk) > 1:
        flag += '  TRAVERSAL DIFFERS (nodes/attempts/evaluated/charged/replayed): ' + \
                '; '.join('/'.join(w) for w in sorted(walk))
    print(f"  {row:16s} " + ' | '.join(cells) + ratio + flag)
if len(bins) == 2 and total[bins[0]] > 0:
    print(f"  {'total':16s} {total[bins[0]]:.1f}s -> {total[bins[1]]:.1f}s "
          f"({total[bins[1]] / total[bins[0]]:.3f}x)")
EOF
}

case ${1:-} in
  run)
    [ $# -eq 4 ] || { echo "usage: $0 run <binary> <b1|b2> <row>" >&2; exit 2; }
    run_one "$2" "$3" "$4" "${TAG:-}" ;;
  ab)
    [ $# -ge 3 ] || { echo "usage: $0 ab <binaryA> <binaryB> [reps]" >&2; exit 2; }
    tag="ab-$(date +%Y%m%dT%H%M%S)"
    for ((rep = 1; rep <= ${4:-2}; rep++)); do
      for row in $B1_ROWS; do
        run_one "$2" b1 "$row" "$tag"
        run_one "$3" b1 "$row" "$tag"
      done
    done
    summary "$tag" ;;
  b2)
    [ $# -eq 2 ] || { echo "usage: $0 b2 <binary>" >&2; exit 2; }
    run_one "$2" b2 "$B2_ROW" "b2-$(date +%Y%m%dT%H%M%S)" ;;
  summary)
    summary "${2:-}" ;;
  *)
    sed -n '2,14p' "$0"; exit 2 ;;
esac
