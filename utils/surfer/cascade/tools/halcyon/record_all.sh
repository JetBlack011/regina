#!/bin/bash
# On halcyon: every CERTIFIED cascade proof so far recorded in a fresh store
# (cascade_record.py), each run under its name, and each run's proof
# directories gathered beside it, ready to copy into the atlas's
# results/cascade/. Log: record_all.log.
B=$HOME/cascade-bench
V=$B/v3
OUT=$B/atlas_store
L=$B/record_all.log
# Incremental: the recorder skips (run, target) pairs already in the store and
# surfaces already recorded, so this may run again as runs finish.
mkdir -p "$OUT"
echo "record_all starts $(date -Is)" >> "$L"
rec() { # <run name> <results dir>
  local run=$1 dir=$B/results/$2
  local targets=()
  for t in "$dir"/*/; do
    [ -f "$t/certificate.json" ] && grep -q 'VERDICT: CERTIFIED' "$t/check.txt" 2>/dev/null && targets+=("${t%/}")
  done
  [ ${#targets[@]} -eq 0 ] && { echo "$run: nothing certified in $dir" >> "$L"; return; }
  "$HOME/.venvs/atlas/bin/python" "$V/utils/surfer/cascade/tools/cascade_record.py" \
    --run "$run" --store "$OUT" --farsidediagram "$V/build/utils/surfer/farsidediagram" \
    --jobs 4 "${targets[@]}" 2>&1 | grep -v tkinter >> "$L"
  mkdir -p "$OUT/$run"
  for t in "${targets[@]}"; do rsync -a "$t" "$OUT/$run/"; done
}
# Runs to record: given as arguments (<run name>=<results dir>), or all.
if [ $# -gt 0 ]; then
  for a in "$@"; do rec "${a%%=*}" "${a#*=}"; done
else
  rec 2026-09-28_knots8 knots8_pure
  rec 2026-09-28_knots9 knots9_master
  rec 2026-09-28_knots10 knots10_pure
  [ -d "$B/results/knots10_master" ] && rec 2026-09-28_knots10_master knots10_master
  rec 2026-09-28_links links_support
  rec 2026-09-28_dg18 dg18
fi
echo "record_all done $(date -Is)" >> "$L"
