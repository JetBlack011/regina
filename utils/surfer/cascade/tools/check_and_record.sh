#!/bin/bash
# Check a run's certificates in parallel (cascade_check.py), then record
# every CERTIFIED proof in the cascade store (cascade_record.py) and copy its
# directory beside it. Incremental: a certificate already CERTIFIED is not
# checked again, and the recorder skips proofs it has recorded.
#
# usage: check_and_record.sh <run name> <results dir> [store] [jobs]
#   results dir  holds one directory per target, each with certificate.json
#   store        default: the atlas's results/cascade/
#   jobs         checks at once (default: all cores less two)
# FSD (farsidediagram) and PY (the atlas venv's python) may be overridden.
set -u
HERE=$(cd "$(dirname "$0")" && pwd)
REPO=$(cd "$HERE/../../../.." && pwd)
RUN=$1 DIR=$2
STORE=${3:-$HOME/Projects/triangles/cobordism-atlas/results/cascade}
JOBS=${4:-$(( $(nproc) - 2 ))}
PY=${PY:-$HOME/.venvs/atlas/bin/python}
FSD=${FSD:-$REPO/build/utils/surfer/farsidediagram}

check_one() {
  local d=$1
  grep -q 'VERDICT: CERTIFIED' "$d/check.txt" 2>/dev/null && return
  "$PY" "$HERE/cascade_check.py" --farsidediagram "$FSD" "$d/certificate.json" \
    > "$d/check.txt" 2>&1
}
export -f check_one
export PY HERE FSD
find "$DIR" -mindepth 2 -maxdepth 2 -name certificate.json -printf '%h\0' |
  xargs -0 -P "$JOBS" -I{} bash -c 'check_one "$1"' _ {}

certified=() failed=()
for d in "$DIR"/*/; do
  d=${d%/}
  [ -f "$d/certificate.json" ] || continue
  if grep -q 'VERDICT: CERTIFIED' "$d/check.txt" 2>/dev/null; then certified+=("$d")
  else failed+=("$(basename "$d")"); fi
done
echo "$RUN: ${#certified[@]} CERTIFIED, ${#failed[@]} not: ${failed[*]}"
[ ${#certified[@]} -eq 0 ] && exit 0
"$PY" "$HERE/cascade_record.py" --run "$RUN" --store "$STORE" --farsidediagram "$FSD" \
  --jobs "$JOBS" "${certified[@]}" 2>&1 | grep -v -i tkinter
mkdir -p "$STORE/$RUN"
for d in "${certified[@]}"; do rsync -a "$d" "$STORE/$RUN/"; done
