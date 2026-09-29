#!/bin/bash
# On halcyon: v9 = v3's sources (everything the ladder runs) + patchJ_lift
# (components lying above or below everything lifted off as split unknots,
# in the cascade and in the checker), in its own tree and build, so the
# running ladder's binary is untouched. 6 jobs, beside the ladder.
# Then the cascade and exactnaming tests and the checker's own tests.
# Log: v9.log.
B=$HOME/cascade-bench
L=$B/v9.log
log() { echo "$*" >> "$L"; }
start=$(date +%s)
rsync -a --exclude build "$B/v3/" "$B/v9/"
cd "$B/v9" || exit 1
# cascade_check_test.py is new since v3's patches: it comes whole (shipped
# beside the patch), and the patch's hunk for it is left out.
cp "$B/cascade_check_test.py" utils/surfer/cascade/tools/
if git apply --exclude=utils/surfer/cascade/tools/cascade_check_test.py "$B/patchJ_lift.patch" \
   && cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DDISABLE_GUI=1 -DDISABLE_PYTHON=1 > configure.v9.log 2>&1 \
   && make -C build -j6 cascadesearch farsidediagram cascade_profile_test cascade_proofgraph_test \
        cascade_diagramiso_test cascade_nodes_test cascade_hopedges_test cascade_hoprunner_test \
        cascade_leaves_test exactnaming_test snappeaisometry_test > build.v9.log 2>&1 \
   && { ( cd build/utils/surfer/cascade && ctest --output-on-failure ) \
        && ( cd build/utils/surfer/exactnaming && ctest --output-on-failure ) \
        && "$HOME/.venvs/atlas/bin/python" utils/surfer/cascade/tools/cascade_check_test.py; } > tests.v9.log 2>&1; then
  log "v9 built and tested in $(( $(date +%s) - start )) s $(date -Is): $(grep -h 'tests passed\|^passed' tests.v9.log | tr '\n' ' ')"
else
  log "v9 FAILED $(date -Is): see $B/v9/build.v9.log and tests.v9.log"
fi
