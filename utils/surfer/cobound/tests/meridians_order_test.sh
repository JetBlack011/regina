#!/bin/sh
# meridians_order_test.sh <cobound> <test data dir>
#
# cobound meridians (peripheral_slopes) answers each input line on its own: the same records in
# the reverse order must give the same output, record by record. Until phase
# 3.0 they did not (10 of 24 `sig` records changed on reversal):
#   - a surface's boundary edges were collected, and split into curves, in
#     the order of their addresses, which follow everything decoded before;
#   - drilling pinched edges in the order of a hash map keyed by addresses;
#   - `sig` simplifies each complement with Regina's random moves, drawn from
#     one process-wide generator, so a record's answer depended on how many
#     draws the records before it had made.
#
# Inputs: the stored cobordisms' pair signatures in data/search_10_3_cobordisms.csv
# (sig, dump) and the small tables' diagrams (dump-link).
set -eu

P=$1
D=$2
T=$(mktemp -d)
trap 'rm -rf "$T"' EXIT

# "<id> <pairsig>": the first 12 witnesses (column 9 is the pair signature).
awk -F, 'NR > 1 && NR <= 13 { print "w" NR, $9 }' "$D/search_10_3_cobordisms.csv" > "$T/pairsigs"
# "<id> <pd>": every row of both small tables.
for f in "$D/knots_to_6.csv" "$D/links_to_6.csv"; do
  awk -F, 'NR > 1 { print $1, $2 }' "$f"
done > "$T/pds"

# One line per record (RECORD..ENDRECORD blocks joined), sorted.
records() {
  awk '/^RECORD / { block = $0; inblock = 1; next }
       inblock { block = block "|" $0; if ($0 ~ /^ENDRECORD/) { print block; inblock = 0 }; next }
       { print }' "$1" | LC_ALL=C sort
}

fail=0
for spec in sig:pairsigs dump:pairsigs dump-link:pds; do
  mode=${spec%%:*}; input=$T/${spec#*:}
  "$P" meridians "$mode" < "$input" > "$T/$mode.fwd"
  awk '{ l[NR] = $0 } END { for (i = NR; i >= 1; i--) print l[i] }' "$input" > "$T/$mode.rin"
  "$P" meridians "$mode" < "$T/$mode.rin" > "$T/$mode.rev"
  records "$T/$mode.fwd" > "$T/$mode.fwd.sorted"
  records "$T/$mode.rev" > "$T/$mode.rev.sorted"
  n=$(wc -l < "$T/$mode.fwd.sorted")
  good=$(grep -vc '^FAILED' "$T/$mode.fwd.sorted" || true)
  if [ "$n" -lt "$(wc -l < "$input")" ] || [ "$good" -eq 0 ]; then
    echo "FAIL: $mode answered $n records for $(wc -l < "$input") inputs, $good not FAILED"; fail=1
  elif cmp -s "$T/$mode.fwd.sorted" "$T/$mode.rev.sorted"; then
    echo "ok: $mode, $n records, the same in both orders"
  else
    echo "FAIL: $mode: $(diff "$T/$mode.fwd.sorted" "$T/$mode.rev.sorted" | grep -c '^<') of $n records change when the input is reversed"
    fail=1
  fi
done
[ $fail = 0 ] && echo "PASS: every record is a function of its own input line"
exit $fail
