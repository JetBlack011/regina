// leaves.h
//
// Which outside facts may become proof-graph leaves. See README.md, "Leaf
// facts".

#pragma once

#include <optional>
#include <string>
#include <utility>

namespace cascade {

/// A table's literature 4-genus, "2" or "[0;1]", as (lo, hi); nullopt if the
/// field is malformed (a malformed value must never become a bound).
std::optional<std::pair<int, int>> parseTableG4(const std::string &s);

/**
 * Whether a node whose exact table class is `nodeClass` may receive that
 * class's literature UPPER bound as a leaf. Never the target's own class:
 * the target's literature value proving the target would be circular, and a
 * duplicate node of the target (a diagram the registry did not recognise)
 * would otherwise do exactly that. Lower bounds are always kept: they only
 * gate contradictions.
 */
bool mayUseLiteratureUpperBound(const std::string &nodeClass,
                                const std::string &targetClass,
                                bool literatureAllowed);

} // namespace cascade
