// leaves.cpp

#include "cobound/bounds/axioms.h"

namespace cascade {

bool mayUseLiteratureUpperBound(const std::string &nodeClass,
                                const std::string &targetClass,
                                bool literatureAllowed) {
  if (!literatureAllowed || nodeClass.empty()) return false;
  return targetClass.empty() || nodeClass != targetClass;
}

} // namespace cascade
