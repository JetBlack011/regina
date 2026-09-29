// leaves_test.cpp: the literature leaf policy (README.md, "Leaf facts").

#include "../leaves.h"
#include "check.h"

using namespace cascade;

int main() {
  // F1: table values parse exactly; malformed ones never become bounds.
  CHECK(parseTableG4("2") == std::make_pair(2, 2), "plain value");
  CHECK(parseTableG4("[0;1]") == std::make_pair(0, 1), "interval");
  CHECK(parseTableG4("[1;2]") == std::make_pair(1, 2), "interval 1..2");
  for (const char *bad : {"", "x", "[1;0]", "[;1]", "[0;]", "[0,1]", "-1", "1.5", "[0;1"})
    CHECK(!parseTableG4(bad).has_value(), std::string("malformed refused: ") + bad);
  // F2: never the target's own class (circular), otherwise when allowed.
  CHECK(!mayUseLiteratureUpperBound("13n_65", "13n_65", true),
        "the target's own literature never proves it");
  CHECK(mayUseLiteratureUpperBound("L7a4{0}", "13n_65", true), "other names may");
  CHECK(!mayUseLiteratureUpperBound("L7a4{0}", "13n_65", false),
        "constructive mode uses no literature");
  CHECK(!mayUseLiteratureUpperBound("", "13n_65", true), "an unnamed node has none");
  CHECK(mayUseLiteratureUpperBound("3_1", "", true),
        "an off-table target excludes nothing");
  return cascadetest::finish("leaves_test");
}
