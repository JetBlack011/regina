// leaves_test.cpp: the literature leaf policy (README.md, "Leaf facts").

#include "cobound/bounds/axioms.h"
#include "linknaming/tests/check.h"

using namespace bounds;

int main() {
  // F1 (table values parse exactly; malformed ones never become bounds) is
  // tables_test's now, with the parser (linknaming/tables.h).
  // F2: never the target's own class (circular), otherwise when allowed.
  CHECK(!mayUseLiteratureUpperBound("13n_65", "13n_65", true),
        "the target's own literature never proves it");
  CHECK(mayUseLiteratureUpperBound("L7a4{0}", "13n_65", true), "other names may");
  CHECK(!mayUseLiteratureUpperBound("L7a4{0}", "13n_65", false),
        "constructive mode uses no literature");
  CHECK(!mayUseLiteratureUpperBound("", "13n_65", true), "an unnamed node has none");
  CHECK(mayUseLiteratureUpperBound("3_1", "", true),
        "an off-table target excludes nothing");
  return checks::finish("axioms_test");
}
