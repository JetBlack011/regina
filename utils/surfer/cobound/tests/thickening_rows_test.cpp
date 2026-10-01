// thickening_rows_test.cpp
//
// buildAmbient()'s component count (ThickenedLink::componentCount, an edge
// walk: edgecycles::countClosedCurves()) on every row of the small test
// tables (data/), against the count the row's name states and Regina's
// count of its diagram's components. Phase 2 checked the same on all 17,153
// table rows once, with a one-off probe; this keeps it under ctest.

#include <string>

#include "diagramtriangulation/thickening/thickening.h"
#include "linknaming/names.h"
#include "linknaming/tables.h"
#include "linknaming/tests/check.h"

int main() {
  size_t rows = 0;
  for (const char *file : {"knots_to_6.csv", "links_to_6.csv"})
    for (const exactnaming::TableRow &row :
         exactnaming::readTableRows(std::string(CASCADE_TEST_DATA) + "/" + file)) {
      ++rows;
      ThickenedLink t;
      buildAmbient(row.pd, 2, 2, /*useCone=*/false, t);
      CHECK_EQ(t.componentCount, cobordismgraph::componentsFromName(row.name),
               row.name + ": the edge walk counts the components the name states");
      CHECK_EQ(static_cast<size_t>(t.componentCount),
               exactnaming::linkFromTablePD(row.pd).countComponents(),
               row.name + ": and Regina's diagram agrees");
    }
  CHECK_EQ(rows, size_t(19), "every row of both small tables");
  return cascadetest::finish("thickening_rows_test");
}
