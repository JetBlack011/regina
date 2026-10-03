//
//  targets.cpp
//

#include "cobound/driver/targets.h"

#include <algorithm>
#include <stdexcept>

#include "diagramtriangulation/pdcode.h"
#include "linknaming/tables.h"

namespace targets {

std::vector<cobordismgraph::InputRow> loadInputCsv(const std::filesystem::path &path) {
  std::vector<cobordismgraph::InputRow> rows;
  for (const exactnaming::TableRow &table : exactnaming::readTableRows(path)) {
    // A row to search must state its literature bounds: a malformed field
    // stops the run (it always did), never becomes a bound.
    const auto g4 = exactnaming::parseTableG4(table.g4);
    if (!g4)
      throw std::runtime_error("malformed 4-genus field '" + table.g4 + "' for " +
                               table.name + " in " + path.string());
    cobordismgraph::InputRow row;
    row.name = table.name;
    row.pdNotation = table.pd;
    row.lo = g4->first;
    row.hi = g4->second;
    const std::string &pd = row.pdNotation;
    // Crossing count is derived from the PD code itself (works uniformly
    // for both knot names like "13n_1109" and link names like "L10a1{0}",
    // which have no leading digit run to parse) rather than from `name`.
    row.crossings = static_cast<int>(knotbuilder::parsePDCode(pd).size());
    rows.push_back(std::move(row));
  }
  return rows;
}

cobordismgraph::InputRow oneDiagram(const std::string &name, const std::string &pd,
                                    const std::vector<std::filesystem::path> &tables) {
  cobordismgraph::InputRow row;
  row.name = name;
  row.pdNotation = pd;
  row.lo = 0;
  row.hi = 99;
  for (const std::filesystem::path &table : tables)
    for (const exactnaming::TableRow &t : exactnaming::readTableRows(table))
      if (t.name == name)
        if (const auto g4 = exactnaming::parseTableG4(t.g4)) {
          row.lo = g4->first;
          row.hi = g4->second;
        }
  row.crossings = static_cast<int>(knotbuilder::parsePDCode(pd).size());
  return row;
}

std::vector<cobordismgraph::InputRow>
searchOrder(const std::vector<cobordismgraph::InputRow> &rows, int maxCrossings,
            std::unordered_map<std::string, verdicts::OutputRow> &outputRows) {
  std::vector<cobordismgraph::InputRow> pending;
  pending.reserve(rows.size());
  for (const auto &row : rows) {
    if (row.crossings > maxCrossings) {
      if (!outputRows.contains(row.name)) {
        verdicts::OutputRow out;
        out.knot = row.name;
        out.status = "skipped";
        out.witnessKind = "none";
        out.literatureLo = row.lo;
        out.literatureHi = row.hi;
        outputRows[row.name] = std::move(out);
      }
      continue;
    }
    pending.push_back(row);
  }
  std::stable_sort(pending.begin(), pending.end(),
                   [](const cobordismgraph::InputRow &a, const cobordismgraph::InputRow &b) {
                     return a.crossings < b.crossings;
                   });
  return pending;
}

std::map<std::string, std::string> tablePDs(const std::vector<std::string> &files) {
  std::map<std::string, std::string> out;
  for (const std::string &file : files)
    for (const exactnaming::TableRow &row : exactnaming::readTableRows(file))
      out[row.name] = row.pd;
  return out;
}

} // namespace targets
