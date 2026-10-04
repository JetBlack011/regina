//
//  targets.cpp
//

#include "cobound/driver/targets.h"

#include <algorithm>
#include <stdexcept>

#include "diagramtriangulation/pdcode.h"
#include "linknaming/tables.h"

namespace targets {

std::vector<solver::InputRow> loadInputCsv(const std::filesystem::path &path) {
  std::vector<solver::InputRow> rows;
  for (const linknaming::TableRow &table : linknaming::readTableRows(path)) {
    // A row to search must state its literature bounds: a malformed field
    // stops the run (it always did), never becomes a bound.
    const auto g4 = linknaming::parseTableG4(table.g4);
    if (!g4)
      throw std::runtime_error("malformed 4-genus field '" + table.g4 + "' for " +
                               table.name + " in " + path.string());
    solver::InputRow row;
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

solver::InputRow oneDiagram(const std::string &name, const std::string &pd,
                                    const std::vector<std::filesystem::path> &tables) {
  solver::InputRow row;
  row.name = name;
  row.pdNotation = pd;
  row.lo = 0;
  row.hi = 99;
  for (const std::filesystem::path &table : tables)
    for (const linknaming::TableRow &t : linknaming::readTableRows(table))
      if (t.name == name)
        if (const auto g4 = linknaming::parseTableG4(t.g4)) {
          row.lo = g4->first;
          row.hi = g4->second;
        }
  row.crossings = static_cast<int>(knotbuilder::parsePDCode(pd).size());
  return row;
}

std::vector<solver::InputRow>
searchOrder(const std::vector<solver::InputRow> &rows, int maxCrossings,
            std::unordered_map<std::string, verdicts::OutputRow> &outputRows) {
  std::vector<solver::InputRow> pending;
  pending.reserve(rows.size());
  for (const auto &row : rows) {
    if (row.crossings > maxCrossings) {
      if (!outputRows.contains(row.name)) {
        verdicts::OutputRow out;
        out.knot = row.name;
        out.status = "skipped";
        out.cobordismKind = "none";
        out.literatureLo = row.lo;
        out.literatureHi = row.hi;
        outputRows[row.name] = std::move(out);
      }
      continue;
    }
    pending.push_back(row);
  }
  std::stable_sort(pending.begin(), pending.end(),
                   [](const solver::InputRow &a, const solver::InputRow &b) {
                     return a.crossings < b.crossings;
                   });
  return pending;
}

std::map<std::string, std::string> tablePDs(const std::vector<std::string> &files) {
  std::map<std::string, std::string> out;
  for (const std::string &file : files)
    for (const linknaming::TableRow &row : linknaming::readTableRows(file))
      out[row.name] = row.pd;
  return out;
}

} // namespace targets
