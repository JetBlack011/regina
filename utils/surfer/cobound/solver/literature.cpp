//
//  literature.cpp
//

#include "cobound/solver/literature.h"

#include "linknaming/names.h"
#include "linknaming/tables.h"

namespace cobordismgraph {

void NameTable::addLiterature(const std::string &name, int lo, int hi) {
    auto [it, inserted] = info_.try_emplace(name);
    if (inserted) {
        it->second.components = componentsFromName(name);
        byBase_[baseName(name)].push_back(name);
    }
    it->second.haveLiterature = true;
    it->second.litLo = lo;
    it->second.litHi = hi;
}

const NameInfo *NameTable::find(const std::string &name) const {
    auto it = info_.find(name);
    return it == info_.end() ? nullptr : &it->second;
}

int NameTable::components(const std::string &name) const {
    auto it = info_.find(name);
    return it == info_.end() ? componentsFromName(name) : it->second.components;
}

std::vector<std::string>
NameTable::candidates(const std::string &name,
                      std::optional<int> observedComponents) const {
    auto it = byBase_.find(baseName(name));
    if (it == byBase_.end() || it->second.empty())
        return {name};

    if (!observedComponents)
        return it->second;

    std::vector<std::string> filtered;
    for (const std::string &candidate : it->second)
        if (components(candidate) == *observedComponents)
            filtered.push_back(candidate);
    // An empty result means every registered variant of this base has the
    // wrong component count, i.e. the identification and the geometry
    // disagree. The far side is none of those variants, so a max/min over
    // them bounds nothing: hand back the name itself, which is unregistered
    // and so bears no bound.
    return filtered.empty() ? std::vector<std::string>{name} : filtered;
}

} // namespace cobordismgraph

namespace witnessstore {

// Loads a literature table for its names and bounds only, skipping the PD
// code entirely. Used for tables that aren't this run's --input: we need
// their names (to expand orientation-blind identifications into candidate
// sets) and their bounds, but never build a triangulation from them, so
// there is no reason to pay parsePDCode()'s cost across 12k+ rows.
size_t loadNameTable(const std::filesystem::path &path,
                     cobordismgraph::NameTable &names) {
  size_t loaded = 0;
  for (const exactnaming::TableRow &row : exactnaming::readTableRows(path))
    if (auto g4 = exactnaming::parseTableG4(row.g4)) {
      names.addLiterature(row.name, g4->first, g4->second);
      ++loaded;
    }
  return loaded;
}

} // namespace witnessstore
