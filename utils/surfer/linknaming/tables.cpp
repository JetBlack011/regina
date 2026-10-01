//
//  exacttables.cpp
//

#include "linknaming/tables.h"

#include <algorithm>
#include <array>
#include <fstream>
#include <stdexcept>

namespace exactnaming {

namespace {

// Rows of a "Name,PD,Genus-4D" table: the name (before the first comma), the
// PD field (between the first and the last) and the 4-genus (after the last).
template <typename Fn>
void eachRow(const std::string &path, Fn &&fn) {
    std::ifstream in(path);
    if (!in) throw regina::InvalidArgument("cannot open " + path);
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        size_t a = line.find(','), b = line.rfind(',');
        if (a == std::string::npos || b <= a) continue;
        fn(line.substr(0, a), line.substr(a + 1, b - a - 1), line.substr(b + 1));
    }
}

Symmetry parseSymmetry(const std::string &s) {
    if (s == "chiral") return Symmetry::chiral;
    if (s == "reversible") return Symmetry::reversible;
    if (s == "positive amphicheiral") return Symmetry::positiveAmphicheiral;
    if (s == "negative amphicheiral") return Symmetry::negativeAmphicheiral;
    if (s == "fully amphicheiral") return Symmetry::fullyAmphicheiral;
    return Symmetry::unknown;
}

} // namespace

regina::Link linkFromTablePD(const std::string &pd) {
    std::vector<std::array<long, 4>> code;
    std::vector<long> n;
    long cur = -1;
    for (char ch : pd + ' ') {
        if (ch >= '0' && ch <= '9') {
            cur = (cur < 0 ? 0 : cur * 10) + (ch - '0');
        } else if (cur >= 0) {
            n.push_back(cur);
            cur = -1;
        }
    }
    if (n.empty() || n.size() % 4 != 0)
        throw regina::InvalidArgument("not a PD code: " + pd);
    for (size_t i = 0; i < n.size(); i += 4)
        code.push_back({n[i], n[i + 1], n[i + 2], n[i + 3]});
    return regina::Link::fromPD(code.begin(), code.end());
}

ExactTables ExactTables::load(const std::string &knotTable, const std::string &linkTable,
                              const std::string &knotSymmetry) {
    ExactTables t;
    auto add = [&](const std::string &name, const std::string &pd, const std::string &g4) {
        TableEntry e;
        e.name = name;
        e.base = name.substr(0, name.find('{'));
        e.diagram = linkFromTablePD(pd);
        e.components = e.diagram.countComponents();
        e.g4 = g4;
        t.byName_.emplace(name, t.entries_.size());
        t.entries_.push_back(std::move(e));
    };
    eachRow(knotTable, add);
    eachRow(linkTable, add);
    if (t.entries_.empty()) throw regina::InvalidArgument("the tables yielded no entries");

    // Indices, only once entries_ no longer moves.
    for (const TableEntry &e : t.entries_) {
        t.byBase_[e.base].push_back(&e);
        const std::string u = e.diagram.sig<2>(true, true, true);
        auto [it, fresh] = t.unoriented_.try_emplace(u, e.base);
        if (!fresh && it->second != e.base)
            throw regina::InvalidArgument("two table bases share one diagram: " + it->second +
                                          " and " + e.base);
        for (int mirror = 0; mirror < 2; ++mirror)
            for (int reverse = 0; reverse < 2; ++reverse) {
                regina::Link v(e.diagram);
                if (mirror) v.reflect();
                if (reverse) v.reverse();
                std::vector<VersionMatch> &slot = t.exact_[v.sig<2>(false, false, true)];
                slot.push_back({&e, mirror == 1, reverse == 1});
            }
    }
    // Classes of entries sharing an exact version: one link up to mirror and
    // global reversal. Union-find over entry indices.
    {
        std::vector<size_t> parent(t.entries_.size());
        for (size_t i = 0; i < parent.size(); ++i) parent[i] = i;
        auto find = [&](size_t x) {
            while (parent[x] != x) x = parent[x] = parent[parent[x]];
            return x;
        };
        for (const auto &[sig, matches] : t.exact_)
            for (const VersionMatch &m : matches)
                parent[find(t.byName_.at(m.entry->name))] = find(t.byName_.at(matches.front().entry->name));
        std::unordered_map<size_t, std::vector<size_t>> classes;
        for (size_t i = 0; i < parent.size(); ++i) classes[find(i)].push_back(i);
        for (const auto &[root, members] : classes) {
            std::string least = t.entries_[members.front()].name;
            bool agree = true;
            for (size_t i : members) {
                least = std::min(least, t.entries_[i].name);
                agree = agree && t.entries_[i].g4 == t.entries_[members.front()].g4;
            }
            for (size_t i : members) t.canonical_[t.entries_[i].name] = least;
            if (!agree) {
                std::string s;
                for (size_t i : members) s += (s.empty() ? "" : " ") + t.entries_[i].name + "=" + t.entries_[i].g4;
                t.inconsistent_.push_back(s);
            }
        }
    }
    if (!knotSymmetry.empty()) {
        std::ifstream in(knotSymmetry);
        if (!in) throw regina::InvalidArgument("cannot open " + knotSymmetry);
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)) {
            size_t a = line.find(','), b = line.find(',', a + 1);
            if (a == std::string::npos) continue;
            t.symmetry_[line.substr(0, a)] =
                parseSymmetry(line.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1));
        }
    }
    return t;
}

const std::vector<VersionMatch> *ExactTables::exact(const std::string &sig) const {
    auto it = exact_.find(sig);
    return it == exact_.end() ? nullptr : &it->second;
}

const std::string *ExactTables::base(const std::string &sig) const {
    auto it = unoriented_.find(sig);
    return it == unoriented_.end() ? nullptr : &it->second;
}

const std::vector<const TableEntry *> &ExactTables::variants(const std::string &base) const {
    static const std::vector<const TableEntry *> none;
    auto it = byBase_.find(base);
    return it == byBase_.end() ? none : it->second;
}

const TableEntry *ExactTables::entry(const std::string &name) const {
    auto it = byName_.find(name);
    return it == byName_.end() ? nullptr : &entries_[it->second];
}

const std::string &ExactTables::canonical(const std::string &name) const {
    auto it = canonical_.find(name);
    return it == canonical_.end() ? name : it->second;
}

Symmetry ExactTables::symmetry(const std::string &knot) const {
    auto it = symmetry_.find(knot);
    return it == symmetry_.end() ? Symmetry::unknown : it->second;
}

} // namespace exactnaming

namespace witnessstore {

// One row of the input knot table (Name,PD Notation,Genus-4D). No RFC-4180
// quoting appears in that file (PD Notation uses ';' internally, never a
// literal comma), so a naive two-comma split suffices -- see the input
// table's own format, confirmed during design.
bool splitInputLine(const std::string &line, std::string &name,
                    std::string &pd, std::string &genusField) {
  size_t c1 = line.find(',');
  if (c1 == std::string::npos)
    return false;
  size_t c2 = line.find(',', c1 + 1);
  if (c2 == std::string::npos)
    return false;
  name = line.substr(0, c1);
  pd = line.substr(c1 + 1, c2 - c1 - 1);
  genusField = line.substr(c2 + 1);
  if (!genusField.empty() && genusField.back() == '\r')
    genusField.pop_back();
  return true;
}

// Parses "N" or "[lo;hi]" into lo/hi (lo == hi in the plain-integer case).
void parseGenusField(const std::string &field, int &lo, int &hi) {
  if (!field.empty() && field.front() == '[') {
    size_t semi = field.find(';');
    lo = std::stoi(field.substr(1, semi - 1));
    hi = std::stoi(field.substr(semi + 1, field.size() - semi - 2));
  } else {
    lo = hi = std::stoi(field);
  }
}

std::vector<LiteratureRow> readLiteratureRows(const std::filesystem::path &path) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open name table: " + path.string());

  std::vector<LiteratureRow> rows;
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    LiteratureRow row;
    std::string genusField;
    if (!splitInputLine(line, row.name, row.pd, genusField))
      continue;
    try {
      parseGenusField(genusField, row.lo, row.hi);
    } catch (const std::exception &) {
      continue;
    }
    rows.push_back(std::move(row));
  }
  return rows;
}

} // namespace witnessstore
