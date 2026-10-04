//
//  tables.cpp
//

#include "linknaming/tables.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <fstream>
#include <map>
#include <stdexcept>
#include <tuple>

#include "linknaming/names.h"

namespace linknaming {

std::vector<TableRow> readTableRows(const std::filesystem::path &path) {
    std::ifstream in(path);
    if (!in) throw regina::InvalidArgument("cannot open " + path.string());
    std::vector<TableRow> rows;
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        if (line.empty()) continue;
        const size_t a = line.find(',');
        if (a == std::string::npos) continue;
        const size_t b = line.find(',', a + 1);
        if (b == std::string::npos) continue;
        rows.push_back({line.substr(0, a), line.substr(a + 1, b - a - 1), line.substr(b + 1)});
    }
    return rows;
}

namespace {
std::optional<int> parseDigits(const std::string &s) {
    if (s.empty()) return std::nullopt;
    for (char c : s)
        if (!std::isdigit(static_cast<unsigned char>(c))) return std::nullopt;
    return std::stoi(s);
}
} // namespace

std::optional<std::pair<int, int>> parseTableG4(const std::string &s) {
    if (s.size() >= 5 && s.front() == '[' && s.back() == ']') {
        const auto semi = s.find(';');
        if (semi == std::string::npos) return std::nullopt;
        auto lo = parseDigits(s.substr(1, semi - 1));
        auto hi = parseDigits(s.substr(semi + 1, s.size() - semi - 2));
        if (!lo || !hi || *lo > *hi) return std::nullopt;
        return std::make_pair(*lo, *hi);
    }
    if (auto v = parseDigits(s)) return std::make_pair(*v, *v);
    return std::nullopt;
}

std::optional<SymmetryType> parseSymmetryType(const std::string &text) {
    if (text == "chiral")
        return SymmetryType::chiral;
    if (text == "reversible")
        return SymmetryType::reversible;
    if (text == "positive amphicheiral")
        return SymmetryType::positiveAmphicheiral;
    if (text == "negative amphicheiral")
        return SymmetryType::negativeAmphicheiral;
    if (text == "fully amphicheiral")
        return SymmetryType::fullyAmphicheiral;
    return std::nullopt;
}

SymmetryTable readSymmetryTable(const std::filesystem::path &path) {
    std::ifstream in(path);
    if (!in) throw regina::InvalidArgument("cannot open " + path.string());
    SymmetryTable table;
    std::string line;
    std::getline(in, line); // header
    while (std::getline(in, line)) {
        if (!line.empty() && line.back() == '\r') line.pop_back();
        const size_t a = line.find(',');
        if (a == std::string::npos) continue;
        const size_t b = line.find(',', a + 1);
        if (auto t = parseSymmetryType(
                line.substr(a + 1, b == std::string::npos ? std::string::npos : b - a - 1)))
            table[line.substr(0, a)] = *t;
    }
    return table;
}

namespace {

/**
 * Composite knots that bound a smooth disk, and so ground a chain exactly as
 * the unknot does.
 *
 * `K # m(K^r)` -- K summed with the reverse of its mirror -- is the identity
 * of the concordance group and bounds an explicit ribbon disk, for every K.
 * Which *spelling* of that qualifies depends on the symmetry of K, so this is
 * an explicit allowlist and not a pattern:
 *
 *   - "3_1#m3_1": 3_1 is invertible, so 3_1^r = 3_1 and m(3_1^r) = m3_1.
 *   - "4_1#4_1":  4_1 is invertible AND amphichiral, so m(4_1^r) = 4_1 and
 *                 the sum of 4_1 with *itself* is already the ribbon case.
 *
 * DO NOT generalize this to the string pattern "A#mA". That is wrong as soon
 * as A is non-invertible (the first such knot is 8_17), where A # mA and
 * A # m(A^r) are different knots and only the latter is slice. Adding an
 * entry means checking the symmetry of the summand first.
 *
 * Names are matched exactly, in the canonical spelling that
 * tools/identify_by_retriangulation.py emits: summands sorted, and the
 * lexicographically smaller of the name and its overall mirror.
 */
bool isSliceComposite(const std::string &name) {
    return name == "3_1#m3_1" || name == "4_1#4_1";
}

} // namespace

bool isElementarySlice(const std::string &name, const SymmetryTable &symmetry) {
    if (isSliceComposite(name))
        return true;
    std::vector<std::string> parts = linknaming::knotSummands(name);
    if (parts.empty())
        return false;
    // Each summand reduced to what its symmetry type leaves meaningful of its
    // marks (m: mirrored, r: reversed), as (knot, m, r) with the meaningless
    // marks cleared; its concordance inverse -K = mrK reduced the same way.
    // Slice when every reduced summand pairs off with its inverse (a
    // self-inverse one with another copy of itself).
    using Key = std::tuple<std::string, bool, bool>;
    auto reduce = [](const std::string &knot, SymmetryType type, bool m, bool r) -> Key {
        switch (type) {
        case SymmetryType::fullyAmphicheiral: return {knot, false, false};
        case SymmetryType::reversible: return {knot, m, false};         // K^r = K
        case SymmetryType::negativeAmphicheiral: return {knot, m != r, false}; // K^r = mK
        case SymmetryType::positiveAmphicheiral: return {knot, false, r};      // mK = K
        default: return {knot, m, r};                                   // chiral
        }
    };
    std::map<Key, int> count;
    std::map<Key, Key> inverse;
    for (const std::string &p : parts) {
        if (p == "Unknot" || p == "mUnknot")
            continue;
        size_t k = 0;
        const bool m = k < p.size() && p[k] == 'm';
        if (m) ++k;
        const bool r = k < p.size() && p[k] == 'r';
        if (r) ++k;
        const std::string knot = p.substr(k);
        auto type = symmetry.find(knot);
        if (type == symmetry.end())
            return false;
        const Key key = reduce(knot, type->second, m, r);
        ++count[key];
        inverse[key] = reduce(knot, type->second, !m, !r);
    }
    for (const auto &[key, n] : count) {
        const Key &inv = inverse.at(key);
        if (inv == key) {
            if (n % 2 != 0)
                return false;
        } else {
            auto it = count.find(inv);
            if (it == count.end() || it->second != n)
                return false;
        }
    }
    return true;
}

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
    // Regina numbers strands from 1; a code holding a 0 is 0-based (a literal
    // 0 can appear in no other), as parsePDCode() also reads it.
    const long shift = std::ranges::find(n, 0L) != n.end() ? 1 : 0;
    for (size_t i = 0; i < n.size(); i += 4)
        code.push_back({n[i] + shift, n[i + 1] + shift, n[i + 2] + shift, n[i + 3] + shift});
    return regina::Link::fromPD(code.begin(), code.end());
}

Tables Tables::load(const std::string &knotTable, const std::string &linkTable,
                              const std::string &knotSymmetry) {
    Tables t;
    auto add = [&](const std::string &name, const std::string &pd, const std::string &g4) {
        TableEntry e;
        e.name = name;
        e.base = linknaming::stripOrientationTag(name);
        e.diagram = linkFromTablePD(pd);
        e.components = e.diagram.countComponents();
        e.g4 = g4;
        t.byName_.emplace(name, t.entries_.size());
        t.entries_.push_back(std::move(e));
    };
    for (const TableRow &row : readTableRows(knotTable)) add(row.name, row.pd, row.g4);
    t.knotEntries_ = t.entries_.size();
    for (const TableRow &row : readTableRows(linkTable)) add(row.name, row.pd, row.g4);
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
    if (!knotSymmetry.empty()) t.symmetry_ = readSymmetryTable(knotSymmetry);
    return t;
}

const std::vector<VersionMatch> *Tables::exact(const std::string &sig) const {
    auto it = exact_.find(sig);
    return it == exact_.end() ? nullptr : &it->second;
}

const std::string *Tables::base(const std::string &sig) const {
    auto it = unoriented_.find(sig);
    return it == unoriented_.end() ? nullptr : &it->second;
}

const std::vector<const TableEntry *> &Tables::variants(const std::string &base) const {
    static const std::vector<const TableEntry *> none;
    auto it = byBase_.find(base);
    return it == byBase_.end() ? none : it->second;
}

const TableEntry *Tables::entry(const std::string &name) const {
    auto it = byName_.find(name);
    return it == byName_.end() ? nullptr : &entries_[it->second];
}

const std::string &Tables::canonical(const std::string &name) const {
    auto it = canonical_.find(name);
    return it == canonical_.end() ? name : it->second;
}

std::optional<SymmetryType> Tables::symmetry(const std::string &knot) const {
    auto it = symmetry_.find(knot);
    if (it == symmetry_.end()) return std::nullopt;
    return it->second;
}

} // namespace linknaming

namespace linknaming {

SignatureTable SignatureTable::fromTables(const std::string &knotTable,
                                          const std::string &linkTable) {
    // A row that does not parse is an error, not a silent gap: an empty or
    // partial table would quietly send every outgoing link back to the
    // complement route (it did, once, when the PD codes were read with
    // 0-based labels).
    // The PD codes as Regina reads them (linknaming::linkFromTablePD()):
    // labels as written, 1..2n. Not knotbuilder::parsePDCode(), which
    // renumbers from 0.
    SignatureTable t;
    if (!knotTable.empty())
        for (const linknaming::TableRow &row : linknaming::readTableRows(knotTable)) {
            t.knots_.try_emplace(linknaming::linkFromTablePD(row.pd).knotSig(true, true),
                                 row.name);
            t.knotNames_.insert(row.name);
        }
    if (!linkTable.empty())
        for (const linknaming::TableRow &row : linknaming::readTableRows(linkTable))
            t.links_.try_emplace(linknaming::linkFromTablePD(row.pd).sig<2>(true, true, true),
                                 linknaming::stripOrientationTag(row.name));
    if ((!knotTable.empty() && t.knots_.empty()) ||
        (!linkTable.empty() && t.links_.empty()))
        throw regina::InvalidArgument("a table yielded no signatures");
    return t;
}

SignatureTable SignatureTable::fromTables(const linknaming::Tables &tables) {
    // Exactly fromTables(knotTable, linkTable) over the same files: each
    // entry's diagram is linkFromTablePD() of its row's PD, in file order,
    // the knot table's first, so every signature and every first-wins name
    // is the same.
    SignatureTable t;
    const std::vector<linknaming::TableEntry> &entries = tables.entries();
    for (size_t i = 0; i < entries.size(); ++i) {
        const linknaming::TableEntry &e = entries[i];
        if (i < tables.knotEntries()) {
            t.knots_.try_emplace(e.diagram.knotSig(true, true), e.name);
            t.knotNames_.insert(e.name);
        } else {
            t.links_.try_emplace(e.diagram.sig<2>(true, true, true), e.base);
        }
    }
    if (t.knots_.empty() || t.links_.empty())
        throw regina::InvalidArgument("a table yielded no signatures");
    return t;
}

const std::string *SignatureTable::knot(const std::string &sig) const {
    auto it = knots_.find(sig);
    return it == knots_.end() ? nullptr : &it->second;
}

const std::string *SignatureTable::link(const std::string &sig) const {
    auto it = links_.find(sig);
    return it == links_.end() ? nullptr : &it->second;
}

} // namespace linknaming
