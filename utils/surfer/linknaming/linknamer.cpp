//
//  exactnamer.cpp
//

#include "linknaming/linknamer.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <functional>
#include <iomanip>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>
#include <thread>

#include "linknaming/census/censusnaming.h"
#include "linknaming/complement/linkcomplement.h"
#include "linknaming/complement/unlinknaming.h"

namespace linknaming {

std::vector<long> linkingNumbers(const GaussDiagram &g) {
    std::vector<long> over(g.crossings(), -1), under(g.crossings(), -1);
    for (size_t c = 0; c < g.components(); ++c)
        for (long v : g.comps[c])
            (v > 0 ? over : under)[static_cast<size_t>(std::abs(v) - 1)] = static_cast<long>(c);
    std::map<std::pair<long, long>, long> twice;
    for (size_t k = 0; k < g.crossings(); ++k)
        if (over[k] != under[k])
            twice[{std::min(over[k], under[k]), std::max(over[k], under[k])}] += g.signs[k];
    std::vector<long> out;
    for (size_t i = 0; i < g.components(); ++i)
        for (size_t j = i + 1; j < g.components(); ++j)
            out.push_back(twice[{static_cast<long>(i), static_cast<long>(j)}] / 2);
    std::sort(out.begin(), out.end());
    return out;
}

std::string PieceName::display() const {
    if (by == By::untabulated) return "diagram:" + sig;
    std::string s;
    for (const std::string &n : names) s += (s.empty() ? "" : "|") + n;
    return s;
}

std::string FarSideName::proof() const {
    size_t exact = 0, isometry = 0, search = 0, untab = 0;
    for (const PieceName &p : pieces)
        (p.by == PieceName::By::exactDiagram           ? exact
         : p.by == PieceName::By::isometry              ? isometry
         : p.by == PieceName::By::searchAndInvariants   ? search
                                                        : untab)++;
    std::ostringstream o;
    o << "diagram " << exact << ", isometry " << isometry << ", search+invariants "
      << search << ", untabulated " << untab;
    return o.str();
}

ExactNamer::ExactNamer(const ExactTables &tables, NamerLimits limits,
                       std::shared_ptr<TableCaches> caches)
    : tables_(tables), limits_(limits),
      caches_(caches ? std::move(caches) : std::make_shared<TableCaches>(tables)) {
    // Every cache is keyed by entries of one ExactTables.
    if (caches_->tables != &tables_)
        throw std::invalid_argument("ExactNamer: table caches built for other tables");
}

const regina::Laurent2<regina::Integer> &ExactNamer::homfly(const TableEntry &e, bool mirror) const {
    std::lock_guard<std::mutex> lock(caches_->cacheMutex);
    auto &cache = caches_->homfly;
    auto it = cache.find({&e, mirror});
    if (it == cache.end()) {
        regina::Link l(e.diagram);
        if (mirror) l.reflect();
        it = cache.emplace(std::make_pair(&e, mirror), l.homfly()).first;
    }
    return it->second;
}

bool ExactNamer::inFlypeOrbit(const TableEntry &e, const std::string &graph) const {
    TableCaches &c = *caches_;
    {
        std::lock_guard<std::mutex> lock(c.cacheMutex);
        if (auto it = c.flypeOrbits.find(&e); it != c.flypeOrbits.end())
            return it->second.contains(graph);
    }
    // Every graph reachable by flypes, compared up to relabelling and
    // reflection of the sphere. Composites of minimal flypes give every
    // flype (ModelLinkGraph::findFlype()), so minimal ones suffice.
    std::set<std::string> orbit;
    std::vector<regina::ModelLinkGraph> todo{e.diagram.graph()};
    orbit.insert(todo.front().canonicalPlantri(true));
    constexpr size_t cap = 100000;
    while (!todo.empty() && orbit.size() < cap) {
        regina::ModelLinkGraph g = std::move(todo.back());
        todo.pop_back();
        for (size_t n = 0; n < g.size(); ++n)
            for (int i = 0; i < 4; ++i) {
                regina::ModelLinkGraphArc from(g.node(n), i);
                auto [left, right] = g.findFlype(from);
                if (!left)
                    continue;
                regina::ModelLinkGraph f = g.flype(from, left, right);
                if (orbit.insert(f.canonicalPlantri(true)).second)
                    todo.push_back(std::move(f));
            }
    }
    std::lock_guard<std::mutex> lock(c.cacheMutex);
    return c.flypeOrbits.emplace(&e, std::move(orbit)).first->second.contains(graph);
}

std::vector<const TableEntry *> ExactNamer::homflyCandidates(
        const regina::Link &l, const regina::Laurent2<regina::Integer> &h) const {
    TableCaches &c = *caches_;
    std::call_once(c.homflyIndexOnce, [this, &c] {
        // Every entry's polynomial and its mirror's (~34,000 through 13
        // crossings) cost ~3.8 s on one thread, paid at the first lookup by
        // each set of table caches. They are independent, so they are found
        // on a pool, cached, and then indexed in table order as before.
        const std::vector<TableEntry> &entries = tables_.entries();
        std::vector<regina::Laurent2<regina::Integer>> found(2 * entries.size());
        std::atomic<size_t> next{0};
        auto work = [&] {
            for (size_t i; (i = next.fetch_add(1)) < found.size();) {
                regina::Link l(entries[i / 2].diagram);
                if (i % 2)
                    l.reflect();
                found[i] = l.homfly();
            }
        };
        std::vector<std::thread> pool;
        for (unsigned t = 1; t < std::max(1u, std::thread::hardware_concurrency()); ++t)
            pool.emplace_back(work);
        work();
        for (std::thread &t : pool)
            t.join();
        {
            std::lock_guard<std::mutex> lock(c.cacheMutex);
            for (size_t i = 0; i < found.size(); ++i)
                c.homfly.try_emplace({&entries[i / 2], i % 2 == 1}, std::move(found[i]));
        }
        for (const TableEntry &e : entries) {
            for (int mirror = 0; mirror < 2; ++mirror) {
                auto &bases = c.homflyIndex[homfly(e, mirror == 1).str()];
                if (std::find(bases.begin(), bases.end(), e.base) == bases.end())
                    bases.push_back(e.base);
            }
        }
    });
    std::vector<const TableEntry *> candidates;
    auto it = c.homflyIndex.find(h.str());
    if (it == c.homflyIndex.end())
        return candidates;
    // A table diagram is minimal, so no bigger than any diagram of its link.
    for (const std::string &b : it->second) {
        const auto &vs = tables_.variants(b);
        if (!vs.empty() && vs.front()->diagram.size() <= l.size())
            candidates.push_back(vs.front());
    }
    // Smaller table diagrams first: cheaper to test, and a HOMFLY coincidence
    // with a bigger entry (8_16 shares 10_156's) is then tried last.
    std::stable_sort(candidates.begin(), candidates.end(),
                     [](const TableEntry *a, const TableEntry *b) {
                         return a->diagram.size() < b->diagram.size();
                     });
    return candidates;
}

const KernelLink &ExactNamer::kernelLinkOf(const TableEntry &e) const {
    TableCaches &c = *caches_;
    {
        std::lock_guard<std::mutex> lock(c.kernelMutex);
        if (auto it = c.kernelLinks.find(&e); it != c.kernelLinks.end())
            return *it->second;
    }
    auto k = std::make_unique<KernelLink>(e.diagram);
    std::lock_guard<std::mutex> lock(c.kernelMutex);
    // Another thread may have built it meanwhile; either copy will do.
    return *c.kernelLinks.try_emplace(&e, std::move(k)).first->second;
}

const std::string &ExactNamer::canonicalName(const TableEntry &e) const {
    TableCaches &c = *caches_;
    {
        std::lock_guard<std::mutex> lock(c.classMutex);
        if (auto it = c.canonicalOf.find(&e); it != c.canonicalOf.end())
            return it->second;
    }
    // The whole base at once: a union-find over its variants, joined when
    // their table diagrams are one link up to mirror and global reversal, or
    // when an isometry carries meridians to meridians with a uniform sign.
    const std::vector<const TableEntry *> &vs = tables_.variants(e.base);
    std::vector<size_t> parent(vs.size());
    std::iota(parent.begin(), parent.end(), size_t(0));
    std::function<size_t(size_t)> find = [&](size_t x) {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
    };
    std::vector<std::string> conflicts;
    for (size_t i = 0; i < vs.size(); ++i)
        for (size_t j = i + 1; j < vs.size(); ++j) {
            if (find(i) == find(j))
                continue;
            bool same = tables_.canonical(vs[i]->name) == tables_.canonical(vs[j]->name);
            if (!same) {
                const KernelLink &a = kernelLinkOf(*vs[i]), &b = kernelLinkOf(*vs[j]);
                same = a.hyperbolic() && b.hyperbolic() && a.sameOrientedLinkAs(b);
                if (same && vs[i]->g4 != vs[j]->g4)
                    conflicts.push_back(vs[i]->name + " (" + vs[i]->g4 + ") = " + vs[j]->name +
                                        " (" + vs[j]->g4 + ")");
            }
            if (same)
                parent[std::max(find(i), find(j))] = std::min(find(i), find(j));
        }
    // Each class named as ExactTables names its first member, so that a
    // class of one keeps the name it always had.
    std::lock_guard<std::mutex> lock(c.classMutex);
    for (size_t i = 0; i < vs.size(); ++i)
        c.canonicalOf.try_emplace(vs[i], tables_.canonical(vs[find(i)]->name));
    c.classConflicts.insert(c.classConflicts.end(), conflicts.begin(), conflicts.end());
    return c.canonicalOf.at(&e);
}

std::vector<std::string> ExactNamer::classConflicts() const {
    std::lock_guard<std::mutex> lock(caches_->classMutex);
    return caches_->classConflicts;
}

std::set<std::string> ExactNamer::invariantSurvivors(const GaussDiagram &piece,
                                                     const regina::Laurent2<regina::Integer> &h,
                                                     const std::string &base,
                                                     std::set<bool> *mirrors) const {
    // Rule out every orientation variant of the base, and its mirror, whose
    // HOMFLY polynomial or linking numbers differ from the piece's (the piece
    // as drawn, so in the surface's orientation). Global reversal changes
    // neither, and a name does not see it.
    const std::vector<long> lk = linkingNumbers(piece);
    std::set<std::string> names;
    for (const TableEntry *e : tables_.variants(base))
        for (int mirror = 0; mirror < 2; ++mirror) {
            if (homfly(*e, mirror == 1) != h) continue;
            std::vector<long> lkE = linkingNumbers(GaussDiagram::of(e->diagram, {}));
            if (mirror)
                for (long &x : lkE) x = -x;
            std::sort(lkE.begin(), lkE.end());
            if (lkE != lk) continue;
            names.insert(canonicalName(*e));
            if (mirrors)
                mirrors->insert(mirror == 1);
        }
    return names;
}

std::optional<ExactNamer::IsometryMatch> ExactNamer::isometryMatch(
        const regina::Link &l, const regina::Laurent2<regina::Integer> &h) const {
    const std::vector<const TableEntry *> candidates = homflyCandidates(l, h);
    if (candidates.empty())
        return std::nullopt; // not a table link: no need for its complement
    // A second and third try on randomised triangulations: canonisation
    // occasionally retriangulates at random, and a miss is not an answer.
    for (int attempt = 0; attempt < 3; ++attempt) {
        const KernelLink piece(l, attempt);
        if (!piece.hyperbolic())
            return std::nullopt; // left to the diagram searches
        for (const TableEntry *e : candidates) {
            if (!piece.sameLinkAs(kernelLinkOf(*e)))
                continue;
            IsometryMatch m;
            m.base = e->base;
            for (const TableEntry *v : tables_.variants(e->base))
                for (const KernelLink::Meridional &iso : piece.meridionalIsometriesTo(kernelLinkOf(*v)))
                    if (iso.uniform()) {
                        m.names.insert(canonicalName(*v));
                        break;
                    }
            if (m.names.size() == 1) {
                // Mirror and reversal relative to the named entry itself.
                if (const TableEntry *named = tables_.entry(*m.names.begin())) {
                    std::set<bool> mir, rev;
                    for (const KernelLink::Meridional &iso :
                         piece.meridionalIsometriesTo(kernelLinkOf(*named)))
                        if (iso.uniform()) {
                            mir.insert(iso.reflects);
                            rev.insert(iso.sign.front() < 0);
                        }
                    if (mir.size() == 1) m.mirror = *mir.begin();
                    if (rev.size() == 1) m.reverse = *rev.begin();
                }
            }
            return m;
        }
    }
    return std::nullopt;
}

std::optional<std::string> ExactNamer::tableSideBase(
        const regina::Link &l, const regina::Laurent2<regina::Integer> &h) const {
    const std::vector<const TableEntry *> candidates = homflyCandidates(l, h);
    if (candidates.empty())
        return std::nullopt;

    // Alternating: the flype orbit of an alternating table diagram of the
    // same size. The graph fixes an alternating diagram up to mirror, and
    // the orientations are pinned later, so a graph match proves the base.
    if (l.isAlternating() && l.size() <= 26) {
        const std::string g = l.graph().canonicalPlantri(true);
        for (const TableEntry *e : candidates)
            if (e->diagram.size() == l.size() && e->diagram.isAlternating() &&
                inFlypeOrbit(*e, g))
                return e->base;
    }

    // Otherwise, outward from the candidates' table diagrams: the piece in
    // every orientation of its components (rewrite() may reverse them all
    // at once and reflect, which sig<2>(true, true, true) forgets too).
    std::set<std::string> want;
    std::vector<regina::Link> orientations{l};
    for (size_t c = 1; c < l.countComponents(); ++c) {
        const size_t n = orientations.size();
        for (size_t k = 0; k < n; ++k) {
            regina::Link r(orientations[k]);
            r.reverse(r.component(c));
            orientations.push_back(std::move(r));
        }
    }
    for (const regina::Link &o : orientations)
        want.insert(o.sig<2>(true, true, true));
    for (int height = 1; height <= limits_.tableSideHeight; ++height)
        for (const TableEntry *e : candidates) {
            // The rewrite() at this height reaches diagrams of at most
            // size + height crossings.
            if (e->diagram.size() + static_cast<size_t>(height) < l.size())
                continue;
            bool hit = false;
            size_t visits = 0;
            regina::Link t(e->diagram);
            t.rewrite(height, 1, nullptr, [&](regina::Link &&d) {
                if (d.size() <= l.size() && want.contains(d.sig<2>(true, true, true))) {
                    hit = true;
                    return true;
                }
                return ++visits >= limits_.tableSideVisits;
            });
            if (hit)
                return e->base;
        }
    return std::nullopt;
}

void ExactNamer::decompose(const GaussDiagram &g, std::vector<GaussDiagram> &primes) const {
    regina::Link l = g.link();
    l.simplify(); // never reflects or reverses
    if (limits_.exhaustiveHeight > 0 && l.size() > 0 && l.size() <= limits_.maxSearchCrossings)
        while (l.size() > 0 && l.simplifyExhaustive(limits_.exhaustiveHeight)) {}
    const GaussDiagram s = GaussDiagram::of(l, g.origin);
    if (s.components() != g.components())
        throw regina::InvalidArgument("simplify() changed the number of components");
    for (const GaussDiagram &part : splitPieces(s)) {
        if (part.crossings() == 0) continue; // unknotted, split: see name()
        if (auto cut = visibleSum(part)) {
            decompose(cut->first, primes);
            decompose(cut->second, primes);
        } else {
            primes.push_back(part);
        }
    }
}

PieceName ExactNamer::identify(const GaussDiagram &piece) const {
    PieceName p;
    p.components = piece.components();
    p.origin = piece.origin;
    const regina::Link l = piece.link();
    p.crossings = l.size();

    // Exact: the diagram IS a version of a table entry.
    auto tryExact = [&](const regina::Link &m) {
        const std::vector<VersionMatch> *ms = tables_.exact(m.sig<2>(false, false, true));
        if (!ms) return false;
        std::set<std::string> names;
        for (const VersionMatch &v : *ms) names.insert(canonicalName(*v.entry));
        p.by = PieceName::By::exactDiagram;
        p.base = ms->front().entry->base;
        p.names.assign(names.begin(), names.end());
        // Mirror and reversal relative to the canonical entry, when the
        // diagram is only one of its versions.
        if (names.size() == 1) {
            std::set<bool> mir, rev;
            for (const VersionMatch &v : *ms)
                if (v.entry->name == *names.begin()) {
                    mir.insert(v.mirror);
                    rev.insert(v.reverse);
                }
            if (mir.size() == 1) p.mirror = *mir.begin();
            if (rev.size() == 1) p.reverse = *rev.begin();
        }
        return true;
    };
    if (limits_.exactDiagram && tryExact(l)) return p;
    for (int t = 0; limits_.exactDiagram && t < limits_.simplifyTries; ++t) {
        regina::Link m(l);
        m.simplify();
        if (tryExact(m)) return p;
    }

    std::optional<std::string> base;
    std::optional<regina::Laurent2<regina::Integer>> homflyOfPiece;
    auto pieceHomfly = [&]() -> const regina::Laurent2<regina::Integer> & {
        if (!homflyOfPiece)
            homflyOfPiece = l.homfly();
        return *homflyOfPiece;
    };

    // Hyperbolic: an isometry of complements carrying meridians to meridians
    // (snappeaisometry.h). Milliseconds, indifferent to which diagram of the
    // piece was drawn -- what the Reidemeister searches below depend on --
    // and its action on the oriented meridians pins the variant too.
    if (limits_.isometry && l.size() > 0) {
        if (auto m = isometryMatch(l, pieceHomfly())) {
            if (limits_.checkIsometryPins && !m->names.empty()) {
                const std::set<std::string> survivors =
                    invariantSurvivors(piece, pieceHomfly(), m->base, nullptr);
                if (!std::includes(survivors.begin(), survivors.end(), m->names.begin(),
                                   m->names.end()))
                    throw std::logic_error("exactnaming: the isometry pins a variant of " +
                                           m->base + " that the invariants rule out");
            }
            if (!m->names.empty()) {
                p.by = PieceName::By::isometry;
                p.base = m->base;
                p.names.assign(m->names.begin(), m->names.end());
                p.mirror = m->mirror;
                p.reverse = m->reverse;
                return p;
            }
            base = m->base; // no variant reached with a uniform sign: never seen
        }
    }

    // Search: a diagram of a table base, up to mirror and orientations.
    if (!base && limits_.searchHeight >= 0 && l.size() <= limits_.maxSearchCrossings) {
        regina::Link m(l);
        if (limits_.exhaustiveHeight > 0)
            while (m.size() > 0 && m.simplifyExhaustive(limits_.exhaustiveHeight)) {}
        if (limits_.exactDiagram && tryExact(m)) return p; // simplifyExhaustive() keeps orientation
        if (const std::string *b = tables_.base(m.sig<2>(true, true, true))) base = *b;
        auto search = [&](int height, size_t cap) {
            size_t visits = 0;
            m.rewrite(height, 1, nullptr, [&](regina::Link &&d) {
                if (const std::string *b = tables_.base(d.sig<2>(true, true, true))) {
                    base = *b;
                    return true;
                }
                return ++visits >= cap;
            });
        };
        if (!base) search(limits_.searchHeight, limits_.searchVisits);
        // From the table's side, before the deep forward search: it reaches
        // what that search cannot (another minimal diagram of the entry, or
        // one stuck above minimal), and costs nothing when no entry shares
        // the piece's HOMFLY polynomial.
        if (!base && limits_.tableSideHeight >= 0 && l.size() > 0 &&
            l.size() <= limits_.maxTableSideCrossings)
            base = tableSideBase(l, pieceHomfly());
        if (!base && limits_.deepHeight > limits_.searchHeight &&
            m.size() <= limits_.maxDeepCrossings)
            search(limits_.deepHeight, limits_.deepVisits);
    } else if (!base && limits_.tableSideHeight >= 0 && l.size() > 0 &&
               l.size() <= limits_.maxTableSideCrossings) {
        base = tableSideBase(l, pieceHomfly());
    }
    if (base) {
        std::set<bool> mir;
        const std::set<std::string> names = invariantSurvivors(piece, pieceHomfly(), *base, &mir);
        if (!names.empty()) {
            p.by = PieceName::By::searchAndInvariants;
            p.base = *base;
            p.names.assign(names.begin(), names.end());
            if (names.size() == 1 && mir.size() == 1) p.mirror = *mir.begin();
            return p;
        }
        // Nothing survives: the search and the invariants contradict each
        // other, which is a bug somewhere. Fall through to untabulated, and
        // say so in the base.
        p.base = "CONTRADICTION:" + *base;
    }
    regina::Link r(l);
    r.reverse();
    p.by = PieceName::By::untabulated;
    p.sig = std::min(l.sig<2>(true, false, true), r.sig<2>(true, false, true));
    return p;
}

namespace {

// A composite knot's name from its summands (every piece a knot). Each
// summand is written with the marks its symmetry type makes meaningful --
// `m` mirror, `r` reversed, relative to the table's diagram -- and the sum
// is spelled with the fewest marks over the four global mirror/reversal
// choices (a name does not see those). An unpinned mark that matters is
// written as every alternative. Returns the alternatives, sorted.
std::vector<std::string> knotSumSpellings(const std::vector<const PieceName *> &summands,
                                          const ExactTables &tables) {
    struct Bits { std::string name; std::optional<SymmetryType> sym; bool m, r; };
    // Every assignment of the unpinned bits.
    std::vector<std::vector<Bits>> assignments{{}};
    for (const PieceName *p : summands) {
        const std::string &k = p->names.front();
        const std::optional<SymmetryType> sym = tables.symmetry(k);
        std::vector<std::vector<Bits>> next;
        for (const auto &a : assignments)
            for (int m = 0; m < 2; ++m) {
                if (p->mirror && *p->mirror != (m == 1)) continue;
                for (int r = 0; r < 2; ++r) {
                    if (p->reverse && *p->reverse != (r == 1)) continue;
                    auto b = a;
                    b.push_back({k, sym, m == 1, r == 1});
                    next.push_back(std::move(b));
                }
            }
        assignments = std::move(next);
    }
    auto token = [](Bits b) {
        switch (b.sym.value_or(SymmetryType::chiral)) {
            case SymmetryType::fullyAmphicheiral: b.m = b.r = false; break;
            case SymmetryType::reversible: b.r = false; break;
            case SymmetryType::positiveAmphicheiral: b.m = false; break;
            case SymmetryType::negativeAmphicheiral: b.m = (b.m != b.r); b.r = false; break;
            default: break; // chiral or unknown: both matter
        }
        return std::string(b.m ? "m" : "") + (b.r ? "r" : "") + b.name;
    };
    std::set<std::string> spellings;
    for (const auto &a : assignments) {
        std::string best;
        size_t bestMarks = SIZE_MAX;
        for (int gm = 0; gm < 2; ++gm)
            for (int gr = 0; gr < 2; ++gr) {
                std::vector<std::string> toks;
                size_t marks = 0;
                for (Bits b : a) {
                    b.m = b.m != (gm == 1);
                    b.r = b.r != (gr == 1);
                    std::string t = token(b);
                    marks += t.size() - b.name.size();
                    toks.push_back(std::move(t));
                }
                std::sort(toks.begin(), toks.end());
                std::string s;
                for (const std::string &t : toks) s += (s.empty() ? "" : "#") + t;
                if (marks < bestMarks || (marks == bestMarks && s < best)) {
                    best = s;
                    bestMarks = marks;
                }
            }
        spellings.insert(best);
    }
    std::vector<std::string> out(spellings.begin(), spellings.end());
    std::sort(out.begin(), out.end(), [](const std::string &a, const std::string &b) {
        auto marks = [](const std::string &s) {
            size_t n = 0;
            for (size_t i = 0; i < s.size(); ++i)
                if ((s[i] == 'm' || s[i] == 'r') && (i == 0 || s[i - 1] == '#' || s[i - 1] == 'm'))
                    ++n;
            return n;
        };
        return std::make_pair(marks(a), a) < std::make_pair(marks(b), b);
    });
    return out;
}

} // namespace

FarSideName ExactNamer::name(const regina::Link &drawn) const {
    const size_t n = drawn.countComponents();
    std::vector<size_t> origin(n);
    std::iota(origin.begin(), origin.end(), 0);
    std::vector<GaussDiagram> primes;
    decompose(GaussDiagram::of(drawn, origin), primes);

    FarSideName out;
    std::vector<bool> covered(n, false);
    for (const GaussDiagram &g : primes) {
        out.pieces.push_back(identify(g));
        for (size_t o : g.origin) covered[o] = true;
    }
    out.splitUnknots = static_cast<size_t>(std::count(covered.begin(), covered.end(), false));
    out.pinned = std::all_of(out.pieces.begin(), out.pieces.end(),
                             [](const PieceName &p) { return p.pinned(); });
    if (primes.empty()) {
        out.name = complement::unlinkName(n);
        out.exact = out.pinned = true;
        out.factors = n;
        return out;
    }

    // Split factors: pieces joined, transitively, by a shared far-side
    // component (the component they are summed along).
    const size_t np = out.pieces.size();
    std::vector<size_t> parent(np);
    std::iota(parent.begin(), parent.end(), 0);
    std::function<size_t(size_t)> find = [&](size_t x) {
        return parent[x] == x ? x : parent[x] = find(parent[x]);
    };
    std::map<size_t, size_t> ownerOf; // far-side component -> a piece using it
    for (size_t i = 0; i < np; ++i)
        for (size_t o : out.pieces[i].origin) {
            auto [it, fresh] = ownerOf.try_emplace(o, i);
            if (!fresh) parent[find(i)] = find(it->second);
        }
    std::map<size_t, std::vector<size_t>> factors;
    for (size_t i = 0; i < np; ++i) factors[find(i)].push_back(i);

    std::vector<std::string> factorNames;
    bool allExact = true;
    for (const auto &[root, members] : factors) {
        std::vector<const PieceName *> ps;
        for (size_t i : members) ps.push_back(&out.pieces[i]);
        if (ps.size() == 1) {
            factorNames.push_back(ps[0]->display());
            allExact = allExact && ps[0]->pinned();
            continue;
        }
        const bool allKnots = std::all_of(ps.begin(), ps.end(),
                                          [](const PieceName *p) { return p->components == 1; });
        const bool allTabulated = std::all_of(ps.begin(), ps.end(), [](const PieceName *p) {
            return p->by != PieceName::By::untabulated && p->names.size() == 1;
        });
        if (allKnots && allTabulated) {
            std::vector<std::string> alts = knotSumSpellings(ps, tables_);
            std::string s;
            for (const std::string &a : alts) s += (s.empty() ? "" : "|") + a;
            factorNames.push_back(s);
            allExact = allExact && alts.size() == 1;
            continue;
        }
        allExact = false;
        if (allKnots) { // a composite knot with an untabulated summand
            std::vector<std::string> toks;
            for (const PieceName *p : ps) toks.push_back(p->display());
            std::sort(toks.begin(), toks.end());
            std::string s;
            for (const std::string &t : toks) s += (s.empty() ? "" : "#") + t;
            factorNames.push_back(s);
            continue;
        }
        // Sums with a link. Sites: far-side components shared by >= 2 pieces.
        std::map<size_t, std::vector<const PieceName *>> sites;
        for (const PieceName *p : ps)
            for (size_t o : p->origin) sites[o].push_back(p);
        std::vector<const PieceName *> links;
        for (const PieceName *p : ps)
            if (p->components > 1) links.push_back(p);
        bool simple = links.size() == 1;
        for (const auto &[o, at] : sites)
            if (at.size() > 2) simple = false;
        if (simple) {
            // Knots summed into components of one link: `K #_? L`.
            std::vector<std::string> terms;
            for (const PieceName *p : ps)
                if (p != links[0]) terms.push_back(p->display() + " #_? ");
            std::sort(terms.begin(), terms.end());
            std::string s;
            for (const std::string &t : terms) s += t;
            factorNames.push_back(s + links[0]->display());
        } else {
            std::vector<std::string> siteStrs;
            for (const auto &[o, at] : sites) {
                if (at.size() < 2) continue;
                std::vector<std::string> toks;
                for (const PieceName *p : at) toks.push_back(p->display() + "[?]");
                std::sort(toks.begin(), toks.end());
                std::string s;
                for (const std::string &t : toks) s += (s.empty() ? "" : " # ") + t;
                siteStrs.push_back(s);
            }
            std::sort(siteStrs.begin(), siteStrs.end());
            std::string s;
            for (const std::string &t : siteStrs) s += (s.empty() ? "" : " ; ") + t;
            factorNames.push_back("#{" + s + "}");
        }
    }
    out.factors = factorNames.size() + out.splitUnknots;
    // A split union is an identity only when all but one factor are split
    // unknots (g_4(L u U) = g_4(L), and L u U is determined by L).
    out.exact = allExact && factorNames.size() == 1;
    for (size_t i = 0; i < out.splitUnknots; ++i) factorNames.push_back("Unknot");
    std::sort(factorNames.begin(), factorNames.end(), [](const std::string &a, const std::string &b) {
        return std::make_pair(a == "Unknot", a) < std::make_pair(b == "Unknot", b);
    });
    for (const std::string &f : factorNames) out.name += (out.name.empty() ? "" : " u ") + f;
    return out;
}

} // namespace linknaming

namespace linknaming {

namespace {

long long microsSince(std::chrono::steady_clock::time_point t) {
    return std::chrono::duration_cast<std::chrono::microseconds>(
               std::chrono::steady_clock::now() - t)
        .count();
}

} // namespace

void NamingStats::noteDuration(long long micros, const char *route,
                               const std::string &name) {
    long long seen = slowestMicros_.load(std::memory_order_relaxed);
    if (micros <= seen) return;
    std::lock_guard<std::mutex> lock(slowestMutex_);
    if (micros <= slowestMicros_.load()) return;
    slowestMicros_.store(micros);
    slowest_ = std::string(route) + ' ' + name;
}

std::string NamingStats::slowest() const {
    std::lock_guard<std::mutex> lock(slowestMutex_);
    return slowest_;
}

std::string NamingStats::summary() const {
    auto secs = [](long long micros) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(1) << micros / 1e6;
        return o.str();
    };
    std::ostringstream o;
    o << calls << " far sides drawn: unknot " << unknots << ", unlink " << unlinks
      << ", table knot " << tableKnots << ", learned knot " << learnedKnots
      << ", table link " << tableLinks << ", other link " << diagramLinks << " (+"
      << jonesLinks << " by Jones), learned link " << learnedLinks
      << "; complement fallbacks " << fallbacks << " (" << learned << " learned, "
      << nonPlanar << " non-planar drawings); exact oriented names " << exactNamed
      << " (+" << exactCacheHits << " cached, " << exactFailed << " failed); diagrams "
      << secs(microsDiagram) << "s, fallbacks " << secs(microsFallback) << "s, exact "
      << secs(microsExact) << "s; slowest " << secs(slowestMicros()) << "s";
    if (const std::string s = slowest(); !s.empty()) o << " (" << s << ")";
    return o.str();
}

LinkNamer::LinkNamer(const SignatureTable &table) : table_(table) {}

std::string LinkNamer::name(const Link &curves,
                            const std::function<DrawnCurves()> &draw) const {
    ++stats_.calls;
    const auto start = std::chrono::steady_clock::now();
    const long long fallbacksBefore = stats_.fallbacks.load();
    std::string out = nameOnce(curves, draw);
    stats_.noteDuration(microsSince(start),
                        stats_.fallbacks.load() != fallbacksBefore ? "complement" : "diagram",
                        out);
    return census::perturbedForTesting(std::move(out));
}

std::string LinkNamer::nameOnce(const Link &curves,
                                const std::function<DrawnCurves()> &draw) const {
    const auto start = std::chrono::steady_clock::now();
    const size_t n = curves.comps_.size();
    std::string key; // what a fallback's answer is remembered under
    try {
        DrawnCurves drawn = draw();
        if (drawn.outcome == DrawnCurves::Outcome::nonPlanar) {
            ++stats_.nonPlanar; // a drawer defect: never name from it, and never learn
        } else if (drawn.outcome == DrawnCurves::Outcome::drawn) {
            const bool someLinking = drawn.someLinking;
            regina::Link simplified = std::move(drawn.diagram);
            simplified.simplify();
            if (simplified.size() == 0) {
                stats_.microsDiagram += microsSince(start);
                if (n == 1) { ++stats_.unknots; return complement::unlinkName(1); }
                ++stats_.unlinks;
                return complement::unlinkName(n);
            }
            if (n == 1) {
                std::string sig = simplified.knotSig(true, true);
                if (const std::string *hit = table_.knot(sig)) {
                    ++stats_.tableKnots;
                    stats_.microsDiagram += microsSince(start);
                    return *hit;
                }
                key = "K" + sig;
            } else {
                std::string sig = simplified.sig<2>(true, true, true);
                if (const std::string *hit = table_.link(sig)) {
                    ++stats_.tableLinks;
                    stats_.microsDiagram += microsSince(start);
                    return *hit;
                }
                key = "L" + sig;
                if (someLinking) {
                    ++stats_.diagramLinks;
                    stats_.microsDiagram += microsSince(start);
                    return "diagram:" + sig;
                }
            }
            {
                std::lock_guard<std::mutex> lock(learnedMutex_);
                if (auto it = learned_.find(key); it != learned_.end()) {
                    ++(n == 1 ? stats_.learnedKnots : stats_.learnedLinks);
                    stats_.microsDiagram += microsSince(start);
                    return it->second;
                }
            }
            if (n > 1) {
                // Every linking number is zero. A Jones polynomial other than the
                // unlink's still proves this is not an unlink.
                regina::Laurent<regina::Integer> unlink;
                {
                    std::lock_guard<std::mutex> lock(learnedMutex_);
                    auto it = unlinkJones_.find(n);
                    if (it == unlinkJones_.end())
                        it = unlinkJones_.emplace(n, regina::Link(n).jones()).first;
                    unlink = it->second;
                }
                if (simplified.jones() != unlink) {
                    ++stats_.jonesLinks;
                    stats_.microsDiagram += microsSince(start);
                    return "diagram:" + key.substr(1);
                }
            }
        }
    } catch (const regina::InvalidArgument &) {
        key.clear();
    }
    stats_.microsDiagram += microsSince(start);

    // The complement route, as before this namer existed; its answer is
    // remembered against the diagram.
    const auto fb = std::chrono::steady_clock::now();
    ++stats_.fallbacks;
    std::string name = census::identify(curves);
    stats_.microsFallback += microsSince(fb);
    if (!key.empty()) {
        std::lock_guard<std::mutex> lock(learnedMutex_);
        if (learned_.try_emplace(key, name).second) ++stats_.learned;
    }
    return name;
}

} // namespace linknaming
