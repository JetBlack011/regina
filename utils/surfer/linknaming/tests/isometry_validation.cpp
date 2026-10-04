//
//  isometry_validation.cpp
//
//  Whole-table validation of the isometry step of linknaming::ExactNamer
//  (../snappeaisometry.h), as the drawer was validated on every table row.
//  Not a ctest target: it takes minutes on the full tables.
//
//  isometry_validation --knots K --links L --symmetry S --out PREFIX
//                      [--threads N] [--stride K]
//                      [--no-search]
//
//  1. census    every table entry's complement in the SnapPea kernel:
//               hyperbolic or not, and its volume (PREFIX.census.tsv, to
//               compare with KnotInfo/LinkInfo).
//  2. table     every entry, under every transform (as written, mirrored,
//               reversed, and each component c >= 1 reversed), scrambled by
//               seeded random Reidemeister moves (orientation kept) into a
//               different diagram of the same oriented link, and named three
//               ways, the last two by identify() on the raw scrambled piece:
//                 ref     the full namer on the UNscrambled transformed
//                         diagram (an exact diagram match: the reference);
//                 iso     the scrambled diagram, isometry step only;
//                 search  the scrambled diagram, Reidemeister searches only.
//               iso and search must each equal ref, or be a proved set of
//               alternatives containing it, or (iso, non-hyperbolic) stay
//               untabulated; anything else is DIFFERENT (PREFIX.table.tsv).
//  3. pairs     negatives: every pair of distinct table links sharing a
//               HOMFLY polynomial (in either mirror), and every pair with
//               isometric complements (found by volume), must NOT be the
//               same link. Positives: every orientation variant of a link
//               (a different table diagram) must be the same link as its
//               first variant (PREFIX.pairs.tsv).
//  4. threads   the iso namings of part 2 redone on one thread, with a fresh
//               namer: identical names.
//

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <mutex>
#include <random>
#include <optional>
#include <set>
#include <sstream>
#include <string>
#include <thread>
#include <vector>

#include <link/link.h>

#include "linknaming/linknamer.h"
#include "linknaming/tables.h"
#include "linknaming/diagrams/gaussdiagram.h"
#include "linknaming/isometry/isometry.h"

using namespace linknaming;

namespace {

double seconds(std::chrono::steady_clock::time_point t0) {
    return std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
}

// A different diagram of the same ORIENTED link: seeded random Reidemeister
// moves, which act on the diagram in place and keep every component's
// orientation (rewrite() would not: it identifies diagrams up to reversing
// any of their components). R2 moves between random strands add crossings
// and break alternation, R3 moves then change the projection, and an R1 kink
// goes in last. Fails soft: whatever moves succeed are kept.
regina::Link scramble(const regina::Link &l, unsigned seed, std::string &how) {
    regina::Link d(l);
    if (d.size() == 0) {
        how = "empty";
        return d;
    }
    std::mt19937 rng(seed);
    std::uniform_int_distribution<int> side(0, 1);
    auto strand = [&]() {
        std::uniform_int_distribution<size_t> c(0, d.size() - 1);
        return regina::StrandRef(d.crossing(c(rng)), side(rng));
    };
    int r1 = 0, r2 = 0, r3 = 0;
    for (int k = 0; k < 2; ++k) // two R2 moves (four crossings)
        for (int tries = 0; tries < 200; ++tries)
            if (d.r2(strand(), side(rng), strand(), side(rng))) {
                ++r2;
                break;
            }
    for (int k = 0; k < 400 && r3 < 8; ++k) // up to eight R3 moves
        if (d.r3(strand(), side(rng)))
            ++r3;
    if (d.r1(strand(), side(rng), side(rng) ? 1 : -1))
        ++r1;
    how = "r1x" + std::to_string(r1) + " r2x" + std::to_string(r2) + " r3x" + std::to_string(r3);
    return d;
}

GaussDiagram gaussOf(const regina::Link &l) {
    std::vector<size_t> origin(l.countComponents());
    for (size_t i = 0; i < origin.size(); ++i)
        origin[i] = i;
    return GaussDiagram::of(l, origin);
}

std::set<std::string> alternatives(const std::string &name) {
    std::set<std::string> out;
    std::stringstream in(name);
    std::string a;
    while (std::getline(in, a, '|'))
        out.insert(a);
    return out;
}

// How `got` compares with the reference name.
std::string verdict(const std::string &ref, const std::string &got) {
    if (got == ref)
        return "same";
    if (got.rfind("diagram:", 0) == 0)
        return "untabulated";
    const auto r = alternatives(ref), g = alternatives(got);
    if (std::includes(g.begin(), g.end(), r.begin(), r.end()))
        return "unpinned";
    return "DIFFERENT";
}

std::string byOf(const PieceName &p) {
    switch (p.by) {
        case PieceName::By::exactDiagram: return "diagram";
        case PieceName::By::isometry: return "isometry";
        case PieceName::By::searchAndInvariants: return "search";
        default: return "untabulated";
    }
}

struct Case {
    size_t entry;
    std::string transform;
    regina::Link diagram;   // transformed, unscrambled
    regina::Link scrambled;
    std::string scrambleHow;
    bool nontrivial = false; // differs from the table diagram, even up to reflection and reversal
    std::string ref, iso, isoBy, search, searchBy;
};

template <typename F>
void parallel(size_t n, unsigned threads, F &&f) {
    std::atomic<size_t> next{0};
    std::vector<std::thread> pool;
    for (unsigned t = 0; t < threads; ++t)
        pool.emplace_back([&] {
            for (size_t i = next++; i < n; i = next++)
                f(i);
        });
    for (auto &th : pool)
        th.join();
}

} // namespace

int main(int argc, char **argv) {
    std::string knots, links, symmetry, out = "isometry_validation";
    unsigned threads = std::max(1u, std::thread::hardware_concurrency());
    size_t stride = 1;
    bool withSearch = true;
    for (int i = 1; i < argc; ++i) {
        const std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) { std::cerr << a << " needs a value\n"; std::exit(2); }
            return argv[++i];
        };
        if (a == "--knots") knots = next();
        else if (a == "--links") links = next();
        else if (a == "--symmetry") symmetry = next();
        else if (a == "--out") out = next();
        else if (a == "--threads") threads = static_cast<unsigned>(std::stoul(next()));
        else if (a == "--stride") stride = std::stoul(next());
        else if (a == "--no-search") withSearch = false;
        else { std::cerr << "unknown option " << a << "\n"; return 2; }
    }
    if (knots.empty() || links.empty()) {
        std::cerr << "usage: isometry_validation --knots K --links L [--symmetry S] --out PREFIX\n"
                     "       [--threads N] [--stride K] [--no-search]\n";
        return 2;
    }

    const auto t0 = std::chrono::steady_clock::now();
    const Tables tables = Tables::load(knots, links, symmetry);
    const std::vector<TableEntry> &entries = tables.entries();
    std::vector<size_t> chosen;
    for (size_t i = 0; i < entries.size(); i += stride)
        chosen.push_back(i);
    std::cout << "tables: " << entries.size() << " entries, " << chosen.size()
              << " chosen (stride " << stride << "), " << threads << " threads\n";

    // ---- 1. census -------------------------------------------------------
    std::vector<std::unique_ptr<KernelLink>> kernel(entries.size());
    {
        const auto tc = std::chrono::steady_clock::now();
        parallel(chosen.size(), threads, [&](size_t k) {
            kernel[chosen[k]] = std::make_unique<KernelLink>(entries[chosen[k]].diagram);
        });
        std::ofstream f(out + ".census.tsv");
        f << "name\tbase\tcrossings\tcomponents\thyperbolic\tvolume\n";
        size_t hyp = 0;
        for (size_t i : chosen) {
            const TableEntry &e = entries[i];
            const KernelLink &k = *kernel[i];
            hyp += k.hyperbolic();
            f << e.name << '\t' << e.base << '\t' << e.diagram.size() << '\t' << e.components
              << '\t' << (k.hyperbolic() ? 1 : 0) << '\t';
            f.precision(12);
            f << k.volume() << '\n';
        }
        std::cout << "census: " << hyp << " of " << chosen.size() << " hyperbolic, "
                  << seconds(tc) << " s\n";
    }

    // ---- 2. table --------------------------------------------------------
    NamerLimits fullL;
    NamerLimits isoL;
    isoL.exactDiagram = false;
    isoL.searchHeight = -1;
    isoL.tableSideHeight = -1;
    isoL.deepHeight = -1;
    NamerLimits searchL;
    searchL.exactDiagram = false;
    searchL.isometry = false;
    const LinkNamer fullN(tables, fullL), isoN(tables, isoL), searchN(tables, searchL);

    std::vector<Case> cases;
    for (size_t i : chosen) {
        const TableEntry &e = entries[i];
        auto add = [&](std::string label, regina::Link l) {
            cases.push_back(Case{i, std::move(label), std::move(l), {}, {}, false, {}, {}, {}, {}, {}});
        };
        add("as-written", e.diagram);
        regina::Link m(e.diagram);
        m.changeAll();
        add("mirror", std::move(m));
        regina::Link r(e.diagram);
        r.reverse();
        add("reversed", std::move(r));
        for (size_t c = 1; c < e.diagram.countComponents(); ++c) {
            regina::Link rc(e.diagram);
            rc.reverse(rc.component(c));
            add("reverse-" + std::to_string(c), std::move(rc));
        }
    }
    {
        const auto tt = std::chrono::steady_clock::now();
        std::atomic<size_t> done{0};
        parallel(cases.size(), threads, [&](size_t k) {
            Case &c = cases[k];
            c.scrambled = scramble(c.diagram, static_cast<unsigned>(k) * 2654435761u + 1u,
                                   c.scrambleHow);
            c.nontrivial = c.scrambled.sig<2>(true, true, true) !=
                           entries[c.entry].diagram.sig<2>(true, true, true);
            c.ref = fullN.name(c.diagram).name;
            // identify() on the scrambled piece itself: no simplify() first,
            // so the complement is built from the scrambled diagram.
            const PieceName iso = isoN.namePiece(gaussOf(c.scrambled));
            c.iso = iso.display();
            c.isoBy = byOf(iso);
            if (withSearch) {
                const PieceName sp = searchN.namePiece(gaussOf(c.scrambled));
                c.search = sp.display();
                c.searchBy = byOf(sp);
            }
            const size_t d = ++done;
            if (d % 5000 == 0)
                std::cerr << "  table: " << d << "/" << cases.size() << ", " << seconds(tt) << " s\n";
        });
        std::ofstream f(out + ".table.tsv");
        f << "entry\ttransform\tscramble\tnontrivial\tcrossings\tref\tiso\tiso_by\tiso_verdict\tsearch"
             "\tsearch_by\tsearch_verdict\thyperbolic\n";
        std::map<std::string, size_t> isoV, searchV, isoByHyp;
        for (const Case &c : cases) {
            const std::string iv = verdict(c.ref, c.iso);
            const std::string sv = withSearch ? verdict(c.ref, c.search) : "-";
            const bool hyp = kernel[c.entry]->hyperbolic();
            ++isoV[iv + (hyp ? " (hyperbolic)" : " (not hyperbolic)")];
            ++searchV[sv];
            f << entries[c.entry].name << '\t' << c.transform << '\t' << c.scrambleHow << '\t'
              << c.nontrivial << '\t' << c.scrambled.size() << '\t' << c.ref << '\t' << c.iso << '\t' << c.isoBy << '\t'
              << iv << '\t' << c.search << '\t' << c.searchBy << '\t' << sv << '\t' << hyp << '\n';
        }
        size_t nontrivial = 0;
        for (const Case &c : cases)
            nontrivial += c.nontrivial;
        std::cout << "table: " << cases.size() << " cases (" << nontrivial
                  << " scrambled to a diagram other than the table's), " << seconds(tt) << " s\n";
        for (const LinkNamer *n : {&fullN, &isoN, &searchN})
            for (const std::string &conflict : n->classConflicts())
                std::cout << "  CLASS CONFLICT (literature g4 differs): " << conflict << "\n";
        for (const auto &[v, n] : isoV)
            std::cout << "  iso    " << v << ": " << n << "\n";
        if (withSearch)
            for (const auto &[v, n] : searchV)
                std::cout << "  search " << v << ": " << n << "\n";
    }

    // ---- 3. pairs --------------------------------------------------------
    {
        const auto tp = std::chrono::steady_clock::now();
        // One representative entry per base (its first variant), and every
        // variant for the positive check.
        std::map<std::string, std::vector<size_t>> byBase;
        for (size_t i : chosen)
            byBase[entries[i].base].push_back(i);
        // HOMFLY (either mirror) of every chosen entry -> bases.
        std::map<std::string, std::set<std::string>> homflyBases;
        for (size_t i : chosen)
            for (int mirror = 0; mirror < 2; ++mirror) {
                regina::Link l(entries[i].diagram);
                if (mirror)
                    l.reflect();
                homflyBases[l.homfly().str()].insert(entries[i].base);
            }
        std::set<std::pair<std::string, std::string>> homflyPairs;
        for (const auto &[h, bases] : homflyBases)
            for (auto a = bases.begin(); a != bases.end(); ++a)
                for (auto b = std::next(a); b != bases.end(); ++b)
                    homflyPairs.insert({*a, *b});
        // Isometric complements: hyperbolic bases grouped by volume.
        std::vector<std::pair<double, std::string>> vols;
        for (const auto &[b, is] : byBase)
            if (kernel[is.front()]->hyperbolic())
                vols.push_back({kernel[is.front()]->volume(), b});
        std::sort(vols.begin(), vols.end());
        std::set<std::pair<std::string, std::string>> volumePairs;
        for (size_t a = 0; a < vols.size(); ++a)
            for (size_t b = a + 1; b < vols.size() && vols[b].first - vols[a].first < 1e-6; ++b)
                volumePairs.insert(std::minmax(vols[a].second, vols[b].second));

        struct Pair { std::string a, b, why; bool complement = false, link = false; };
        std::vector<Pair> pairs;
        for (const auto &[a, b] : homflyPairs)
            pairs.push_back({a, b, volumePairs.contains({a, b}) ? "homfly+volume" : "homfly"});
        for (const auto &[a, b] : volumePairs)
            if (!homflyPairs.contains({a, b}))
                pairs.push_back({a, b, "volume"});
        // Positives: every variant against the base's first variant.
        std::vector<std::pair<size_t, size_t>> variantPairs;
        for (const auto &[b, is] : byBase)
            for (size_t k = 1; k < is.size(); ++k)
                variantPairs.push_back({is.front(), is[k]});

        parallel(pairs.size(), threads, [&](size_t k) {
            Pair &p = pairs[k];
            const KernelLink &x = *kernel[byBase[p.a].front()], &y = *kernel[byBase[p.b].front()];
            p.complement = x.sameComplementAs(y);
            p.link = x.sameLinkAs(y);
        });
        std::vector<char> variantSame(variantPairs.size(), 0), variantHyp(variantPairs.size(), 0),
            variantOriented(variantPairs.size(), 0);
        parallel(variantPairs.size(), threads, [&](size_t k) {
            const KernelLink &x = *kernel[variantPairs[k].first], &y = *kernel[variantPairs[k].second];
            variantHyp[k] = x.hyperbolic() && y.hyperbolic();
            variantSame[k] = x.sameLinkAs(y);
            variantOriented[k] = x.sameOrientedLinkAs(y);
        });

        std::ofstream f(out + ".pairs.tsv");
        f << "a\tb\twhy\tboth_hyperbolic\tsame_complement\tsame_link\n";
        size_t hyperbolicPairs = 0, complementTwins = 0, falseLinks = 0;
        for (const Pair &p : pairs) {
            const bool hyp = kernel[byBase[p.a].front()]->hyperbolic() &&
                             kernel[byBase[p.b].front()]->hyperbolic();
            hyperbolicPairs += hyp;
            complementTwins += p.complement;
            falseLinks += p.link;
            f << p.a << '\t' << p.b << '\t' << p.why << '\t' << hyp << '\t' << p.complement << '\t'
              << p.link << '\n';
        }
        size_t vHyp = 0, vSame = 0, vOriented = 0, vOrientedG4 = 0;
        for (size_t k = 0; k < variantPairs.size(); ++k) {
            vHyp += variantHyp[k];
            vSame += variantSame[k] && variantHyp[k];
            vOriented += variantOriented[k] && variantHyp[k];
            const TableEntry &a = entries[variantPairs[k].first], &b = entries[variantPairs[k].second];
            if (variantOriented[k] && variantHyp[k] && a.g4 != b.g4)
                ++vOrientedG4;
            f << a.name << '\t' << b.name << "\tvariant\t" << int(variantHyp[k]) << "\t"
              << int(variantOriented[k]) << '\t' << int(variantSame[k]) << '\n';
        }
        std::cout << "pairs: " << homflyPairs.size() << " HOMFLY-twin and " << volumePairs.size()
                  << " equal-volume pairs of distinct bases (" << pairs.size() << " in all, "
                  << hyperbolicPairs << " both hyperbolic): " << complementTwins
                  << " with isometric complements, " << falseLinks
                  << " called the same link (must be 0)\n";
        std::cout << "  variants: " << vSame << " of " << vHyp
                  << " hyperbolic variant pairs are the same link (must be all); " << vOriented
                  << " are the same ORIENTED link up to mirror and global reversal (one name), "
                  << vOrientedG4 << " of them with different literature g4 (must be 0), "
                  << seconds(tp) << " s\n";
    }

    // ---- 4. threads ------------------------------------------------------
    {
        const auto tt = std::chrono::steady_clock::now();
        const LinkNamer fresh(tables, isoL);
        size_t differ = 0;
        std::ofstream f(out + ".threads.tsv");
        f << "entry\ttransform\tthreaded\tsingle\n";
        for (const Case &c : cases) {
            const std::string again = fresh.namePiece(gaussOf(c.scrambled)).display();
            if (again != c.iso) {
                ++differ;
                f << entries[c.entry].name << '\t' << c.transform << '\t' << c.iso << '\t' << again << '\n';
            }
        }
        std::cout << "threads: " << differ << " of " << cases.size()
                  << " iso names differ on one thread (must be 0), " << seconds(tt) << " s\n";
    }
    std::cout << "total " << seconds(t0) << " s\n";
    return 0;
}
