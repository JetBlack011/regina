//
//  linkingnumber_test.cpp
//
//  linkingnumber::linkingNumber() (see ../linkingnumber.h) against two
//  independent references:
//
//    1. Diagrams. buildLink() builds a 2-component link into a triangulated
//       S^3; its components' linking number must equal the diagram's, which
//       Regina reads straight off the PD code's crossing signs
//       (Link::linking()). Covers |lk| = 0, 1, 2, 3.
//    2. The old route. Random vertex-disjoint cycles in those triangulations,
//       compared with EdgeComplement::linkingNumberWith() (drilling plus
//       homology), which verifyslicegenus used until now.
//
//  A null answer (the method declining; see linkingnumber.h) is allowed and
//  counted, never a wrong one. The old route costs ~1 s a call, so ctest runs
//  a sample of test 2 and checks rerouted curves against the diagrams only;
//  --full compares with the old route everywhere (~10 minutes). An optional
//  link table CSV of Name,PD,... rows sweeps every 2-component row of it
//  with test 1 (not part of ctest):
//
//    ./linkingnumber_test [--full] [../../../../cobordism-atlas/data/links_..._pd_codes.csv]
//
//  The heavy check on real vertex links is `cobound run` with audit_linking = 1.
//

#include <algorithm>
#include <chrono>
#include <fstream>
#include <optional>
#include <iostream>
#include <random>
#include <set>
#include <sstream>
#include <string>
#include <unistd.h>
#include <vector>

#include <link/link.h>
#include <triangulation/dim3.h>

#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/linkcomplement.h"
#include "surfer/submanifold/linkingnumber.h"

static int passed = 0, failed_count = 0;

namespace {
bool colorEnabled() {
    static bool enabled = isatty(fileno(stdout));
    return enabled;
}
std::ostream &green(std::ostream &os) {
    return colorEnabled() ? os << "\033[32m" : os;
}
std::ostream &red(std::ostream &os) {
    return colorEnabled() ? os << "\033[31m" : os;
}
std::ostream &bold(std::ostream &os) {
    return colorEnabled() ? os << "\033[1m" : os;
}
std::ostream &resetColor(std::ostream &os) {
    return colorEnabled() ? os << "\033[0m" : os;
}
} // namespace

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n";  \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << red << "  FAIL: " << (desc) << "\n"                   \
                      << "        expected " << _e << ", got " << _a           \
                      << resetColor << "\n";                                   \
            ++failed_count;                                                    \
        }                                                                      \
    } while (0)

namespace {

struct Case {
    std::string name;
    std::string pd;
    long expected; // |lk| as tabulated, checked against Regina too
};

// The diagram as Regina's own parser reads the table's PD string, so the
// expected answer never passes through buildLink(). That parser takes the
// table's 1-based labels as they are, but not its ';' separators.
regina::Link reginaLink(std::string pd) {
    std::replace(pd.begin(), pd.end(), ';', ',');
    return regina::Link::fromPD(pd);
}

// |lk| of a 2-component link, from the diagram alone.
long diagramLinking(const std::string &pd) {
    long lk = reginaLink(pd).linking();
    return lk < 0 ? -lk : lk;
}

// Both components of the drawn link, or nothing unless there are exactly 2.
struct Drawn {
    diagramtriangulation::TriangulationWithLink built;
    std::vector<const regina::Edge<3> *> a, b;
};
std::optional<Drawn> draw(const diagramtriangulation::PDCode &pd) {
    Drawn d{diagramtriangulation::buildLink(pd), {}, {}};
    Link link(d.built.tri, d.built.edges);
    if (link.countComponents() != 2)
        return std::nullopt;
    d.a = link.comps_[0].edges();
    d.b = link.comps_[1].edges();
    return d;
}

int declined = 0; // null answers, allowed but counted

// --full: every comparison with the old (drilling) route. It costs ~1 s a
// call, ~10 minutes in all -- too slow for ctest, which runs a sample.
bool fullComparison = false;

// Test 1 on one diagram; returns false on a wrong answer.
bool checkDiagram(const std::string &name, const std::string &pdString,
                  std::optional<long> tabulated, bool verbose) {
    const long expected = diagramLinking(pdString);
    const diagramtriangulation::PDCode pd = diagramtriangulation::parsePDCode(pdString);
    if (tabulated)
        EXPECT_EQ(expected, *tabulated,
                  name + ": Regina's diagram linking number is the tabulated one");
    auto d = draw(pd);
    if (!d) {
        if (verbose)
            std::cout << "  (" << name << ": not 2 components as drawn)\n";
        return true;
    }
    linkingnumber::Complex cx(d->built.tri);
    EXPECT_EQ(cx.valid(), true, name + ": the drawn S^3 is closed and oriented");
    std::optional<long> ab = linkingnumber::linkingNumber(cx, d->a, d->b);
    std::optional<long> ba = linkingnumber::linkingNumber(cx, d->b, d->a);
    for (const auto &got : {ab, ba}) {
        if (!got) {
            ++declined;
            continue;
        }
        if (*got != expected) {
            EXPECT_EQ(*got, expected, name + ": |lk| from the triangulation");
            return false;
        }
        ++passed;
    }
    if (verbose)
        std::cout << "  " << name << ": |lk| " << expected << ", computed "
                  << (ab ? std::to_string(*ab) : "declined") << " / "
                  << (ba ? std::to_string(*ba) : "declined") << " ("
                  << d->built.tri.size() << " tetrahedra)\n";
    return true;
}

const std::vector<Case> CASES = {
    {"L2a1{0} (Hopf)", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]", 1},
    {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]",
     2},
    {"L5a1{0} (Whitehead)",
     "PD[X[6; 1; 7; 2]; X[10; 7; 5; 8]; X[4; 5; 1; 6]; X[2; 10; 3; 9]; "
     "X[8; 4; 9; 3]]",
     0},
    {"L6a3{0}",
     "PD[X[8; 1; 9; 2]; X[2; 9; 3; 10]; X[10; 3; 11; 4]; X[12; 5; 7; 6]; "
     "X[6; 7; 1; 8]; X[4; 11; 5; 12]]",
     3},
    {"L7n1{0}",
     "PD[X[6; 1; 7; 2]; X[12; 7; 13; 8]; X[4; 13; 1; 14]; X[5; 10; 6; 11]; "
     "X[3; 8; 4; 9]; X[9; 14; 10; 5]; X[11; 2; 12; 3]]",
     2},
};

void test_diagrams() {
    for (const Case &c : CASES)
        checkDiagram(c.name, c.pd, c.expected, true);
}

// A random simple cycle in tri's 1-skeleton avoiding `avoid` (vertex
// indices): a random walk, loop-erased at its first return to a vertex it
// has visited.
std::vector<const regina::Edge<3> *>
randomCycle(const regina::Triangulation<3> &tri, const std::set<size_t> &avoid,
            std::mt19937 &rng) {
    std::vector<std::vector<const regina::Edge<3> *>> at(tri.countVertices());
    for (const regina::Edge<3> *e : tri.edges()) {
        if (e->vertex(0) == e->vertex(1))
            continue; // loops make degenerate cycles; not what petals are
        at[e->vertex(0)->index()].push_back(e);
        at[e->vertex(1)->index()].push_back(e);
    }
    for (int attempt = 0; attempt < 200; ++attempt) {
        size_t v = std::uniform_int_distribution<size_t>(
            0, tri.countVertices() - 1)(rng);
        if (avoid.count(v))
            continue;
        std::vector<size_t> path{v};
        std::vector<const regina::Edge<3> *> edges;
        const regina::Edge<3> *last = nullptr;
        for (int step = 0; step < 400; ++step) {
            std::vector<const regina::Edge<3> *> options;
            for (const regina::Edge<3> *e : at[path.back()]) {
                size_t w = e->vertex(0)->index() == path.back()
                               ? e->vertex(1)->index()
                               : e->vertex(0)->index();
                if (e != last && !avoid.count(w))
                    options.push_back(e);
            }
            if (options.empty())
                break;
            const regina::Edge<3> *e = options[std::uniform_int_distribution<
                size_t>(0, options.size() - 1)(rng)];
            size_t w = e->vertex(0)->index() == path.back()
                           ? e->vertex(1)->index()
                           : e->vertex(0)->index();
            auto seen = std::find(path.begin(), path.end(), w);
            if (seen != path.end()) {
                size_t from = seen - path.begin();
                std::vector<const regina::Edge<3> *> cycle(
                    edges.begin() + from, edges.end());
                cycle.push_back(e);
                if (cycle.size() >= 3)
                    return cycle;
                break;
            }
            path.push_back(w);
            edges.push_back(e);
            last = e;
        }
    }
    return {};
}

void test_random_against_old_route() {
    std::mt19937 rng(20260928);
    int compared = 0, nonzero = 0;
    // The sample: the three smallest triangulations, 8 pairs each.
    const size_t cases = fullComparison ? CASES.size() : 3;
    const int trials = fullComparison ? 12 : 8;
    for (size_t k = 0; k < cases; ++k) {
        const Case &c = CASES[k];
        auto d = draw(diagramtriangulation::parsePDCode(c.pd));
        if (!d)
            continue;
        const regina::Triangulation<3> &tri = d->built.tri;
        linkingnumber::Complex cx(tri);
        for (int trial = 0; trial < trials; ++trial) {
            auto a = randomCycle(tri, {}, rng);
            if (a.empty())
                continue;
            std::set<size_t> used;
            for (const regina::Edge<3> *e : a) {
                used.insert(e->vertex(0)->index());
                used.insert(e->vertex(1)->index());
            }
            auto b = randomCycle(tri, used, rng);
            if (b.empty())
                continue;
            const long old = Knot(tri, a).linkingNumberWith(Knot(tri, b));
            const std::optional<long> fast = linkingnumber::linkingNumber(cx, a, b);
            if (!fast) {
                ++declined;
                continue;
            }
            ++compared;
            if (old != 0)
                ++nonzero;
            if (*fast != old)
                EXPECT_EQ(*fast, old,
                          c.name + ": random pair, new route agrees with old");
            else
                ++passed;
        }
    }
    std::cout << "  compared " << compared << " random pairs (" << nonzero
              << " linked)\n";
    EXPECT_EQ(compared >= 20, true, "enough random pairs were compared");
}

// Linked pairs, the case random cycles almost never produce: a table link's
// component B, rerouted again and again across triangles, against its
// partner A. Replacing an edge (u, w) of B by the other two sides of a
// triangle (u, x, w), with x on neither curve, cannot change lk(A, B) --
// A lies in the 1-skeleton and cannot meet the triangle -- so every
// rerouting must still give the diagram's value, by both routes.
void test_rerouted_components() {
    std::mt19937 rng(9281);
    int compared = 0;
    for (const Case &c : CASES) {
        auto d = draw(diagramtriangulation::parsePDCode(c.pd));
        if (!d)
            continue;
        const regina::Triangulation<3> &tri = d->built.tri;
        linkingnumber::Complex cx(tri);
        std::set<size_t> onA;
        for (const regina::Edge<3> *e : d->a) {
            onA.insert(e->vertex(0)->index());
            onA.insert(e->vertex(1)->index());
        }
        std::vector<const regina::Edge<3> *> b = d->b;
        for (int round = 0; round < 6; ++round) {
            for (int move = 0; move < 5; ++move) {
                std::set<size_t> onB;
                for (const regina::Edge<3> *e : b) {
                    onB.insert(e->vertex(0)->index());
                    onB.insert(e->vertex(1)->index());
                }
                const size_t slot =
                    std::uniform_int_distribution<size_t>(0, b.size() - 1)(rng);
                const regina::Edge<3> *e = b[slot];
                bool rerouted = false;
                for (const auto &emb : e->embeddings()) {
                    const regina::Tetrahedron<3> *tet = emb.simplex();
                    const regina::Perm<4> v = emb.vertices();
                    for (int k : {v[2], v[3]}) {
                        const size_t x = tet->vertex(k)->index();
                        if (onA.count(x) || onB.count(x))
                            continue;
                        const regina::Edge<3> *ux = tet->edge(v[0], k);
                        const regina::Edge<3> *xw = tet->edge(k, v[1]);
                        if (ux == xw)
                            continue;
                        b.erase(b.begin() + static_cast<std::ptrdiff_t>(slot));
                        b.push_back(ux);
                        b.push_back(xw);
                        rerouted = true;
                        break;
                    }
                    if (rerouted)
                        break;
                }
            }
            const long expected = diagramLinking(c.pd);
            const std::optional<long> fast = linkingnumber::linkingNumber(cx, d->a, b);
            if (fullComparison) {
                const long old =
                    Knot(tri, d->a).linkingNumberWith(Knot(tri, b));
                EXPECT_EQ(old, expected, c.name + ": rerouted B, old route "
                                                  "keeps the diagram's |lk|");
            }
            if (!fast) {
                ++declined;
                continue;
            }
            ++compared;
            EXPECT_EQ(*fast, expected,
                      c.name + ": rerouted B (" + std::to_string(b.size()) +
                          " edges), new route keeps the diagram's |lk|");
        }
    }
    EXPECT_EQ(compared > 20, true, "enough rerouted pairs were compared");
}

void test_refuses_non_cycles() {
    auto d = draw(diagramtriangulation::parsePDCode(CASES[1].pd));
    linkingnumber::Complex cx(d->built.tri);
    std::vector<const regina::Edge<3> *> open(d->a.begin(), d->a.end() - 1);
    EXPECT_EQ(linkingnumber::linkingNumber(cx, open, d->b).has_value(), false,
              "an open arc is refused, not given a number");
    EXPECT_EQ(linkingnumber::linkingNumber(cx, d->a, d->a).has_value(), false,
              "curves sharing vertices are refused");
}

void sweep(const std::string &path) {
    std::ifstream in(path);
    std::string line;
    std::getline(in, line); // header
    int rows = 0, twoComponent = 0, wrong = 0;
    const auto start = std::chrono::steady_clock::now();
    while (std::getline(in, line)) {
        const size_t comma = line.find(',');
        const size_t next = line.find(',', comma + 1);
        if (comma == std::string::npos || next == std::string::npos)
            continue;
        const std::string name = line.substr(0, comma);
        // Every orientation variant is the same unoriented link; one will do.
        if (name.find("{0}") == std::string::npos)
            continue;
        ++rows;
        const std::string pd = line.substr(comma + 1, next - comma - 1);
        if (reginaLink(pd).countComponents() != 2)
            continue;
        ++twoComponent;
        if (!checkDiagram(name, pd, std::nullopt, false))
            ++wrong;
    }
    std::cout << "  swept " << twoComponent << " two-component links of "
              << rows << " rows in "
              << std::chrono::duration<double>(std::chrono::steady_clock::now() -
                                               start)
                     .count()
              << " s: " << wrong << " wrong\n";
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main(int argc, char **argv) {
    std::string table;
    for (int i = 1; i < argc; ++i) {
        if (std::string(argv[i]) == "--full")
            fullComparison = true;
        else
            table = argv[i];
    }
    run("diagrams", test_diagrams);
    run("random_against_old_route", test_random_against_old_route);
    run("rerouted_components", test_rerouted_components);
    run("refuses_non_cycles", test_refuses_non_cycles);
    if (!table.empty()) {
        std::cout << bold << "\n=== sweep " << table << " ===" << resetColor
                  << "\n";
        sweep(table);
    }
    std::cout << "  declined (null answers): " << declined << "\n";
    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
