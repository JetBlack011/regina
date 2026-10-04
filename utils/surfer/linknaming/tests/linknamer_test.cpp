//
//  linknamer_test.cpp
//
//  linknaming::LinkNamer on diagrams built by hand from table PD codes:
//
//    1. every table entry drawn as itself is named as itself (canonically),
//       exactly;
//    2. reversing one component of L7n1{0} (g_4 2) gives L7n1{1} (g_4 0),
//       and back: the orientation variant is pinned, not guessed;
//    3. composite knots keep their relative chirality: the granny
//       3_1#3_1 and the square 3_1#m3_1 get different names, and a sum and
//       its global reversal get the same one;
//    4. splits, split unknots, unlinks, knots summed into a link and links
//       summed along components get the atlas syntax, exact only when the
//       name is an identity;
//    5. the Gauss-diagram cut finds exactly the visible summands;
//    6. namers sharing table caches name as namers with their own do.
//

#include <algorithm>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

#include <link/link.h>

#include "linknaming/linknamer.h"
#include "linknaming/tables.h"
#include "linknaming/diagrams/gaussdiagram.h"

using namespace linknaming;

static int passed = 0;
static int failed_count = 0;

#define EXPECT_EQ(actual, expected, desc)                                     \
    do {                                                                      \
        auto _a = (actual);                                                   \
        auto _e = (expected);                                                 \
        if (_a == _e) {                                                       \
            ++passed;                                                         \
        } else {                                                              \
            std::cout << "  FAIL: " << (desc) << "\n        expected " << _e \
                      << ", got " << _a << "\n";                              \
            ++failed_count;                                                   \
        }                                                                     \
    } while (0)

namespace {

const char *KNOTS = "exactnaming_test_knots.csv";
const char *LINKS = "exactnaming_test_links.csv";
const char *SYMMETRY = "exactnaming_test_symmetry.csv";

const std::vector<std::pair<std::string, std::string>> KNOT_ROWS = {
    {"3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]],1"},
    {"4_1", "[[4;2;5;1];[8;6;1;5];[6;3;7;4];[2;7;3;8]],1"},
    {"5_2", "[[1;5;2;4];[3;9;4;8];[5;1;6;10];[7;3;8;2];[9;7;10;6]],1"},
    {"8_17", "[[6;2;7;1];[14;8;15;7];[8;3;9;4];[2;13;3;14];[12;5;13;6];[4;9;5;10];[16;12;1;11];[10;16;11;15]],1"},
    {"10_151", "[[2;15;3;16];[4;19;5;20];[6;3;7;4];[8;14;9;13];[10;8;11;7];[11;18;12;19];[14;10;15;9];[16;1;17;2];[17;12;18;13];[20;5;1;6]],1"},
    {"11a_18", "[[4;2;5;1];[8;4;9;3];[12;6;13;5];[18;16;19;15];[16;10;17;9];[10;18;11;17];[22;19;1;20];[20;13;21;14];[14;21;15;22];[6;12;7;11];[2;8;3;7]],1"},
    {"8_16", "[[2;7;3;8];[4;10;5;9];[6;1;7;2];[8;14;9;13];[10;15;11;16];[12;6;13;5];[14;3;15;4];[16;11;1;12]],1"},
    {"10_156", "[[1;13;2;12];[3;8;4;9];[5;14;6;15];[7;18;8;19];[10;15;11;16];[11;1;12;20];[13;6;14;7];[16;9;17;10];[17;4;18;5];[19;3;20;2]],1"},
};
const std::vector<std::pair<std::string, std::string>> LINK_ROWS = {
    {"L2a1{0}", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]],0"},
    {"L2a1{1}", "PD[X[4; 2; 3; 1]; X[2; 4; 1; 3]],0"},
    {"L4a1{0}", "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]],0"},
    {"L4a1{1}", "PD[X[6; 2; 7; 1]; X[8; 4; 5; 3]; X[2; 8; 3; 7]; X[4; 6; 1; 5]],1"},
    {"L5a1{0}", "PD[X[6; 1; 7; 2]; X[10; 7; 5; 8]; X[4; 5; 1; 6]; X[2; 10; 3; 9]; X[8; 4; 9; 3]],1"},
    {"L5a1{1}", "PD[X[8; 2; 9; 1]; X[10; 7; 5; 8]; X[4; 10; 1; 9]; X[2; 5; 3; 6]; X[6; 3; 7; 4]],1"},
    {"L7n1{0}", "PD[X[6; 1; 7; 2]; X[12; 7; 13; 8]; X[4; 13; 1; 14]; X[5; 10; 6; 11]; X[3; 8; 4; 9]; X[9; 14; 10; 5]; X[11; 2; 12; 3]],2"},
    {"L7n1{1}", "PD[X[12; 2; 13; 1]; X[6; 11; 7; 12]; X[4; 6; 1; 5]; X[13; 8; 14; 9]; X[3; 11; 4; 10]; X[9; 14; 10; 5]; X[7; 3; 8; 2]],0"},
};

void writeTables() {
    std::ofstream k(KNOTS);
    k << "Name,PD Notation,Genus-4D\n";
    for (const auto &[n, rest] : KNOT_ROWS) k << n << "," << rest << "\n";
    std::ofstream l(LINKS);
    l << "Name,PD Notation (KnotTheory),Genus-4D\n";
    for (const auto &[n, rest] : LINK_ROWS) l << n << "," << rest << "\n";
    std::ofstream s(SYMMETRY);
    s << "name,symmetry_type,basis\n3_1,reversible,x\n4_1,fully amphicheiral,x\n"
         "5_2,reversible,x\n8_17,negative amphicheiral,x\n10_151,chiral,x\n11a_18,chiral,x\n"
         "8_16,reversible,x\n10_156,reversible,x\n";
}

regina::Link table(const std::string &name) {
    for (const auto &rows : {KNOT_ROWS, LINK_ROWS})
        for (const auto &[n, rest] : rows)
            if (n == name) return linkFromTablePD(rest.substr(0, rest.rfind(',')));
    throw std::runtime_error("no row " + name);
}

GaussDiagram gauss(const regina::Link &l) {
    std::vector<size_t> origin(l.countComponents());
    for (size_t i = 0; i < origin.size(); ++i) origin[i] = i;
    return GaussDiagram::of(l, origin);
}

// The disjoint union of two diagrams (a split link).
GaussDiagram unite(const GaussDiagram &a, const GaussDiagram &b) {
    GaussDiagram u = a;
    const long off = static_cast<long>(a.crossings());
    u.signs.insert(u.signs.end(), b.signs.begin(), b.signs.end());
    for (const auto &c : b.comps) {
        std::vector<long> &o = u.comps.emplace_back();
        for (long v : c) o.push_back(v > 0 ? v + off : v - off);
    }
    u.origin.clear();
    for (size_t i = 0; i < u.comps.size(); ++i) u.origin.push_back(i);
    return u;
}

// b's component 0 summed into a's component `into`, at both base points.
GaussDiagram sum(const GaussDiagram &a, size_t into, const GaussDiagram &b) {
    GaussDiagram u = unite(a, b);
    const size_t bc = a.comps.size(); // b's component 0, now in u
    u.comps[into].insert(u.comps[into].end(), u.comps[bc].begin(), u.comps[bc].end());
    u.comps.erase(u.comps.begin() + static_cast<long>(bc));
    u.origin.pop_back();
    return u;
}

regina::Link mirrored(regina::Link l) { l.reflect(); return l; }
regina::Link reversedAll(regina::Link l) { l.reverse(); return l; }
regina::Link reversedComponent(regina::Link l, size_t c) { l.reverse(l.component(c)); return l; }

void test_tables(const Tables &t) {
    EXPECT_EQ(t.size(), KNOT_ROWS.size() + LINK_ROWS.size(), "every row loaded");
    EXPECT_EQ(t.canonical("L2a1{1}"), t.canonical("L2a1{0}"),
              "the Hopf link's two orientations are mirror images: one class");
    EXPECT_EQ(t.canonical("L7n1{0}") == t.canonical("L7n1{1}"), false,
              "L7n1's two orientations are different links");
    EXPECT_EQ(t.inconsistentClasses().size(), size_t(0), "no class has two literature values");
}

void test_each_entry_names_itself(const LinkNamer &namer, const Tables &t) {
    for (const auto &rows : {KNOT_ROWS, LINK_ROWS})
        for (const auto &[name, rest] : rows) {
            for (const regina::Link &l : {table(name), mirrored(table(name)), reversedAll(table(name))}) {
                LinkName n = namer.name(l);
                EXPECT_EQ(n.name, t.canonical(name), name + " (or its mirror or reverse) is named as itself");
                EXPECT_EQ(n.isName, true, name + " is an exact name");
            }
        }
}

void test_orientation_variant_is_pinned(const LinkNamer &namer, const Tables &t) {
    EXPECT_EQ(namer.name(reversedComponent(table("L7n1{0}"), 1)).name, t.canonical("L7n1{1}"),
              "L7n1{0} with one component reversed is L7n1{1} (g4 0, not 2)");
    EXPECT_EQ(namer.name(reversedComponent(table("L7n1{1}"), 0)).name, t.canonical("L7n1{0}"),
              "and L7n1{1} with one component reversed is L7n1{0}");
    EXPECT_EQ(namer.name(reversedComponent(table("L4a1{0}"), 1)).name, t.canonical("L4a1{1}"),
              "L4a1{0} with one component reversed is L4a1{1}");
}

void test_composite_knots(const LinkNamer &namer) {
    const GaussDiagram t = gauss(table("3_1")), mt = gauss(mirrored(table("3_1")));
    LinkName granny = namer.name(sum(t, 0, t).link());
    LinkName square = namer.name(sum(t, 0, mt).link());
    EXPECT_EQ(granny.name, std::string("3_1#3_1"), "the granny knot");
    EXPECT_EQ(square.name, std::string("3_1#m3_1"), "the square knot");
    EXPECT_EQ(granny.isName && square.isName, true, "both are exact names");
    EXPECT_EQ(namer.name(sum(mt, 0, mt).link()).name, std::string("3_1#3_1"),
              "the mirror of the granny is named as the granny");
    EXPECT_EQ(namer.name(sum(gauss(table("4_1")), 0, t).link()).name, std::string("3_1#4_1"),
              "an amphicheiral summand carries no mark");
    // 8_17 is negative amphicheiral (r8_17 = m8_17): 3_1#m8_17 is the global
    // reverse of 3_1#8_17 (3_1 is reversible), so the same name.
    const GaussDiagram k = gauss(table("8_17")), mk = gauss(mirrored(table("8_17")));
    EXPECT_EQ(namer.name(sum(t, 0, k).link()).name, namer.name(sum(t, 0, mk).link()).name,
              "3_1#8_17 and its global reverse 3_1#m8_17 get one name");
    LinkName three = namer.name(sum(sum(t, 0, t), 0, gauss(table("5_2"))).link());
    EXPECT_EQ(three.name, std::string("3_1#3_1#5_2"), "three summands");
}

void test_splits_and_sums(const LinkNamer &namer, const Tables &t) {
    const GaussDiagram t31 = gauss(table("3_1")), f41 = gauss(table("4_1"));
    GaussDiagram unknot;
    unknot.comps = {{}};
    unknot.origin = {0};
    LinkName split = namer.name(unite(t31, f41).link());
    EXPECT_EQ(split.name, std::string("3_1 u 4_1"), "a split of two knots");
    EXPECT_EQ(split.isName, false, "which is a description, not an identity");
    LinkName withUnknot = namer.name(unite(t31, unknot).link());
    EXPECT_EQ(withUnknot.name, std::string("3_1 u Unknot"), "a knot and a split unknot");
    EXPECT_EQ(withUnknot.isName, true, "which is an identity");
    EXPECT_EQ(namer.name(unite(unknot, unknot).link()).name, std::string("2-component unlink"),
              "two split unknots");
    const GaussDiagram hopf = gauss(table("L2a1{0}"));
    LinkName knotIntoLink = namer.name(sum(hopf, 0, t31).link());
    EXPECT_EQ(knotIntoLink.name, "3_1 #_? " + t.canonical("L2a1{0}"),
              "a trefoil summed into a component of the Hopf link");
    EXPECT_EQ(knotIntoLink.isName, false, "which is a description");
    // Two Hopf links summed along one component each: the 3-chain.
    GaussDiagram chain = unite(hopf, hopf);
    chain.comps[1].insert(chain.comps[1].end(), chain.comps[2].begin(), chain.comps[2].end());
    chain.comps.erase(chain.comps.begin() + 2);
    chain.origin = {0, 1, 2};
    LinkName three = namer.name(chain.link());
    const std::string h = t.canonical("L2a1{0}");
    EXPECT_EQ(three.name, "#{" + h + "[?] # " + h + "[?]}", "two Hopf links summed along a component");
    EXPECT_EQ(three.pinned, true, "each piece pinned");
}

// The search path, forced: no simplification, so only a Reidemeister search
// can reach the table diagram -- which proves the link only up to mirror and
// orientations -- and the orientation variant must then come from the
// invariants of the piece as drawn.
void test_search_then_invariants(const Tables &t) {
    NamerLimits limits;
    limits.isometry = false; // it would answer first (test_isometry)
    limits.simplifyTries = 0;
    limits.exhaustiveHeight = 0;
    limits.searchHeight = 1;
    LinkNamer searchOnly(t, limits);
    for (const std::string &name : {std::string("L7n1{1}"), std::string("L7n1{0}"), std::string("5_2")}) {
        regina::Link l = table(name);
        l.r1(l.component(0), 0, 1); // a kink: no longer the table diagram
        PieceName p = searchOnly.namePiece(gauss(l));
        EXPECT_EQ(p.by == PieceName::By::searchAndInvariants, true,
                  name + " with a kink is found by search, not by its diagram");
        EXPECT_EQ(p.display(), t.canonical(name),
                  name + " with a kink: the invariants pin its variant");
    }
    // Reversing one component after the search: the other variant.
    regina::Link l = reversedComponent(table("L7n1{0}"), 1);
    l.r1(l.component(0), 0, 1);
    EXPECT_EQ(searchOnly.namePiece(gauss(l)).display(), t.canonical("L7n1{1}"),
              "L7n1{0}, one component reversed, with a kink: L7n1{1}");
}

// Two outgoing links from the master (2026-09-27) that the forward search left
// untabulated and the SnapPy pipeline named: each is another minimal diagram
// of its table knot. The table-side step names them, and without it (or with
// its forward-search-only limits) they stay untabulated.
void test_table_side(const Tables &t) {
    // 10_151 (non-alternating): reached by rewrite() outward from the table
    // diagram at height 2. 11a_18 (alternating): in the table diagram's
    // flype orbit.
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"k-LbSTpoqLnsCyc", "10_151"}, {"l3-aqR6Fa4qPEHyf", "11a_18"}};
    NamerLimits on;
    on.isometry = false;  // it would answer first (test_isometry)
    on.searchHeight = -1; // the forward search cannot reach these; skip it
    on.tableSideHeight = 2;
    NamerLimits off = on;
    off.tableSideHeight = -1;
    LinkNamer withTableSide(t, on), without(t, off);
    for (const auto &[sig, knot] : cases) {
        const regina::Link l = regina::Link::fromSig(sig);
        const PieceName p = withTableSide.namePiece(gauss(l));
        EXPECT_EQ(p.display(), knot, "a gap diagram of " + knot + " is named from the table side");
        EXPECT_EQ(p.by == PieceName::By::searchAndInvariants, true,
                  knot + ": proved by search, the variant by invariants");
        EXPECT_EQ(without.namePiece(gauss(l)).by == PieceName::By::untabulated, true,
                  knot + ": untabulated without the table-side step");
    }
}

// The same outgoing links, and a 10-crossing drawing of 8_16 (two crossings above
// minimal; a table-side rewrite needs height 4 and ~30 s), named by an
// isometry of complements carrying meridians to meridians, with every
// Reidemeister search off. 8_16 shares its HOMFLY polynomial with 10_156,
// which is in the test tables: the isometry must pick 8_16.
void test_isometry(const Tables &t) {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"k-LbSTpoqLnsCyc", "10_151"}, {"l3-aqR6Fa4qPEHyf", "11a_18"},
        {"k-ygSLpCidvGDyc", "8_16"}};
    NamerLimits on;
    on.searchHeight = -1;
    on.tableSideHeight = -1;
    NamerLimits off = on;
    off.isometry = false;
    LinkNamer withIsometry(t, on), without(t, off);
    for (const auto &[sig, knot] : cases) {
        const regina::Link l = regina::Link::fromSig(sig);
        const PieceName p = withIsometry.namePiece(gauss(l));
        EXPECT_EQ(p.display(), knot, "a gap diagram of " + knot + " is named by isometry");
        EXPECT_EQ(p.by == PieceName::By::isometry, true,
                  knot + ": proved by isometry, the variant by invariants");
        EXPECT_EQ(without.namePiece(gauss(l)).by == PieceName::By::untabulated, true,
                  knot + ": untabulated with every step off");
    }
    // A torus knot's complement is not hyperbolic: its kinked diagram is left
    // to the searches, and with them off it stays untabulated.
    regina::Link kinked = table("3_1");
    kinked.r1(kinked.component(0), 0, 1);
    NamerLimits noSimplify = on;
    noSimplify.simplifyTries = 0;
    noSimplify.exhaustiveHeight = 0;
    EXPECT_EQ(LinkNamer(t, noSimplify).namePiece(gauss(kinked)).by == PieceName::By::untabulated,
              true, "a kinked 3_1 is not named by isometry");
}

// Namers over the same tables can share what they learn about them (a
// goal run's link namer and every search's outgoing namer do): a second namer
// names through the first's HOMFLY index and kernel complements exactly as a
// namer with its own would. Caches built for other tables -- even loaded from
// the same files -- are refused, since they are keyed by those tables'
// entries.
void test_shared_caches(const Tables &t) {
    const std::vector<std::pair<std::string, std::string>> cases = {
        {"k-LbSTpoqLnsCyc", "10_151"}, {"l3-aqR6Fa4qPEHyf", "11a_18"},
        {"k-ygSLpCidvGDyc", "8_16"}};
    NamerLimits on;
    on.searchHeight = -1;
    on.tableSideHeight = -1;
    LinkNamer first(t, on), alone(t, on);
    LinkNamer second(t, on, first.caches());
    EXPECT_EQ(second.caches() == first.caches(), true, "the second namer holds the first's caches");
    EXPECT_EQ(alone.caches() == first.caches(), false, "a namer given none makes its own");
    for (const auto &[sig, knot] : cases) {
        const regina::Link l = regina::Link::fromSig(sig);
        const PieceName a = first.namePiece(gauss(l));
        const PieceName b = second.namePiece(gauss(l));
        const PieceName c = alone.namePiece(gauss(l));
        EXPECT_EQ(a.display(), knot, knot + ": named by the first namer");
        EXPECT_EQ(b.display(), a.display(), knot + ": the same name through shared caches");
        EXPECT_EQ(b.by == PieceName::By::isometry, true, knot + ": by isometry through shared caches");
        EXPECT_EQ(c.display(), a.display(), knot + ": the same name through its own caches");
    }
    const Tables other = Tables::load(KNOTS, LINKS, SYMMETRY);
    bool refused = false;
    try {
        LinkNamer wrong(other, on, first.caches());
    } catch (const std::invalid_argument &) {
        refused = true;
    }
    EXPECT_EQ(refused, true, "caches built for other tables are refused");
}

void test_visible_sum(void) {
    GaussDiagram t = gauss(table("3_1")), f = gauss(table("4_1"));
    GaussDiagram s = sum(t, 0, f);
    auto cut = visibleSum(s);
    EXPECT_EQ(cut.has_value(), true, "3_1 # 4_1 shows a visible sum");
    if (cut) {
        std::vector<size_t> sizes{cut->first.crossings(), cut->second.crossings()};
        std::sort(sizes.begin(), sizes.end());
        EXPECT_EQ(sizes[0] == 3 && sizes[1] == 4, true, "cut into 3 + 4 crossings");
        EXPECT_EQ(cut->first.link().isClassical() && cut->second.link().isClassical(), true,
                  "both summands planar");
    }
    EXPECT_EQ(visibleSum(t).has_value(), false, "the trefoil shows none");
}

} // namespace

int main() {
    writeTables();
    Tables tables = Tables::load(KNOTS, LINKS, SYMMETRY);
    LinkNamer namer(tables);
    test_tables(tables);
    test_each_entry_names_itself(namer, tables);
    test_orientation_variant_is_pinned(namer, tables);
    test_composite_knots(namer);
    test_splits_and_sums(namer, tables);
    test_search_then_invariants(tables);
    test_table_side(tables);
    test_isometry(tables);
    test_shared_caches(tables);
    test_visible_sum();
    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}
