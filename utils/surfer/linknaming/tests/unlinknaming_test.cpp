// unlinknaming_test.cpp
//
// Tests for identify::certifiesUnlink() and identify::capInCone() (see
// ../complement/unlinknaming.h): the census-free recognition of unknots and
// unlinks that KnottedSurface's local flatness and resolution checks rest
// on.
//
// Split from identifycomplement_test.cpp (refactor phase 2).

#include <iostream>
#include <string>
#include <unistd.h>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/example3.h>

#include "linknaming/complement/unlinknaming.h"
#include "diagramtriangulation/fromdiagram.h"
#include "linknaming/complement/linkcomplement.h"

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

// ─────────────────────────────────────────────────────────────────────────
// identify::certifiesUnlink(): the certificate behind --resolve-unlinked.
// It is one-sided, so the tests that matter most are the negatives: a
// false "unlink" would admit a genuinely self-intersecting surface.
// ─────────────────────────────────────────────────────────────────────────

// PD codes from Regina (Link::pdData()), each checked there, independently of
// anything in utils/surfer, via Link::complement().group():
//  - kUnlink3PD: Link(3) with two classical R2 moves (component 0 over 1, then
//    1 over 2), so a CONNECTED 4-crossing diagram -- knotbuilder needs one --
//    of the 3-component unlink. Its group simplifies to Free(3).
//  - Whitehead, Borromean: ExampleLink::whitehead()/borromean(). Groups with
//    2 generators/1 relation and 3 generators/2 relations. Both have
//    pairwise linking number 0, which is exactly why they are the controls:
//    P_transverse cannot see them, so only the certificate stands in the way.
//
// Passed through parsePDCode(), never as PDCode literals: buildLink() wants
// 0-indexed labels, and parsePDCode() is what normalizes Regina's 1-indexed
// ones.
const knotbuilder::PDCode kUnlink2PD = {{0, 3, 1, 2}, {1, 3, 0, 2}};
const knotbuilder::PDCode kUnlink3PD =
    knotbuilder::parsePDCode("3 1 4 2 6 1 3 2 7 5 8 6 8 5 7 4");
const knotbuilder::PDCode kWhiteheadPD =
    knotbuilder::parsePDCode("1 7 2 8 3 10 4 9 5 2 6 3 8 5 9 4 10 6 7 1");
const knotbuilder::PDCode kBorromeanPD = knotbuilder::parsePDCode(
    "1 10 2 11 3 9 4 12 5 2 6 3 7 1 8 4 9 6 10 7 11 5 12 8");

void test_certifies_unlink_positive() {
    {
        auto [tri, edges, reversed] = knotbuilder::buildLink(kUnlink2PD);
        EXPECT_EQ(identify::certifiesUnlink(tri, edges, 2), true,
                  "the 2-component unlink is certified");
        EXPECT_EQ(identify::certifiesUnlink(tri, edges, 1), false,
                  "...but not as a 1-component unlink (wrong m)");
        EXPECT_EQ(identify::certifiesUnlink(tri, edges, 3), false,
                  "...nor as a 3-component one");

        Link link(tri, edges);
        std::vector<const regina::Edge<3> *> one = link.comps_[0].edges();
        EXPECT_EQ(identify::certifiesUnlink(tri, one, 1), true,
                  "one component of it alone is certified the unknot");
    }
    {
        auto [tri, edges, reversed] = knotbuilder::buildLink(kUnlink3PD);
        EXPECT_EQ(Link(tri, edges).countComponents(), 3,
                  "fixture sanity: the R2-tangled unlink has 3 components");
        EXPECT_EQ(identify::certifiesUnlink(tri, edges, 3), true,
                  "the 3-component unlink is certified");
    }
}

void test_certifies_unlink_negative() {
    auto check = [](const char *label, const knotbuilder::PDCode &pd,
                    size_t m) {
        auto [tri, edges, reversed] = knotbuilder::buildLink(pd);
        EXPECT_EQ(Link(tri, edges).countComponents(), m,
                  std::string("fixture sanity: ") + label + " has " +
                      std::to_string(m) + " components");
        EXPECT_EQ(identify::certifiesUnlink(tri, edges, m), false,
                  std::string(label) + " is not certified an unlink");
        // ...and for the right reason: the complement is built, and its
        // group keeps a relation, rather than certifiesUnlink() having bailed
        // out on an exception or on malformed input.
        auto complement = EdgeComplement(tri, edges).buildComplement();
        EXPECT_EQ(complement.group().countRelations() >= 1, true,
                  std::string(label) +
                      ": its complement's group keeps at least one relation");
    };
    check("Hopf link", knotbuilder::parsePDCode("1 4 2 3 3 2 4 1"), 2);
    check("Whitehead link (lk = 0)", kWhiteheadPD, 2);
    check("Borromean rings (pairwise lk = 0)", kBorromeanPD, 3);
    check("trefoil", knotbuilder::parsePDCode("1 4 2 5 3 6 4 1 5 2 6 3"), 1);
}

void test_certifies_unlink_malformed() {
    auto [tri, edges, reversed] = knotbuilder::buildLink(kUnlink2PD);
    std::vector<const regina::Edge<3> *> dropped(edges.begin() + 1,
                                                 edges.end());
    EXPECT_EQ(identify::certifiesUnlink(tri, dropped, 2), false,
              "an edge set with an arc in it is refused, whatever its "
              "complement");
    std::vector<const regina::Edge<3> *> doubled = edges;
    doubled.push_back(edges.front());
    EXPECT_EQ(identify::certifiesUnlink(tri, doubled, 2), false,
              "a repeated edge is refused");
    EXPECT_EQ(identify::certifiesUnlink(tri, {}, 0), false,
              "the empty link is refused");
}

// A 3-ball containing a proper arc: `sphere` (from knotbuilder) with one
// tetrahedron containing a knot edge removed. That edge now lies in the
// boundary sphere, so the rest of the knot is an arc whose capping -- by
// that very edge -- recovers the original knot. So the arc is knotted
// exactly when the knot is, which is what capInCone() + isUnknot() must see.
struct BallWithArc {
    regina::Triangulation<3> ball;
    std::vector<const regina::Edge<3> *> arc; // points into `ball`
};

bool buildBallWithArc(const regina::Triangulation<3> &sphere,
                      const std::vector<const regina::Edge<3> *> &knot,
                      BallWithArc &out) {
    for (const auto *cut : knot) {
        for (const auto &emb : *cut) {
            const auto *tet = emb.tetrahedron();
            bool distinct = true;
            for (int i = 0; i < 4 && distinct; ++i)
                for (int j = i + 1; j < 4 && distinct; ++j)
                    distinct = tet->vertex(i) != tet->vertex(j);
            if (!distinct)
                continue;

            out.ball = sphere;
            // Simplex pointers survive removeSimplex() of a different
            // simplex; face pointers do not, so re-resolve afterwards.
            std::vector<std::pair<regina::Tetrahedron<3> *, int>> desc;
            bool ok = true;
            for (const auto *k : knot) {
                if (k == cut)
                    continue;
                bool found = false;
                for (const auto &e2 : *k)
                    if (e2.tetrahedron()->index() != tet->index()) {
                        desc.emplace_back(
                            out.ball.tetrahedron(e2.tetrahedron()->index()),
                            e2.edge());
                        found = true;
                        break;
                    }
                if (!found) {
                    ok = false;
                    break;
                }
            }
            if (!ok)
                continue;
            out.ball.removeSimplex(out.ball.tetrahedron(tet->index()));
            if (!out.ball.isBall())
                continue;
            out.arc.clear();
            for (auto [t, i] : desc)
                out.arc.push_back(t->edge(i));
            return true;
        }
    }
    return false;
}

void test_cap_in_cone_single_tetrahedron() {
    // One tetrahedron: a 3-ball whose every edge is on the boundary. An
    // edge is an arc with both ends on the boundary; capped through the
    // apex it is a triangle -- trivially an unknot in the coned-off S^3.
    regina::Triangulation<3> ball;
    ball.newTetrahedron();
    identify::CappedCurves capped;
    EXPECT_EQ(identify::capInCone(ball, {ball.edge(0)}, capped), true,
              "a boundary edge of a ball is capped");
    EXPECT_EQ(capped.edges.size(), size_t{3},
              "...into a 3-edge loop through the apex");
    EXPECT_EQ(capped.components, size_t{1}, "...forming one closed curve");
    EXPECT_EQ(capped.tri.isSphere(), true,
              "the coned-off ball is a 3-sphere");
    EXPECT_EQ(identify::certifiesUnlink(capped.tri, capped.edges, 1), true,
              "the capped edge is certified unknotted");

    identify::CappedCurves none;
    regina::Triangulation<3> closed = regina::Example<3>::threeSphere();
    EXPECT_EQ(identify::capInCone(closed, {closed.edge(0)}, none), false,
              "a closed triangulation has no boundary to cone off");
}

void test_cap_in_cone_detects_knotted_arc() {
    struct Case {
        const char *label;
        const char *pd;
        bool unknotted;
    };
    for (const Case &c :
         {Case{"R1-kinked unknot", "1 2 2 1", true},
          Case{"trefoil", "1 4 2 5 3 6 4 1 5 2 6 3", false}}) {
        auto result = knotbuilder::buildLink(knotbuilder::parsePDCode(c.pd));
        BallWithArc fixture;
        bool built = buildBallWithArc(result.tri, result.edges, fixture);
        EXPECT_EQ(built, true,
                  std::string(c.label) +
                      ": some tetrahedron's removal leaves a ball");
        if (!built)
            continue;
        identify::CappedCurves capped;
        EXPECT_EQ(identify::capInCone(fixture.ball, fixture.arc, capped), true,
                  std::string(c.label) + ": the arc is capped");
        EXPECT_EQ(capped.components, size_t{1},
                  std::string(c.label) + ": into one closed curve");
        EXPECT_EQ(identify::isUnknot(Knot(capped.tri, capped.edges)),
                  c.unknotted,
                  std::string(c.label) + ": capped arc is " +
                      (c.unknotted ? "unknotted" : "knotted"));
    }
}

} // namespace

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("certifies_unlink_positive", test_certifies_unlink_positive);
    run("certifies_unlink_negative", test_certifies_unlink_negative);
    run("certifies_unlink_malformed", test_certifies_unlink_malformed);
    run("cap_in_cone_single_tetrahedron", test_cap_in_cone_single_tetrahedron);
    run("cap_in_cone_detects_knotted_arc",
        test_cap_in_cone_detects_knotted_arc);


    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
