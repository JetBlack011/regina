// cobordismbuilder_test.cpp
// Tests for CobordismBuilder::thicken() — verifies SimplicialPrism
// construction and gluing produce a valid product cobordism. See
// knotbuilder_test.cpp for integration with knotbuilder's PD-code ->
// triangulated-S³ pipeline.

#include <iostream>
#include <string>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>
#include <unistd.h>

#include "diagramtriangulation/thickening/thickening.h"

static int passed = 0, failed_count = 0;

// ANSI color manipulators for PASS/FAIL output. Auto-disabled when stdout
// isn't a terminal (e.g. piped to a file), so logs don't fill up with
// escape codes.
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

// ─────────────────────────────────────────────────────────────────────────────
// Single tetrahedron (B³), one layer.
//
// No glue() calls (all four faces are boundary).  Tests the prism constructor
// alone.  B³ × I ≅ B⁴:
//   - 1 tet × 4 prism pentachora = 4 pentachora
//   - ∂(B³ × I) ≅ S³ → 1 boundary component
// ─────────────────────────────────────────────────────────────────────────────
void test_ball_thicken_one_layer() {
    std::cout << "\n--- B³ (single tet), thicken(1): no glue() calls ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();

    CobordismBuilder<3> cob(ball);
    auto &result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 4,
              "1 tet × 4 prism simplices = 4 pentachora");
    EXPECT_EQ(result.isValid(), true, "B⁴ triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "B⁴ triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              "∂(B³×I) ≅ S³: 1 boundary component");
}

// ─────────────────────────────────────────────────────────────────────────────
// regina::Example<3>::threeSphere() is *not* two tetrahedra glued to each
// other: it's a single tetrahedron with two pairs of its own faces
// self-identified (verified directly: tri.size() == 1). It is also not
// ordered as constructed, so this also exercises CobordismBuilder's
// constructor-time call to Triangulation<3>::order(). This is a good
// self-gluing stress case for glue(), since the two self-identifications
// use distinct, non-matching facet/otherFacet pairs once ordered.
// ─────────────────────────────────────────────────────────────────────────────
void test_s3_selfglued_thicken_one_layer() {
    std::cout << "\n--- S³ as a single self-glued tet, thicken(1) ---\n";

    auto s3 = regina::Example<3>::threeSphere();
    EXPECT_EQ((int)s3.size(), 1, "threeSphere() is a single self-glued tet");

    CobordismBuilder<3> cob(s3);
    auto &result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 4,
              "1 tet × 4 prism simplices = 4 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×I triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "S³×I triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(0) ≅ S³");
    EXPECT_EQ(
        result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(1) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// Same self-glued S³, two layers: exercises stitchTop()'s inter-layer seam
// on top of glue()'s self-gluing case.
// ─────────────────────────────────────────────────────────────────────────────
void test_s3_selfglued_thicken_two_layers() {
    std::cout << "\n--- S³ as a single self-glued tet, thicken(2) ---\n";

    auto s3 = regina::Example<3>::threeSphere();

    CobordismBuilder<3> cob(s3);
    auto &result = cob.thicken(2);

    EXPECT_EQ((int)result.size(), 8, "2 layers × 4 pentachora = 8 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×[0,2] triangulation is valid");
    EXPECT_EQ(result.isConnected(), true,
              "S³×[0,2] triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(0) ≅ S³");
    EXPECT_EQ(
        result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(1) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// A genuine two-*distinct*-tetrahedra triangulation of S³: two tetrahedra
// doubled along their entire boundary via the identity map on each face
// (i.e. attaching two 3-balls along the identity of their boundary S²).
// No self-gluing anywhere, and glue() is called once per facet (4 times).
// ─────────────────────────────────────────────────────────────────────────────
regina::Triangulation<3> doubledTetrahedra() {
    regina::Triangulation<3> tri;
    auto a = tri.newTetrahedron();
    auto b = tri.newTetrahedron();
    for (int f = 0; f < 4; ++f)
        a->join(f, b, regina::Perm<4>());
    return tri;
}

void test_doubled_tetrahedra_thicken_one_layer() {
    std::cout << "\n--- Doubled tetrahedra (2 distinct tets), thicken(1): "
                 "4 glue() calls ---\n";

    auto s3 = doubledTetrahedra();
    EXPECT_EQ((int)s3.size(), 2, "doubled tetrahedra has 2 distinct tets");
    EXPECT_EQ(s3.isSphere(), true, "doubled tetrahedra triangulates S³");

    CobordismBuilder<3> cob(s3);
    auto &result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 8,
              "2 tets × 4 prism simplices = 8 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×I triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "S³×I triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(0) ≅ S³");
    EXPECT_EQ(
        result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(1) ≅ S³");
}

void test_doubled_tetrahedra_thicken_two_layers() {
    std::cout << "\n--- Doubled tetrahedra, thicken(2): inter-layer gluing "
                 "with 2 distinct tets ---\n";

    auto s3 = doubledTetrahedra();

    CobordismBuilder<3> cob(s3);
    auto &result = cob.thicken(2);

    EXPECT_EQ((int)result.size(), 16,
              "2 layers × 8 pentachora = 16 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×[0,2] triangulation is valid");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(0) ≅ S³");
    EXPECT_EQ(
        result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(1) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// ∂Δ⁴ as S³ (5 tetrahedra), one layer.
//
// In this triangulation every tet is adjacent to every other, so glue() is
// called 10 times, with a mix of matching and non-matching facet/otherFacet
// pairs. This stresses the glue() permutation logic more than the 2-tet
// case. The extracted boundary triangulation is not ordered as built, so
// this also exercises the CobordismBuilder constructor's ordering step.
// S³ × I:
//   - 5 tets × 4 pentachora = 20 pentachora
//   - 2 boundary components, each ≅ ∂Δ⁴ S³
// ─────────────────────────────────────────────────────────────────────────────
void test_bdy_delta4_thicken_one_layer() {
    std::cout << "\n--- ∂Δ⁴ as S³ (5 tets), thicken(1): 10 glue() calls ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    auto s3 = fourBall.boundaryComponent(0)->build();

    std::cout << "    " << s3.size() << " tetrahedra\n";

    CobordismBuilder<3> cob(s3);
    auto &result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 20,
              "5 tets × 4 prism simplices = 20 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×I triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "S³×I triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(0) ≅ S³");
    EXPECT_EQ(
        result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(1) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// Exhaustive regression test for SimplicialPrism::glue()'s vertex
// permutation. For 4-vertex tetrahedra there are exactly 16 (facet,
// otherFacet) pairs; for each one we glue two tetrahedra along exactly that
// one facet pair (via the unique order-preserving gluing permutation), so
// the base triangulation is always a ball (two tetrahedra glued along a
// single common face). Thickening this ball must always give a valid
// pentachoron triangulation with exactly one boundary component, and that
// boundary component must be S³ (the double of a ball). This is exactly
// the class of case that exposed the original bug, where glue()'s inner
// vertex permutation was computed as "order-preserving on raw local
// indices" instead of decoding through the actual (base vertex, top/bottom)
// meaning of each index — that bug only manifested when facet != otherFacet,
// so every combination needs to be checked, not just a couple of examples.
// ─────────────────────────────────────────────────────────────────────────────
void test_all_facet_pairs_two_tets() {
    std::cout << "\n--- All 16 (facet, otherFacet) pairs, 2 tets glued at "
                 "one face ---\n";

    int numOk = 0;
    for (int facet = 0; facet < 4; ++facet) {
        for (int otherFacet = 0; otherFacet < 4; ++otherFacet) {
            regina::Triangulation<3> tri;
            auto a = tri.newTetrahedron();
            auto b = tri.newTetrahedron();

            // The unique order-preserving bijection {0,..,3}\{facet} ->
            // {0,..,3}\{otherFacet}, with facet -> otherFacet.
            std::array<int, 4> image;
            for (int i = 0; i < 4; ++i) {
                if (i == facet) {
                    image[i] = otherFacet;
                    continue;
                }
                int rank = (i < facet) ? i : i - 1;
                image[i] = (rank < otherFacet) ? rank : rank + 1;
            }
            a->join(facet, b, regina::Perm<4>(image));

            CobordismBuilder<3> cob(tri);
            auto &result = cob.thicken(1);

            bool ok = result.isValid() && result.size() == 8 &&
                      result.countBoundaryComponents() == 1 &&
                      result.boundaryComponent(0)->build().isSphere();
            if (ok) {
                ++numOk;
            } else {
                std::cout << "  FAIL at facet=" << facet
                          << " otherFacet=" << otherFacet
                          << ": valid=" << result.isValid()
                          << " size=" << result.size()
                          << " bcs=" << result.countBoundaryComponents()
                          << "\n";
            }
        }
    }

    EXPECT_EQ(
        numOk, 16,
        "all 16 (facet, otherFacet) pairs give a valid B⁴ with S³ boundary");
}

// ─────────────────────────────────────────────────────────────────────────────
// Regression test: CobordismBuilder::isOrdered() previously looped `for (int
// f = 0; f < d; ++f)`, which only checks facets 0..d-1 of a d-simplex —
// but a d-simplex has d+1 facets (0..d), so the last facet's gluing was
// never inspected. This let non-ordered triangulations whose only bad
// gluing was on a last facet slip through as "ordered", which in turn
// meant CobordismBuilder's constructor would skip calling order()
// (dim == 3) or skip throwing (dim != 3) when it should not have.
// ─────────────────────────────────────────────────────────────────────────────
void test_isOrdered_checks_every_facet() {
    std::cout << "\n--- isOrdered() must check every facet, including the "
                 "last one ---\n";

    // dim = 3: two tetrahedra glued only at facet 3 (the last facet), with
    // a non-order-preserving permutation (swap 0,1) on the remaining
    // vertices.
    {
        regina::Triangulation<3> tri;
        auto a = tri.newTetrahedron();
        auto b = tri.newTetrahedron();
        a->join(3, b, regina::Perm<4>(std::array<int, 4>{1, 0, 2, 3}));

        EXPECT_EQ(CobordismBuilder<3>::isOrdered(tri), tri.isOrdered(),
                  "dim=3: our isOrdered() agrees with Regina's native check "
                  "for a bad gluing on the last facet");
        EXPECT_EQ(CobordismBuilder<3>::isOrdered(tri), false,
                  "dim=3: gluing swapping vertices 0,1 on facet 3 is "
                  "correctly detected as not ordered");
    }

    // dim = 2: same idea, but dim=2 has no order() fallback, so an
    // undetected bad gluing means the constructor wrongly proceeds instead
    // of throwing.
    {
        regina::Triangulation<2> tri;
        auto a = tri.newTriangle();
        auto b = tri.newTriangle();
        a->join(2, b, regina::Perm<3>(1, 0));

        EXPECT_EQ(CobordismBuilder<2>::isOrdered(tri), false,
                  "dim=2: gluing swapping vertices 0,1 on facet 2 (the last "
                  "facet of a triangle) is correctly detected as not "
                  "ordered");

        bool threw = false;
        try {
            CobordismBuilder<2> cob(tri);
        } catch (const std::exception &) {
            threw = true;
        }
        EXPECT_EQ(threw, true,
                  "dim=2: constructor throws on this genuinely unordered "
                  "triangulation instead of silently proceeding");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Doubled triangle (S²) via CobordismBuilder<2>, exercising SimplicialPrism<3>
// instead of SimplicialPrism<4>. This checks that glue()'s formulas
// (derived generically in terms of the template parameter `dim`, not
// hardcoded to any particular dimension) actually hold across dimensions,
// not just for the dim=4 prisms used elsewhere in this file.
// ─────────────────────────────────────────────────────────────────────────────
void test_dim2_doubled_triangle() {
    std::cout << "\n--- Doubled triangle (S²), thicken(2): "
                 "CobordismBuilder<2> / SimplicialPrism<3> ---\n";

    regina::Triangulation<2> tri;
    auto a = tri.newTriangle();
    auto b = tri.newTriangle();
    for (int f = 0; f < 3; ++f)
        a->join(f, b, regina::Perm<3>());

    EXPECT_EQ(CobordismBuilder<2>::isOrdered(tri), true,
              "doubled triangle is already ordered");
    EXPECT_EQ(tri.isSphere(), true, "doubled triangle triangulates S²");

    CobordismBuilder<2> cob(tri);
    auto &result = cob.thicken(2);

    EXPECT_EQ((int)result.size(), 12,
              "2 layers × (2 triangles × 3 prism tets) = 12 tetrahedra");
    EXPECT_EQ(result.isValid(), true, "S²×[0,2] triangulation is valid");
    EXPECT_EQ(result.isConnected(), true,
              "S²×[0,2] triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S²)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(tri).has_value(),
        true, "bc(0) ≅ S²");
    EXPECT_EQ(
        result.boundaryComponent(1)->build().isIsomorphicTo(tri).has_value(),
        true, "bc(1) ≅ S²");
}

// ─────────────────────────────────────────────────────────────────────────────
// baseBoundaryComponent(): identifies the bottom (untouched) boundary
// component, not the top -- the hard case, since thicken() literally builds
// (base) x [0,1], so bottom and top are always combinatorially *identical*
// (isomorphic, in fact often isoSig-identical), meaning no external,
// isomorphism-based test can ever tell "correctly identified the bottom"
// apart from "happened to return a component isomorphic to the base" (both
// would pass). This instead cross-checks against an ambient VERTEX
// independently known to lie on the bottom: SimplicialPrism's own encoding
// (see decode_()/encode_() in simplicialprism.cpp) places the bottom copy
// of any base vertex v at local index v within simplex(v) of a layer's own
// prism, and thicken_()'s wall-gluing never touches that specific local
// vertex (only facets) -- so the *first* layer's own simplex(0)'s local
// vertex 0, captured via currentTopSimplex() immediately after the first
// thicken() call (before any later call overwrites it), is the base
// tetrahedron's own unmoved vertex 0, on the permanently untouched bottom.
// This uses a different access path (.vertex(0) on a directly-captured
// Simplex<4>*) than baseBoundaryComponent()'s own implementation (.face<
// dim>(dim+1) on a facet), so an indexing bug in either would show up as a
// mismatch here.
// ─────────────────────────────────────────────────────────────────────────────
void test_base_boundary_component_is_the_bottom_not_top() {
    std::cout << "\n--- baseBoundaryComponent(): identifies bottom, not top, "
                 "across multiple thickenings ---\n";

    auto s3 = doubledTetrahedra();
    CobordismBuilder<3> cob(s3);

    bool threwBeforeAnyThicken = false;
    try {
        cob.baseBoundaryComponent();
    } catch (const std::exception &) {
        threwBeforeAnyThicken = true;
    }
    EXPECT_EQ(threwBeforeAnyThicken, true,
              "throws if called before any thicken() call");

    const regina::Simplex<3> *tet0 = cob.baseTriangulation().simplex(0);
    cob.thicken(); // layer 1
    regina::Simplex<4> *firstLayerSimplex0 = cob.currentTopSimplex(tet0, 0);

    cob.thicken(2); // 2 more layers -- topPrisms_ now reflects layer 3, not layer 1
    auto &result = cob.getCobordism();
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "3 layers, no cone(): 2 distinct (but isomorphic) boundary "
              "components");

    regina::BoundaryComponent<4> *base = cob.baseBoundaryComponent();
    regina::Vertex<4> *bottomVertex = firstLayerSimplex0->vertex(0);

    EXPECT_EQ(bottomVertex->boundaryComponent() == base, true,
              "an ambient vertex independently known to lie on the "
              "permanently untouched bottom is reported as being on "
              "baseBoundaryComponent()");

    // Sanity check in the opposite direction: a vertex from the *most
    // recent* layer's own top-facing side should NOT be on the base
    // boundary component (it's on the *other* one instead), confirming
    // this isn't trivially true for every vertex.
    const regina::Simplex<3> *tet1 = cob.baseTriangulation().simplex(1);
    regina::Simplex<4> *lastLayerSimplex0 = cob.currentTopSimplex(tet1, 0);
    // local vertex `dim+1` (= 4) is the TOP copy of base vertex 0 -- see
    // decode_()'s m==dim case.
    regina::Vertex<4> *topVertex = lastLayerSimplex0->vertex(4);
    EXPECT_EQ(topVertex->boundaryComponent() == base, false,
              "a vertex independently known to lie on the most recent "
              "layer's top is NOT reported as being on "
              "baseBoundaryComponent()");
}

// ─────────────────────────────────────────────────────────────────────────────
// buildAmbient()'s collar seed comes in index order (phase 3.0), never in the
// order of CollarBuilder::resolve()'s pointer set: a search numbers its
// surfaces' triangles from the seed's, so each component's orientation, and
// so where each outgoing curve starts, followed the heap layout.
// ─────────────────────────────────────────────────────────────────────────────
void test_build_ambient_seed_in_index_order() {
    std::cout << "\n--- buildAmbient(): the collar seed in index order ---\n";
    const char *pds[] = {"[[1;5;2;4];[3;1;4;6];[5;3;6;2]]",
                         "PD[X[6; 1; 7; 2]; X[8; 3; 5; 4]; X[2; 5; 3; 6]; X[4; 7; 1; 8]]"};
    for (const char *pd : pds) {
        ThickenedLink thickened;
        buildAmbient(pd, 2, 2, thickened);
        bool increasing = !thickened.seedFaces.empty();
        for (size_t i = 1; i < thickened.seedFaces.size(); ++i)
            if (thickened.seedFaces[i - 1] >= thickened.seedFaces[i]) increasing = false;
        EXPECT_EQ(increasing, true,
                  std::string(pd) + ": seed faces non-empty and strictly increasing");
    }
}

template <typename F>
void run(const char *name, F fn) {
    std::cout << "\nRunning " << name << "...\n";
    try {
        fn();
    } catch (const std::exception &e) {
        std::cout << red << "  EXCEPTION: " << e.what() << resetColor << "\n";
        ++failed_count;
    }
}

int main() {
    run("test_ball_thicken_one_layer", test_ball_thicken_one_layer);
    run("test_s3_selfglued_thicken_one_layer",
        test_s3_selfglued_thicken_one_layer);
    run("test_s3_selfglued_thicken_two_layers",
        test_s3_selfglued_thicken_two_layers);
    run("test_doubled_tetrahedra_thicken_one_layer",
        test_doubled_tetrahedra_thicken_one_layer);
    run("test_doubled_tetrahedra_thicken_two_layers",
        test_doubled_tetrahedra_thicken_two_layers);
    run("test_bdy_delta4_thicken_one_layer", test_bdy_delta4_thicken_one_layer);
    run("test_all_facet_pairs_two_tets", test_all_facet_pairs_two_tets);
    run("test_isOrdered_checks_every_facet", test_isOrdered_checks_every_facet);
    run("test_dim2_doubled_triangle", test_dim2_doubled_triangle);
    run("test_base_boundary_component_identifies_bottom_not_top",
        test_base_boundary_component_is_the_bottom_not_top);
    run("test_build_ambient_seed_in_index_order",
        test_build_ambient_seed_in_index_order);

    std::cout << "\n"
              << bold << (failed_count > 0 ? red : green) << "=== " << passed
              << " passed, " << failed_count << " failed ===" << resetColor
              << "\n";
    return failed_count > 0 ? 1 : 0;
}
