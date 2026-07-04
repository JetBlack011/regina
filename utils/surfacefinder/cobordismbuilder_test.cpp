// cobordismbuilder_test.cpp
// Tests for CobordismBuilder::thicken() and cone() — verifies SimplicialPrism
// construction and gluing produce a valid product cobordism and that it can
// be capped into a ball, including integration with knotbuilder's PD-code
// -> triangulated-S³ pipeline.

#include <iostream>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>

#include "cobordismbuilder.h"
#include "knotbuilder.h"
#include "knottedsurfaces.h"

static int passed = 0, failed_count = 0;

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << "  PASS: " << (desc) << "\n";                         \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << "  FAIL: " << (desc) << "\n"                          \
                      << "        expected " << _e << ", got " << _a << "\n";  \
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

    EXPECT_EQ((int)result.size(), 4, "1 tet × 4 prism simplices = 4 pentachora");
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

    EXPECT_EQ((int)result.size(), 8, "2 tets × 4 prism simplices = 8 pentachora");
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

    EXPECT_EQ((int)result.size(), 16, "2 layers × 8 pentachora = 16 pentachora");
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
// cone() on its own (no thicken() calls first): Cone(S³) should be a valid
// B⁴ with a single boundary component ≅ the original S³, and one pentachoron
// per base tetrahedron (2 tets -> 2 pentachora). eulerCharManifold() == 1
// is an independent check that this is genuinely contractible, i.e.
// actually a ball and not just "some valid closed-up thing".
// ─────────────────────────────────────────────────────────────────────────────
void test_cone_alone() {
    std::cout << "\n--- cone() alone, no thicken(): Cone(S³) ---\n";

    auto s3 = doubledTetrahedra();

    CobordismBuilder<3> cob(s3);
    auto &result = cob.cone();

    EXPECT_EQ((int)result.size(), 2, "1 cone pentachoron per base tet = 2");
    EXPECT_EQ(result.isValid(), true, "Cone(S³) triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "Cone(S³) triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              "Cone(S³) has exactly 1 boundary component");
    EXPECT_EQ(result.eulerCharManifold(), 1L,
              "Cone(S³) is contractible (Euler characteristic 1)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "bc(0) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// Regression test: cone() previously decided whether to "start fresh" or
// "glue onto existing content" by checking cob_.isConnected() — but a
// cobordism that already has one or more thickened layers is *also*
// connected, so this check incorrectly took the "start fresh" branch and
// silently discarded everything thicken() had built (cone() would just
// overwrite cob_ with a bare 2-pentachoron cone, losing the 8 pentachora
// from thicken(1)). The fix uses capTop(), which glues the cone directly
// onto the most recent layer the same way stitchTop() glues layers to each
// other. Checking the *size* here is essential: a size of 2 instead of 10
// is exactly the failure signature of the original bug.
// ─────────────────────────────────────────────────────────────────────────────
void test_cone_after_one_thicken() {
    std::cout << "\n--- thicken(1) then cone(): must not discard the "
                 "thickened layer ---\n";

    auto s3 = doubledTetrahedra();

    CobordismBuilder<3> cob(s3);
    cob.thicken(1);
    auto &result = cob.cone();

    EXPECT_EQ((int)result.size(), 10,
              "8 pentachora (thicken) + 2 (cone) = 10, not just the cone's 2");
    EXPECT_EQ(result.isValid(), true, "thicken(1)+cone() triangulation is valid");
    EXPECT_EQ(result.isConnected(), true,
              "thicken(1)+cone() triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              "capping the top leaves exactly 1 boundary component (the "
              "original bottom S³)");
    EXPECT_EQ(result.eulerCharManifold(), 1L,
              "thicken(1)+cone() is contractible (Euler characteristic 1)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "the remaining boundary component ≅ the original S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// Multiple thickenings capped by a single cone() — this is the actual shape
// of the pipeline the knotbuilder integration tests below use: thicken(n)
// to get some combinatorial room, then cone() once at the end to cap the
// cobordism into a ball whose boundary is the untouched bottom layer.
// ─────────────────────────────────────────────────────────────────────────────
void test_cone_after_multiple_thickens() {
    std::cout << "\n--- thicken(3) then cone() ---\n";

    auto s3 = doubledTetrahedra();

    CobordismBuilder<3> cob(s3);
    cob.thicken(3);
    auto &result = cob.cone();

    EXPECT_EQ((int)result.size(), 26,
              "3 layers × 8 pentachora + 2 (cone) = 26");
    EXPECT_EQ(result.isValid(), true, "thicken(3)+cone() triangulation is valid");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              "thicken(3)+cone() has exactly 1 boundary component");
    EXPECT_EQ(result.eulerCharManifold(), 1L,
              "thicken(3)+cone() is contractible (Euler characteristic 1)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
        true, "the remaining boundary component ≅ the original S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// thicken()+cone() on the self-glued-tetrahedron S³ and on the 5-tet ∂Δ⁴ S³,
// to make sure capTop() (like glue() and stitchTop() before it) holds up on
// triangulations that aren't the simple 2-distinct-tets case: self-gluing,
// and a mix of matching/non-matching facet pairs.
// ─────────────────────────────────────────────────────────────────────────────
void test_cone_after_thicken_selfglued_and_bdy_delta4() {
    std::cout
        << "\n--- thicken()+cone() on the self-glued tet and ∂Δ⁴ S³s ---\n";

    {
        auto s3 = regina::Example<3>::threeSphere();
        CobordismBuilder<3> cob(s3);
        cob.thicken(2);
        auto &result = cob.cone();

        EXPECT_EQ((int)result.size(), 9,
                  "self-glued tet: 2 layers × 4 + 1 (cone) = 9 pentachora");
        EXPECT_EQ(result.isValid(), true,
                  "self-glued tet: thicken(2)+cone() is valid");
        EXPECT_EQ((int)result.countBoundaryComponents(), 1,
                  "self-glued tet: thicken(2)+cone() has 1 boundary component");
        EXPECT_EQ(result.eulerCharManifold(), 1L,
                  "self-glued tet: thicken(2)+cone() is contractible");
        EXPECT_EQ(
            result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
            true, "self-glued tet: remaining boundary ≅ S³");
    }

    {
        regina::Triangulation<4> fourBall;
        fourBall.newSimplex();
        auto s3 = fourBall.boundaryComponent(0)->build();

        CobordismBuilder<3> cob(s3);
        cob.thicken(1);
        auto &result = cob.cone();

        EXPECT_EQ((int)result.size(), 25,
                  "∂Δ⁴: 1 layer × 20 + 5 (cone) = 25 pentachora");
        EXPECT_EQ(result.isValid(), true, "∂Δ⁴: thicken(1)+cone() is valid");
        EXPECT_EQ((int)result.countBoundaryComponents(), 1,
                  "∂Δ⁴: thicken(1)+cone() has 1 boundary component");
        EXPECT_EQ(result.eulerCharManifold(), 1L,
                  "∂Δ⁴: thicken(1)+cone() is contractible");
        EXPECT_EQ(
            result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
            true, "∂Δ⁴: remaining boundary ≅ S³");
    }
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

    EXPECT_EQ(numOk, 16,
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
    EXPECT_EQ(result.isConnected(), true, "S²×[0,2] triangulation is connected");
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
// cone() in dim=2 as well, so capTop() (like glue() and stitchTop()) is
// checked across dimensions rather than only for the dim=4 pentachora used
// everywhere else in this file.
// ─────────────────────────────────────────────────────────────────────────────
void test_dim2_cone_after_thicken() {
    std::cout << "\n--- Doubled triangle (S²), thicken(2)+cone(): "
                 "CobordismBuilder<2> capTop() ---\n";

    regina::Triangulation<2> tri;
    auto a = tri.newTriangle();
    auto b = tri.newTriangle();
    for (int f = 0; f < 3; ++f)
        a->join(f, b, regina::Perm<3>());

    CobordismBuilder<2> cob(tri);
    cob.thicken(2);
    auto &result = cob.cone();

    EXPECT_EQ((int)result.size(), 14,
              "2 layers × 6 tets + 2 (cone) = 14 tetrahedra");
    EXPECT_EQ(result.isValid(), true, "S²: thicken(2)+cone() is valid");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              "S²: thicken(2)+cone() has exactly 1 boundary component");
    EXPECT_EQ(result.eulerCharManifold(), 1L,
              "S²: thicken(2)+cone() is contractible (Euler characteristic 1)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(tri).has_value(),
        true, "S²: remaining boundary component ≅ S²");
}

// ═════════════════════════════════════════════════════════════════════════
// Knotbuilder integration.
//
// PD codes below use the standard KnotInfo/Knot Atlas convention: one
// X(a,b,c,d) per crossing (1-indexed), flattened into a single list.
// They are not just trusted: each is checked to actually produce a valid,
// closed triangulation of S³ (isSphere()), and Link is checked to find the
// expected number of components from the returned link edges, before being
// used in any CobordismBuilder pipeline test.
// ═════════════════════════════════════════════════════════════════════════
const char *TREFOIL_PD = "1 4 2 5 3 6 4 1 5 2 6 3";              // 3_1
const char *HOPF_LINK_PD = "1 4 2 3 3 2 4 1";                    // 2-component link
const char *FIGURE_EIGHT_PD = "4 2 5 1 8 6 1 5 6 3 7 4 2 7 3 8"; // 4_1

regina::Triangulation<3> buildFromPD(const char *pd,
                                     std::vector<const regina::Edge<3> *> &edges) {
    regina::Triangulation<3> tri;
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pd);
    if (!knotbuilder::buildLink(tri, pdcode, edges))
        throw regina::InvalidArgument("buildFromPD: buildLink() failed");
    return tri;
}

// ─────────────────────────────────────────────────────────────────────────────
// Layer 1: knotbuilder's raw output, on its own, independent of
// CobordismBuilder. If buildLink() ever stops producing a valid closed S³
// (e.g. a knotbuilder change breaks the block gluing pattern), this should
// fail here rather than surfacing as a confusing CobordismBuilder/glue()
// failure downstream.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_trefoil_is_valid_s3() {
    std::cout << "\n--- knotbuilder: trefoil PD code -> closed S³ ---\n";

    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(TREFOIL_PD, edges);

    EXPECT_EQ(tri.isValid(), true, "trefoil triangulation is valid");
    EXPECT_EQ(tri.isClosed(), true, "trefoil triangulation is closed");
    EXPECT_EQ(tri.isIdeal(), false, "trefoil triangulation has no ideal vertices");
    EXPECT_EQ(tri.isSphere(), true, "trefoil triangulation is S³");
    EXPECT_EQ((int)edges.size(), 9, "3 crossings × 3 link edges per block = 9");

    Link link(tri, edges);
    EXPECT_EQ((int)link.comps_.size(), 1, "trefoil is a 1-component knot");
}

void test_knotbuilder_hopf_link_two_components() {
    std::cout << "\n--- knotbuilder: Hopf link PD code -> 2-component link ---\n";

    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(HOPF_LINK_PD, edges);

    EXPECT_EQ(tri.isValid(), true, "Hopf link triangulation is valid");
    EXPECT_EQ(tri.isSphere(), true, "Hopf link triangulation is S³");

    Link link(tri, edges);
    EXPECT_EQ((int)link.comps_.size(), 2, "Hopf link has 2 components");
}

// ─────────────────────────────────────────────────────────────────────────────
// Regression test: this specific diagram is what originally exposed the bug
// in Block::glue(). It's the Hopf link's shadow with crossing 1's tuple
// rotated by one position, which flips that crossing's over/under role and
// makes the diagram non-alternating (both transitions become same-parity
// instead of alternating). The rotated diagram is the 2-component unlink,
// not the Hopf link.
//
// The old parity-based glue() produced a triangulation that looked fine
// (isValid() == true, no non-manifold points) but was structurally wrong:
// the two S²×I boundary components came out as tori instead of spheres, so
// finiteToIdeal() turned them into genuine cusps (isIdeal() == true) instead
// of closing up to S³. isSphere() == false was the tell. The current
// glue() (matching walls by role -- own/anchor/private -- rather than by
// wall-index parity) has no such branch to get wrong, but this is kept as a
// dedicated, minimal, always-run check specifically for this failure mode.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_nonalternating_regression() {
    std::cout << "\n--- knotbuilder: non-alternating regression (flipped "
                 "Hopf shadow -> unlink) ---\n";

    // Same 4-valent shadow as the Hopf link (crossing 0 unchanged), but
    // crossing 1's tuple [2,1,3,0] is cyclically rotated by one position to
    // [1,3,0,2], flipping which pair of arms is the under-strand.
    knotbuilder::PDCode pd = {{0, 3, 1, 2}, {1, 3, 0, 2}};

    regina::Triangulation<3> tri;
    std::vector<const regina::Edge<3> *> edges;
    bool ok = knotbuilder::buildLink(tri, pd, edges);

    EXPECT_EQ(ok, true, "non-alternating shadow: buildLink() succeeds");
    EXPECT_EQ(tri.isValid(), true,
              "non-alternating shadow: triangulation is valid");
    EXPECT_EQ(tri.isClosed(), true,
              "non-alternating shadow: triangulation is closed");
    EXPECT_EQ(tri.isIdeal(), false,
              "non-alternating shadow: no ideal vertices (an ideal vertex "
              "here would mean a torus boundary wasn't closed into a "
              "sphere before capping)");
    EXPECT_EQ(tri.isSphere(), true,
              "non-alternating shadow: triangulation is S³");

    Link link(tri, edges);
    EXPECT_EQ((int)link.comps_.size(), 2,
              "non-alternating shadow: 2 components (the unlink)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Broader robustness sweep: real PD codes for a range of named knots (3 to
// 13 crossings, both alternating and non-alternating -- the "a"/"n" suffix
// in the 11+ crossing names is KnotInfo's own alternating/non-alternating
// classification), from KnotInfo (https://knotinfo.org/). Each is checked
// for a valid, closed S³ triangulation with the correct single-component
// count, and its complement's homology is cross-checked against Regina's
// own independent Link::fromPD() construction of the same PD code -- a
// cheap (no exponential search needed) but genuine correctness signal
// beyond "didn't crash", especially for the larger examples where an exact
// isomorphism/retriangulation check would be too slow to run routinely.
// ─────────────────────────────────────────────────────────────────────────────
struct NamedKnot {
    const char *name;
    const char *pd;
};

const NamedKnot NAMED_KNOTS[] = {
    {"3_1", "1,5,2,4 3,1,4,6 5,3,6,2"},
    {"4_1", "4,2,5,1 8,6,1,5 6,3,7,4 2,7,3,8"},
    {"5_1", "2,8,3,7 4,10,5,9 6,2,7,1 8,4,9,3 10,6,1,5"},
    {"5_2", "1,5,2,4 3,9,4,8 5,1,6,10 7,3,8,2 9,7,10,6"},
    {"6_1", "1,7,2,6 3,10,4,11 5,3,6,2 7,1,8,12 9,4,10,5 11,9,12,8"},
    {"6_2", "1,8,2,9 3,11,4,10 5,1,6,12 7,2,8,3 9,7,10,6 11,5,12,4"},
    {"6_3", "4,2,5,1 8,4,9,3 12,9,1,10 10,5,11,6 6,11,7,12 2,8,3,7"},
    {"7_4",
     "2,10,3,9 4,12,5,11 6,14,7,13 8,4,9,3 10,2,11,1 12,8,13,7 14,6,1,5"},
    {"8_19 (non-alternating)",
     "2,14,3,13 5,11,6,10 7,15,8,14 9,5,10,4 11,7,12,6 12,2,13,1 15,9,16,8 "
     "16,4,1,3"},
    {"8_20 (non-alternating)",
     "1,7,2,6 4,13,5,14 5,9,6,8 7,3,8,2 10,15,11,16 12,9,13,10 14,3,15,4 "
     "16,11,1,12"},
    {"8_21 (non-alternating)",
     "1,7,2,6 4,13,5,14 5,9,6,8 7,3,8,2 9,13,10,12 11,1,12,16 14,3,15,4 "
     "15,11,16,10"},
    {"10_165",
     "1,8,2,9 3,12,4,13 4,17,5,18 7,2,8,3 9,14,10,15 11,17,12,16 13,6,14,7 "
     "15,20,16,1 18,5,19,6 19,11,20,10"},
    {"11a_1",
     "4,2,5,1 10,6,11,5 8,3,9,4 2,9,3,10 16,12,17,11 14,7,15,8 6,15,7,16 "
     "20,14,21,13 22,18,1,17 18,22,19,21 12,20,13,19"},
    {"11n_1 (non-alternating)",
     "4,2,5,1 10,6,11,5 8,3,9,4 2,9,3,10 11,16,12,17 14,7,15,8 6,15,7,16 "
     "13,20,14,21 17,22,18,1 21,18,22,19 19,12,20,13"},
    {"12a_1",
     "1,5,2,4 3,8,4,9 5,11,6,10 7,14,8,15 9,2,10,3 11,17,12,16 13,21,14,20 "
     "15,6,16,7 17,22,18,23 19,13,20,12 21,24,22,1 23,18,24,19"},
    {"12n_1 (non-alternating)",
     "1,5,2,4 3,8,4,9 5,11,6,10 14,8,15,7 9,2,10,3 16,11,17,12 20,13,21,14 "
     "6,16,7,15 22,18,23,17 12,19,13,20 24,22,1,21 18,24,19,23"},
    {"13a_1",
     "3,1,4,26 1,8,2,9 7,2,8,3 9,5,10,4 5,14,6,15 13,6,14,7 15,11,16,10 "
     "11,20,12,21 19,12,20,13 21,17,22,16 17,24,18,25 23,18,24,19 "
     "25,23,26,22"},
    {"13n_1 (non-alternating)",
     "3,1,4,26 1,8,2,9 7,2,8,3 9,5,10,4 14,6,15,5 6,14,7,13 10,15,11,16 "
     "20,12,21,11 12,20,13,19 16,21,17,22 24,18,25,17 18,24,19,23 "
     "22,25,23,26"},
};

void test_knotbuilder_many_named_knots() {
    std::cout << "\n--- knotbuilder: sweep over named knots (3-13 "
                 "crossings, alternating and non-alternating), from "
                 "KnotInfo ---\n";

    for (const auto &knot : NAMED_KNOTS) {
        std::vector<const regina::Edge<3> *> edges;
        auto tri = buildFromPD(knot.pd, edges);

        EXPECT_EQ(tri.isValid(), true,
                  std::string(knot.name) + ": triangulation is valid");
        EXPECT_EQ(tri.isClosed(), true,
                  std::string(knot.name) + ": triangulation is closed");
        EXPECT_EQ(tri.isSphere(), true,
                  std::string(knot.name) + ": triangulation is S³");

        Link link(tri, edges);
        EXPECT_EQ((int)link.comps_.size(), 1,
                  std::string(knot.name) + ": is a 1-component knot");

        // Cross-check against Regina's own independent Link::fromPD()
        // pipeline for the same PD code: the two complements' homology
        // (H_1 of a knot complement is always Z, so this mainly catches
        // gross errors, but it's cheap enough to run on every case,
        // including the larger ones where an exact isomorphism check via
        // retriangulation would be too slow).
        auto mine = link.buildComplement();
        regina::Link reginaLink = regina::Link::fromPD(knot.pd);
        auto reginaComp = reginaLink.complement();
        EXPECT_EQ(mine.homology() == reginaComp.homology(), true,
                  std::string(knot.name) +
                      ": complement homology matches Regina's independent "
                      "Link::fromPD() construction");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Layer 2: the specific property the user asked to check. knotbuilder's raw
// output is *not* ordered (block construction doesn't aim for that), but it
// must be orderable, since that's what CobordismBuilder<3>'s constructor
// relies on (it calls Triangulation<3>::order() itself when needed — see
// cobordismbuilder.h). Recording the raw isOrdered() value directly here,
// rather than assuming it, means a future knotbuilder change that breaks
// orderability entirely (as opposed to just not ordering by default) gets
// caught right here instead of as a mysterious CobordismBuilder exception.
// ─────────────────────────────────────────────────────────────────────────────
void test_knotbuilder_output_is_orderable() {
    std::cout
        << "\n--- knotbuilder output: raw ordering + orderability ---\n";

    const std::pair<const char *, const char *> knots[] = {
        {"trefoil", TREFOIL_PD},
        {"Hopf link", HOPF_LINK_PD},
        {"figure-8", FIGURE_EIGHT_PD},
    };

    for (const auto &[name, pd] : knots) {
        std::vector<const regina::Edge<3> *> edges;
        auto tri = buildFromPD(pd, edges);

        std::cout << "  " << name
                   << ": raw isOrdered() = " << tri.isOrdered() << "\n";

        regina::Triangulation<3> ordered(tri);
        EXPECT_EQ(ordered.order(), true,
                  std::string(name) + ": knotbuilder's output is orderable");
        EXPECT_EQ(ordered.isOrdered(), true,
                  std::string(name) +
                      ": order() actually produces an ordered triangulation");

        // The property that matters end to end: CobordismBuilder<3>'s
        // constructor requires isOrdered() and calls order() itself, so it
        // must accept knotbuilder's raw (unordered) output without throwing.
        bool constructedOk = true;
        try {
            CobordismBuilder<3> cob(tri);
        } catch (const std::exception &) {
            constructedOk = false;
        }
        EXPECT_EQ(constructedOk, true,
                  std::string(name) +
                      ": CobordismBuilder<3> accepts knotbuilder's raw output");
    }
}

// ─────────────────────────────────────────────────────────────────────────────
// Layer 3: the full pipeline requested — knotbuilder builds S³ with the
// link as an edge circuit, CobordismBuilder thickens it (giving room for
// surfacefinder to later search in), then cone() caps the top into a
// genuine triangulated 4-ball. The one boundary component left afterward
// must be combinatorially identical to knotbuilder's original S³ — that's
// what will let a future surfacefinder test locate the link's edges inside
// the boundary of this B⁴.
// ─────────────────────────────────────────────────────────────────────────────
void checkKnotToBallPipeline(const char *name, const char *pd, int layers) {
    std::vector<const regina::Edge<3> *> edges;
    auto tri = buildFromPD(pd, edges);

    CobordismBuilder<3> cob(tri);
    if (layers > 0)
        cob.thicken(layers);
    auto &result = cob.cone();

    EXPECT_EQ(result.isValid(),
              true, std::string(name) + ": thicken(" + std::to_string(layers) +
                        ")+cone() triangulation is valid");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              std::string(name) + ": thicken(" + std::to_string(layers) +
                  ")+cone() has exactly 1 boundary component");
    EXPECT_EQ(result.eulerCharManifold(), 1L,
              std::string(name) + ": thicken(" + std::to_string(layers) +
                  ")+cone() is contractible (Euler characteristic 1)");
    EXPECT_EQ(
        result.boundaryComponent(0)->build().isIsomorphicTo(tri).has_value(),
        true,
        std::string(name) + ": remaining boundary ≅ knotbuilder's original S³");
}

void test_knotbuilder_trefoil_to_ball_pipeline() {
    std::cout << "\n--- knotbuilder trefoil -> thicken(n)+cone() -> B⁴ ---\n";
    checkKnotToBallPipeline("trefoil (0 layers)", TREFOIL_PD, 0);
    checkKnotToBallPipeline("trefoil (1 layer)", TREFOIL_PD, 1);
    checkKnotToBallPipeline("trefoil (3 layers)", TREFOIL_PD, 3);
}

void test_knotbuilder_hopf_link_to_ball_pipeline() {
    std::cout << "\n--- knotbuilder Hopf link -> thicken(n)+cone() -> B⁴ ---\n";
    checkKnotToBallPipeline("Hopf link (1 layer)", HOPF_LINK_PD, 1);
    checkKnotToBallPipeline("Hopf link (2 layers)", HOPF_LINK_PD, 2);
}

void test_knotbuilder_figure_eight_to_ball_pipeline() {
    std::cout << "\n--- knotbuilder figure-8 -> thicken(n)+cone() -> B⁴ ---\n";
    checkKnotToBallPipeline("figure-8 (1 layer)", FIGURE_EIGHT_PD, 1);
    checkKnotToBallPipeline("figure-8 (2 layers)", FIGURE_EIGHT_PD, 2);
}

template <typename F>
void run(const char *name, F fn) {
    std::cout << "\nRunning " << name << "...\n";
    try {
        fn();
    } catch (const std::exception &e) {
        std::cout << "  EXCEPTION: " << e.what() << "\n";
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
    run("test_cone_alone", test_cone_alone);
    run("test_cone_after_one_thicken", test_cone_after_one_thicken);
    run("test_cone_after_multiple_thickens", test_cone_after_multiple_thickens);
    run("test_cone_after_thicken_selfglued_and_bdy_delta4",
        test_cone_after_thicken_selfglued_and_bdy_delta4);
    run("test_bdy_delta4_thicken_one_layer", test_bdy_delta4_thicken_one_layer);
    run("test_all_facet_pairs_two_tets", test_all_facet_pairs_two_tets);
    run("test_isOrdered_checks_every_facet", test_isOrdered_checks_every_facet);
    run("test_dim2_doubled_triangle", test_dim2_doubled_triangle);
    run("test_dim2_cone_after_thicken", test_dim2_cone_after_thicken);
    run("test_knotbuilder_trefoil_is_valid_s3",
        test_knotbuilder_trefoil_is_valid_s3);
    run("test_knotbuilder_hopf_link_two_components",
        test_knotbuilder_hopf_link_two_components);
    run("test_knotbuilder_nonalternating_regression",
        test_knotbuilder_nonalternating_regression);
    run("test_knotbuilder_many_named_knots", test_knotbuilder_many_named_knots);
    run("test_knotbuilder_output_is_orderable",
        test_knotbuilder_output_is_orderable);
    run("test_knotbuilder_trefoil_to_ball_pipeline",
        test_knotbuilder_trefoil_to_ball_pipeline);
    run("test_knotbuilder_hopf_link_to_ball_pipeline",
        test_knotbuilder_hopf_link_to_ball_pipeline);
    run("test_knotbuilder_figure_eight_to_ball_pipeline",
        test_knotbuilder_figure_eight_to_ball_pipeline);

    std::cout << "\n=== " << passed << " passed, " << failed_count
              << " failed ===\n";
    return failed_count > 0 ? 1 : 0;
}
