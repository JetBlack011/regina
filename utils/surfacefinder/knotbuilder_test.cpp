// knotbuilder_test.cpp
// Tests for CobordismBuilder::thicken() — verifies SimplicialPrism construction
// and gluing produce a valid product cobordism.

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <triangulation/example3.h>
#include <triangulation/example4.h>
#include <iostream>

#include "knotbuilder.h"

static int passed = 0, failed_count = 0;

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << "  PASS: " << (desc) << "\n";                        \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << "  FAIL: " << (desc) << "\n"                         \
                      << "        expected " << _e << ", got " << _a << "\n"; \
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

    knotbuilder::CobordismBuilder<3> cob(ball);
    auto& result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 4,
              "1 tet × 4 prism simplices = 4 pentachora");
    EXPECT_EQ(result.isValid(), true, "B⁴ triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "B⁴ triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 1,
              "∂(B³×I) ≅ S³: 1 boundary component");
}

// ─────────────────────────────────────────────────────────────────────────────
// Minimal S³ (2 tetrahedra sharing all 4 faces), one layer.
//
// glue() is called 4 times (once per shared face-pair).  S³ × I:
//   - 2 tets × 4 pentachora = 8 pentachora
//   - 2 boundary components, each ≅ S³
// ─────────────────────────────────────────────────────────────────────────────
void test_s3_thicken_one_layer() {
    std::cout << "\n--- Minimal S³ (2 tets), thicken(1): 4 glue() calls ---\n";

    auto s3 = regina::Example<3>::threeSphere();

    knotbuilder::CobordismBuilder<3> cob(s3);
    auto& result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 8,
              "2 tets × 4 prism simplices = 8 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×I triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "S³×I triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
              true, "bc(0) ≅ S³");
    EXPECT_EQ(result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
              true, "bc(1) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// Minimal S³, two layers.
//
// Second call to thicken_() adds a new prism layer and must glue it to the
// first layer along their shared boundary (bc 1 ↔ bc 2).  S³ × [0,2]:
//   - 16 pentachora
//   - 2 boundary components, each ≅ S³
// ─────────────────────────────────────────────────────────────────────────────
void test_s3_thicken_two_layers() {
    std::cout << "\n--- Minimal S³ (2 tets), thicken(2): inter-layer gluing ---\n";

    auto s3 = regina::Example<3>::threeSphere();

    knotbuilder::CobordismBuilder<3> cob(s3);
    auto& result = cob.thicken(2);

    EXPECT_EQ((int)result.size(), 16,
              "2 layers × 8 pentachora = 16 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×[0,2] triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "S³×[0,2] triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
              true, "bc(0) ≅ S³");
    EXPECT_EQ(result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
              true, "bc(1) ≅ S³");
}

// ─────────────────────────────────────────────────────────────────────────────
// ∂Δ⁴ as S³ (5 tetrahedra), one layer.
//
// In this triangulation every tet is adjacent to every other, so glue() is
// called 10 times.  This stresses the glue() permutation logic more than the
// 2-tet case.  S³ × I:
//   - 5 tets × 4 pentachora = 20 pentachora
//   - 2 boundary components, each ≅ ∂Δ⁴ S³
// ─────────────────────────────────────────────────────────────────────────────
void test_bdy_delta4_thicken_one_layer() {
    std::cout << "\n--- ∂Δ⁴ as S³ (5 tets), thicken(1): 10 glue() calls ---\n";

    regina::Triangulation<4> fourBall;
    fourBall.newSimplex();
    auto s3 = fourBall.boundaryComponent(0)->build();

    std::cout << "    " << s3.size() << " tetrahedra\n";

    knotbuilder::CobordismBuilder<3> cob(s3);
    auto& result = cob.thicken(1);

    EXPECT_EQ((int)result.size(), 20,
              "5 tets × 4 prism simplices = 20 pentachora");
    EXPECT_EQ(result.isValid(), true, "S³×I triangulation is valid");
    EXPECT_EQ(result.isConnected(), true, "S³×I triangulation is connected");
    EXPECT_EQ((int)result.countBoundaryComponents(), 2,
              "2 boundary components (top and bottom S³)");
    EXPECT_EQ(result.boundaryComponent(0)->build().isIsomorphicTo(s3).has_value(),
              true, "bc(0) ≅ S³");
    EXPECT_EQ(result.boundaryComponent(1)->build().isIsomorphicTo(s3).has_value(),
              true, "bc(1) ≅ S³");
}

template <typename F>
void run(const char* name, F fn) {
    std::cout << "\nRunning " << name << "...\n";
    try {
        fn();
    } catch (const std::exception& e) {
        std::cout << "  EXCEPTION: " << e.what() << "\n";
        ++failed_count;
    }
}

int main() {
    run("test_ball_thicken_one_layer",    test_ball_thicken_one_layer);
    run("test_s3_thicken_one_layer",      test_s3_thicken_one_layer);
    run("test_s3_thicken_two_layers",     test_s3_thicken_two_layers);
    run("test_bdy_delta4_thicken_one_layer", test_bdy_delta4_thicken_one_layer);

    std::cout << "\n=== " << passed << " passed, " << failed_count
              << " failed ===\n";
    return failed_count > 0 ? 1 : 0;
}
