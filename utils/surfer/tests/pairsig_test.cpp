// pairsig_test.cpp
//
// Tests for pairSig()/fromPairSig() (see ../pairsig.h): a signature that
// encodes both an ambient triangulation and a marked subcomplex of its
// subdim-skeleton, generalizing Regina's isomorphism signature to pairs.
//
// The property that matters most here -- and the reason a naive "isoSig(T)
// + isoSigDetail()'s arbitrarily-relabeled face list" approach was rejected
// during design -- is that pairSig() must be a genuine isomorphism
// invariant: two isomorphic (ambient, marked subcomplex) pairs must produce
// byte-identical signatures, even when the ambient triangulation has
// automorphisms that move the marked subcomplex around. Tests 2 and 3 below
// exercise exactly that.
//
// Note: several tests below embed the literal delimiter character '_'
// (pairSig()'s internal Base64Encoder::spare[0]) directly into hand-built
// malformed signature strings. This couples the tests to that specific
// choice of delimiter -- acceptable here since the point of those tests is
// to pin down fromPairSig()'s error handling given an already-known-bad
// string, not to re-derive the delimiter from scratch.

#include <iostream>
#include <string>
#include <vector>
#include <unistd.h>

#include <maths/perm.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "../pairsig.h"

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

// Maps ambient face `f` of `src` to its image face index in `dst`, under
// the isomorphism `iso : src -> dst`. This is a standalone re-derivation of
// the same "front().vertices() + FaceNumbering::faceNumber" technique
// pairSig() itself uses internally -- written independently here so that
// tests 2 and 3 can construct a *known* relabeled marked set without
// depending on pairSig()'s own automorphism-minimization machinery.
template <int dim, int subdim>
size_t mapFace(const regina::Triangulation<dim> &src,
                const regina::Triangulation<dim> &dst,
                const regina::Isomorphism<dim> &iso, int f) {
    const auto &emb = src.template face<subdim>(f)->front();
    size_t srcSimplex = emb.simplex()->index();
    regina::Perm<dim + 1> p = emb.vertices();
    auto destSimplex = static_cast<size_t>(iso.simpImage(srcSimplex));
    regina::Perm<dim + 1> q = iso.facetPerm(srcSimplex) * p;
    int localFace = regina::FaceNumbering<dim, subdim>::faceNumber(q);
    return dst.simplex(destSimplex)->template face<subdim>(localFace)->index();
}

} // namespace

// ─────────────────────────────────────────────────────────────────────────────
// Test 1: encode -> decode -> re-encode reproduces the same signature.
// ─────────────────────────────────────────────────────────────────────────────
void test_round_trip_identity() {
    std::cout << "\n--- round-trip identity ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();
    Skeleton<3, 2> skel(ball);
    EmbeddedSubmanifold<3, 2> sub(skel, {0, 1});

    std::string sig1 = pairSig<3, 2>(skel, sub);
    auto decoded = fromPairSig<3, 2>(sig1);
    std::string sig2 = pairSig<3, 2>(*decoded.skeleton, *decoded.submanifold);

    EXPECT_EQ(sig1, sig2,
              "re-encoding a decoded pair signature reproduces the same "
              "string");
    EXPECT_EQ(decoded.submanifold->markedFaces().size(), (size_t)2,
              "decoded submanifold has 2 marked faces");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 2: pairSig is invariant under an arbitrary relabelling of the
// ambient triangulation (and a correspondingly relabelled marked set).
// This is the property a naive "isoSigDetail() + relabeled indices" scheme
// would already satisfy trivially for a *single* relabelling -- the real
// test is Test 3 below, which additionally varies over automorphisms of a
// single fixed triangulation.
// ─────────────────────────────────────────────────────────────────────────────
void test_isomorphism_invariance_under_relabeling() {
    std::cout << "\n--- isomorphism invariance under relabeling ---\n";

    regina::Triangulation<3> tri;
    auto *t0 = tri.newTetrahedron();
    auto *t1 = tri.newTetrahedron();
    t0->join(0, t1, regina::Perm<4>());

    std::vector<int> marked = {
        static_cast<int>(t0->triangle(1)->index()),
        static_cast<int>(t0->triangle(2)->index())};

    auto randomIso = regina::Isomorphism<3>::random(tri.size());
    regina::Triangulation<3> tri2 = randomIso(tri);

    std::vector<int> marked2;
    for (int f : marked)
        marked2.push_back(
            static_cast<int>(mapFace<3, 2>(tri, tri2, randomIso, f)));

    std::string sig1 = pairSig<3, 2>(tri, marked);
    std::string sig2 = pairSig<3, 2>(tri2, marked2);
    EXPECT_EQ(sig1, sig2,
              "a randomly relabeled copy of the same pair produces an "
              "identical signature");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 3: the real test of the automorphism-minimization step. A single,
// ungued tetrahedron has a full S4 automorphism group (any permutation of
// its 4 vertices is an automorphism), so marking one boundary triangle and
// then applying a nontrivial automorphism moves the mark to a *different*
// ambient face index -- yet the pair is unmistakably isomorphic to itself.
// Without minimizing over Aut(canon), pairSig() could pick a different
// arbitrary relabeling for each and produce two different signatures.
// ─────────────────────────────────────────────────────────────────────────────
void test_non_automorphism_invariant_marked_set() {
    std::cout << "\n--- marked set moved by ambient's own automorphism ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();

    regina::Isomorphism<3> beta(1);
    bool found = false;
    ball.findAllIsomorphisms(ball,
        [&](const regina::Isomorphism<3> &alpha) {
            if (mapFace<3, 2>(ball, ball, alpha, 0) != 0) {
                beta = alpha;
                found = true;
                return true; // stop as soon as one is found
            }
            return false;
        });
    EXPECT_EQ(found, true,
              "found a nontrivial automorphism moving marked face 0");
    if (!found)
        return;

    std::vector<int> M = {0};
    auto f2 = static_cast<int>(mapFace<3, 2>(ball, ball, beta, 0));
    std::vector<int> M2 = {f2};
    EXPECT_EQ(f2 != 0, true, "the automorphism actually moves the marked face");

    std::string sig1 = pairSig<3, 2>(ball, M);
    std::string sig2 = pairSig<3, 2>(ball, M2);
    EXPECT_EQ(sig1, sig2,
              "pairSig is invariant under the ambient triangulation's own "
              "automorphisms, not just under an arbitrary single relabeling");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 4: KnottedSurface-specific round trip, reusing the two-pentachora
// triangulation from embeddedsubmanifold_test.cpp's known-good {6,7} pair
// (see test_buildgraph_known_incompleteness there).
// ─────────────────────────────────────────────────────────────────────────────
void test_knotted_surface_round_trip() {
    std::cout << "\n--- KnottedSurface round trip ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    Skeleton<4, 2> skeleton(tri);
    KnottedSurface surface(skeleton);
    bool added6 = surface.addFace(6);
    bool added7 = surface.addFace(7);
    EXPECT_EQ(added6 && added7, true,
              "the known-good {6,7} pair is accepted");
    if (!(added6 && added7))
        return;

    std::string sig = pairSig<4, 2>(skeleton, surface);
    auto decoded = fromKnottedSurfaceSig(sig);

    // decoded.surface's marked faces are indices into fromSig(sig)'s
    // canonical reconstruction, not into `tri` -- generally a different
    // (though isomorphic) numbering, so the round-trip check here is
    // re-encoding (as in test_round_trip_identity), not comparing raw
    // index lists directly.
    std::string sig2 = pairSig<4, 2>(*decoded.skeleton, *decoded.surface);
    EXPECT_EQ(sig, sig2,
              "re-encoding the decoded KnottedSurface reproduces the same "
              "signature");
    EXPECT_EQ(KnottedSurface::formatSurfaceType(decoded.surface->surfaceType()),
              KnottedSurface::formatSurfaceType(surface.surfaceType()),
              "decoded KnottedSurface has the same surface type");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 5: the empty marked set is a valid degenerate case, and its
// signature's prefix (before the delimiter) equals a plain isoSig().
// ─────────────────────────────────────────────────────────────────────────────
void test_empty_marked_set() {
    std::cout << "\n--- empty marked set ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();
    Skeleton<3, 2> skel(ball);
    EmbeddedSubmanifold<3, 2> sub(skel);

    std::string sig = pairSig<3, 2>(skel, sub);
    std::string plainIsoSig = ball.isoSig();

    auto pos = sig.find('_');
    EXPECT_EQ(pos != std::string::npos, true,
              "empty-marked-set signature contains the delimiter");
    EXPECT_EQ(sig.substr(0, pos), plainIsoSig,
              "empty marked set: signature prefix equals plain isoSig()");

    auto decoded = fromPairSig<3, 2>(sig);
    EXPECT_EQ(decoded.submanifold->markedFaces().empty(), true,
              "decoded submanifold has no marked faces");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 6: malformed signature strings throw regina::InvalidArgument rather
// than crashing or silently misbehaving.
// ─────────────────────────────────────────────────────────────────────────────
void test_malformed_input() {
    std::cout << "\n--- malformed input ---\n";

    regina::Triangulation<3> ball;
    ball.newTetrahedron();
    std::string plainSig = ball.isoSig();

    bool threw = false;
    try {
        fromPairSig<3, 2>(plainSig); // no delimiter at all
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "missing delimiter throws InvalidArgument");

    threw = false;
    try {
        fromPairSig<3, 2>(plainSig + "_abc"); // non-numeric index
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true, "non-numeric index throws InvalidArgument");

    // A triple self-fold (all 3 edges of a triangle identified together):
    // both of its 2 triangles are irreparably self-folded (see
    // embeddingsearch_test.cpp's test_triple_self_fold_excluded), so
    // whichever one fromSig()'s reconstruction happens to number 0, marking
    // it alone is not addable at all, let alone jointly.
    regina::Triangulation<3> foldTri;
    auto *t = foldTri.newTetrahedron();
    t->join(0, t, regina::Perm<4>(1, 2, 0, 3));
    t->join(2, t, regina::Perm<4>(0, 2, 3, 1));
    std::string foldSig = foldTri.isoSig();

    threw = false;
    try {
        fromPairSig<3, 2>(foldSig + "_0");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true,
              "a self-folded singleton face throws InvalidArgument (not "
              "jointly addable)");
}

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("round_trip_identity", test_round_trip_identity);
    run("isomorphism_invariance_under_relabeling",
        test_isomorphism_invariance_under_relabeling);
    run("non_automorphism_invariant_marked_set",
        test_non_automorphism_invariant_marked_set);
    run("knotted_surface_round_trip", test_knotted_surface_round_trip);
    run("empty_marked_set", test_empty_marked_set);
    run("malformed_input", test_malformed_input);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
