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
//
// The marked-face suffix itself is fixed-width base64 (the same alphabet
// and IsoSigPrintable-style scheme isoSig() uses for its own gluing data --
// see pairSig()'s implementation), not decimal digits. Test 6's hand-built
// suffixes are built accordingly: a base64 character does *not* mean what
// its face value would in decimal (e.g. the character '0' decodes to 52,
// not 0 -- 'a' is base64 for 0), and a suffix's length must be an exact
// multiple of the per-index width w = Base64Encoder::integerWidth(M - 1).

#include <iostream>
#include <atomic>
#include <string>
#include <thread>
#include <vector>
#include <unistd.h>

#include <maths/perm.h>
#include <triangulation/dim3.h>
#include <utilities/sigutils.h>
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

    // A character genuinely outside the base64 alphabet (a..zA..Z0..9+-):
    // Base64Encoder::spare[1] ('.') is documented to never be one of those
    // 64 characters, so this is guaranteed invalid regardless of `ball`'s
    // own per-index width w.
    threw = false;
    try {
        fromPairSig<3, 2>(
            plainSig + "_" + std::string(1, regina::Base64Encoder::spare[1]));
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true,
              "a non-base64 character in the suffix throws InvalidArgument");

    // A suffix whose length isn't a multiple of the per-index width w: for
    // `ball` (only 4 triangles) w == 1, and every length is trivially a
    // multiple of 1, so this failure mode needs a triangulation large
    // enough that w >= 2. Chain enough tetrahedra together (each glued to
    // the next along one facet) to push the triangle count past 64.
    regina::Triangulation<3> big;
    regina::Perm<4> swap01(1, 0, 2, 3);
    auto *prev = big.newTetrahedron();
    for (int i = 1; i < 25; ++i) {
        // Facet 0 of `prev` (outgoing) glues to facet 1 of `next`
        // (incoming), via a permutation swapping 0 and 1 -- keeping
        // "incoming" and "outgoing" on distinct facet numbers so each
        // interior tetrahedron's facet 0 is used exactly once, as the
        // *next* link's outgoing facet.
        auto *next = big.newTetrahedron();
        prev->join(0, next, swap01);
        prev = next;
    }
    std::string bigSig = big.isoSig();
    size_t bigM = big.countFaces<2>();
    int bigW = regina::Base64Encoder::integerWidth(bigM == 0 ? 0 : bigM - 1);
    EXPECT_EQ(bigW >= 2, true,
              "the chained triangulation has enough triangles that its "
              "per-index width is at least 2 (a precondition for the next "
              "check to be meaningful)");

    threw = false;
    try {
        // A single base64 character: length 1, not a multiple of bigW.
        fromPairSig<3, 2>(bigSig + "_a");
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true,
              "a suffix length that isn't a multiple of the per-index "
              "width throws InvalidArgument");

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

    // Build a properly-encoded "index 0" suffix via Base64Encoder directly,
    // rather than hand-guessing a literal character: this is exactly how
    // pairSig() itself would encode index 0 for foldTri's own per-index
    // width, whatever that happens to be.
    size_t foldM = foldTri.countFaces<2>();
    int foldW = regina::Base64Encoder::integerWidth(foldM == 0 ? 0 : foldM - 1);
    regina::Base64Encoder foldEnc;
    foldEnc.encodeInt(0, foldW);
    std::string zeroSuffix = std::move(foldEnc).str();

    threw = false;
    try {
        fromPairSig<3, 2>(foldSig + "_" + zeroSuffix);
    } catch (const regina::InvalidArgument &) {
        threw = true;
    }
    EXPECT_EQ(threw, true,
              "a self-folded singleton face throws InvalidArgument (not "
              "jointly addable)");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 7: PairSigContext produces byte-identical output to the free
// pairSig(), across many distinct marked sets over one ambient.
//
// This is the entire contract of the context: it is a cheaper ROUTE to the
// same string, never a different encoding. results/cobordisms.csv stores
// these strings and peripheral_slopes reconstructs surfaces from them, so
// any divergence here would silently invalidate recorded witnesses.
//
// Reusing one context across every marked set is deliberate -- it is how the
// search uses it, and it is what would expose state leaking between calls.
// ─────────────────────────────────────────────────────────────────────────────
void test_context_matches_free_function() {
    std::cout << "\n--- PairSigContext matches pairSig() ---\n";

    // 3-D: two glued tetrahedra, every marked subset of a few triangles.
    regina::Triangulation<3> tri;
    auto *t0 = tri.newTetrahedron();
    auto *t1 = tri.newTetrahedron();
    t0->join(0, t1, regina::Perm<4>());

    PairSigContext<3, 2> ctx3(tri);
    int mismatches3 = 0;
    auto nTriangles = static_cast<int>(tri.countFaces<2>());
    for (int a = 0; a < nTriangles; ++a) {
        for (int b = a + 1; b < nTriangles; ++b) {
            std::vector<int> marked = {a, b};
            if (ctx3.sig(marked) != pairSig<3, 2>(tri, marked))
                ++mismatches3;
        }
        std::vector<int> single = {a};
        if (ctx3.sig(single) != pairSig<3, 2>(tri, single))
            ++mismatches3;
    }
    EXPECT_EQ(mismatches3, 0,
              "3-D: context matches pairSig() for every 1- and 2-face "
              "marked set");

    // The empty branch short-circuits before any automorphism work, so it
    // is the one path that never touches autos_ -- check it explicitly.
    // (Hoisted into locals: the comma in `pairSig<3, 2>` would otherwise be
    // read by the preprocessor as an extra macro argument.)
    std::string ctx3Empty = ctx3.sig({});
    std::string free3Empty = pairSig<3, 2>(tri, {});
    EXPECT_EQ(ctx3Empty, free3Empty,
              "3-D: context matches pairSig() for the empty marked set");

    // 4-D: the dimension that actually matters for the drain.
    regina::Triangulation<4> tri4;
    auto *p = tri4.newPentachoron();
    auto *q = tri4.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    PairSigContext<4, 2> ctx4(tri4);
    int mismatches4 = 0;
    auto nTri4 = static_cast<int>(tri4.countFaces<2>());
    for (int a = 0; a < nTri4; ++a) {
        std::vector<int> single = {a};
        if (ctx4.sig(single) != pairSig<4, 2>(tri4, single))
            ++mismatches4;
        for (int b = a + 1; b < nTri4; ++b) {
            std::vector<int> marked = {a, b};
            if (ctx4.sig(marked) != pairSig<4, 2>(tri4, marked))
                ++mismatches4;
        }
    }
    EXPECT_EQ(mismatches4, 0,
              "4-D: context matches pairSig() for every 1- and 2-face "
              "marked set");
    std::string ctx4Empty = ctx4.sig({});
    std::string free4Empty = pairSig<4, 2>(tri4, {});
    EXPECT_EQ(ctx4Empty, free4Empty,
              "4-D: context matches pairSig() for the empty marked set");

    // The shared prefix is the thing the context exists to compute once.
    EXPECT_EQ(ctx4.ambientSig(), tri4.isoSig(),
              "context's ambient prefix equals the ambient's own isoSig()");
    EXPECT_EQ(ctx4.automorphismCount() >= 1, true,
              "at least the identity automorphism was collected");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 8: signatures a context produces still DECODE.
//
// Byte equality with the old path (test 7) proves the encoding did not
// change; it does not prove the result is still a usable signature. This is
// the property peripheral_slopes.cpp relies on when it rebuilds a surface
// from a recorded pairsig, so exercise the real decoders.
// ─────────────────────────────────────────────────────────────────────────────
void test_context_output_decodes() {
    std::cout << "\n--- PairSigContext output decodes ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    Skeleton<4, 2> skeleton(tri);
    KnottedSurface surface(skeleton);
    bool added = surface.addFace(6) && surface.addFace(7);
    EXPECT_EQ(added, true, "the known-good {6,7} pair is accepted");
    if (!added)
        return;

    PairSigContext<4, 2> ctx(tri);
    std::string sig = ctx.sig(surface.markedFaces());

    // fromKnottedSurfaceSig(): reconstruct, then re-encode. The decoded
    // marked-face indices are into fromSig()'s canonical reconstruction, not
    // into `tri`, so idempotence (not raw index equality) is the check --
    // exactly as test_knotted_surface_round_trip does for the free function.
    auto decoded = fromKnottedSurfaceSig(sig);
    PairSigContext<4, 2> decodedCtx(*decoded.ambient);
    EXPECT_EQ(decodedCtx.sig(decoded.surface->markedFaces()), sig,
              "re-encoding a context-produced signature through a fresh "
              "context reproduces it exactly");
    EXPECT_EQ(decoded.ambient->size(), tri.size(),
              "decoded ambient has the same pentachoron count");
    EXPECT_EQ(decoded.ambient->countFaces<2>(), tri.countFaces<2>(),
              "decoded ambient has the same triangle count");
    EXPECT_EQ(decoded.surface->markedFaces().size(),
              surface.markedFaces().size(),
              "decoded surface has the same number of marked faces");

    // fromPairSig(): the same, through the base-class decoder.
    auto decodedBase = fromPairSig<4, 2>(sig);
    std::string reencodedBase =
        pairSig<4, 2>(*decodedBase.skeleton, *decodedBase.submanifold);
    EXPECT_EQ(reencodedBase, sig,
              "fromPairSig() round-trips a context-produced signature");
}

// ─────────────────────────────────────────────────────────────────────────────
// Test 9: many threads sharing one LazyPairSigContext, from a standing start.
//
// Two distinct races are in scope here, and only the first is obvious:
//
//  1. The lazy build itself. Every drain thread reaches get() at once on the
//     first witness; std::call_once must serialise the build and publish a
//     fully-constructed context to all of them.
//
//  2. The SHARED canonical triangulation. This is the subtler one, and it is
//     new: the old code built a local `canon` inside every pairSig() call, so
//     each thread had its own. A context holds ONE canon_ that every thread
//     then reads concurrently. Regina computes a triangulation's skeleton
//     lazily on first access, so if canon_ reached the threads with its
//     skeleton uncomputed, the first concurrent face<subdim>() calls would
//     race. What prevents that is the constructor forcing the skeleton
//     (countFaces<subdim>(), then findAllIsomorphisms()) BEFORE call_once
//     publishes the object -- load-bearing, and the reason this test exists.
//
// Correct output under contention is the assertion; run under ThreadSanitizer
// to check the memory model itself.
// ─────────────────────────────────────────────────────────────────────────────
void test_context_concurrent_first_use() {
    std::cout << "\n--- LazyPairSigContext under concurrent first use ---\n";

    regina::Triangulation<4> tri;
    auto *p = tri.newPentachoron();
    auto *q = tri.newPentachoron();
    p->join(4, q, regina::Perm<5>());
    q->join(0, q, regina::Perm<5>(1, 0, 2, 3, 4));

    // Reference values, computed single-threaded through the OLD path.
    auto nTri = static_cast<int>(tri.countFaces<2>());
    std::vector<std::vector<int>> markedSets;
    for (int a = 0; a < nTri; ++a) {
        markedSets.push_back({a});
        for (int b = a + 1; b < nTri; ++b)
            markedSets.push_back({a, b});
    }
    std::vector<std::string> reference;
    reference.reserve(markedSets.size());
    for (const auto &m : markedSets)
        reference.push_back(pairSig<4, 2>(tri, m));

    // One shared context, N threads, all starting together.
    LazyPairSigContext<4, 2> lazy(tri);
    constexpr int THREADS = 8; // oversubscribed on purpose
    std::atomic<bool> go{false};
    std::vector<std::thread> threads;
    std::vector<int> mismatches(THREADS, 0);

    for (int t = 0; t < THREADS; ++t) {
        threads.emplace_back([&, t] {
            while (!go.load(std::memory_order_acquire))
                ; // spin, so every thread hits get() as close to together
                  // as possible -- a staggered start would let one thread
                  // finish the build before the others ever contend.
            for (size_t i = 0; i < markedSets.size(); ++i) {
                if (lazy.get().sig(markedSets[i]) != reference[i])
                    ++mismatches[t];
            }
        });
    }
    go.store(true, std::memory_order_release);
    for (auto &th : threads)
        th.join();

    int total = 0;
    for (int m : mismatches)
        total += m;
    EXPECT_EQ(total, 0,
              "8 threads sharing one lazy context all reproduce the "
              "single-threaded signatures exactly");
    // Safe to read only now that every thread has joined.
    EXPECT_EQ(lazy.built(), true, "the shared context was actually built");
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
    run("context_matches_free_function", test_context_matches_free_function);
    run("context_output_decodes", test_context_output_decodes);
    run("context_concurrent_first_use", test_context_concurrent_first_use);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
