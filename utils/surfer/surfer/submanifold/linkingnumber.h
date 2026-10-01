//
//  linkingnumber.h
//

#ifndef LINKINGNUMBER_H

#define LINKINGNUMBER_H

#include <array>
#include <atomic>
#include <cstdint>
#include <optional>
#include <vector>

#include <triangulation/dim3.h>

/*! \file utils/surfer/linkingnumber.h
 *  \brief The linking number of two disjoint edge cycles in a triangulated
 *  3-sphere, by a sparse cochain computation.
 *
 * KnottedSurface::addFace() needs, for two closed petals at an interior
 * vertex v, whether their traces A and B link in Lk(v) ~ S^3. The old route,
 * EdgeComplement::linkingNumberWith(), drills A out of a copy of Lk(v) and
 * reads B's class off a MarkedAbelianGroup: dense bignum Smith normal form,
 * ~1-6 s a call. This computes the same number from the cochain complex of
 * Lk(v) itself, without drilling or subdividing.
 *
 * THE METHOD. Let T be a closed, oriented triangulation of S^3 and A, B
 * vertex-disjoint cycles in its 1-skeleton.
 *
 * 1. Push B off into the dual cell structure. Orient B as b_0 -e_0-> b_1
 *    -e_1-> ... For each edge e_k pick a tetrahedron t_k containing it, and
 *    for each vertex b_k a path of tetrahedra from t_{k-1} to t_k, each
 *    containing b_k, consecutive ones sharing a triangle that contains b_k.
 *    The dual edges of those triangles form a dual 1-cycle B*.
 *
 *    B* is homologous to B in S^3 - A. Near b_k, B (from the midpoint of
 *    e_{k-1} through b_k to the midpoint of e_k) and the corresponding piece
 *    of B* (through barycentres of the tetrahedra and triangles on the path)
 *    both lie in the open star of b_k. When every tetrahedron containing b_k
 *    meets it at a single corner, that open star is an open cone on the link
 *    of b_k, hence contractible; and it misses A, since A is a union of
 *    closed edges none of which has b_k as a vertex. So the two pieces are
 *    homotopic rel endpoints in S^3 - A. The pieces joining consecutive
 *    vertices are traversed twice in opposite directions and cancel.
 *
 * 2. Take beta = PD(B*), the 2-cochain counting B*'s signed crossings of
 *    each triangle. Since H^2(S^3) = 0 there is a 1-cochain x with
 *    delta x = beta; dually, x is a 2-chain c* of dual cells with boundary
 *    B*. Then lk(A, B) = lk(A, B*) = c* . A = +-x(A): a dual 2-cell meets
 *    exactly one primal edge, once, transversally. Any solution x works --
 *    two differ by a coboundary delta f (H^1(S^3) = 0), and
 *    (delta f)(A) = f(dA) = 0.
 *
 * Solved over GF(p), p = 2^61 - 1, the value mod p is exact because
 * |lk(A, B)| is at most B's edge count, far below p/2.
 *
 * SELF-CHECKS. Every answer is checked, not trusted: beta must be a cocycle
 * (delta beta = 0 on every tetrahedron) and x must satisfy delta x = beta on
 * every triangle. The precondition of step 1 is checked too. A failed check
 * returns nullopt, never a number, and the caller falls back to the old
 * route.
 */

namespace linkingnumber {

/**
 * The oriented cell structure of a closed, orientable Triangulation<3>,
 * flattened: each simplex by index, with the signs of the simplicial
 * coboundary maps relative to Regina's own vertex orderings. Built once per
 * triangulation (see PetalCache::linkComplex()) and read-only afterwards,
 * so one instance may be shared across threads.
 */
class Complex {
  public:
    /** Flattens `tri`, which must outlive this object. */
    explicit Complex(const regina::Triangulation<3> &tri);

    /**
     * Whether `tri` was closed and its tetrahedra consistently oriented.
     * linkingNumber() refuses (returns nullopt) otherwise.
     */
    bool valid() const { return valid_; }

    /** The triangulation this was built from. */
    const regina::Triangulation<3> &triangulation() const { return *tri_; }

  private:
    friend std::optional<long>
    linkingNumber(const Complex &, const std::vector<const regina::Edge<3> *> &,
                  const std::vector<const regina::Edge<3> *> &);

    const regina::Triangulation<3> *tri_;
    bool valid_ = true;
    int nV_ = 0, nE_ = 0, nF_ = 0, nT_ = 0;
    std::vector<std::array<int, 2>> edgeEnds_; /**< Edge -> its vertex(0), vertex(1). */
    std::vector<std::array<int, 3>> faceEdge_; /**< Triangle -> edge opposite each vertex. */
    std::vector<std::array<int8_t, 3>> faceEdgeSign_;
        /**< (delta x)(f) = sum_k faceEdgeSign_[f][k] * x(faceEdge_[f][k]). */
    std::vector<std::array<int, 4>> tetFace_; /**< Tetrahedron -> triangle opposite each corner. */
    std::vector<std::array<int8_t, 4>> tetSign_;
        /**< The coefficient of each face in the boundary of the oriented
             fundamental cycle: s(t) (-1)^i o_i(t). A closed, oriented T has
             the two sides of every triangle opposite. */
    std::vector<std::array<int, 4>> tetAdj_; /**< Adjacent tetrahedron through each face. */
    std::vector<std::array<regina::Perm<4>, 4>> tetGluing_; /**< The gluing through each face. */
    std::vector<std::array<int, 4>> tetVertex_; /**< Tetrahedron -> vertex at each corner. */
};

/**
 * |lk(A, B)| for vertex-disjoint cycles A and B in the 1-skeleton of the
 * 3-sphere `complex` was built from; see this file's documentation.
 *
 * \return the magnitude, or nullopt when the method cannot vouch for an
 * answer: `complex` is not valid(), A or B is not a single closed curve,
 * some tetrahedron meets a vertex of B at more than one corner, or a
 * self-check failed. Never a wrong number.
 */
std::optional<long>
linkingNumber(const Complex &complex,
              const std::vector<const regina::Edge<3> *> &A,
              const std::vector<const regina::Edge<3> *> &B);

/**
 * When set, KnottedSurface::addFace() also computes every linking number it
 * needs by the old route (EdgeComplement::linkingNumberWith()) and counts
 * any disagreement (see PetalCache::Stats). Validation only: slow. Set once,
 * before any search thread starts.
 */
extern std::atomic<bool> auditLinkingNumbers;

} // namespace linkingnumber

#endif // LINKINGNUMBER_H
