//
//  pairsig.h
//
//  Created by John Teague on 07/26/2026.
//

#ifndef PAIRSIG_H

#define PAIRSIG_H

#include <memory>
#include <string>
#include <vector>

#include "embeddedsubmanifold.h"
#include "skeleton.h"

/*! \file utils/surfer/pairsig.h
 *  \brief Isomorphism signatures of (ambient triangulation, embedded
 *  subcomplex) pairs.
 *
 *  A pair signature generalizes Regina's isomorphism signature (see
 *  regina::Triangulation::isoSig()) to also encode a marked subcomplex of
 *  the ambient triangulation's subdim-skeleton -- the same subcomplex
 *  EmbeddedSubmanifold<dim,subdim>/KnottedSurface track. Like a plain
 *  isoSig, it is a true isomorphism invariant: two isomorphic (ambient,
 *  marked subcomplex) pairs produce byte-identical signatures, even when
 *  the ambient triangulation has automorphisms that move the marked
 *  subcomplex around. See pairSig()'s implementation for the
 *  canonicalization argument (it minimizes over every automorphism of the
 *  isoSig-canonical reconstruction, not just one arbitrary relabelling).
 *
 *  All of this is built entirely from Regina's existing public API
 *  (isoSigDetail(), fromSig(), findAllIsomorphisms(), Isomorphism, and
 *  FaceEmbedding/FaceNumbering) -- no engine changes are required.
 */

/**
 * Computes a pair signature for (`ambient`, `markedFaces`).
 *
 * `markedFaces` is a list of indices into `ambient`'s subdim-faces (the
 * same representation EmbeddedSubmanifold::markedFaces() returns/consumes).
 *
 * \pre `ambient` is non-empty and connected (the same precondition
 * isoSigDetail() itself carries).
 * \pre every entry of `markedFaces` is a valid index into `ambient`'s
 * subdim-faces.
 *
 * \exception regina::FailedPrecondition `ambient` is empty or
 * disconnected.
 *
 * \tparam dim the dimension of the ambient triangulation.
 * \tparam subdim the dimension of the marked subcomplex's faces.
 */
template <int dim, int subdim>
std::string pairSig(const regina::Triangulation<dim> &ambient,
                     const std::vector<int> &markedFaces);

/**
 * As above, computing a pair signature directly from an already-built
 * EmbeddedSubmanifold (or KnottedSurface, via the base-class reference).
 */
template <int dim, int subdim>
std::string pairSig(const Skeleton<dim, subdim> &skeleton,
                     const EmbeddedSubmanifold<dim, subdim> &s) {
    return pairSig<dim, subdim>(skeleton.triangulation(), s.markedFaces());
}

/**
 * Everything fromPairSig() reconstructs from a pair signature.
 *
 * \warning `skeleton` holds a non-owning pointer into `*ambient` (see
 * Skeleton), and `submanifold` holds a non-owning reference into
 * `*skeleton` (see EmbeddedSubmanifold). This struct owns all three
 * indirectly via unique_ptr specifically so that moving/returning a
 * DecodedPairSig by value never relocates the pointed-to objects, only the
 * pointers to them.
 */
template <int dim, int subdim>
struct DecodedPairSig {
    std::unique_ptr<regina::Triangulation<dim>> ambient;
    std::unique_ptr<Skeleton<dim, subdim>> skeleton;
    std::unique_ptr<EmbeddedSubmanifold<dim, subdim>> submanifold;
};

/**
 * Reconstructs the (ambient triangulation, embedded subcomplex) pair
 * encoded by `sig`.
 *
 * \exception regina::InvalidArgument `sig` is malformed (missing
 * delimiter, an unparseable index list, or an index list that isn't
 * jointly addable -- see EmbeddedSubmanifold's seeded constructor).
 */
template <int dim, int subdim>
DecodedPairSig<dim, subdim> fromPairSig(const std::string &sig);

/** As DecodedPairSig<4,2>, but reconstructing a KnottedSurface instead of the base class. */
struct DecodedKnottedSurfaceSig {
    std::unique_ptr<regina::Triangulation<4>> ambient;
    std::unique_ptr<Skeleton<4, 2>> skeleton;
    std::unique_ptr<KnottedSurface> surface;
};

/** As fromPairSig<4,2>(), but reconstructing a KnottedSurface instead of the base class. */
DecodedKnottedSurfaceSig fromKnottedSurfaceSig(const std::string &sig);

extern template std::string pairSig<3, 2>(
    const regina::Triangulation<3> &, const std::vector<int> &);
extern template std::string pairSig<4, 2>(
    const regina::Triangulation<4> &, const std::vector<int> &);

extern template DecodedPairSig<3, 2> fromPairSig<3, 2>(const std::string &);
extern template DecodedPairSig<4, 2> fromPairSig<4, 2>(const std::string &);

#endif // PAIRSIG_H
