//
//  pairsig.h
//
//  Created by John Teague on 07/26/2026.
//

#ifndef PAIRSIG_H

#define PAIRSIG_H

#include <memory>
#include <mutex>
#include <string>
#include <vector>

#include <triangulation/isomorphism.h>

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
 *  (isoSigDetail(), fromSig(), findAllIsomorphisms(), Isomorphism,
 *  FaceEmbedding/FaceNumbering, and Base64Encoder/Base64Decoder for the
 *  marked-face suffix -- see pairSig()'s implementation) -- no engine
 *  changes are required.
 */

/**
 * A subdim-face's image under a chain of isomorphisms, expressed as a
 * (destination simplex, vertex permutation) pair rather than a resolved
 * face index. This lets a face be pushed through several isomorphisms in
 * sequence (e.g. ambient -> canon, then canon -> canon for each of canon's
 * automorphisms) without resolving to a concrete face object until the
 * final step -- what pairSig()'s per-automorphism inner loop relies on to
 * stay cheap, and equally useful for canonicalizing a marked face set
 * against a single fixed triangulation's own automorphism group (see
 * identify::BoundarySignatureCache in identifycomplement.h).
 */
template <int dim>
struct FaceDescriptor {
    size_t simplex;
    regina::Perm<dim + 1> vertexPerm;
};

/**
 * The descriptor for `ambient`'s own subdim-face `f` (index into
 * `ambient`'s subdim-faces), before any isomorphism has been applied.
 *
 * Defined here (not in pairsig.cpp) so that any translation unit --
 * including identifycomplement.cpp's BoundarySignatureCache, which has no
 * other reason to link against pairsig.cpp's own explicit instantiations
 * (pairSig()/fromPairSig() and everything they in turn pull in, e.g.
 * EmbeddedSubmanifold/KnottedSurface/Skeleton's .cpp files) -- can
 * instantiate it locally with no added link-time dependency.
 */
template <int dim, int subdim>
FaceDescriptor<dim> faceDescriptor(const regina::Triangulation<dim> &ambient,
                                    int f) {
    const auto &emb = ambient.template face<subdim>(f)->front();
    return {static_cast<size_t>(emb.simplex()->index()), emb.vertices()};
}

/** The descriptor of `desc`'s image under `iso`. See faceDescriptor() for why this is header-defined. */
template <int dim>
FaceDescriptor<dim> applyIsomorphism(const regina::Isomorphism<dim> &iso,
                                      const FaceDescriptor<dim> &desc) {
    return {static_cast<size_t>(iso.simpImage(desc.simplex)),
            iso.facetPerm(desc.simplex) * desc.vertexPerm};
}

/**
 * Resolves `desc` -- already expressed in `codomain`'s own numbering, e.g.
 * as returned by applyIsomorphism() -- to a concrete subdim-face index of
 * `codomain`. See faceDescriptor() for why this is header-defined.
 */
template <int dim, int subdim>
size_t resolveFaceIndex(const regina::Triangulation<dim> &codomain,
                         const FaceDescriptor<dim> &desc) {
    int localFace =
        regina::FaceNumbering<dim, subdim>::faceNumber(desc.vertexPerm);
    return codomain.simplex(desc.simplex)
        ->template face<subdim>(localFace)
        ->index();
}

/**
 * Everything in a pair signature that depends on the AMBIENT triangulation
 * alone, computed once so that many surfaces over the same ambient can be
 * signed without repeating it.
 *
 * This exists because the ambient part is overwhelmingly the expensive part,
 * and in a search it is the same every time. Measured on a 1,728-pentachoron
 * cobordism (a 9-crossing link, 2 layers):
 *
 *     pairSig<4,2>()             32,868 ms
 *       isoSigDetail(ambient)    32,809 ms   (99.8%)  <- ambient only
 *       findAllIsomorphisms         265 ms            <- ambient only
 *       fromSig(sig)                  0.5 ms          <- ambient only
 *
 * `verifyslicegenus` builds one cobordism per row and then signs every
 * witness found in it, so recomputing the above per witness made the drain's
 * cost proportional to the number of witnesses rather than the amount of
 * work: `perf` put 79% of the whole run in IsoSigData<1,4>::fillFrom plus
 * IsoSigPrintable::encode<4>, against 0.46% in the surface-dependent part of
 * pairSig() itself.
 *
 * sig() is byte-for-byte identical to pairSig<dim,subdim>(ambient, ...) --
 * it is the same code reading precomputed members. That equality is the
 * whole contract: `results/cobordisms.csv` stores these strings and
 * peripheral_slopes reconstructs surfaces from them, so a context must never
 * be a different encoding, only a cheaper route to the same one.
 *
 * Not copyable: `autos_` are isomorphisms *of* `canon_`, so copying the
 * members independently would leave them describing a different object.
 *
 * Thread-safe for concurrent sig() calls once constructed: every member is
 * read-only thereafter.
 */
template <int dim, int subdim>
class PairSigContext {
  public:
    /**
     * \pre `ambient` is non-empty and connected (as isoSigDetail() requires).
     * \exception regina::FailedPrecondition `ambient` is empty or
     * disconnected.
     */
    explicit PairSigContext(const regina::Triangulation<dim> &ambient)
        : PairSigContext(ambient, detailFor(ambient)) {}

    PairSigContext(const PairSigContext &) = delete;
    PairSigContext &operator=(const PairSigContext &) = delete;

    /** Identical to pairSig<dim,subdim>(ambient, markedFaces). */
    std::string sig(const std::vector<int> &markedFaces) const;

    /** The ambient's own isoSig, i.e. every signature's shared prefix. */
    const std::string &ambientSig() const { return sig_; }

    /** How many automorphisms of the canonical ambient sig() minimises over. */
    size_t automorphismCount() const { return autos_.size(); }

  private:
    using Detail = std::pair<std::string, regina::Isomorphism<dim>>;

    static Detail detailFor(const regina::Triangulation<dim> &ambient);

    PairSigContext(const regina::Triangulation<dim> &ambient, Detail detail);

    const regina::Triangulation<dim> *ambient_;
    std::string sig_;
    regina::Isomorphism<dim> psi0_;
    regina::Triangulation<dim> canon_;
    std::vector<regina::Isomorphism<dim>> autos_;
    size_t numFaces_;
    int width_;
};

/**
 * A PairSigContext built on first use and shared thereafter.
 *
 * Laziness is the point, not an implementation detail. Building a context
 * costs an isoSigDetail() of the whole ambient -- ~33 s on a
 * 1,728-pentachoron cobordism -- and MOST SEARCH ROWS NEVER SIGN ANYTHING,
 * because they find no witness at all. Constructing eagerly (say, alongside
 * the per-row Skeleton) would hand every barren row a large bill it does not
 * currently pay, turning a win on productive rows into a loss overall. So the
 * cost is deferred to the first sig() that actually happens.
 *
 * get() is safe to call concurrently: std::call_once both serialises the
 * build and establishes happens-before for every other caller, so the
 * returned context is fully constructed before any thread observes it. This
 * is the same pattern BoundarySignatureCache::ensureAutomorphismGroup_() and
 * SurfaceSearch::ensureBoundarySigCaches_() already use.
 *
 * Holds a non-owning pointer: `ambient` must outlive this object.
 */
template <int dim, int subdim>
class LazyPairSigContext {
  public:
    explicit LazyPairSigContext(const regina::Triangulation<dim> &ambient)
        : ambient_(&ambient) {}

    LazyPairSigContext(const LazyPairSigContext &) = delete;
    LazyPairSigContext &operator=(const LazyPairSigContext &) = delete;

    const PairSigContext<dim, subdim> &get() const {
        std::call_once(once_, [this] {
            ctx_ = std::make_unique<PairSigContext<dim, subdim>>(*ambient_);
        });
        return *ctx_;
    }

    /** Whether the context has actually been built yet (for tests/reporting). */
    bool built() const { return static_cast<bool>(ctx_); }

  private:
    const regina::Triangulation<dim> *ambient_;
    mutable std::once_flag once_;
    mutable std::unique_ptr<PairSigContext<dim, subdim>> ctx_;
};

/**
 * Computes a pair signature for (`ambient`, `markedFaces`).
 *
 * `markedFaces` is a list of indices into `ambient`'s subdim-faces (the
 * same representation EmbeddedSubmanifold::markedFaces() returns/consumes).
 *
 * Builds a throwaway PairSigContext, so it pays the full ambient cost on
 * every call. Signing many surfaces over one ambient (which is what a search
 * does) should build one PairSigContext and reuse it.
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
 * delimiter, a marked-face suffix whose length isn't a multiple of the
 * per-index base64 width, a suffix containing a character outside
 * Regina's base64 alphabet, or an index list that isn't jointly addable
 * -- see EmbeddedSubmanifold's seeded constructor).
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

extern template class PairSigContext<3, 2>;
extern template class PairSigContext<4, 2>;

extern template std::string pairSig<3, 2>(
    const regina::Triangulation<3> &, const std::vector<int> &);
extern template std::string pairSig<4, 2>(
    const regina::Triangulation<4> &, const std::vector<int> &);

extern template DecodedPairSig<3, 2> fromPairSig<3, 2>(const std::string &);
extern template DecodedPairSig<4, 2> fromPairSig<4, 2>(const std::string &);

#endif // PAIRSIG_H
