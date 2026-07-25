//
//  embeddedsubmanifold.h
//
//  Created by John Teague on 07/21/2026.
//

#ifndef EMBEDDEDSUBMANIFOLD_H

#define EMBEDDEDSUBMANIFOLD_H

#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <utilities/typeutils.h>

#include "linkcomplement.h"
#include "rollbackunionfind.h"
#include "skeleton.h"

/*! \file utils/surfer/embeddedsubmanifold.h
 *  \brief Incrementally tracks an embedded subcomplex of a triangulation.
 */

/**
 * A restriction on which boundary components a tracked subcomplex must
 * have; see EmbeddedSubmanifold::satisfies().
 */
enum class BoundaryCondition : uint8_t {
    all,       /**< No restriction. */
    closed,    /**< The subcomplex has no boundary. */
    proper,    /**< Every boundary facet of the subcomplex lies on the ambient boundary. */
    connected, /**< Each boundary component of the subcomplex maps to a distinct ambient boundary component. */
};

/**
 * Incrementally tracks a subcomplex of an ambient triangulation as
 * subdim-faces are added and removed one at a time, maintaining whether
 * the result is embedded, closed, and/or proper.
 *
 * \warning As of the codimension >= 2 embeddedness check being disabled
 * (see addFace()), this only guarantees genuine embeddedness at the facet
 * (codimension-1) level. Faces may be accepted whose addition produces a
 * subcomplex that is not injective at codimension >= 2 -- use
 * isEmbedded() to check the current state's actual codimension >= 2
 * injectivity.
 *
 * \tparam dim the dimension of the ambient triangulation.
 * \tparam subdim the dimension of the tracked subcomplex.
 */
template <int dim, int subdim> class EmbeddedSubmanifold {
protected:
  using Face = const typename Skeleton<dim, subdim>::Face;
      /**< A face of the ambient skeleton being tracked. */
  template <int k> using LowDimVec = std::vector<int>;

  const Skeleton<dim, subdim> &skeleton_; /**< The ambient skeleton. */
  regina::Triangulation<subdim> subtri_;
      /**< The triangulation of the tracked subcomplex. */
  std::vector<regina::Simplex<subdim> *> faces_;
      /**< faces_[f]: the subtri_ simplex realizing ambient face f, or
           null if f is not currently in the subcomplex. */

private:
  // For k in [0, subdim-1]: number of submanifold simplices containing each
  // k-face of tri_. For k == subdim-1 (facets): 0 = absent, 1 = boundary,
  // ≥2 = interior. For k < subdim-1: positive iff the face is in the
  // submanifold.
  regina::TupleOverRange<0, subdim, LowDimVec> faceCount_;

  /**
   * Tracking for isEmbedded(): whether the realization map subtri_ ->
   * the ambient triangulation is injective at codimension >= 2 (k in
   * [0, subdim-2]); codimension-1 injectivity is already guaranteed
   * structurally by addFace()'s facet-level check.
   *
   * A "slot" is (ambient face index f, dimension k, local k-face index i),
   * identified by f * FaceNumbering<subdim,k>::nFaces + i. dsu_[k] unions
   * slots that addFace()'s realized gluings have identified in subtri_.
   * classRoots_[k][v] holds the currently-distinct DSU root ids, among
   * slots mapping to ambient k-face v, that are presently registered --
   * |classRoots_[k][v]| >= 2 means v is a singularity. singularCount_ is
   * the number of (k, v) pairs currently singular; isEmbedded_ :=
   * (singularCount_ == 0).
   *
   * Deliberately not derived from subtri_'s own Regina-computed skeleton:
   * every join()/newSimplex()/removeSimplex() call invalidates it
   * entirely, so re-deriving it would cost a full rebuild on every single
   * addFace()/removeFace() call.
   */
  std::vector<std::vector<std::vector<int>>> classRoots_;
  std::vector<RollbackUnionFind> dsu_;

  /**
   * Tracking for isProper(): facetIsAmbientBoundary_ records, once,
   * whether each ambient facet lies on the ambient triangulation's
   * boundary. badProperCount_ counts ambient facets currently at
   * boundary-facet-count 1 (see faceCount_) that are not on the ambient
   * boundary; isProper() := (badProperCount_ == 0). Updated alongside
   * faceCount_ in addFace()/removeFace(), so it needs no rollback log of
   * its own.
   */
  std::vector<bool> facetIsAmbientBoundary_;
  int badProperCount_ = 0;

  struct RegistryUndoEntry {
    size_t v;
    int root;
    bool wasInsert;
        /**< true: undo by erasing `root` from classRoots_[k][v]; false:
             undo by re-inserting it. */
    int singularDelta;
        /**< The exact singularCount_ change this event caused; undone
             via singularCount_ -= singularDelta. */
  };
  std::vector<std::vector<RegistryUndoEntry>> registryUndoLog_; /**< Indexed by k. */

  struct Checkpoint {
    std::vector<size_t> dsuMark; /**< dsu_[k].checkpoint(), per k. */
    std::vector<size_t> registryMark; /**< registryUndoLog_[k].size(), per k. */
  };
  std::vector<Checkpoint> checkpoints_; /**< Indexed by ambient face index f. */

  int singularCount_ = 0;
  bool isEmbedded_ = true; /**< Kept in lockstep with (singularCount_ == 0). */

  /**
   * Unions the DSU slot for (implicit ambient face, k, `slotA`) with the
   * slot for (implicit destination face, k, `slotB`), both mapping to
   * ambient k-face `v`. If this merges two previously-independent
   * classes at `v`, updates classRoots_[k][v]/singularCount_ and logs the
   * change for removeFace() to undo.
   */
  void unite_(int k, size_t v, int slotA, int slotB);

  /**
   * Registers slot `r`'s current DSU root against ambient k-face `v`'s
   * set of known classes, if not already present. If this is the second
   * distinct class registered at `v`, marks `v` as a new singularity.
   */
  void registerRoot_(int k, size_t v, int r);

protected:
  /**
   * Read-only access to isEmbedded()-tracking state, for subclasses
   * layering additional codimension-2 checks (e.g. KnottedSurface's
   * transverse-self-intersection/local-flatness checks) on top of this
   * combinatorial machinery.
   *
   * \pre Ambient face `f` is currently present in the tracked subcomplex.
   */
  int vertexClassRoot(int f, int localVertex) const {
    return dsu_[0].find(f * regina::FaceNumbering<subdim, 0>::nFaces +
                        localVertex);
  }

  /** Returns the currently-registered DSU class roots at ambient vertex `v`. */
  const std::vector<int> &registeredClassRoots(size_t v) const {
    return classRoots_[0][v];
  }

  /** Returns the current submanifold-membership count of ambient facet `ambientFacetIdx`. */
  int facetCount(size_t ambientFacetIdx) const {
    return std::get<subdim - 1>(faceCount_)[ambientFacetIdx];
  }

public:
  /** Creates an empty tracked subcomplex over `skeleton`. */
  EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton);

  /**
   * As above, then adds every face in `seedFaces` via addFace(), in the
   * given order.
   *
   * \pre `seedFaces` is a connected, jointly addable set of faces.
   *
   * \throws regina::InvalidArgument if any addFace() call fails, rather
   * than silently leaving a partial result.
   */
  EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton,
                      const std::vector<int> &seedFaces);

  /**
   * Attempts to add face `f` to the tracked subcomplex.
   *
   * \return \c false (with no change made) if adding `f` would violate
   * the facet-level (codimension-1) embedding condition; otherwise
   * commits the addition and returns \c true.
   *
   * \warning See the class documentation: this only enforces
   * codimension-1 embeddedness. Use isEmbedded() to check the current
   * state's codimension >= 2 injectivity.
   */
  bool addFace(int f);

  /**
   * Attempts to add every face in `faces`, in order, via addFace().
   *
   * Transactional like addFace(): on success every face is committed and
   * this returns \c true; on failure, everything added during this call
   * is rolled back and this returns \c false.
   */
  bool addFaces(const std::vector<int> &faces);

  /**
   * Removes face `f` from the tracked subcomplex.
   *
   * \pre `f` is the most recently added face not yet removed (removals
   * must undo additions in LIFO order).
   */
  void removeFace(int f);

  /** Returns the triangulation of the tracked subcomplex. */
  const regina::Triangulation<subdim> &triangulation() const { return subtri_; }

  /**
   * Returns whether the tracked subcomplex has no boundary.
   *
   * \note Implemented directly via hasBoundaryFacets() rather than
   * Triangulation::isClosed(), which is only defined for dimensions 2, 3,
   * and 4.
   */
  bool isClosed() const { return !subtri_.hasBoundaryFacets(); }

  /**
   * Returns whether every boundary facet of the tracked subcomplex maps
   * to a facet that is itself on the ambient triangulation's boundary.
   */
  bool isProper() const;

  /**
   * Returns whether each boundary component of the tracked subcomplex
   * maps into a distinct boundary component of the ambient
   * triangulation.
   */
  bool boundaryComponentsMapInjectively() const;

  /**
   * Returns whether the realization map is injective at every dimension
   * for every face currently present, not just at the facet level.
   *
   * O(1): incrementally maintained by addFace()/removeFace().
   */
  bool isEmbedded() const { return isEmbedded_; }

  /** Returns whether the tracked subcomplex satisfies `cond`. */
  bool satisfies(BoundaryCondition cond) const;

  /**
   * Returns whether `gluings` contains a self-gluing that could never be
   * realized: a local facet glued to two or more distinct partner facets
   * of the same simplex.
   */
  static bool hasIrreparableSelfGluing(
      const std::vector<typename Skeleton<dim, subdim>::Gluing> &gluings);
};

extern template class EmbeddedSubmanifold<3, 2>;
extern template class EmbeddedSubmanifold<4, 2>;

/**
 * A surface embedded in a 4-manifold, in the sense of classical knotted
 * surface theory: EmbeddedSubmanifold<4, 2>, plus boundaryLinks() (which
 * identifies the knots/links bounded by the surface's boundary) and
 * surfaceTypeKey()/formatSurfaceType() (which classify the surface's own
 * topology -- genus, orientability, punctures).
 */
class KnottedSurface : public EmbeddedSubmanifold<4, 2> {
public:
  using SurfaceTypeKey = std::tuple<bool, int, int>;

private:
  std::vector<regina::Triangulation<3>> bdryComponents_;
      /**< The ambient triangulation's boundary components. */

  /**
   * (ambient face index, local vertex index) pairs of the surface
   * currently touching ambient vertex `v` (indexed by `v`), in insertion
   * order.
   *
   * Stored as a pair rather than just the face index: a face can touch
   * `v` via more than one of its own local vertices (a legitimate case
   * arising from the face's own self-gluings, which are committed as
   * self-joins when the face is added), and each such occurrence has its
   * own distinct opposite edge in the vertex link, so the local index
   * can't be re-derived from the face alone.
   *
   * Maintained by this class's own addFace()/removeFace() overrides via
   * push_back()/pop_back() -- safe with no extra bookkeeping, since
   * removeFace() only ever undoes the most-recently-added face.
   */
  std::vector<std::vector<std::pair<int, int>>> facesAtVertex_;

  /**
   * Per-ambient-vertex, lazily built: (pentachoron index * 5 + local
   * vertex index) -> index into that vertex's own
   * embeddings()/buildLink() tetrahedron numbering. A vertex's link never
   * changes shape across a search, so this is safe to memoize for the
   * instance's lifetime.
   */
  mutable std::vector<std::optional<std::unordered_map<uint64_t, size_t>>>
      vertexEmbedIndexCache_;

  /**
   * Translates ambient triangle `f`'s local vertex `localVertex` (which
   * must equal `ambientVertex`) into the corresponding edge of
   * `ambientVertex`'s link -- the edge traced by `f` in the vertex link,
   * i.e. `f`'s facet opposite `localVertex`.
   */
  const regina::Edge<3> *linkEdgeForTriangle_(
      const regina::Vertex<4> *ambientVertex, int f, int localVertex) const;

  /**
   * Returns whether the DSU class `root` at ambient vertex `v` is
   * currently closed: every spoke edge among its member triangles is at
   * facetCount() == 2.
   */
  bool isPetalClosed_(size_t v, int root) const;

  /** Returns the (unordered) edge list of the closed petal `root` at `v`. */
  std::vector<const regina::Edge<3> *>
  closedPetalCurve_(const regina::Vertex<4> *ambientVertex, size_t v,
                    int root) const;

public:
  /** Creates an empty tracked surface over `skeleton`. */
  KnottedSurface(const Skeleton<4, 2> &skeleton);

  /** As above, then adds every face in `seedFaces`; see the base class's seeded constructor. */
  KnottedSurface(const Skeleton<4, 2> &skeleton,
                 const std::vector<int> &seedFaces);

  /**
   * As EmbeddedSubmanifold::addFace(), plus additional checks: if adding
   * `f` closes one of its touched vertices' petals, rejects the addition
   * (rolling back via removeFace(f)) when the closed curve is knotted
   * (non-local-flatness) or nonzero-linked with another already-closed
   * petal at the same vertex (transverse self-intersection). Both
   * conditions are hereditary under removeFace(), so this is safe to
   * enforce as a hard, permanent rejection during a search.
   */
  bool addFace(int f);

  bool addFaces(const std::vector<int> &faces);

  void removeFace(int f);

  /** Returns the knot/link bounded by each boundary component of the surface, keyed by boundary component index. */
  std::vector<std::pair<size_t, Link>> boundaryLinks() const;

  /** Classifies `surface`'s topology as (orientable, genus, number of punctures). */
  static SurfaceTypeKey surfaceTypeKey(const regina::Triangulation<2> &surface);

  /** Formats a SurfaceTypeKey as a human-readable description. */
  static std::string formatSurfaceType(const SurfaceTypeKey &key);

  /** Returns this surface's own SurfaceTypeKey. */
  SurfaceTypeKey surfaceType() const { return surfaceTypeKey(triangulation()); }
};

#endif // EMBEDDEDSUBMANIFOLD_H
