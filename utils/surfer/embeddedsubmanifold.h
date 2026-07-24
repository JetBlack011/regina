#ifndef EMBEDDEDSUBMANIFOLD_H

#define EMBEDDEDSUBMANIFOLD_H

#include <string>
#include <tuple>
#include <vector>

#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <utilities/typeutils.h>

#include "linkcomplement.h"
#include "rollbackunionfind.h"
#include "skeleton.h"

enum class BoundaryCondition : uint8_t { all, closed, proper, connected };

template <int dim, int subdim> class EmbeddedSubmanifold {
protected:
  using Face = const typename Skeleton<dim, subdim>::Face;
  template <int k> using LowDimVec = std::vector<int>;

  const Skeleton<dim, subdim> &skeleton_;
  regina::Triangulation<subdim> subtri_;
  std::vector<regina::Simplex<subdim> *> faces_;

private:
  // For k in [0, subdim-1]: number of submanifold simplices containing each
  // k-face of tri_. For k == subdim-1 (facets): 0 = absent, 1 = boundary,
  // ≥2 = interior. For k < subdim-1: positive iff the face is in the
  // submanifold.
  regina::TupleOverRange<0, subdim, LowDimVec> faceCount_;

  // Tracking for isEmbedded(): whether the realization map subtri_ -> tri_
  // is injective at codimension >= 2 (k in [0, subdim-2]). Facet-level
  // (codimension-1) injectivity is already guaranteed structurally by
  // addFace()'s Condition 1 above and needs no separate tracking.
  //
  // A "slot" is (ambient face index f, dimension k, local k-face index i),
  // identified by the integer f * FaceNumbering<subdim,k>::nFaces + i.
  // dsu_[k] unions slots that addFace()'s realized gluings have actually
  // identified in subtri_. classRoots_[k][v] holds the currently-distinct
  // DSU root ids, among slots mapping to ambient k-face v, that are
  // presently registered -- |classRoots_[k][v]| >= 2 means v is a
  // singularity (two or more genuinely distinct parts of the submanifold
  // map to the same ambient k-face). singularCount_ is the number of (k, v)
  // pairs currently in that state; isEmbedded_ := (singularCount_ == 0).
  //
  // Deliberately NOT derived from subtri_'s own Regina-computed skeleton:
  // every join()/newSimplex()/removeSimplex() call destroys and marks
  // uncomputed subtri_'s ENTIRE cached skeleton (see addFace()'s
  // definition), so any Face<subdim,k>* obtained before such a call is
  // unsafe to compare afterward, and re-deriving it costs a full rebuild
  // proportional to subtri_'s current size -- unusable for a check meant to
  // run once per addFace()/removeFace() call.
  std::vector<std::vector<std::vector<int>>> classRoots_;
  std::vector<RollbackUnionFind> dsu_;

  struct RegistryUndoEntry {
    size_t v;
    int root;
    bool wasInsert;    // true: undo by erasing root from classRoots_[k][v];
                        // false: undo by re-inserting it.
    int singularDelta; // exact singularCount_ change this event caused;
                        // undone via singularCount_ -= singularDelta.
  };
  std::vector<std::vector<RegistryUndoEntry>> registryUndoLog_; // indexed [k]

  struct Checkpoint {
    std::vector<size_t> dsuMark;      // dsu_[k].checkpoint(), per k
    std::vector<size_t> registryMark; // registryUndoLog_[k].size(), per k
  };
  std::vector<Checkpoint> checkpoints_; // indexed by ambient face index f

  int singularCount_ = 0;
  bool isEmbedded_ = true; // kept in lockstep with (singularCount_ == 0)

  // Unions slot (f-implicit via caller, k, local index encoded in slotA)
  // with slot (dst-implicit, k, local index encoded in slotB), both mapping
  // to ambient k-face v. If this merges two previously-independently-
  // registered classes at v, updates classRoots_[k][v]/singularCount_ and
  // logs the change for removeFace() to undo.
  void unite_(int k, size_t v, int slotA, int slotB);

  // Registers slot r's current DSU root against ambient k-face v's set of
  // known classes, if not already present. If this is the second distinct
  // class registered at v, marks v as a new singularity.
  void registerRoot_(int k, size_t v, int r);

public:
  EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton);

  // Seeded constructor: starts from the same empty state as above, then adds
  // every face in seedFaces via addFace(), in the given order. seedFaces is
  // a precondition assumed valid by the caller (a connected, jointly
  // addable set of faces): if any addFace() call fails, throws
  // regina::InvalidArgument rather than silently leaving a partial result.
  EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton,
                      const std::vector<int> &seedFaces);

  // Attempts to add face f to the current submanifold. This fails if adding the
  // face results in something which is not embedded, and returns false.
  // Otherwise, returns true and updates the state of the current submanifold
  // triangulation.
  //
  // NOTE: as of the Phase 2 check being disabled (see addFace()'s
  // definition), this only guarantees genuine embeddedness at the
  // facet (codimension-1) level -- it may accept faces whose addition
  // produces a submanifold that is not injective at codimension >= 2
  // (self-touching vertices/edges below the facet level, always cusp
  // rather than transversal intersections). Use isEmbedded() to check
  // codimension >= 2 injectivity of the current state.
  bool addFace(int f);

  // Attempts to add every face in `faces`, in the given order, via
  // addFace(). Transactional like addFace(): on success, every face in
  // `faces` is committed and this returns true; on failure, rolls back
  // everything added during this call and returns false.
  bool addFaces(const std::vector<int> &faces);

  void removeFace(int f);

  const regina::Triangulation<subdim> &triangulation() const { return subtri_; }

  // Unfortunately, Regina::Triangulation::isClosed() is only valid for
  // triangulations of dimensions 2, 3, and 4
  bool isClosed() const { return !subtri_.hasBoundaryFacets(); }

  bool isProper() const;

  bool boundaryComponentsMapInjectively() const;

  // True iff the realization map subtri_ -> tri_ is injective on every face
  // currently present, at every dimension -- not just facets (Phase 1 above
  // already guarantees that structurally). Incrementally maintained by
  // addFace()/removeFace(); O(1).
  bool isEmbedded() const { return isEmbedded_; }

  bool satisfies(BoundaryCondition cond) const;

  static bool hasIrreparableSelfGluing(
      const std::vector<typename Skeleton<dim, subdim>::Gluing> &gluings);

  // True iff `face` has two of its own local k-subfaces (k <= subdim-2, i.e.
  // vertices when subdim == 2) coinciding as the same ambient object WITHOUT
  // this being forced by one of face's own facet-level self-gluing entries
  // in `gluings`. A collision that IS explained by a self-gluing entry (e.g.
  // a one-triangle Mobius band self-fold) is legitimate: addFace()'s Phase 3
  // self-joins mirror it inside subtri_ too, keeping the realization map
  // injective. An unexplained one means subtri_ would keep the two abstract
  // objects distinct while the ambient triangulation has already identified
  // them -- not a valid embedding on its own, regardless of the rest of the
  // subcomplex (though it may become valid in combination with other faces
  // not present here; this check is deliberately a sound-but-incomplete,
  // single-face-only approximation, not a full hereditary characterization).
  static bool hasUnexplainedSelfCollision(
      const typename Skeleton<dim, subdim>::Face *face,
      const std::vector<typename Skeleton<dim, subdim>::Gluing> &gluings);
};

extern template class EmbeddedSubmanifold<3, 2>;
extern template class EmbeddedSubmanifold<4, 2>;

// A submanifold embedding specialized to dim == 4, subdim == 2: a surface
// embedded in a 4-manifold, in the sense of classical knotted surface theory.
// Everything beyond EmbeddedSubmanifold<4, 2> itself is surface-specific:
// boundaryLinks() identifies the links bounded by the surface's boundary
// (which needs 3-manifold boundary components to make sense of), and
// surfaceTypeKey()/formatSurfaceType() classify the surface's own topology
// (genus, orientability, punctures).
class KnottedSurface : public EmbeddedSubmanifold<4, 2> {
public:
  using SurfaceTypeKey = std::tuple<bool, int, int>;

private:
  std::vector<regina::Triangulation<3>> bdryComponents_;

public:
  KnottedSurface(const Skeleton<4, 2> &skeleton);

  KnottedSurface(const Skeleton<4, 2> &skeleton,
                 const std::vector<int> &seedFaces);

  std::vector<std::pair<size_t, Link>> boundaryLinks() const;

  static SurfaceTypeKey surfaceTypeKey(const regina::Triangulation<2> &surface);

  static std::string formatSurfaceType(const SurfaceTypeKey &key);

  SurfaceTypeKey surfaceType() const { return surfaceTypeKey(triangulation()); }
};

#endif // EMBEDDEDSUBMANIFOLD_H
