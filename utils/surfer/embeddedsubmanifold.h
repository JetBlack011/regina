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
  // rather than transversal intersections). See
  // utils/surfer/ADDFACE_VERTEX_COLLISION_BUG.md.
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
