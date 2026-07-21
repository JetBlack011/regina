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

  // Attempts to add face f to the current submanifold. This fails if adding the
  // face results in something which is not embedded, and returns false.
  // Otherwise, returns true and updates the state of the current submanifold
  // triangulation.
  bool addFace(int f);

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

  std::vector<std::pair<size_t, Link>> boundaryLinks() const;

  static SurfaceTypeKey surfaceTypeKey(const regina::Triangulation<2> &surface);

  static std::string formatSurfaceType(const SurfaceTypeKey &key);

  SurfaceTypeKey surfaceType() const { return surfaceTypeKey(triangulation()); }
};

#endif // EMBEDDEDSUBMANIFOLD_H
