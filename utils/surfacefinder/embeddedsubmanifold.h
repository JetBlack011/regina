#ifndef EMBEDDEDSUBMANIFOLD_H

#define EMBEDDEDSUBMANIFOLD_H

#include <vector>

#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <utilities/typeutils.h>

#include "linkcomplement.h"
#include "skeleton.h"

enum class BoundaryCondition : uint8_t { all, closed, proper, connected };

template <int dim, int subdim> class EmbeddedSubmanifold {
private:
  using Face = const typename Skeleton<dim, subdim>::Face;
  template <int k> using LowDimVec = std::vector<int>;

  const Skeleton<dim, subdim> &skeleton_;
  regina::Triangulation<subdim> subtri_;
  std::vector<regina::Simplex<subdim> *> faces_;
  // For k in [0, subdim-1]: number of submanifold simplices containing each
  // k-face of tri_. For k == subdim-1 (facets): 0 = absent, 1 = boundary,
  // ≥2 = interior. For k < subdim-1: positive iff the face is in the
  // submanifold.
  regina::TupleOverRange<0, subdim, LowDimVec> faceCount_;

  std::vector<regina::Triangulation<dim - 1>> bdryComponents_;

public:
  EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton);

  bool addFace(int f);

  void removeFace(int f);

  const regina::Triangulation<subdim> &triangulation() const { return subtri_; }

  // Unfortunately, Regina::Triangulation::isClosed() is only valid for
  // triangulations of dimensions 2, 3, and 4
  bool isClosed() const { return !subtri_.hasBoundaryFacets(); }

  bool isProper() const;

  bool boundaryComponentsMapInjectively() const;

  std::vector<std::pair<size_t, Link>> boundaryLinks() const;

  bool satisfies(BoundaryCondition cond) const;

  static bool hasIrreparableSelfGluing(
      const std::vector<typename Skeleton<dim, subdim>::Gluing> &gluings);
};

extern template class EmbeddedSubmanifold<3, 2>;
extern template class EmbeddedSubmanifold<4, 2>;

#endif // EMBEDDEDSUBMANIFOLD_H
