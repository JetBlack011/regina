//
//  embeddedsubmanifold.cpp
//

#include "embeddedsubmanifold.h"

#include <array>
#include <cassert>
#include <unordered_map>
#include <unordered_set>

template <int dim, int subdim>
EmbeddedSubmanifold<dim, subdim>::EmbeddedSubmanifold(
    const Skeleton<dim, subdim> &skeleton)
    : skeleton_(skeleton), faces_(skeleton.numFaces()) {
  const auto &tri = skeleton_.triangulation();
  regina::for_constexpr<0, subdim>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    std::get<k>(faceCount_).assign(tri.template countFaces<k>(), 0);
  });
  if constexpr (dim == 4 && subdim == 2) {
    bdryComponents_.reserve(tri.countBoundaryComponents());
    for (size_t c = 0; c < tri.countBoundaryComponents(); ++c)
      bdryComponents_.push_back(tri.boundaryComponent(c)->build());
  }
}

template <int dim, int subdim>
bool EmbeddedSubmanifold<dim, subdim>::addFace(int f) {
  const auto &node = skeleton_.getNodes()[f];

  // Phase 1: Check Condition 1 (facet condition) and build F_bitmask.
  // F_bitmask has bit j set iff local facet j of the new simplex is being
  // glued to a boundary facet of the current submanifold.
  std::array<bool, subdim + 1> facetIsUsed = {};
  int F_bitmask = 0;

  for (const auto &g : node.gluings) {
    const int dstFacet = g.gluing[g.srcFacet];

    if (g.dstIndex == static_cast<size_t>(f)) {
      // Self-gluing: facets g.srcFacet and dstFacet of this simplex are
      // the same (subdim-1)-face in tri_. Each pair is listed twice in
      // gluings; process only srcFacet < dstFacet to avoid double-work.
      if (g.srcFacet >= dstFacet)
        continue;

      if (facetIsUsed[g.srcFacet] || facetIsUsed[dstFacet])
        return false;

      // Adding this simplex will increment the count for the shared face by
      // 2 (once per local facet). It must currently be absent.
      if (std::get<subdim - 1>(faceCount_)[node.face->template face<subdim - 1>(g.srcFacet)->index()] != 0)
        return false;

      facetIsUsed[g.srcFacet] = true;
      facetIsUsed[dstFacet] = true;
      // Self-gluings are NOT added to F_bitmask: they are internal to the
      // new simplex, not gluings to the existing submanifold.

    } else if (faces_[g.dstIndex] != nullptr) {
      // External gluing to an existing simplex in the submanifold.
      if (facetIsUsed[g.srcFacet])
        return false;

      // The destination facet must be a boundary facet (count == 1).
      if (std::get<subdim - 1>(faceCount_)[node.face->template face<subdim - 1>(g.srcFacet)->index()] != 1)
        return false;

      facetIsUsed[g.srcFacet] = true;
      F_bitmask |= (1 << g.srcFacet);
    }
    // Gluings to faces not in the submanifold are ignored.
  }

  // Phase 2: Check Condition 2 (higher-codimension condition).
  // Every k-face (k <= subdim-2) of the new simplex that is already in the
  // submanifold must be a subface of some facet in F.
  bool valid = true;
  regina::for_constexpr<0, subdim - 1>([&](auto kW) {
    if (!valid)
      return;
    constexpr int k = decltype(kW)::value;
    const auto &counts = std::get<k>(faceCount_);
    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i) {
      if (counts[node.face->template face<k>(i)->index()] == 0)
        continue;

      // Bitmask of local vertex indices (0..subdim) of the i-th k-face.
      const auto perm = regina::FaceNumbering<subdim, k>::ordering(i);
      int S = 0;
      for (int j = 0; j <= k; ++j)
        S |= (1 << perm[j]);

      // Facet j (opposite vertex j) contains this k-face iff j NOT in S.
      // The face is covered by some facet in F iff (F_bitmask & ~S) != 0.
      if ((F_bitmask & ~S) == 0) {
        valid = false;
        return;
      }
    }
  });
  if (!valid)
    return false;

  // Phase 3: All conditions pass — commit state.
  // For k == subdim-1, iterating all subdim+1 local facet indices handles
  // self-gluings automatically: if two local facets share a tri_-face, its
  // count increments twice, making it interior immediately.
  regina::for_constexpr<0, subdim>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    auto &counts = std::get<k>(faceCount_);
    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i)
      counts[node.face->template face<k>(i)->index()]++;
  });

  auto *src = subtri_.newSimplex();
  faces_[f] = src;

  for (const auto &g : node.gluings) {
    if (faces_[g.dstIndex] == nullptr)
      continue;
    if (src->adjacentSimplex(g.srcFacet) != nullptr)
      continue; // already joined (second direction of a self-gluing)
    src->join(g.srcFacet, faces_[g.dstIndex], g.gluing);
  }

  return true;
}

template <int dim, int subdim>
void EmbeddedSubmanifold<dim, subdim>::removeFace(int f) {
  auto *face = faces_[f];
  if (face == nullptr)
    throw regina::InvalidArgument(
        "EmbeddedSubmanifold::removeFace(): Was asked to remove a face "
        "which is not in the embedded submanifold.");

  const auto &node = skeleton_.getNodes()[f];

  regina::for_constexpr<0, subdim>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    auto &counts = std::get<k>(faceCount_);
    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i)
      counts[node.face->template face<k>(i)->index()]--;
  });

  faces_[f] = nullptr;
  subtri_.removeSimplex(face);
}

template <int dim, int subdim>
bool EmbeddedSubmanifold<dim, subdim>::isProper() const {
  for (size_t f = 0; f < faces_.size(); ++f) {
    const auto *simplex = faces_[f];
    if (simplex == nullptr)
      continue; // face f is not part of this embedding

    for (int i = 0; i <= subdim; ++i) {
      if (simplex->adjacentSimplex(i) != nullptr)
        continue; // internal facet of subtri_

      const auto *ambientFacet =
          skeleton_.getNodes()[f].face->template face<subdim - 1>(i);
      if (!ambientFacet->isBoundary())
        return false;
    }
  }
  return true;
}

template <int dim, int subdim>
bool EmbeddedSubmanifold<dim, subdim>::boundaryComponentsMapInjectively()
    const {
  std::unordered_map<size_t, const regina::BoundaryComponent<dim> *>
      subToAmbient;
  std::unordered_set<const regina::BoundaryComponent<dim> *> usedAmbient;

  for (size_t f = 0; f < faces_.size(); ++f) {
    const auto *simplex = faces_[f];
    if (simplex == nullptr)
      continue;

    for (int i = 0; i <= subdim; ++i) {
      if (simplex->adjacentSimplex(i) != nullptr)
        continue;

      const auto *ambientFacet =
          skeleton_.getNodes()[f].face->template face<subdim - 1>(i);
      const auto *ambientBC = ambientFacet->boundaryComponent();
      if (ambientBC == nullptr)
        return false;

      const auto *subFacet = simplex->template face<subdim - 1>(i);
      const auto *subBC = subFacet->boundaryComponent();
      assert(subBC != nullptr);

      size_t subIdx = subBC->index();
      auto it = subToAmbient.find(subIdx);
      if (it == subToAmbient.end()) {
        if (usedAmbient.contains(ambientBC))
          return false;
        subToAmbient.emplace(subIdx, ambientBC);
        usedAmbient.insert(ambientBC);
      } else {
        assert(it->second == ambientBC);
        if (it->second != ambientBC)
          return false;
      }
    }
  }
  return true;
}

template <int dim, int subdim>
std::vector<std::pair<size_t, Link>>
EmbeddedSubmanifold<dim, subdim>::boundaryLinks() const {
  static_assert(dim == 4 && subdim == 2,
                "boundaryLinks() only makes sense for surfaces (dim == "
                "4, subdim == 2)");

  std::vector<std::unordered_set<const regina::Edge<3> *>> edgesByComponent(
      bdryComponents_.size());

  for (size_t f = 0; f < faces_.size(); ++f) {
    const auto *simplex = faces_[f];
    if (simplex == nullptr)
      continue;

    for (int i = 0; i <= subdim; ++i) {
      if (simplex->adjacentSimplex(i) != nullptr)
        continue; // internal facet of subtri_

      const auto *ambientFacet =
          skeleton_.getNodes()[f].face->template face<subdim - 1>(i);
      const auto *ambientBC = ambientFacet->boundaryComponent();
      if (ambientBC == nullptr)
        continue; // shouldn't happen when cond is proper/connected

      size_t c = ambientBC->index();

      for (int k = 0; k < ambientBC->countEdges(); ++k) {
        if (ambientBC->edge(k) == ambientFacet) {
          edgesByComponent[c].insert(bdryComponents_[c].edge(k));
          break;
        }
      }
    }
  }

  std::vector<std::pair<size_t, Link>> result;
  for (size_t c = 0; c < edgesByComponent.size(); ++c) {
    if (edgesByComponent[c].empty())
      continue;

    std::vector<const regina::Edge<3> *> edgeList(edgesByComponent[c].begin(),
                                                  edgesByComponent[c].end());
    result.emplace_back(c, Link(bdryComponents_[c], edgeList));
  }
  return result;
}

template <int dim, int subdim>
bool EmbeddedSubmanifold<dim, subdim>::satisfies(
    BoundaryCondition cond) const {
  switch (cond) {
  case BoundaryCondition::all:
    return true;
  case BoundaryCondition::closed:
    return isClosed();
  case BoundaryCondition::proper:
    return isProper();
  case BoundaryCondition::connected:
    return boundaryComponentsMapInjectively();
  }

  throw regina::InvalidArgument(
      "EmbeddedSubmanifold::satisfies(): Invalid BoundaryCondition");
}

template <int dim, int subdim>
bool hasIrreparableSelfGluing(
    const std::vector<typename Skeleton<dim, subdim>::Gluing> &gluings) {
  std::array<std::array<bool, subdim + 1>, subdim + 1> partnerSeen{};
  for (const auto &g : gluings)
    if (g.srcIndex == g.dstIndex)
      partnerSeen[g.srcFacet][g.gluing[g.srcFacet]] = true;

  for (int i = 0; i <= subdim; ++i) {
    int distinctPartners = 0;
    for (int j = 0; j <= subdim; ++j)
      distinctPartners += partnerSeen[i][j];
    if (distinctPartners >= 2)
      return true;
  }
  return false;
}

// EmbeddedSubmanifold<4, 2> instantiates cleanly as a whole class: its
// boundaryLinks() static_assert(dim == 4 && subdim == 2) holds.
template class EmbeddedSubmanifold<4, 2>;

// EmbeddedSubmanifold<3, 2>, on the other hand, cannot use the blanket
// `template class` form: that would force boundaryLinks() to instantiate
// too, and its static_assert fails for dim == 3. Nothing calls
// boundaryLinks() on a dim=3 submanifold (it's surface-specific), so
// explicitly instantiate every other member individually instead.
template EmbeddedSubmanifold<3, 2>::EmbeddedSubmanifold(
    const Skeleton<3, 2> &);
template bool EmbeddedSubmanifold<3, 2>::addFace(int);
template void EmbeddedSubmanifold<3, 2>::removeFace(int);
template bool EmbeddedSubmanifold<3, 2>::isProper() const;
template bool
EmbeddedSubmanifold<3, 2>::boundaryComponentsMapInjectively() const;
template bool EmbeddedSubmanifold<3, 2>::satisfies(BoundaryCondition) const;

template bool hasIrreparableSelfGluing<3, 2>(
    const std::vector<typename Skeleton<3, 2>::Gluing> &);
template bool hasIrreparableSelfGluing<4, 2>(
    const std::vector<typename Skeleton<4, 2>::Gluing> &);
