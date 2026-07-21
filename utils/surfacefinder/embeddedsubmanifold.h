#ifndef EMBEDDEDSUBMANIFOLD_H

#define EMBEDDEDSUBMANIFOLD_H

#include <array>
#include <unordered_map>
#include <unordered_set>
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
  EmbeddedSubmanifold(const Skeleton<dim, subdim> &skeleton)
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

  bool addFace(int f) {
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

  void removeFace(int f) {
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

  const regina::Triangulation<subdim> &triangulation() const { return subtri_; }

  // Unfortunately, Regina::Triangulation::isClosed() is only valid for
  // triangulations of dimensions 2, 3, and 4
  bool isClosed() const { return !subtri_.hasBoundaryFacets(); }

  bool isProper() const {
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

  bool boundaryComponentsMapInjectively() const {
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

  std::vector<std::pair<size_t, Link>> boundaryLinks() const {
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

  bool satisfies(BoundaryCondition cond) const {
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
};

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

#endif // EMBEDDEDSUBMANIFOLD_H
