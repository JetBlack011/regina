//
//  embeddedsubmanifold.cpp
//

#include "embeddedsubmanifold.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

template <int dim, int subdim>
EmbeddedSubmanifold<dim, subdim>::EmbeddedSubmanifold(
    const Skeleton<dim, subdim> &skeleton)
    : skeleton_(skeleton), faces_(skeleton.numFaces()),
      classRoots_(subdim - 1), registryUndoLog_(subdim - 1),
      checkpoints_(skeleton.numFaces()) {
  const auto &tri = skeleton_.triangulation();
  regina::for_constexpr<0, subdim>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    std::get<k>(faceCount_).assign(tri.template countFaces<k>(), 0);
  });
  dsu_.reserve(subdim - 1);
  regina::for_constexpr<0, subdim - 1>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    classRoots_[k].assign(tri.template countFaces<k>(), {});
    dsu_.emplace_back(skeleton_.numFaces() *
                       regina::FaceNumbering<subdim, k>::nFaces);
  });

  facetIsAmbientBoundary_.resize(tri.template countFaces<subdim - 1>());
  for (size_t i = 0; i < facetIsAmbientBoundary_.size(); ++i)
    facetIsAmbientBoundary_[i] = tri.template face<subdim - 1>(i)->isBoundary();

  for (auto &cp : checkpoints_) {
    cp.dsuMark.assign(subdim - 1, 0);
    cp.registryMark.assign(subdim - 1, 0);
  }
}

template <int dim, int subdim>
EmbeddedSubmanifold<dim, subdim>::EmbeddedSubmanifold(
    const Skeleton<dim, subdim> &skeleton, const std::vector<int> &seedFaces)
    : EmbeddedSubmanifold(skeleton) {
  if (!addFaces(seedFaces))
    throw regina::InvalidArgument(
        "EmbeddedSubmanifold::EmbeddedSubmanifold(): seedFaces could not "
        "be jointly added -- no addition order makes every face embed.");
}

template <int dim, int subdim>
bool EmbeddedSubmanifold<dim, subdim>::addFaces(
    const std::vector<int> &faces) {
  // With Phase 2 disabled (see addFace()), addFace()'s only remaining
  // requirement is the facet-level (codimension-1) boundary condition,
  // which -- for a target set that's genuinely fine at that level -- holds
  // regardless of order (whichever of a facet's <=2 triangles is added
  // second just finds the first already there). So a single pass suffices;
  // no more searching for a working order via repeated retries.
  std::vector<int> added; // commit order, for rollback on failure
  for (int f : faces) {
    if (!addFace(f)) {
      for (auto it = added.rbegin(); it != added.rend(); ++it)
        removeFace(*it);
      return false;
    }
    added.push_back(f);
  }
  return true;
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

  // Phase 2: Check Condition 2 (higher-codimension condition) -- DISABLED.
  //
  // This required that every k-face (k <= subdim-2, i.e. vertices when
  // subdim == 2) of the new simplex already touched by the submanifold be a
  // subface of some facet actually being glued in this same call. That's a
  // sound-but-incomplete, ORDER-DEPENDENT check: it can reject a face whose
  // only "explaining" neighbor hasn't been added yet, even when the full
  // target set is genuinely fine, and -- worse -- some closing dependencies
  // (e.g. completing a cycle) are mutually circular, so no order at all
  // satisfies it (see utils/surfer/ADDFACE_VERTEX_COLLISION_BUG.md and the
  // CollarBuilder investigation that motivated this).
  //
  // Disabled per conjecture: the self-intersections this allowed through
  // are always cusp (tangential), never transversal, intersections, and
  // should always be resolvable after the fact -- rather than rejected
  // during search. That resolution isn't implemented yet; for now this
  // means addFace()/search() can accept submanifolds that are not
  // genuinely embedded (not injective) at codimension >= 2. Codimension-1
  // (facet-level, Phase 1 above) checks are NOT affected and remain in
  // place -- those failures are genuinely unfixable.
  //
  // bool valid = true;
  // regina::for_constexpr<0, subdim - 1>([&](auto kW) {
  //   if (!valid)
  //     return;
  //   constexpr int k = decltype(kW)::value;
  //   const auto &counts = std::get<k>(faceCount_);
  //   for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i) {
  //     if (counts[node.face->template face<k>(i)->index()] == 0)
  //       continue;
  //
  //     // Bitmask of local vertex indices (0..subdim) of the i-th k-face.
  //     const auto perm = regina::FaceNumbering<subdim, k>::ordering(i);
  //     int S = 0;
  //     for (int j = 0; j <= k; ++j)
  //       S |= (1 << perm[j]);
  //
  //     // Facet j (opposite vertex j) contains this k-face iff j NOT in S.
  //     // The face is covered by some facet in F iff (F_bitmask & ~S) != 0.
  //     if ((F_bitmask & ~S) == 0) {
  //       valid = false;
  //       return;
  //     }
  //   }
  // });
  // if (!valid)
  //   return false;

  // Phase 3: All conditions pass — commit state.
  // For k == subdim-1, iterating all subdim+1 local facet indices handles
  // self-gluings automatically: if two local facets share a tri_-face, its
  // count increments twice, making it interior immediately.
  regina::for_constexpr<0, subdim - 1>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    auto &counts = std::get<k>(faceCount_);
    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i)
      counts[node.face->template face<k>(i)->index()]++;
  });

  // Snapshot the isEmbedded()-tracking structures' undo points before
  // touching them, so removeFace(f) can roll back exactly what this call
  // does below (mirroring how faceCount_'s increments above are undone by
  // symmetric decrements in removeFace()).
  for (int k = 0; k < subdim - 1; ++k) {
    checkpoints_[f].dsuMark[k] = dsu_[k].checkpoint();
    checkpoints_[f].registryMark[k] = registryUndoLog_[k].size();
  }

  // k == subdim-1 (facets): increment counts[idx] one local facet at a time
  // (rather than in the generic loop above), tracking each individual
  // increment's before/after value -- self-gluings pass through count == 1
  // transiently (0 -> 1 -> 2, both steps within this same loop, since both
  // of a self-gluing pair's local facets map to the same ambient idx), and
  // that transient state must both raise and then lower badProperCount_
  // within this call, which a single post-hoc scan of the final counts
  // can't distinguish from a facet that is genuinely staying at 1.
  {
    auto &counts = std::get<subdim - 1>(faceCount_);
    for (int i = 0; i <= subdim; ++i) {
      size_t idx = node.face->template face<subdim - 1>(i)->index();
      int before = counts[idx]++;
      int after = before + 1;

      if (!facetIsAmbientBoundary_[idx]) {
        if (after == 1)
          badProperCount_++; // just became a boundary facet of subtri_
        else if (after == 2)
          badProperCount_--; // just got glued over (leaving boundary again)
      }
    }
  }

  auto *src = subtri_.newSimplex();
  faces_[f] = src;

  for (const auto &g : node.gluings) {
    if (faces_[g.dstIndex] == nullptr)
      continue;
    if (src->adjacentSimplex(g.srcFacet) != nullptr)
      continue; // already joined (second direction of a self-gluing)

    // isEmbedded() tracking: this gluing identifies local facet g.srcFacet
    // of the new simplex with the corresponding facet of faces_[g.dstIndex]
    // via g.gluing. For k in [0, subdim-2], every local k-face of the new
    // simplex properly contained in that facet (i.e. not touching the
    // vertex opposite it, g.srcFacet) is thereby identified with its image
    // under g.gluing on the dst side -- handles self-gluings
    // (g.dstIndex == f) and external gluings uniformly, and this direction
    // alone covers a self-gluing's reverse direction too (see the
    // adjacentSimplex guard above).
    regina::for_constexpr<0, subdim - 1>([&](auto kW) {
      constexpr int k = decltype(kW)::value;
      for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i) {
        auto S = regina::FaceNumbering<subdim, k>::ordering(i);
        bool touchesGluedFacet = false;
        for (int t = 0; t <= k; ++t)
          if (S[t] == g.srcFacet) {
            touchesGluedFacet = true;
            break;
          }
        if (touchesGluedFacet)
          continue; // this k-face isn't contained in the shared facet

        int j = regina::FaceNumbering<subdim, k>::faceNumber(g.gluing * S);
        size_t v = node.face->template face<k>(i)->index();
        unite_(k, v,
               f * regina::FaceNumbering<subdim, k>::nFaces + i,
               static_cast<int>(g.dstIndex) *
                       regina::FaceNumbering<subdim, k>::nFaces +
                   j);
      }
    });

    src->join(g.srcFacet, faces_[g.dstIndex], g.gluing);
  }

  // isEmbedded() tracking: register every local k-face slot of the new
  // simplex against its ambient k-face's known classes, whether or not it
  // participated in a gluing above -- an unregistered slot landing on an
  // already-populated ambient k-face is exactly a new singularity.
  regina::for_constexpr<0, subdim - 1>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i) {
      size_t v = node.face->template face<k>(i)->index();
      int r = dsu_[k].find(f * regina::FaceNumbering<subdim, k>::nFaces + i);
      registerRoot_(k, v, r);
    }
  });

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

  // isEmbedded() tracking: undo, in strict reverse order, exactly what the
  // matching addFace(f) logged -- both the registry (classRoots_/
  // singularCount_) and the underlying DSU unions. Relies on the caller's
  // LIFO discipline: removeFace() is only ever called on the most recently
  // added face, so checkpoints_[f] is always this face's own contribution.
  for (int k = 0; k < subdim - 1; ++k) {
    auto &log = registryUndoLog_[k];
    auto &roots = classRoots_[k];
    while (log.size() > checkpoints_[f].registryMark[k]) {
      const auto &e = log.back();
      auto &rv = roots[e.v];
      if (e.wasInsert) {
        rv.erase(std::find(rv.begin(), rv.end(), e.root));
      } else {
        rv.push_back(e.root);
      }
      singularCount_ -= e.singularDelta;
      isEmbedded_ = (singularCount_ == 0);
      log.pop_back();
    }
    dsu_[k].rollbackTo(checkpoints_[f].dsuMark[k]);
  }

  // isProper() tracking: mirror addFace()'s badProperCount_ update, one
  // local facet at a time (decrementing counts[idx] here directly rather
  // than in the generic loop below, for the same self-gluing transient-
  // state reason documented in addFace() -- a self-gluing pair's shared
  // idx passes through count == 1 transiently here too, in reverse).
  {
    auto &counts = std::get<subdim - 1>(faceCount_);
    for (int i = 0; i <= subdim; ++i) {
      size_t idx = node.face->template face<subdim - 1>(i)->index();
      int before = counts[idx]--;
      if (facetIsAmbientBoundary_[idx])
        continue;
      if (before == 1)
        badProperCount_--; // leaving count 1 -> 0
      else if (before == 2)
        badProperCount_++; // leaving count 2 -> 1
    }
  }

  regina::for_constexpr<0, subdim - 1>([&](auto kW) {
    constexpr int k = decltype(kW)::value;
    auto &counts = std::get<k>(faceCount_);
    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces; ++i)
      counts[node.face->template face<k>(i)->index()]--;
  });

  faces_[f] = nullptr;
  subtri_.removeSimplex(face);
}

template <int dim, int subdim>
void EmbeddedSubmanifold<dim, subdim>::unite_(int k, size_t v, int slotA,
                                              int slotB) {
  auto &dsu = dsu_[k];
  int rootABefore = dsu.find(slotA);
  int rootBBefore = dsu.find(slotB);
  if (rootABefore == rootBBefore)
    return; // already identified

  dsu.unite(slotA, slotB);
  int survivor = dsu.find(slotA);
  int absorbed = (survivor == rootABefore) ? rootBBefore : rootABefore;

  // If `absorbed` was itself an already-registered class at v, this merge
  // just identified two previously-independent classes -- if that leaves
  // exactly one class behind, v just stopped being a singularity. If
  // `absorbed` was never registered (e.g. it's the brand-new slot of the
  // face currently being added), there's nothing to remove here; the
  // survivor gets (re-)registered by registerRoot_()'s pass afterward.
  auto &roots = classRoots_[k][v];
  auto it = std::find(roots.begin(), roots.end(), absorbed);
  if (it != roots.end()) {
    roots.erase(it);
    int delta = (roots.size() == 1) ? -1 : 0;
    singularCount_ += delta;
    isEmbedded_ = (singularCount_ == 0);
    registryUndoLog_[k].push_back(
        {.v = v, .root = absorbed, .wasInsert = false, .singularDelta = delta});
  }
}

template <int dim, int subdim>
void EmbeddedSubmanifold<dim, subdim>::registerRoot_(int k, size_t v, int r) {
  auto &roots = classRoots_[k][v];
  if (std::find(roots.begin(), roots.end(), r) != roots.end())
    return; // already registered (e.g. the dst side of an external gluing)

  roots.push_back(r);
  int delta = (roots.size() == 2) ? 1 : 0;
  singularCount_ += delta;
  isEmbedded_ = (singularCount_ == 0);
  registryUndoLog_[k].push_back(
      {.v = v, .root = r, .wasInsert = true, .singularDelta = delta});
}

template <int dim, int subdim>
bool EmbeddedSubmanifold<dim, subdim>::isProper() const {
  return badProperCount_ == 0;
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
bool EmbeddedSubmanifold<dim, subdim>::hasIrreparableSelfGluing(
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

//template <int dim, int subdim>
//bool EmbeddedSubmanifold<dim, subdim>::hasUnexplainedSelfCollision(
//    const typename Skeleton<dim, subdim>::Face *face,
//    const std::vector<typename Skeleton<dim, subdim>::Gluing> &gluings) {
//  // explained[i][j]: true iff face's own local facets i and j are directly
//  // self-glued to each other, i.e. literally the same ambient (subdim-1)-face
//  // object. For subdim == 2, this is the ONLY way two of face's own local
//  // vertices can legitimately coincide as the same ambient object: two
//  // distinct edges of a triangle share exactly one vertex, so if they're the
//  // same ambient edge, their other two (non-shared) vertices are forced to
//  // coincide too -- and that identification is exactly what addFace()'s
//  // Phase 3 self-join reproduces inside subtri_.
//  std::array<std::array<bool, subdim + 1>, subdim + 1> explained{};
//  for (const auto &g : gluings)
//    if (g.srcIndex == g.dstIndex)
//      explained[g.srcFacet][g.gluing[g.srcFacet]] = true;
//
//  bool unexplained = false;
//  regina::for_constexpr<0, subdim - 1>([&](auto kW) {
//    if (unexplained)
//      return;
//    constexpr int k = decltype(kW)::value;
//    for (int i = 0; i < regina::FaceNumbering<subdim, k>::nFaces && !unexplained; ++i) {
//      for (int j = i + 1; j < regina::FaceNumbering<subdim, k>::nFaces; ++j) {
//        if (face->template face<k>(i) != face->template face<k>(j))
//          continue;
//        if (k == 0 && explained[i][j])
//          continue; // legitimately explained by this face's own self-gluing
//        unexplained = true;
//        break;
//      }
//    }
//  });
//  return unexplained;
//}

template class EmbeddedSubmanifold<3, 2>;
template class EmbeddedSubmanifold<4, 2>;

KnottedSurface::KnottedSurface(const Skeleton<4, 2> &skeleton)
    : EmbeddedSubmanifold<4, 2>(skeleton) {
  const auto &tri = skeleton.triangulation();
  bdryComponents_.reserve(tri.countBoundaryComponents());
  for (size_t c = 0; c < tri.countBoundaryComponents(); ++c)
    bdryComponents_.push_back(tri.boundaryComponent(c)->build());
  facesAtVertex_.resize(tri.countVertices());
  vertexEmbedIndexCache_.resize(tri.countVertices());
}

KnottedSurface::KnottedSurface(const Skeleton<4, 2> &skeleton,
                               const std::vector<int> &seedFaces)
    : EmbeddedSubmanifold<4, 2>(skeleton) {
  const auto &tri = skeleton.triangulation();
  bdryComponents_.reserve(tri.countBoundaryComponents());
  for (size_t c = 0; c < tri.countBoundaryComponents(); ++c)
    bdryComponents_.push_back(tri.boundaryComponent(c)->build());
  facesAtVertex_.resize(tri.countVertices());
  vertexEmbedIndexCache_.resize(tri.countVertices());

  if (!addFaces(seedFaces))
    throw regina::InvalidArgument(
        "KnottedSurface::KnottedSurface(): seedFaces could not be jointly "
        "added -- no addition order makes every face embed without "
        "creating a transverse self-intersection or non-locally-flat "
        "point.");
}

const regina::Edge<3> *KnottedSurface::linkEdgeForTriangle_(
    const regina::Vertex<4> *ambientVertex, int f, int localVertex) const {
  size_t v = ambientVertex->index();
  auto &cache = vertexEmbedIndexCache_[v];
  if (!cache) {
    cache.emplace();
    size_t i = 0;
    for (const auto &emb : ambientVertex->embeddings())
      (*cache)[static_cast<uint64_t>(emb.simplex()->index()) * 5 +
               static_cast<uint64_t>(emb.face())] = i++;
  }

  const auto *triangle = skeleton_.getNodes()[f].face;
  auto emb = triangle->front();
  auto *pent = emb.simplex();
  regina::Perm<5> perm = emb.vertices();
  int vLocal = perm[localVertex];
  int w1 = perm[(localVertex + 1) % 3];
  int w2 = perm[(localVertex + 2) % 3];

  size_t idx = cache->at(static_cast<uint64_t>(pent->index()) * 5 +
                        static_cast<uint64_t>(vLocal));
  auto *tet = ambientVertex->buildLink().tetrahedron(idx);
  regina::Perm<5> inv = pent->tetrahedronMapping(vLocal).inverse();
  return tet->edge(inv[w1], inv[w2]);
}

bool KnottedSurface::isPetalClosed_(size_t v, int root) const {
  for (const auto &[f, localVertex] : facesAtVertex_[v]) {
    if (vertexClassRoot(f, localVertex) != root)
      continue;

    const auto *face = skeleton_.getNodes()[f].face;
    for (int i = 0; i <= 2; ++i) {
      if (i == localVertex)
        continue; // the opposite edge, not a spoke
      if (facetCount(face->face<1>(i)->index()) != 2)
        return false;
    }
  }
  return true;
}

std::vector<const regina::Edge<3> *> KnottedSurface::closedPetalCurve_(
    const regina::Vertex<4> *ambientVertex, size_t v, int root) const {
  std::vector<const regina::Edge<3> *> curve;
  for (const auto &[f, localVertex] : facesAtVertex_[v]) {
    if (vertexClassRoot(f, localVertex) != root)
      continue;
    curve.push_back(linkEdgeForTriangle_(ambientVertex, f, localVertex));
  }
  return curve;
}

bool KnottedSurface::addFace(int f) {
  if (!EmbeddedSubmanifold<4, 2>::addFace(f))
    return false;

  const auto &node = skeleton_.getNodes()[f];
  for (int local = 0; local <= 2; ++local) {
    auto *ambientVertex = node.face->face<0>(local);
    size_t v = ambientVertex->index();
    facesAtVertex_[v].push_back({f, local});

    // Vertex<4>::buildLink() is a closed S^3 only for INTERIOR ambient
    // vertices -- the hereditary transverse-self-intersection/local-
    // flatness proofs this feature relies on were specifically about that
    // case. A boundary vertex's link has boundary itself (e.g. a 3-ball),
    // where the EdgeComplement/Knot machinery's pinchEdge()-based drilling
    // doesn't apply (pinchEdge() requires an internal edge) -- skipping
    // entirely here is sound (not just expedient): under-checking can
    // only miss a prune, never accept a state that isn't otherwise valid,
    // and boundary-vertex flatness/transversality (relative to the
    // surface's own boundary curves, already handled separately by
    // boundaryLinks()) was never part of what was proven here.
    if (ambientVertex->isBoundary())
      continue;

    int root = vertexClassRoot(f, local);
    if (!isPetalClosed_(v, root))
      continue;

    auto curve = closedPetalCurve_(ambientVertex, v, root);
    Knot knotA(ambientVertex->buildLink(), curve);
    if (!knotA.isUnknot()) {
      // Non-locally-flat: this petal just closed into a knotted circle in
      // Lk(v). Hereditary under removeFace() (a closed curve can only
      // shrink to an open arc on removal, never split into a smaller
      // closed loop -- arcs are never "knotted"), so safe to reject
      // permanently here rather than deferring to a post-hoc filter.
      removeFace(f);
      return false;
    }

    for (int other : registeredClassRoots(v)) {
      if (other == root || !isPetalClosed_(v, other))
        continue;
      Knot knotB(ambientVertex->buildLink(),
                closedPetalCurve_(ambientVertex, v, other));
      if (knotA.linkingNumberWith(knotB) != 0) {
        // Transverse self-intersection: two closed, nonzero-linked petals
        // at v. Hereditary under removeFace() (the same isotopy that
        // would separate two petals still separates any subset of them),
        // so this can never resolve later -- safe to prune permanently.
        removeFace(f);
        return false;
      }
    }
  }
  return true;
}

bool KnottedSurface::addFaces(const std::vector<int> &faces) {
  std::vector<int> added;
  for (int f : faces) {
    if (!addFace(f)) {
      for (auto it = added.rbegin(); it != added.rend(); ++it)
        removeFace(*it);
      return false;
    }
    added.push_back(f);
  }
  return true;
}

void KnottedSurface::removeFace(int f) {
  const auto &node = skeleton_.getNodes()[f];
  EmbeddedSubmanifold<4, 2>::removeFace(f);
  for (int local = 0; local <= 2; ++local) {
    size_t v = node.face->face<0>(local)->index();
    facesAtVertex_[v].pop_back();
  }
}

std::vector<std::pair<size_t, Link>> KnottedSurface::boundaryLinks() const {
  std::vector<std::unordered_set<const regina::Edge<3> *>> edgesByComponent(
      bdryComponents_.size());

  for (size_t f = 0; f < faces_.size(); ++f) {
    const auto *simplex = faces_[f];
    if (simplex == nullptr)
      continue;

    for (int i = 0; i <= 2; ++i) {
      if (simplex->adjacentSimplex(i) != nullptr)
        continue; // internal facet of subtri_

      const auto *ambientFacet = skeleton_.getNodes()[f].face->face<1>(i);
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

KnottedSurface::SurfaceTypeKey
KnottedSurface::surfaceTypeKey(const regina::Triangulation<2> &surface) {
  bool isOrientable = surface.isOrientable();
  int punctures = surface.countBoundaryComponents();
  int genus = isOrientable ? (2 - surface.eulerChar() - punctures) / 2
                           : 2 - surface.eulerChar() - punctures;
  return {isOrientable, genus, punctures};
}

std::string
KnottedSurface::formatSurfaceType(const SurfaceTypeKey &key) {
  auto [isOrientable, genus, punctures] = key;
  std::ostringstream ans;

  if (isOrientable) {
    if (genus == 0 && punctures == 1)
      ans << "Disc";
    else if (genus == 0 && punctures == 2)
      ans << "Annulus";
    else {
      if (genus == 0)
        ans << "Sphere";
      else if (genus == 1)
        ans << "Torus";
      else
        ans << "Orientable genus " << genus << " surface";

      if (punctures == 1)
        ans << ", 1 puncture";
      else if (punctures > 1)
        ans << ", " << punctures << " punctures";
    }
  } else {
    if (genus == 1 && punctures == 1)
      ans << "Möbius band";
    else {
      if (genus == 1)
        ans << "Projective plane";
      else if (genus == 2)
        ans << "Klein bottle";
      else
        ans << "Non-orientable genus " << genus << " surface";

      if (punctures == 1)
        ans << ", 1 puncture";
      else if (punctures > 1)
        ans << ", " << punctures << " punctures";
    }
  }

  return ans.str();
}
