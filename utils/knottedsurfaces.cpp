//
//  knottedsurfaces.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include "knottedsurfaces.h"

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <string>
#include <utility>
#include <vector>

#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/faceembedding.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

/** KnotEmbedding Implementation */

/** Knot Implementation */

Knot::Knot(const regina::Triangulation<3> &tri, const Curve<3> &edges)
    : tri_(tri), edges_(edges) {
    // Precompute the set of edges and the count of edges in each tetrahedron
    for (const regina::Edge<3> *edge : edges) {
        edgeSet_.insert(edge);
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            auto search = tetEdgeCount_.find(emb.tetrahedron());
            if (search == tetEdgeCount_.end()) {
                tetEdgeCount_.insert({emb.tetrahedron(), 1});
            } else {
                ++search->second;
            }
        }
    }
}

void Knot::simplify() {
    if (isUnknot_) {
        return;
    }
}

/** Link Implementation */

// std::pair<EdgeMap, regina::Triangulation<3>> drillDualCurve(
//     const regina::Triangulation<3> &tri, const DualCurve<3> &curve) {
//
// }

Link::Link(const regina::Triangulation<3> &tri,
           const std::vector<Curve<3>> &components)
    : tri_(tri) {
    for (const Curve<3> &edges : components) {
        components_.emplace_back(tri, edges);
    }
}

// regina::Triangulation<3> Link::linkComplement(const regina::Triangulation<3>
// &t,
//                                         const Link &l) {
//     regina::Triangulation<3> tri = t;
//     Link link = l;
//
//     while (!link.empty()) {
//         const Knot &knot = link.front();
//         // 1) Push off the current link component
//         DualCurve<3> curve = knotToDualCurve(knot);
//
//         /// 2) Drill, keep track of where the remaining link components end
//
//         const auto &[edgeMap, newTri] = drillDualCurve(tri, curve);
//
//         Link newLink;
//         newLink.reserve(link.size() - 1);
//         for (const Knot &knot2 : link) {
//             if (knot2 == knot) continue;  // Ignore the now drilled-out
//     ot
//
//             Knot newKnot;
//             newKnot.reserve(knot.size());
//             for (const regina::Edge<3> *edge : knot2) {
//                 auto search = edgeMap.find(edge->index());
//                 if (search == edgeMap.end()) {
//                     throw regina::InvalidArgument(
//                         "Could not find corresponding edge in edgeMap!");
//                 }
//                 newKnot.push_back(newTri.edge(search->second));
//             }
//             newLink.push_back(newKnot);
//         }
//
//         link = newLink;
//         tri = newTri;
//         // Repeat with newLink in the drilled triangulation
//     }
//
//     // Now tri is a triangulation of the complement of the original link
//     std::cout << tri.detail() << "\n";
//
//     return tri;
// }

/** Gluing Implementation */

template <int n>
std::ostream &operator<<(std::ostream &os, const Gluing<n> &g) {
    return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
              << g.dst->index() << ", " << g.gluing << ")";
}

/** TriangulationEmbedding Implementation */

template <int dim, int subdim>
TriangulationEmbedding<dim, subdim>::TriangulationEmbedding(
    const regina::Triangulation<dim> &tri)
    : tri_(tri) {
    static_assert(4 <= dim, "Must have dim >= 4");
    static_assert(2 <= subdim <= dim, "Must have 2 <= subdim <= dim");
}

template <int dim, int subdim>
bool TriangulationEmbedding<dim, subdim>::isProper() const {
    return isProper_;
}

template <int dim, int subdim>
void TriangulationEmbedding<dim, subdim>::updateIsProper() {
    for (const regina::BoundaryComponent<subdim> *comp :
         sub_.boundaryComponents()) {
        for (const regina::Face<subdim, subdim - 1> *f : comp->facets()) {
            regina::Face<dim, subdim - 1> *im = image(f);
            if (!im->isBoundary()) {
                isProper_ = false;
                return;
            }
        }
    }

    isProper_ = true;
}

template <int dim, int subdim>
template <int facedim>
regina::Face<dim, facedim> *TriangulationEmbedding<dim, subdim>::image(
    const regina::Face<subdim, facedim> *f) const {
    static_assert(0 <= facedim <= subdim, "Must have 0 <= facedim <= subdim");

    regina::FaceEmbedding<subdim, facedim> &fEmb = f->front();
    return emb_.at(fEmb.simplex()).template face<facedim>(fEmb.face());
}

template <int dim, int subdim>
regina::Face<dim, subdim> *TriangulationEmbedding<dim, subdim>::image(
    const regina::Simplex<subdim> *s) const {
    return emb_.at(s);
}

template <int dim, int subdim>
regina::Face<subdim, subdim> *TriangulationEmbedding<dim, subdim>::preimage(
    const regina::Face<dim, subdim> *s) const {
    return inv_.at(s);
}

template <int dim, int subdim>
regina::Simplex<subdim> *TriangulationEmbedding<dim, subdim>::addFace(
    const regina::Face<dim, subdim> *f) {
    auto search = inv_.find(f);

    if (search != inv_.end()) {
        return search->second;
    }

    regina::Simplex<subdim> *s = sub_.newSimplex();
    emb_[s] = f;
    inv_[f] = s;

    return s;
}

template <int dim, int subdim>
bool TriangulationEmbedding<dim, subdim>::addGluing(Gluing<dim> g) {
    regina::Simplex<subdim> *src = inv_.at(g.src);
    regina::Simplex<subdim> *dst = addFace(g.dst);

    try {
        src->join(g.srcFacet, dst, g.gluing);
    } catch (regina::InvalidArgument &e) {
        removeGluing(g);
        return false;
    }

    if (isProper_) {
        // Quick optimization, boundary containment is entirely determined by
        // the new facets we've made
        for (int i = 0; i < subdim + 1; ++i) {
            if (i == g.gluing[g.srcFacet]) continue;

            const regina::Face<subdim, subdim - 1> *facet =
                dst->template face<subdim - 1>(i);
            const regina::Face<dim, subdim - 1> *im = image(facet);
            if (facet->isBoundary() && !im->isBoundary()) {
                isProper_ = false;
            }
        }
    } else {
        updateIsProper();
    }

    return true;
}

template <int dim, int subdim>
void TriangulationEmbedding<dim, subdim>::removeGluing(Gluing<dim> g) {
    regina::Simplex<subdim> *src = inv_.at(g.src);

    src->unjoin(g.srcFacet);

    if (!sub_.isConnected()) {
        regina::Simplex<subdim> *dst = inv_.at(g.dst);
        sub_.removeSimplex(dst);
        emb_.erase(dst);
        inv_.erase(g.dst);
    }
}

template <int dim, int subdim>
bool operator==(const TriangulationEmbedding<dim, subdim> &lhs,
                const TriangulationEmbedding<dim, subdim> &rhs) {
    return lhs.cmp_ == rhs.cmp_;
}

template <int dim, int subdim>
bool operator<(const TriangulationEmbedding<dim, subdim> &lhs,
               const TriangulationEmbedding<dim, subdim> &rhs) {
    return lhs.cmp_ < rhs.cmp_;
}

template <int dim, int subdim>
std::ostream &operator<<(std::ostream &os,
                         const TriangulationEmbedding<dim, subdim> &emb) {
    if (emb.emb_.empty()) {
        return os << "{}";
    }

    os << "{ ";
    auto it = emb.emb_.begin();
    while (it + 1 != emb.emb_.end()) {
        os << "(" << it->first << "," << it->second << "), ";
        ++it;
    }

    return os << "( " << it->first << "," << it->second << ")";
}

/** KnottedSurface Implementation */

KnottedSurface::KnottedSurface(const regina::Triangulation<4> &tri)
    : TriangulationEmbedding(tri) {
    computeInvariants_();
}

const regina::Triangulation<2> &KnottedSurface::surface() const { return sub_; }

std::string KnottedSurface::detail() const {
    // if (!detail_.empty()) return detail_;

    /* Generate surface details */
    bool isOrientable = sub_.isOrientable();
    int punctures = sub_.countBoundaryComponents();
    int genus = isOrientable ? (2 - sub_.eulerChar() - punctures) / 2
                             : 2 - sub_.eulerChar() - punctures;
    std::ostringstream ans;

    if (!sub_.isConnected()) {
        throw regina::InvalidArgument("Surface must be connceted!");
    }

    if (isOrientable) {
        // Special names for surface_s with boundary:
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
                ans << "genus " << genus << " torus";

            if (punctures == 1)
                ans << ", 1 puncture";
            else if (punctures > 1)
                ans << ", " << punctures << " punctures";
        }
    } else {
        // Special names for surface_s with boundary:
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
    // return detail_ = ans.str();
}

void KnottedSurface::computeInvariants_() {
    bool isOrientable = sub_.isOrientable();
    int punctures = sub_.countBoundaryComponents();
    int genus = isOrientable ? (2 - sub_.eulerChar() - punctures) / 2
                             : 2 - sub_.eulerChar() - punctures;

    invariants_ = {isOrientable, punctures, genus};
}
