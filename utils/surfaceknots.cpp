//
//  embeddedsurfaces.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include "surfaceknots.h"

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
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

std::ostream &operator<<(std::ostream &os, const TriangleMap &triMap) {
    os << "{ ";
    for (const auto &[key, val] : triMap) {
        os << "(" << key << "," << val << ") ";
    }
    return os << "}";
}

Knot::Knot(const regina::Triangulation<3> &tri, Curve<3> &edgeList)
    : tri_(tri), edgeList_(edgeList) {
    for (const regina::Edge<3> *edge : edgeList) {
        edgeSet_.insert(edge);
    }
}

Link::Link(const regina::Triangulation<3> &tri,
           std::vector<std::vector<regina::Edge<3> *>>)
    : tri_(tri) {}

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

std::pair<std::vector<int>, std::vector<int>> getTriIndexList(
    const SurfaceKnot &lhs, const SurfaceKnot &rhs) {
    std::vector<int> lhsKeys, rhsKeys;

    for (const auto &[key, _] : lhs.triMap_) {
        lhsKeys.push_back(key);
    }
    for (const auto &[key, _] : rhs.triMap_) {
        rhsKeys.push_back(key);
    }

    std::sort(lhsKeys.begin(), lhsKeys.end());
    std::sort(rhsKeys.begin(), rhsKeys.end());

    return {lhsKeys, rhsKeys};
}

SurfaceKnot::SurfaceKnot(const regina::Triangulation<4> *tri,
                         const regina::Triangulation<2> &surface,
                         const TriangleMap &triMap, SurfaceCondition cond)
    : tri_(tri),
      surface_(surface),
      triMap_(triMap),
      cond_(cond),
      isOrientable_(surface.isOrientable()),
      punctures_(surface.countBoundaryComponents()),
      genus_(isOrientable_ ? (2 - surface.eulerChar() - punctures_) / 2
                           : 2 - surface.eulerChar() - punctures_) {
    std::ostringstream ans;
    if (!surface_.isConnected()) {
        throw regina::InvalidArgument("Surface must be connceted!");
    }

    if (isOrientable_) {
        // Special names for surface_s with boundary:
        if (genus_ == 0 && punctures_ == 1)
            ans << "Disc";
        else if (genus_ == 0 && punctures_ == 2)
            ans << "Annulus";
        else {
            if (genus_ == 0)
                ans << "Sphere";
            else if (genus_ == 1)
                ans << "Torus";
            else
                ans << "genus_ " << genus_ << " torus";

            if (punctures_ == 1)
                ans << ", 1 puncture";
            else if (punctures_ > 1)
                ans << ", " << punctures_ << " punctures";
        }
    } else {
        // Special names for surface_s with boundary:
        if (genus_ == 1 && punctures_ == 1)
            ans << "Möbius band";
        else {
            if (genus_ == 1)
                ans << "Projective plane";
            else if (genus_ == 2)
                ans << "Klein bottle";
            else
                ans << "Non-orientable genus " << genus_ << " surface";

            if (punctures_ == 1)
                ans << ", 1 puncture";
            else if (punctures_ > 1)
                ans << ", " << punctures_ << " punctures";
        }
    }

    detail_ = ans.str();
}

bool SurfaceKnot::isBoundaryContainment() {
    // I hate this function
    // Any way to do this that isn't O(n^2)?
    for (const regina::BoundaryComponent<2> *comp :
         surface_.boundaryComponents()) {
        for (const regina::Edge<2> *edge : comp->edges()) {
            assert(edge->embeddings().size() == 1);  // Sanity check

            const regina::Edge<4> *triEdge = findEdgeInTri(edge);
            if (!triEdge->isBoundary()) {
                return false;
            }
        }
    }

    return true;
}

// TODO: if this is actually slow in practice make a second map
const regina::Edge<4> *SurfaceKnot::findEdgeInTri(
    const regina::Edge<2> *edge) const {
    for (const regina::EdgeEmbedding<2> &embedding : edge->embeddings()) {
        int triIndex = embedding.triangle()->index();
        const regina::Edge<4> *triEdge;
        for (const auto &[key, val] : triMap_) {
            if (val == triIndex) {
                // TODO: Check, hmm...
                return tri_->triangle(key)->edge(edge->index());
            }
        }
    }

    throw regina::InvalidArgument("Given edge doesn't exist in triangle map!");
}

bool SurfaceKnot::isAdmissible() {
    if (cond_ == closed)
        return isOrientable_ && surface_.isClosed();
    else
        return isOrientable_ && isBoundaryContainment();
}

const regina::Triangulation<2> &SurfaceKnot::surface() const {
    return surface_;
}

bool SurfaceKnot::isOrientable() const { return isOrientable_; }

int SurfaceKnot::punctures() const { return punctures_; }

int SurfaceKnot::genus() const { return genus_; }

std::string SurfaceKnot::detail() const { return detail_; }

bool operator==(const SurfaceKnot &lhs, const SurfaceKnot &rhs) {
    auto [lhsKeys, rhsKeys] = getTriIndexList(lhs, rhs);
    return lhsKeys == rhsKeys && lhs.isOrientable_ == rhs.isOrientable_ &&
           lhs.punctures_ == rhs.punctures_ && lhs.genus_ == rhs.genus_;
}

bool operator<(const SurfaceKnot &lhs, const SurfaceKnot &rhs) {
    auto [lhsKeys, rhsKeys] = getTriIndexList(lhs, rhs);
    lhsKeys.push_back(lhs.isOrientable_);
    lhsKeys.push_back(lhs.punctures_);
    lhsKeys.push_back(lhs.genus_);
    rhsKeys.push_back(rhs.isOrientable_);
    rhsKeys.push_back(rhs.punctures_);
    rhsKeys.push_back(rhs.genus_);
    return lhsKeys < rhsKeys;
}

std::ostream &operator<<(std::ostream &os, const SurfaceKnot &surfaceKnot) {
    return os << surfaceKnot.detail_;
}