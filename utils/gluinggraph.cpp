//
//  gluinggraph.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include "gluinggraph.h"

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <set>
#include <vector>

#include "maths/perm.h"
#include "knottedsurfaces.h"
#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

template <int n>
GluingGraph<n>::GluingNode::GluingNode(const regina::Triangle<n> *src,
                                       int srcFacet,
                                       const regina::Triangle<n> *dst,
                                       regina::Perm<3> gluing)
    : gluing_(src, srcFacet, dst, gluing) {}

template <int n>
bool operator==(const typename GluingGraph<n>::GluingNode &lhs,
                const typename GluingGraph<n>::GluingNode &rhs) {
    // We never want two gluing nodes with the same gluing, use this to
    // distinguish the nodes
    return lhs.gluing_.src == rhs.gluing_.src &&
           lhs.gluing_.srcFacet == rhs.gluing_.srcFacet &&
           lhs.gluing_.dst == rhs.gluing_.dst &&
           lhs.gluing_.gluing == rhs.gluing_.gluing;
}

template <int n>
bool operator<(const typename GluingGraph<n>::GluingNode &lhs,
               const typename GluingGraph<n>::GluingNode &rhs) {
    // Arbitrary tie breaker. If by some miracle you manage to
    // stick in a triangulation with over LONG_MAX triangles, God help
    // you.
    std::vector<long> lhsLex = {
        lhs.gluing_.src == nullptr
            ? -1
            : static_cast<long>(lhs.gluing_.src->index()),
        static_cast<long>(lhs.gluing_.dst->index()),
        lhs.gluing_.srcFacet,
        lhs.gluing_.gluing[0],
        lhs.gluing_.gluing[1],
        lhs.gluing_.gluing[2]};
    std::vector<long> rhsLex = {
        rhs.gluing_.src == nullptr
            ? -1
            : static_cast<long>(rhs.gluing_.src->index()),
        static_cast<long>(rhs.gluing_.dst->index()),
        rhs.gluing_.srcFacet,
        rhs.gluing_.gluing[0],
        rhs.gluing_.gluing[1],
        rhs.gluing_.gluing[2]};

    return lhsLex < rhsLex;
}

template <int n>
std::ostream &operator<<(std::ostream &os,
                         const typename GluingGraph<n>::GluingNode &node) {
    os << "(" << node.gluing_ << " { ";
    for (const typename GluingGraph<n>::GluingNode *adj : node.adjList_) {
        os << adj->gluing_ << " ";
    }

    os << "})";

    return os;
}

template <int n>
std::set<typename GluingGraph<n>::GluingNode> GluingGraph<n>::addTriangle_(
    const regina::Triangle<n> *triangle) {
    std::set<GluingGraph<n>::GluingNode> nodes;

    // For each edge of our triangle, find all of the ways adjacent
    // triangles are glued to that edge in the given triangulation
    for (int facet = 0; facet < 3; ++facet) {
        const regina::Edge<n> *edge = triangle->edge(facet);
        regina::Perm<n + 1> edgeToTriangle = triangle->edgeMapping(facet);

        for (const regina::Triangle<n> *other : tri_->triangles()) {
            for (int i = 0; i < 3; ++i) {
                // Only glue identified edges and ignore identity mappings
                if (other->edge(i) != edge ||
                    (other->edge(i) == edge && triangle == other))
                    continue;

                regina::Perm<n + 1> edgeToOther = other->edgeMapping(i);
                // regina::Perm is immutable, use std::array instead
                std::array<int, 3> p;

                for (int i = 0; i < 3; ++i) {
                    p[edgeToTriangle[i]] = edgeToOther[i];
                }

                GluingGraph<n>::GluingNode node = {
                    triangle, facet, other, {p[0], p[1], p[2]}};
                // std::cout << "NEW NODE " << node;
                nodes.insert(node);
            }
        }
    }

    // std::cout << "\n";
    return nodes;
}

template <int n>
void GluingGraph<n>::buildGluingNodes_(const regina::Triangulation<n> *tri) {
    std::set<GluingGraph<n>::GluingNode> nodes;

    // TODO: Can we do better than O(|tri.triangles()|^2)?
    for (const regina::Triangle<n> *triangle : tri->triangles()) {
        std::set<GluingGraph<n>::GluingNode> newNodes = addTriangle_(triangle);

        for (const GluingGraph<n>::GluingNode &node : newNodes) {
            nodes.insert(node);
        }
    }

    for (const GluingGraph<n>::GluingNode &node : nodes) {
        nodes_.push_back(node);
    }
}

template <int n>
void GluingGraph<n>::buildGluingEdges_() {
    for (GluingGraph<n>::GluingNode &node : nodes_) {
        for (GluingGraph<n>::GluingNode &other : nodes_) {
            // std::cout << node << ", " << other;
            if (node.gluing_.dst == other.gluing_.src) {
                node.adjList_.insert(&other);
            }
        }
    }
}

template <int n>
void GluingGraph<n>::dfs_(std::set<SurfaceKnot> &surfaceKnots,
                          regina::Triangle<2> *srcTriangle,
                          GluingGraph<n>::GluingNode *node, int layer) {
    if (node->visited_ || !node->valid_) {
        // Don't perform the same gluing twice, and stop if adding this
        // gluing would result in an invalid triangulation
        return;
    }

    ++calls_;

    node->visited_ = true;

    regina::Triangle<2> *dstTriangle;

    auto search = triangleMap_.find(node->gluing_.dst->index());
    bool isNewTriangle = search == triangleMap_.end();

    if (isNewTriangle) {
        dstTriangle = surface_.newTriangle();
        triangleMap_.insert({node->gluing_.dst->index(), dstTriangle->index()});
    } else {
        dstTriangle = surface_.triangle(search->second);
    }

    // TODO: Is there some way to move this root check outside the
    // recursion? Probably doesn't matter
    if (srcTriangle != nullptr) {
        try {
            srcTriangle->join(node->gluing_.srcFacet, dstTriangle,
                              node->gluing_.gluing);
        } catch (regina::InvalidArgument &e) {
            if (isNewTriangle) {
                // If this is the first time we've seen this triangle and
                // the gluing is invalid, we're safe to remove it from
                // surface_
                surface_.removeTriangle(dstTriangle);
                // Note that find is guaranteed to return a valid iterator
                triangleMap_.erase(
                    triangleMap_.find(node->gluing_.dst->index()));
            }

            return;
        }
    }

    SurfaceKnot surfaceKnot(tri_, surface_, triangleMap_, cond_);
    if (surfaceKnot.isAdmissible()) {  // TODO: Don't check this here
        surfaceKnots.insert(surfaceKnot);
    }
    // nspaces_(layer);
    // std::cout << layer << ": adding " << surfaceKnot << "\n";

    // Recursive step: walk along each possible next gluing
    for (GluingGraph<n>::GluingNode *nextNode : node->adjList_) {
        dfs_(surfaceKnots, dstTriangle, nextNode, layer + 1);
    }

    // Remember to remove triangles as we move back up the recursion and
    // mark this node as not visited
    if (isNewTriangle) {
        surface_.removeTriangle(dstTriangle);
        triangleMap_.erase(triangleMap_.find(node->gluing_.dst->index()));
    }

    node->visited_ = false;
}

template <int n>
void GluingGraph<n>::reset_() {
    surface_ = {};
    triangleMap_.clear();

    for (GluingGraph<n>::GluingNode &node : nodes_) {
        node.valid_ = true;
        node.visited_ = false;
    }
}

template <int n>
std::set<SurfaceKnot> GluingGraph<n>::findSurfaceKnots_(
    const regina::Triangle<n> *triangle) {
    std::set<SurfaceKnot> ans;

    reset_();
    GluingGraph<n>::GluingNode root = {nullptr, -1, triangle, {}};
    for (auto &node : nodes_) {
        if (node.gluing_.src == triangle) {
            root.adjList_.insert(&node);
        }
    }
    // Call DFS with dummy root node
    dfs_(ans, nullptr, &root, 0);

    std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";

    // std::cout << "\n";

    return ans;
}

template <int n>
GluingGraph<n>::GluingGraph(const regina::Triangulation<n> *tri,
                            SurfaceCondition cond)
    : tri_(tri), cond_(cond) {
    buildGluingNodes_(tri);
    buildGluingEdges_();
    // buildGluingInvalids_();
}

template <int n>
std::set<SurfaceKnot> GluingGraph<n>::findSurfaceKnots() {
    // TODO: Account for surface condition in choosing the triangles we
    // iterate over
    std::set<SurfaceKnot> ans;

    for (const regina::Triangle<n> *triangle : tri_->triangles()) {
        auto surfaceKnots = findSurfaceKnots(triangle);
        for (const SurfaceKnot &surface : surfaceKnots) {
            ans.insert(surface);
        }
    }

    std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";

    return ans;
}
