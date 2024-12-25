//
//  gluinggraph.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef GLUING_GRAPH_H

#define GLUING_GRAPH_H

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cstddef>
#include <cstring>
#include <map>

#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"

#include "knottedsurfaces.h"

template <int dim>
class GluingGraph {
  private:
    size_t calls_ = 0;

    SurfaceCondition cond_;

    const regina::Triangulation<dim> &tri_;
    /**< The given triangulation of a dim-manifold */

    std::map<regina::Triangle<dim> *, GluingNode<dim>> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    std::set<KnottedSurface<dim>> surfaces_;

    KnottedSurface<dim> surface_;

  public:
    GluingGraph(const regina::Triangulation<dim> &tri, SurfaceCondition cond)
        : tri_(tri), surface_(&tri), cond_(cond) {
        buildGluingNodes_();
        std::cout << "[+] Built gluing nodes\n";
        buildGluingEdges_();
        std::cout << "[+] Built gluing edges\n";
    }

    void pruneSurfaces() {
        // TODO: Do this during DFS
        for (auto it = surfaces_.begin(); it != surfaces_.end();) {
            if ((*it).hasSelfIntersection()) {
                it = surfaces_.erase(it);
                continue;
            }
            ++it;
        }
    }

    std::set<KnottedSurface<dim>> &findSurfaces() {
        if (!surfaces_.empty())
            return surfaces_;

        std::cout << "\n";
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            reset_();
            dfs_(&nodes_.at(triangle), 0);
            // std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";
        }
        std::cout << "\n\n";
        std::cout << "[+] Total search calls = " << calls_ << "\n\n";

        pruneSurfaces();

        return surfaces_;
    }

    int countNodes() const { return nodes_.size(); }

    int countEdges() const {
        int count = 0;
        for (auto &[_, node] : nodes_) {
            count += node.adjList.size();
        }
        return count;
    }

  private:
    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    /**
     * Start building our gluing graph. Associate a node to each triangle in the
     * triangulation.
     */
    void buildGluingNodes_() {
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            nodes_.insert({triangle, triangle});
        }
    }

    /**
     * Populate the node's adjacency lists.
     */
    void buildGluingEdges_() {
        for (auto &[triangle, node] : nodes_) {
            std::unordered_set<int> selfGluingEdges;
            for (int i = 0; i < 3; ++i) {
                regina::Edge<dim> *edge = triangle->edge(i);
                regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(i);

                // TODO: Would be nice for this to be O(1), not sure if possible
                for (auto &[otherTriangle, otherNode] : nodes_) {
                    for (int j = 0; j < 3; ++j) {
                        // Only glue identified edges and ignore identity mappings/opposite self
                        // gluings
                        if (otherTriangle->edge(j) != edge ||
                            (triangle == otherTriangle &&
                             (i == j || selfGluingEdges.find(i) != selfGluingEdges.end())))
                            continue;

                        if (triangle == otherTriangle) {
                            selfGluingEdges.insert(j);
                        }

                        regina::Perm<dim + 1> edgeToOther = otherTriangle->edgeMapping(j);
                        // regina::Perm is immutable, use std::array instead
                        std::array<int, 3> p;
                        for (int k = 0; k < 3; ++k) {
                            p[edgeToTriangle[k]] = edgeToOther[k];
                        }

                        node.adjList[&otherNode] = {triangle, i, otherTriangle, p};
                    }
                }
            }
        }
    }

    void dfs_(GluingNode<dim> *node, int layer) {
        ++calls_;
        if (calls_ % 100000 == 0)
            std::cout << "\r" << std::flush << "[*] Call " << calls_ << ", Layer " << layer
                      << ", Surfaces " << surfaces_.size() << "          ";

        // Avoid cycles
        if (node->visited) {
            return;
        }

        // Attempt to add this node's face to the surface triangulation
        if (!surface_.addTriangle(node->f, node->adjList)) {
            return;
        }

        node->visited = true;

        // Win condition: we've found a surface that meets the given boundary requirement
        if (surface_.isProper()) {
            surfaces_.insert(surface_);
        }

        // Recursive step: search over all adjacent triangles
        for (auto &[nextNode, g] : node->adjList) {
            dfs_(nextNode, layer + 1);
        }

        // Don't forget to backtrack
        surface_.removeTriangle(node->f);
    }

    void reset_() {
        surface_ = {&tri_};

        for (auto &[_, node] : nodes_) {
            node.visited = false;
        }
    }
};

#endif
