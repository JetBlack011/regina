//
//  gluinggraph.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef GLUING_GRAPH_H

#define GLUING_GRAPH_H

#include <gmpxx.h>
#include <link/link.h>
#include <set>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cstddef>
#include <cstring>

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

    std::vector<GluingNode<dim>> nodes_;
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

    std::set<KnottedSurface<dim>> &surfaces() { return surfaces_; }

    std::set<KnottedSurface<dim>> &findSurfaces() {
        surfaces_.clear();

        std::cout << "\n";
        for (auto &node : nodes_) {
            // Reset the surface and visited flags
            surface_ = {&tri_};
            for (auto &node : nodes_) {
                node.visited = false;
            }

            dfs_(&node, 0);
        }
        std::cout << "\n\n";
        std::cout << "[+] Total search calls = " << calls_ << "\n\n";

        return surfaces_;
    }

    std::set<KnottedSurface<dim>> &findSurfaces(
        std::unordered_set<regina::Triangle<dim> *> &startingTriangles) {
        surfaces_.clear();

        std::unordered_set<GluingNode<dim> *> startingNodes;
        // Find the nodes corresponding to the starting triangles (can be made
        // O(1) if we use a map, but this is fine for now)
        for (auto &node : nodes_) {
            if (startingTriangles.find(node.f) != startingTriangles.end()) {
                startingNodes.insert(&node);
            }
        }
        if (startingNodes.size() != startingTriangles.size()) {
            throw regina::InvalidArgument("GluingGraph::findSurfaces: Starting "
                                          "triangle(s) not found in the graph");
        }

        std::cout << "\n";
        for (auto &startingNode : startingNodes) {
            // Reset the surface and visited flags
            surface_ = {&tri_};
            for (auto &node : nodes_) {
                if (node != startingNode &&
                    startingNodes.find(&node) != startingNodes.end()) {
                    node.visited = true; // Don't revisit starting nodes
                } else {
                    node.visited = false;
                }
            }

            dfs_(&startingNode, 0);
        }
        std::cout << "\n\n";
        std::cout << "[+] Total search calls = " << calls_ << "\n\n";

        return surfaces_;
    }

    int countNodes() const { return nodes_.size(); }

    int countEdges() const {
        int count = 0;
        for (const auto &node : nodes_) {
            count += node.adjList.size();
        }
        return count;
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const GluingGraph<dim> &graph) {
        os << "GluingGraph<" << dim << ">(nodes = " << graph.countNodes()
           << ", edges = " << graph.countEdges() << ")\n";
        for (const auto &node : graph.nodes_) {
            os << "  Triangle: " << node.f->index() << "\n";
            for (const auto &[adjNode, gluing] : node.adjList) {
                os << "    Adjacent: " << adjNode->f->index() << ", " << gluing
                   << "\n";
            }
        }
        return os;
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
            nodes_.emplace_back(triangle);
        }
    }

    /**
     * Populate the node's adjacency lists.
     */
    void buildGluingEdges_() {
        for (auto &node : nodes_) {
            std::unordered_set<int> selfGluingEdges;
            for (int i = 0; i < 3; ++i) {
                regina::Triangle<dim> *triangle = node.f;
                regina::Edge<dim> *edge = triangle->edge(i);
                regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(i);

                // TODO: Would be nice for this to be O(1), not sure if possible
                for (auto &otherNode : nodes_) {
                    regina::Triangle<dim> *otherTriangle = otherNode.f;
                    for (int j = 0; j < 3; ++j) {
                        // Only glue identified edges and ignore identity
                        // mappings/opposite self gluings
                        if (otherTriangle->edge(j) != edge ||
                            (triangle == otherTriangle &&
                             (i == j || selfGluingEdges.find(i) !=
                                            selfGluingEdges.end())))
                            continue;

                        if (triangle == otherTriangle) {
                            selfGluingEdges.insert(j);
                        }

                        regina::Perm<dim + 1> edgeToOther =
                            otherTriangle->edgeMapping(j);
                        // regina::Perm is immutable, use std::array instead
                        std::array<int, 3> p;
                        for (int k = 0; k < 3; ++k) {
                            p[edgeToTriangle[k]] = edgeToOther[k];
                        }

                        node.adjList[&otherNode] = {triangle, i, otherTriangle,
                                                    p};
                    }
                }
            }
        }
    }

    int maxLayer_ = 0;

    void dfs_(GluingNode<dim> *node, int layer) {
        maxLayer_ = std::max(layer, maxLayer_);
        ++calls_;
        if (calls_ % 100000 == 0)
            std::cout << "\r" << std::flush << "[*] Call " << calls_
                      << ", Layer " << layer << ", Max layer = " << maxLayer_
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

        // Win condition: we've found a surface that meets the given boundary
        // requirement and has no self-intersection
        if ((cond_ == SurfaceCondition::all) ||
            (cond_ == SurfaceCondition::closed && surface_.isClosed()) ||
            (cond_ == SurfaceCondition::boundary && surface_.isProper())) {
            surfaces_.insert(surface_);
            // if (surface_.isClosed()) {
            //     std::cout << "[+] Found closed surface: " <<
            //     surface_.detail()
            //               << " with " << surface_.surface().size() << "\n";
            // } else if (cond_ == SurfaceCondition::boundary &&
            //            !surface_.isClosed()) {
            //     Link bdry = surface_.boundary();

            //    std::cout << surface_.detail() << " has " <<
            //    surface_.hasSelfIntersection() << " and  boundary " << bdry
            //              << "\n";
            //    bdry.recognizeComplement();
            //}
        }

        // Recursive step: search over all adjacent triangles
        for (auto &[nextNode, g] : node->adjList) {
            dfs_(nextNode, layer + 1);
        }

        // Don't forget to backtrack
        surface_.removeTriangle(node->f);
    }
};

#endif
