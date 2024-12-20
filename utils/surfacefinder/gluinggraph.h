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
        // buildGluingInvalids_();
        // std::cout << "Built gluing invalids\n";
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
        if (!surfaces_.empty()) return surfaces_;

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

   private:
    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    void buildGluingNodes_() {
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            nodes_.insert({triangle, triangle});
        }
    }

    void buildGluingEdges_() {
        for (auto &[triangle, node] : nodes_) {
            std::unordered_set<int> selfGluingEdges;
            for (int i = 0; i < 3; ++i) {
                regina::Edge<dim> *edge = triangle->edge(i);
                regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(i);

                // for (const regina::TriangleEmbedding<dim> &emb :
                //      triangle->embeddings()) {
                //     // TODO: GAH
                //     const regina::Simplex<dim> *s = emb.simplex();
                //     for (int j = 0; j < 10; ++j) {
                //         regina::Triangle<dim> *other = s->triangle(j);
                //         GluingNode<dim> &otherNode = nodes_.at(other);

                //        for (int k = 0; k < 3; ++k) {
                //            // Only glue identified edges and ignore identity
                //            // mappings/opposite self gluings
                //            if (other->edge(k) != edge ||
                //                (triangle == other &&
                //                 (i == k || selfGluingEdges.find(i) !=
                //                                selfGluingEdges.end())))
                //                continue;

                //            if (triangle == other) {
                //                selfGluingEdges.insert(k);
                //            }

                //            regina::Perm<dim + 1> edgeToOther =
                //                other->edgeMapping(k);
                //            // regina::Perm is immutable, use std::array
                //            instead std::array<int, 3> p; for (int k = 0; k <
                //            3; ++k) {
                //                p[edgeToTriangle[k]] = edgeToOther[k];
                //            }

                //            node.adjList[&otherNode] = {triangle, i, other,
                //            p};
                //        }
                //    }
                //}

                for (auto &[other, otherNode] : nodes_) {
                    for (int j = 0; j < 3; ++j) {
                        // Only glue identified edges and ignore identity
                        // mappings/opposite self gluings
                        if (other->edge(j) != edge ||
                            (triangle == other &&
                             (i == j || selfGluingEdges.find(i) !=
                                            selfGluingEdges.end())))
                            continue;

                        if (triangle == other) {
                            selfGluingEdges.insert(j);
                        }

                        regina::Perm<dim + 1> edgeToOther =
                            other->edgeMapping(j);
                        // regina::Perm is immutable, use std::array instead
                        std::array<int, 3> p;
                        for (int k = 0; k < 3; ++k) {
                            p[edgeToTriangle[k]] = edgeToOther[k];
                        }

                        node.adjList[&otherNode] = {triangle, i, other, p};
                    }
                }
            }
        }
    }

    // void buildGluingInvalids_() {
    //     for (auto &[f, node] : nodes_) {
    //         std::unordered_set<const regina::Triangle<dim> *> adjTriangles;
    //         std::unordered_set<const regina::Vertex<dim> *> vertices;

    //        for (const auto [adj, _] : node.adjList) {
    //            adjTriangles.insert(adj->f);
    //        }

    //        for (int i = 0; i < 3; ++i) {
    //            vertices.insert(f->vertex(i));
    //        }

    //        // TODO: Optimize
    //        for (auto &[adjTriangle, adjNode] : nodes_) {

    //            if (adjTriangles.find(adjTriangle) != adjTriangles.end())
    //                continue;

    //            bool containsVertex = false;
    //            for (int i = 0; i < 3; ++i) {
    //                if (vertices.find(adjTriangle->vertex(i)) !=
    //                vertices.end())
    //                    containsVertex = true;
    //            }

    //            if (containsVertex) node.invalids.push_back(&adjNode);
    //        }
    //    }
    //}

    void dfs_(GluingNode<dim> *node, int layer) {
        ++calls_;
        if (calls_ % 100000 == 0)
            std::cout << "\r" << std::flush << "[*] Call " << calls_
                      << ", Layer " << layer << ", Surfaces "
                      << surfaces_.size() << "          ";

        if (node->visited) {
            return;
        }

        // nspaces_(layer);
        // std::cout << layer << ": ";
        if (!surface_.addTriangle(node->f, node->adjList)) {
            // std::cout << "Failed to add " << *node << "\n";
            return;
        }

        // std::cout << "Added " << *node << "\n";

        node->visited = true;

        if (surface_.isProper()) {
            //  nspaces_(layer);
            //  std::cout << "IMPROPER EDGES: { ";
            //  for (const regina::Edge<dim> *edge : surface_.improperEdges_) {
            //     std::cout << edge->index() << " ";
            // }
            //  std::cout << "}\n";
            surfaces_.insert(surface_);
        }

        for (auto &[nextNode, g] : node->adjList) {
            ////nspaces_(layer);
            // std::cout << layer << ": Gluing " << node->f->index() << " to "
            //           << nextNode->f->index() << " along (" << g.srcFacet
            //           << ", " << g.gluing[g.srcFacet] << ")\n";
            dfs_(nextNode, layer + 1);
        }

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
