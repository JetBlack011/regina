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

#include <cstring>
#include <iostream>
#include <ostream>
#include <set>
#include <unordered_set>
#include <vector>

#include "maths/perm.h"
#include "surfaceknots.h"
#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"

template <int n>
class GluingGraph {
   private:
    struct Gluing {
        // Only included to make edge building easier
        const regina::Triangle<n> *src;
        const int srcFacet;
        const regina::Triangle<n> *dst;
        const regina::Perm<3> gluing;

        template <int>
        friend std::ostream &operator<<(std::ostream &os, const Gluing &gluing);
    };

    class GluingNode {
       public:
        const Gluing gluing_;
        std::unordered_set<GluingNode *> adjList_;
        std::unordered_set<GluingNode *> invalids_;  // Make sure to precompute
        bool valid_ = true;
        bool visited_ = false;

        // TODO: figure out C++
        GluingNode(const regina::Triangle<n> *src, int srcFacet,
                   const regina::Triangle<n> *dst, regina::Perm<3> gluing);

        /**
         * Comparison operators let us guarantee gluing uniqueness using
         * std::set, since I really don't want to come up with a hash function.
         */
        template <int>
        friend bool operator==(const GluingNode &lhs, const GluingNode &rhs);

        template <int>
        friend bool operator<(const GluingNode &lhs, const GluingNode &rhs);

        template <int>
        friend std::ostream &operator<<(std::ostream &os,
                                        const GluingNode &node);
    };

    size_t calls_ = 0;

    SurfaceCondition cond_;

    const regina::Triangulation<n> *tri_;
    /**< The given triangulation of an n-manifold */

    std::vector<GluingNode> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    regina::Triangulation<2> surface_;
    /**< To keep track of surfaces we built during DFS */

    TriangleMap triangleMap_;
    /**< A mapping between the triangles in the given triangulation and those in
     * the surface we build during DFS */

    std::set<GluingNode> addTriangle_(const regina::Triangle<n> *triangle);

    void buildGluingNodes_(const regina::Triangulation<n> *tri);

    void buildGluingEdges_();

    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    void dfs_(std::set<SurfaceKnot> &surfaceKnots,
              regina::Triangle<2> *srcTriangle, GluingNode *node, int layer);

    void reset_();

    std::set<SurfaceKnot> findSurfaceKnots_(
        const regina::Triangle<n> *triangle);

   public:
    GluingGraph(const regina::Triangulation<n> *tri, SurfaceCondition cond);

    std::set<SurfaceKnot> findSurfaceKnots();
};

#endif