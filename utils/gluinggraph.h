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
#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"

#include "knottedsurfaces.h"

template <int n>
class GluingGraph {
   private:
    class GluingNode {
       public:
        const Gluing<n> gluing_;
        std::unordered_set<GluingNode *> adjList_;
        std::unordered_set<GluingNode *> invalids_;  // Make sure to precompute
        bool valid_ = true;
        bool visited_ = false;

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

    const regina::Triangulation<n> &tri_;
    /**< The given triangulation of an n-manifold */

    std::vector<GluingNode> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    KnottedSurface surface_;

   public:
    GluingGraph(const regina::Triangulation<n> &tri, SurfaceCondition cond);

    std::set<KnottedSurface> findSurfaces();

   private:
    std::set<GluingNode> addTriangle_(const regina::Triangle<n> *triangle);

    void buildGluingNodes_(const regina::Triangulation<n> &tri);

    void buildGluingEdges_();

    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    void dfs_(GluingNode *node, std::set<KnottedSurface> &surfaces, int layer);

    void reset_();

    std::set<KnottedSurface> findSurfaces_(const regina::Triangle<n> *triangle);
};

#endif
