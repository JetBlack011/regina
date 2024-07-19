//
//  embeddedsurfaces.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef SURFACE_KNOTS_H

#define SURFACE_KNOTS_H

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <iostream>
#include <ostream>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/triangulation.h"

/* Custom types */
using FaceMap = std::unordered_map<int, int>;
using EdgeMap = FaceMap;
using TriangleMap = FaceMap;

enum SurfaceCondition { closed, boundary, links };

template <int n>
using Curve = std::vector<regina::Edge<n> *>;

template <int n>
using DualCurve = std::vector<regina::Triangle<n> *>;

class Knot {
   private:
    using EdgeSet = std::unordered_set<const regina::Edge<3> *>;

    const regina::Triangulation<3> &tri_;
    Curve<3> edgeList_;
    EdgeSet edgeSet_;

    bool isUnknot_ = false;

    /**
     * Finds a push off of the given knot onto the dual 1-skeleton of the
     * triangulation it's sitting in. Note that this may change the underlying
     * triangulation.
     */
    DualCurve<3> knotToDualCurve_(const Knot &knot);

    std::pair<EdgeMap, regina::Triangulation<3>> drillDualCurve(
        const regina::Triangulation<3> &tri, const DualCurve<3> &curve);

   public:
    Knot(const regina::Triangulation<3> &tri, Curve<3> &edgeList);

    void simplify() {}

    regina::Triangulation<3> complement();
};

class Link {
   private:
    const regina::Triangulation<3> tri_;
    std::vector<Knot> components_;

   public:
    Link(const regina::Triangulation<3> &tri,
         std::vector<std::vector<regina::Edge<3> *>>);

    /**
     * Identifies a link in a given 3-manifold triangulation.
     *
     * Outline:
     * For each component of the link,
     * 1) Push off the given link into the dual 1-skeleton
     * 2) Drill the dual curve to obtain the link complement, keeping track of
     * where the other components of the link end up in the complement.
     * 3) Identify the link using the complement. SnapPy has a routine for this,
     * might end up running Python code in here.
     */
    regina::Triangulation<3> linkComplement(const regina::Triangulation<3> &t,
                                            const Link &l);
};

/* Useful classes */
class SurfaceKnot {
   private:
    const regina::Triangulation<4> *tri_;
    const regina::Triangulation<2> surface_;
    const TriangleMap triMap_;

    const SurfaceCondition cond_;
    const bool isOrientable_;
    const int punctures_;
    const int genus_;
    std::string detail_;

    friend std::pair<std::vector<int>, std::vector<int>> getTriIndexList(
        const SurfaceKnot &lhs, const SurfaceKnot &rhs);

   public:
    SurfaceKnot(const regina::Triangulation<4> *tri,
                const regina::Triangulation<2> &surface,
                const TriangleMap &triMap, SurfaceCondition cond);

    /**
     * Checks whether the boundary of the embedded surface is contained entirely
     * within the boundary of the triangulation it's embedded in.
     *
     * Basic idea: for each edge in the boundary of surface, find the
     * corresponding edge in tri and check if it's in the boundary.
     */
    bool isBoundaryContainment();

    const regina::Edge<4> *findEdgeInTri(const regina::Edge<2> *edge) const;

    /**
     * Based on the provided boundary condition, check whether we should keep
     * this surface or not
     */
    bool isAdmissible();

    const regina::Triangulation<2> &surface() const;

    bool isOrientable() const;

    int punctures() const;

    int genus() const;

    Link boundary();

    std::string detail() const;

    /**
     * Two embedded surfaces are considered equal if they use exactly the same
     * triangles in tri_ (the keys of a TriangleMap)
     */
    friend bool operator==(const SurfaceKnot &lhs, const SurfaceKnot &rhs);

    /**
     * Arbitrary tie breaker for sorting in a set by lexicographic ordering on
     * the key set's Keys. This may be O(nlog(n)), but, ideally, it happens
     * so infrequently relative to the recursive step of the program that it's
     * ok
     */
    friend bool operator<(const SurfaceKnot &lhs, const SurfaceKnot &rhs);

    friend std::ostream &operator<<(std::ostream &os,
                                    const SurfaceKnot &surfaceKnot);
};

#endif