//
//  linkcomplement.h
//
//  Created by John Teague on 07/16/2026.
//

#ifndef LINKCOMPLEMENT_H

#define LINKCOMPLEMENT_H

#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include <algebra/markedabeliangroup.h>
#include <census/census.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>

/*! \file utils/surfer/linkcomplement.h
 *  \brief Builds and recognizes the complement of a set of edges inside a
 *  triangulation.
 */

/**
 * A mutex guarding regina::Census::lookup() and the memoized cache of its
 * results (see linkcomplement.cpp): Regina's census databases are not safe
 * to query from more than one thread at a time -- concurrent calls have been
 * observed to crash outright (Tokyo Cabinet reports a "threading error"),
 * not just contend. Building or simplifying a complement is unaffected, and
 * stays parallel across callers.
 */
extern std::mutex censusLookupMutex;

/**
 * A set of edges in a Triangulation<3>, together with the ability to build
 * and recognize its complement.
 *
 * This is a general-purpose base for anything encoded as a list of edges
 * whose complement is of interest -- it does not by itself assume the
 * edges form a knot or link.
 */
class EdgeComplement {
  private:
    const regina::Triangulation<3> *tri_; /**< The ambient triangulation. */
    std::vector<const regina::Edge<3> *> edges_; /**< The tracked edges. */
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        tetEdges_; /**< edges_, indexed by incident tetrahedron and local edge number. */

  public:
    /**
     * Tracks `edges` inside `tri`.
     */
    EdgeComplement(const regina::Triangulation<3> &tri,
                   const std::vector<const regina::Edge<3> *> &edges);

    /** Sets this to a copy of `other`. */
    EdgeComplement &operator=(const EdgeComplement &other);

    /**
     * Builds the complement of the tracked edges: `tri_` with edges_
     * collapsed away and simplified.
     */
    regina::Triangulation<3> buildComplement() const;

    /**
     * Prints and returns whether the complement is recognized: either as
     * a genus-1 handlebody (the complement of a single unknotted
     * component), or as a census hit.
     */
    bool recognizeComplement() const;

    /**
     * Non-printing counterpart to recognizeComplement(): the census name
     * if recognized, "Unknot" if a genus-1 handlebody, or else the bare
     * isoSig as a fallback identifier.
     */
    std::string identify() const;

    /**
     * Cheap test for whether the complement is a genus-1 handlebody (a
     * solid torus) -- unlike recognizeComplement()/identify(), this never
     * falls back to the slower Census::lookup().
     *
     * \return \c true if and only if the complement is a genus-1
     * handlebody.
     */
    bool isUnknot() const;

    /**
     * The linking number of this edge set (curve A) with `other` (curve
     * B).
     *
     * \note Knot/link-specific, and arguably out of place on this
     * otherwise general-purpose base class.
     *
     * \pre This and `other` are disjoint simple closed curves in the same
     * ambient Triangulation<3>.
     *
     * \exception InvalidArgument The precondition was violated in a way
     * that is detectable here: either `other`'s edges do not form a single
     * closed curve, or they meet the ideal vertex left behind by drilling
     * this curve (which means the two curves are not disjoint).
     *
     * \return |lk(A,B)|; only the magnitude is meaningful, since no
     * canonical orientation is imposed on either curve.
     */
    long linkingNumberWith(const EdgeComplement &other) const;

    /** Orders by edge count, for use in ordered containers. */
    friend bool operator<(const EdgeComplement &e1, const EdgeComplement &e2);

    /** Writes the tracked edges' indices as `{i, j, ...}`. */
    friend std::ostream &operator<<(std::ostream &os, const EdgeComplement &e);

  private:
    /**
     * Clones `*tri_` and drills this object's own edges_, as
     * buildComplement() does but without simplify(), so `trackEdges`
     * stays trackable through it.
     *
     * \pre `trackEdges` are edges of the original `*tri_`, disjoint from
     * edges_.
     *
     * \return the drilled complement, together with `trackEdges`' images
     * in it (same order).
     */
    std::pair<regina::Triangulation<3>, std::vector<const regina::Edge<3> *>>
    drillTrackingEdges_(
        const std::vector<const regina::Edge<3> *> &trackEdges) const;
};

/**
 * An EdgeComplement whose edges are assumed to form a single knotted
 * component.
 */
class Knot : public EdgeComplement {
  public:
    /** Tracks the knot formed by `edges` inside `tri`. */
    Knot(const regina::Triangulation<3> &tri,
         const std::vector<const regina::Edge<3> *> &edges)
        : EdgeComplement(tri, edges) {}

    /** Sets this to a copy of `other`. */
    Knot &operator=(const Knot &other) {
        if (this != &other) {
            EdgeComplement::operator=(other);
        }
        return *this;
    }
};

/**
 * A (possibly multi-component) link: an EdgeComplement whose edges are
 * split into connected components on construction, each tracked as a
 * Knot.
 */
class Link : public EdgeComplement {
  public:
    std::vector<Knot> comps_; /**< This link's components. */

  public:
    /**
     * Splits `edges` into connected components and tracks each as a Knot.
     */
    Link(const regina::Triangulation<3> &tri,
         const std::vector<const regina::Edge<3> *> &edges);

    /** Sets this to a copy of `other`. */
    Link &operator=(const Link &other);

    /** Builds the complement of a single component. */
    regina::Triangulation<3> buildComplement(int component) const {
        return comps_[component].buildComplement();
    }

    /** Builds the complement of the whole link. */
    regina::Triangulation<3> buildComplement() const {
        return EdgeComplement::buildComplement();
    }

    /** Returns the number of components. */
    int countComponents() const { return comps_.size(); }

    /** Prints whether each component's complement is recognized; see EdgeComplement::recognizeComplement(). */
    void recognizeComplement() const;

    /** Orders by component count, then by edge count, for use in ordered containers. */
    friend bool operator<(const Link &l1, const Link &l2);

    /** Writes each component's tracked edges as `[{...}, {...}, ...]`. */
    friend std::ostream &operator<<(std::ostream &os, const Link &l);
};

#endif // LINKCOMPLEMENT_H
