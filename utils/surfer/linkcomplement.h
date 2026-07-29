//
//  linkcomplement.h
//
//  Created by John Teague on 07/16/2026.
//

#ifndef LINKCOMPLEMENT_H

#define LINKCOMPLEMENT_H

#include <atomic>
#include <functional>
#include <mutex>
#include <optional>
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
 * A mutex guarding only regina::Census::lookup() itself (see
 * linkcomplement.cpp): Regina's census databases are not safe to query from
 * more than one thread at a time -- concurrent calls have been observed to
 * crash outright (Tokyo Cabinet reports a "threading error"), not just
 * contend. Building or simplifying a complement is unaffected, and stays
 * parallel across callers. The memoized cache of recognition results is
 * guarded separately (see recognitionCacheMutex in linkcomplement.cpp), so
 * a cache hit never blocks behind an in-flight lookup on another thread.
 */
extern std::mutex censusLookupMutex;

/**
 * Whether EdgeComplement::buildComplement() should call
 * Triangulation<3>::simplify() on the drilled complement before returning
 * it. Defaults to \c true (existing behavior); set to \c false only for
 * profiling -- skipping simplify() yields much larger triangulations whose
 * isomorphism signatures rarely recur, so both the recognition cache below
 * and Census::lookup() itself become far less effective, in exchange for a
 * possibly much cheaper buildComplement().
 *
 * Set once, before any search worker thread is spawned, and never written
 * again -- std::thread's constructor already establishes happens-before to
 * every new thread, so a relaxed load here is safe.
 */
extern std::atomic<bool> simplifyComplements;

/**
 * The outcome of recognizing a complement, memoized by isomorphism
 * signature (see linkcomplement.cpp's recognitionCache).
 */
struct RecognitionResult {
    /**
     * Triangulation<3>::recogniseHandlebody()'s result, if it has been
     * computed for this isoSig: 1 means the unknot (a genus-1 handlebody),
     * -1 means "not any handlebody" (the only case eligible for a census
     * lookup), and any other value is the pre-existing "almost definitely a
     * bug" case (see Link::recognizeComplement()). nullopt means not yet
     * computed.
     */
    std::optional<ssize_t> genus;

    /**
     * Whether Census::lookup() has been run for this isoSig. Only ever
     * attempted when genus == -1.
     */
    bool censusChecked = false;

    /** The census hit's name, if censusChecked and a hit was found. */
    std::optional<std::string> censusName;
};

/**
 * Counters for how much recomputation the recognition cache is actually
 * avoiding. See EdgeComplement::recognitionCacheStats().
 */
struct RecognitionCacheStats {
    long long genusChecks = 0;
    long long genusCacheHits = 0;
    long long censusChecks = 0;
    long long censusCacheHits = 0;

    /**
     * On a genus cache miss, how the genus was actually resolved: via the
     * fast fundamental-group-is-Z check (proves genus == 1), the fast
     * SnapPea hyperbolicity check (proves genus == -1), or the full
     * recogniseHandlebody() fallback (neither fast check applied, e.g. a
     * torus/satellite knot). These three always sum to
     * genusChecks - genusCacheHits.
     */
    long long groupFastPathHits = 0;
    long long snapPeaFastPathHits = 0;
    long long recogniseHandlebodyFallbacks = 0;
};

/** A snapshot of the recognition cache's current hit/miss counters. */
RecognitionCacheStats recognitionCacheStats();

/** The recognition cache's current entry count (distinct isoSigs seen). */
size_t recognitionCacheSize();

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

    /**
     * The tracked edges' indices in the ambient triangulation, sorted.
     *
     * Cheap to compute (no triangulation/simplification involved) --
     * intended as a pre-triangulation cache key; see BoundarySignatureCache.
     */
    std::vector<size_t> edgeIndices() const;

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

/** Counters for how much recomputation BoundarySignatureCache is actually avoiding. */
struct BoundarySignatureCacheStats {
    long long checks = 0;
    long long hits = 0;
};

/**
 * Memoizes EdgeComplement::identify() results for marked edge sets in one
 * fixed ambient boundary-component triangulation, canonicalized against
 * that triangulation's own automorphism group -- so a boundary curve
 * already seen (exactly, or up to a symmetry of the boundary component)
 * short-circuits before buildComplement()/simplify()/isoSig() ever runs,
 * rather than only being caught by the isoSig-keyed recognitionCache
 * *after* paying for that triangulation (see linkcomplement.cpp).
 *
 * recognitionCache is still needed on top of this: two curves can be the
 * same knot/link type without being combinatorially related by any
 * automorphism of the boundary component (e.g. two non-isomorphic edge
 * paths that happen to drill out to the same manifold), in which case only
 * recognitionCache -- keyed on the post-simplify isoSig -- catches the
 * duplicate. This cache is a cheaper pre-filter in front of that one, not a
 * replacement for it.
 *
 * The ambient boundary-component triangulation this tracks must stay fixed
 * for this cache's lifetime, and must be the exact triangulation
 * KnottedSurface::boundaryLinks() draws its edges from (i.e. built the same
 * way, via BoundaryComponent::build()) -- see SurfaceSearch, which owns one
 * instance per boundary component, shared across every worker thread's
 * KnottedSurface, the same way it shares one PetalCache.
 *
 * Thread-safe: the (lazily built, then read-only) automorphism group is
 * guarded by a std::once_flag, and the memoization table by its own mutex
 * -- so a cache hit never blocks behind another thread still computing the
 * group, and computing the group happens exactly once regardless of how
 * many threads race to trigger it.
 */
class BoundarySignatureCache {
  public:
    /** Tracks `boundary`, which must outlive this cache. */
    explicit BoundarySignatureCache(const regina::Triangulation<3> &boundary);

    /**
     * Returns the identify() result for `edgeIndices` -- this boundary
     * component's marked edges, sorted, e.g. from
     * EdgeComplement::edgeIndices() -- computing it via `compute` on a
     * cache miss and memoizing the result. `compute` is only invoked on a
     * miss; typically `[&]{ return curve.identify(); }`.
     */
    std::string identifyCached(const std::vector<size_t> &edgeIndices,
                                const std::function<std::string()> &compute);

    /**
     * Returns a snapshot of this cache's current hit/miss counters; see
     * BoundarySignatureCacheStats. Returned by value (not by reference)
     * since another thread may be concurrently updating the live counters.
     */
    BoundarySignatureCacheStats stats() const;

    /** This cache's current entry count (distinct canonical signatures seen). */
    size_t size() const;

  private:
    const regina::Triangulation<3> *boundary_;

    std::once_flag groupOnce_;
    std::vector<std::vector<size_t>> automorphismEdgeImages_;
        /**< Lazily built by ensureAutomorphismGroup_():
             automorphismEdgeImages_[g][e] is the image of edge `e` of
             *boundary_ under the g-th automorphism found. Read-only once
             built, so safe to read from any thread without further
             synchronization (std::once_flag already establishes
             happens-before for every caller, not just the one that ran
             the initializer). */

    mutable std::mutex cacheMutex_; /**< Guards cache_ and stats_. */
    std::unordered_map<std::string, std::string> cache_;
    mutable BoundarySignatureCacheStats stats_;

    /** Builds automorphismEdgeImages_, exactly once, on first use. */
    void ensureAutomorphismGroup_();

    /**
     * Canonicalizes `edgeIndices` by minimizing its sorted image over
     * every automorphism in automorphismEdgeImages_ (triggering
     * ensureAutomorphismGroup_() first), returning the result as a
     * comma-joined string suitable as a cache_ key.
     */
    std::string canonicalKey_(const std::vector<size_t> &edgeIndices);
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
