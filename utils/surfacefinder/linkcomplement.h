#ifndef LINKCOMPLEMENT_H

#define LINKCOMPLEMENT_H

#include <mutex>
#include <string>
#include <unordered_map>
#include <unordered_set>

#include <census/census.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>

// regina::Census::lookup()'s Tokyo Cabinet-backed database is reopened from
// scratch, with no internal locking, on every single call (see
// CensusDB::lookupKey() in census-impl.h) -- so unlike most of Regina's
// engine, it cannot safely be called from more than one thread at a time.
// This holds regardless of whether the census databases have already been
// loaded once: concurrent calls fail unpredictably with "Could not open
// Tokyo Cabinet database ... -> threading error" (observed directly).
// Every regina::Census::lookup() call below serializes through this mutex;
// buildComplement()'s pinch + simplify work is unaffected and stays
// parallel across callers.
extern std::mutex censusLookupMutex;

/** Knot/Link implementation, specialized for taking complements */
class EdgeComplement {
  private:
    const regina::Triangulation<3> *tri_;
    std::vector<const regina::Edge<3> *> edges_;
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        tetEdges_;

  public:
    EdgeComplement(const regina::Triangulation<3> &tri,
                   const std::vector<const regina::Edge<3> *> &edges);

    EdgeComplement &operator=(const EdgeComplement &other);

    regina::Triangulation<3> buildComplement() const;

    // Checks the cheap, per-object recogniseHandlebody() genus test *before*
    // falling back to the serialized Census::lookup() below: census entries
    // are closed manifolds or genuinely knotted/hyperbolic complements, so a
    // genus-1 handlebody (an unknotted component's complement -- a solid
    // torus) is never actually going to be a census hit. Skipping straight
    // to "Unknot" for that (extremely common in practice) case avoids
    // contending on censusLookupMutex entirely -- profiling a real search
    // where every boundary happened to be an unknot showed >99% of wall
    // time in the old census-first ordering was spent blocked waiting on
    // that mutex, not doing useful work.
    bool recognizeComplement() const;

    // Non-printing counterpart to recognizeComplement(), for programmatic
    // use: the census name if Regina's census recognizes the complement,
    // "Unknot" if it's a genus-1 handlebody (the complement of an unknotted
    // component), or else the bare isoSig as a fallback identifier. See
    // recognizeComplement() above for why the genus check comes first.
    std::string identify() const;

    friend bool operator<(const EdgeComplement &e1, const EdgeComplement &e2);

    friend std::ostream &operator<<(std::ostream &os, const EdgeComplement &e);
};

class Knot : public EdgeComplement {
  public:
    Knot(const regina::Triangulation<3> &tri,
         const std::vector<const regina::Edge<3> *> &edges)
        : EdgeComplement(tri, edges) {}

    Knot &operator=(const Knot &other) {
        if (this != &other) {
            EdgeComplement::operator=(other);
        }
        return *this;
    }
};

class Link : public EdgeComplement {
  public:
    std::vector<Knot> comps_;

  public:
    Link(const regina::Triangulation<3> &tri,
         const std::vector<const regina::Edge<3> *> &edges);

    Link &operator=(const Link &other);

    regina::Triangulation<3> buildComplement(int component) const {
        return comps_[component].buildComplement();
    }

    regina::Triangulation<3> buildComplement() const {
        return EdgeComplement::buildComplement();
    }

    int countComponents() const { return comps_.size(); }

    void recognizeComplement() const;

    friend bool operator<(const Link &l1, const Link &l2);

    friend std::ostream &operator<<(std::ostream &os, const Link &l);
};

#endif // LINKCOMPLEMENT_H
