//
//  linkcomplement.cpp
//
//  Created by John Teague on 07/21/2026.
//

#include "linkcomplement.h"

#include <iostream>
#include <list>
#include <optional>

#include <triangulation/dim3/homologicaldata.h>
#include <snappea/snappeatriangulation.h>

std::mutex censusLookupMutex;
std::atomic<bool> simplifyComplements{true};

namespace {

// Memoizes recognition results (both recogniseHandlebody()'s genus and, for
// non-handlebody complements, Census::lookup()'s name) by isomorphism
// signature. The same boundary complement (e.g. a hyperbolic knot
// guaranteed to be a census hit) tends to recur across many found surfaces
// -- every real Census::lookup() reopens six on-disk census databases from
// scratch under censusLookupMutex, and recogniseHandlebody() itself is not
// free either, so caching turns "one recognition per surface" into "one
// recognition per distinct complement". Guarded by its own
// recognitionCacheMutex, separate from censusLookupMutex, so that a cache
// hit -- the common case once a search has been running a while, and the
// only case EdgeComplement::isUnknot()'s hot path ever takes -- never
// blocks behind a slow in-flight Census::lookup() on another thread.
std::mutex recognitionCacheMutex;
std::unordered_map<std::string, RecognitionResult> recognitionCache;
RecognitionCacheStats recognitionStats;

// Returns a snapshot of sig's current cache entry, or nullopt if unseen.
// By value, not by reference: another thread may concurrently complete
// this same entry (e.g. filling in the census fields after this snapshot
// was taken), so a reference into the map would be a data race on the
// struct's fields even though the map itself never erases (and hence never
// invalidates references to existing elements).
std::optional<RecognitionResult> lookupRecognition(const std::string &sig) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    auto it = recognitionCache.find(sig);
    if (it == recognitionCache.end())
        return std::nullopt;
    return it->second;
}

// Merges `update` into sig's entry monotonically -- genus, once computed,
// is a deterministic function of sig, so first-write-wins is safe there;
// censusChecked only ever moves false -> true. This is what keeps a thread
// racing to complete an entry from clobbering another thread's already-
// finished result. Returns the post-merge snapshot.
RecognitionResult storeRecognition(const std::string &sig,
                                    const RecognitionResult &update) {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    RecognitionResult &entry = recognitionCache[sig];
    if (update.genus && !entry.genus)
        entry.genus = update.genus;
    if (update.censusChecked && !entry.censusChecked) {
        entry.censusChecked = true;
        entry.censusName = update.censusName;
    }
    return entry;
}

// The actual (uncached) Census::lookup() call, serialized under
// censusLookupMutex per its documented contract. No memoization here --
// that's entirely recognitionCache's job now.
std::optional<std::string> censusLookupName(
    const regina::Triangulation<3> &complement) {
    std::lock_guard<std::mutex> lock(censusLookupMutex);
    std::list<regina::CensusHit> hits = regina::Census::lookup(complement);
    return hits.empty() ? std::nullopt
                        : std::make_optional(hits.front().name());
}

// Fast, sound, one-sided proof that `t`'s genus is 1 (the unknot): its
// fundamental group is Z if and only if it's the unknot (Dehn's lemma).
// group() already tries to simplify the presentation internally (same
// idiom as engine/triangulation/dim3/knot.cpp's Poincare-conjecture
// fast path), so a presentation of exactly one generator and no relations
// -- i.e. <a|>, which *is* Z by construction -- is conclusive. A "false"
// here is only ever inconclusive (simplify() didn't collapse it that far),
// never a wrong answer: this can only shorten the path to genus == 1, so
// it's always safe to fall back to recogniseHandlebody() when it fails.
bool groupProvesUnknot(const regina::Triangulation<3> &t) {
    const regina::GroupPresentation &g = t.group();
    return g.countGenerators() == 1 && g.countRelations() == 0;
}

// Fast, sound, one-sided proof that `t`'s genus is -1 (not any
// handlebody, of any genus): a genuine hyperbolic structure rules out
// every handlebody by geometrization (no handlebody -- solid tori
// included -- admits a complete hyperbolic structure, since its boundary
// is compressible). Symmetric to groupProvesUnknot() above: a "false"
// here is inconclusive, never wrong, so it's always safe to fall back to
// recogniseHandlebody(). SnapPeaTriangulation's construction and
// solutionType() are not documented thread-safe, but each caller here
// builds one from its own local, independently-owned Triangulation<3> --
// no state is shared across threads -- and this was stress-tested under
// concurrent load (8 threads, thousands of calls) with no crashes or
// inconsistent results.
bool hyperbolicityProvesNotHandlebody(const regina::Triangulation<3> &t) {
    regina::SnapPeaTriangulation snappea(t);
    return snappea.solutionType() ==
           regina::SnapPeaTriangulation::Solution::Geometric;
}

// Cached recogniseHandlebody(): never touches censusLookupMutex, so this is
// safe to call from EdgeComplement::isUnknot()'s hot, highly-parallel path
// without risking contention with an in-flight Census::lookup(). Tries the
// two fast, sound one-sided checks above before falling back to the
// expensive normal-surface-theory path -- both are cheap regardless of
// outcome, and either resolving conclusively skips recogniseHandlebody()
// entirely. Neither helps for knots that are neither the unknot nor
// hyperbolic (torus/satellite knots), which still fall through to
// recogniseHandlebody() same as before.
ssize_t cachedGenus(const regina::Triangulation<3> &complement,
                    const std::string &sig) {
    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        ++recognitionStats.genusChecks;
        auto it = recognitionCache.find(sig);
        if (it != recognitionCache.end() && it->second.genus) {
            ++recognitionStats.genusCacheHits;
            return *it->second.genus;
        }
    }

    ssize_t genus;
    enum class Path { Group, SnapPea, Fallback } path;
    if (groupProvesUnknot(complement)) {
        genus = 1;
        path = Path::Group;
    } else if (hyperbolicityProvesNotHandlebody(complement)) {
        genus = -1;
        path = Path::SnapPea;
    } else {
        genus = complement.recogniseHandlebody();
        path = Path::Fallback;
    }

    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        if (path == Path::Group)
            ++recognitionStats.groupFastPathHits;
        else if (path == Path::SnapPea)
            ++recognitionStats.snapPeaFastPathHits;
        else
            ++recognitionStats.recogniseHandlebodyFallbacks;
    }
    return *storeRecognition(sig, RecognitionResult{.genus = genus}).genus;
}

// Full resolution: genus, and (only if genus == -1) a census check.
RecognitionResult resolveRecognition(const regina::Triangulation<3> &complement,
                                     const std::string &sig) {
    ssize_t genus = cachedGenus(complement, sig);
    if (genus != -1)
        return *lookupRecognition(sig); // fully resolved; census never applies

    {
        std::lock_guard<std::mutex> lock(recognitionCacheMutex);
        ++recognitionStats.censusChecks;
        auto it = recognitionCache.find(sig);
        if (it != recognitionCache.end() && it->second.censusChecked) {
            ++recognitionStats.censusCacheHits;
            return it->second;
        }
    }

    auto name = censusLookupName(complement);
    return storeRecognition(
        sig, RecognitionResult{.genus = -1, .censusChecked = true,
                               .censusName = name});
}

} // namespace

RecognitionCacheStats recognitionCacheStats() {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    return recognitionStats;
}

size_t recognitionCacheSize() {
    std::lock_guard<std::mutex> lock(recognitionCacheMutex);
    return recognitionCache.size();
}

EdgeComplement::EdgeComplement(
    const regina::Triangulation<3> &tri,
    const std::vector<const regina::Edge<3> *> &edges)
    : tri_(&tri), edges_(edges) {
    // Add the edge to the correct tetrahedra
    for (const regina::Edge<3> *edge : edges) {
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            tetEdges_[emb.tetrahedron()].insert(emb.face());
        }
    }
}

EdgeComplement &EdgeComplement::operator=(const EdgeComplement &other) {
    if (this != &other) {
        tri_ = other.tri_;
        edges_ = other.edges_;
        tetEdges_ = other.tetEdges_;
    }
    return *this;
}

regina::Triangulation<3> EdgeComplement::buildComplement() const {
    regina::Triangulation<3> complement(*tri_);
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        complementTetEdges;

    for (const auto &[tet, edges] : tetEdges_) {
        complementTetEdges.emplace(complement.tetrahedron(tet->index()),
                                   edges);
    }

    while (!complementTetEdges.empty()) {
        auto &[tet, edges] = *complementTetEdges.begin();
        regina::Edge<3> *e = tet->edge(*edges.begin());

        for (const regina::EdgeEmbedding<3> &emb : e->embeddings()) {
            regina::Tetrahedron<3> *embTet = emb.tetrahedron();
            complementTetEdges[embTet].erase(emb.face());
            if (complementTetEdges[embTet].empty()) {
                complementTetEdges.erase(embTet);
            }
        }
        complement.pinchEdge(e);
    }

    // complement.idealToFinite();
    if (simplifyComplements.load(std::memory_order_relaxed))
        complement.simplify();

    return complement;
}

bool EdgeComplement::recognizeComplement() const {
    auto complement = buildComplement();
    std::string sig = complement.isoSig();
    RecognitionResult result = resolveRecognition(complement, sig);

    if (result.genus == 1) {
        std::cout << "      unknot, " << sig << "\n";
        return true;
    }
    if (result.censusName) {
        std::cout << "      recognized as " << *result.censusName << ", "
                  << sig << "\n";
        return true;
    }
    return false;
}

std::string EdgeComplement::identify() const {
    auto complement = buildComplement();
    std::string sig = complement.isoSig();
    RecognitionResult result = resolveRecognition(complement, sig);

    if (result.genus == 1)
        return "Unknot";
    if (result.censusName)
        return *result.censusName;
    return sig;
}

bool EdgeComplement::isUnknot() const {
    auto complement = buildComplement();
    return cachedGenus(complement, complement.isoSig()) == 1;
}

std::pair<regina::Triangulation<3>, std::vector<const regina::Edge<3> *>>
EdgeComplement::drillTrackingEdges_(
    const std::vector<const regina::Edge<3> *> &trackEdges) const {
    regina::Triangulation<3> complement(*tri_);

    // trackEdges is tracked via (Tetrahedron*, localEdge) descriptors,
    // which stay valid across pinchEdge() (pinchEdge() only ever appends
    // new tetrahedra, never removes or renumbers existing ones).
    std::vector<std::pair<regina::Tetrahedron<3> *, int>> trackedDescs;
    trackedDescs.reserve(trackEdges.size());
    for (const regina::Edge<3> *e : trackEdges) {
        auto emb = e->front();
        trackedDescs.emplace_back(
            complement.tetrahedron(emb.tetrahedron()->index()), emb.edge());
    }

    // Drill this object's own edges_ -- identical loop to
    // buildComplement(), just without simplify(): homology doesn't need
    // simplification, and simplifying would make tracking trackedDescs
    // through it intractable.
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        complementTetEdges;
    for (const auto &[tet, edges] : tetEdges_)
        complementTetEdges.emplace(complement.tetrahedron(tet->index()),
                                   edges);

    while (!complementTetEdges.empty()) {
        auto &[tet, edges] = *complementTetEdges.begin();
        regina::Edge<3> *e = tet->edge(*edges.begin());

        for (const regina::EdgeEmbedding<3> &emb : e->embeddings()) {
            regina::Tetrahedron<3> *embTet = emb.tetrahedron();
            complementTetEdges[embTet].erase(emb.face());
            if (complementTetEdges[embTet].empty()) {
                complementTetEdges.erase(embTet);
            }
        }
        complement.pinchEdge(e);
    }

    std::vector<const regina::Edge<3> *> tracked;
    tracked.reserve(trackedDescs.size());
    for (const auto &[tet, localEdge] : trackedDescs)
        tracked.push_back(tet->edge(localEdge));

    return {std::move(complement), std::move(tracked)};
}

long EdgeComplement::linkingNumberWith(const EdgeComplement &other) const {
    auto [complement, trackedOther] = drillTrackingEdges_(other.edges_);

    // Drilling curve A collapses it to a single ideal vertex -- a torus
    // cusp. H_1 must therefore be computed with that vertex *truncated*
    // (a genuine torus boundary), not treated as an ordinary 0-cell:
    // coning the cusp torus to a point kills every loop on it, including
    // A's meridian, which is precisely the class the linking number is
    // read off. A hand-rolled (vertices x edges, edges x triangles) pair
    // of boundary maps does exactly that wrong thing, and used to make
    // this routine return 0 for every link it was given.
    //
    // HomologicalData computes the correctly-truncated group and, unlike
    // Triangulation<3>::homology(), keeps it marked so a specific cycle
    // can be located inside it. Its standard cellular coordinates order
    // the 1-cells as "edges.begin() to edges.end(), followed by the ideal
    // edges of faces" (see its class documentation), so the triangulation's
    // own edges are a *prefix* of the coordinate space -- which is what
    // lets curve B's edge-indexed cycle vector below be used as-is, padded
    // with zeros over the trailing ideal-truncation cells.
    regina::HomologicalData hd(complement);
    const regina::MarkedAbelianGroup &h1 = hd.homology(1);

    // H_1 of a knot complement in S^3 is always Z, so this should now be
    // unreachable; kept as a defensive check rather than an assertion,
    // since reporting "no detected linking" stays sound (this routine
    // never falsely reports linking) if the precondition is ever violated.
    if (!h1.isZ())
        return 0;

    // Curve B must avoid the drilled cusp for the zero-padding described
    // above to be valid: an edge ending at the ideal vertex terminates on
    // one of the *new* truncation 0-cells instead of an ordinary vertex,
    // so such a cycle would no longer close up in these coordinates and
    // snfRep() would throw rather than return a wrong answer.
    // Geometrically this cannot happen -- an edge of B running into the
    // ideal vertex would mean B meets the drilled-out curve A,
    // contradicting this method's disjointness precondition -- so it is
    // checked, not assumed.
    for (const regina::Edge<3> *e : trackedOther)
        if (e->vertex(0)->isIdeal() || e->vertex(1)->isIdeal())
            throw regina::InvalidArgument(
                "EdgeComplement::linkingNumberWith(): other's curve meets "
                "the drilled ideal vertex, so the two curves are not "
                "disjoint");

    // Walk trackedOther's edges into an oriented cycle (the same
    // shared-vertex walk Link::Link() uses to split a multi-component edge
    // set into per-component knots), then read off its class in H_1. Only
    // the magnitude is meaningful, so no canonical orientation needs to be
    // imposed: any consistent walk direction works. The walk must be
    // *oriented*: an unsigned sum of B's edges is not a cycle in these
    // coordinates and snfRep() would reject it.
    regina::Vector<regina::Integer> cycle(hd.countStandardCells(1));
    if (!trackedOther.empty()) {
        std::vector<const regina::Edge<3> *> remaining(trackedOther.begin(),
                                                        trackedOther.end());
        const regina::Edge<3> *first = remaining.front();
        const regina::Vertex<3> *currVert = first->vertex(1);
        cycle[first->index()] += 1;
        remaining.erase(remaining.begin());

        while (!remaining.empty()) {
            bool found = false;
            for (size_t i = 0; i < remaining.size(); ++i) {
                const regina::Edge<3> *e = remaining[i];
                if (e->vertex(0) == currVert) {
                    cycle[e->index()] += 1;
                    currVert = e->vertex(1);
                } else if (e->vertex(1) == currVert) {
                    cycle[e->index()] -= 1;
                    currVert = e->vertex(0);
                } else {
                    continue;
                }
                remaining.erase(remaining.begin() +
                                static_cast<ptrdiff_t>(i));
                found = true;
                break;
            }
            if (!found)
                throw regina::InvalidArgument(
                    "EdgeComplement::linkingNumberWith(): other's edges do "
                    "not form a single closed curve");
        }
    }

    return h1.snfRep(cycle)[0].abs().safeValue<long>();
}

bool operator<(const EdgeComplement &e1, const EdgeComplement &e2) {
    return e1.edges_.size() < e2.edges_.size();
}

std::ostream &operator<<(std::ostream &os, const EdgeComplement &e) {
    os << "{";
    for (int i = 0; i < e.edges_.size(); ++i) {
        os << e.edges_[i]->index();
        if (i != e.edges_.size() - 1) {
            os << ", ";
        }
    }
    return os << "}";
}

Link::Link(const regina::Triangulation<3> &tri,
          const std::vector<const regina::Edge<3> *> &edges)
    : EdgeComplement(tri, edges) {
    // Add the edge to the correct component
    std::vector<std::vector<const regina::Edge<3> *>> edgesByComp;
    std::unordered_set<const regina::Edge<3> *> edgeSet(edges.begin(),
                                                        edges.end());
    if (edgeSet.size() != edges.size()) {
        throw regina::InvalidArgument(
            "Link::Link: Duplicate edges in link");
    }

    const regina::Vertex<3> *currVert;
    bool newComponent = true;
    while (!edgeSet.empty()) {
        if (newComponent) {
            const auto edge = *edgeSet.begin();
            edgesByComp.push_back({edge});
            currVert = edge->vertex(1);
            edgeSet.erase(edge);
        }

        newComponent = true;
        for (const regina::Edge<3> *edge : edgeSet) {
            if (edge->vertex(0) == currVert ||
                edge->vertex(1) == currVert) {
                edgesByComp.back().push_back(edge);
                currVert = edge->vertex(0) == currVert ? edge->vertex(1)
                                                       : edge->vertex(0);
                edgeSet.erase(edge);
                newComponent = false;
                break;
            }
        }
    }

    if (edgesByComp.size() == 1) {
        comps_.emplace_back(tri, edges);
        return;
    }

    for (const auto &compEdges : edgesByComp) {
        comps_.emplace_back(tri, compEdges);
    }
}

Link &Link::operator=(const Link &other) {
    if (this != &other) {
        EdgeComplement::operator=(other);
        comps_ = other.comps_;
    }
    return *this;
}

void Link::recognizeComplement() const {
    if (EdgeComplement::recognizeComplement()) {
        return;
    }

    // buildComplement()/simplify() runs again here (EdgeComplement::
    // recognizeComplement() above already built the same complement) --
    // a pre-existing inefficiency this cache doesn't address, since
    // simplify() has no isoSig to key off of until after it's run. The
    // cache does at least make this second recogniseHandlebody() free.
    auto complement = buildComplement();
    std::string sig = complement.isoSig();
    ssize_t genus = cachedGenus(complement, sig);
    int numComponents = countComponents();

    if (genus != -1 && numComponents == 1) {
        std::cout << "[!] WARNING! Recognized as a genus " << genus
                  << " handlebody, " << sig << "\n";
        std::cout << "[!] This is almost definitely a bug, please "
                     "report it!\n";
    } else if (numComponents == 1) {
        std::cout << "      NOT unknot, " << sig << "\n";
    } else if (numComponents > 1) {
        for (int i = 0; i < numComponents; ++i) {
            regina::Triangulation<3> compI = buildComplement(i);
            std::string sigI = compI.isoSig();
            RecognitionResult result = resolveRecognition(compI, sigI);

            std::cout << "    Component " << i + 1 << ": ";
            if (result.genus == 1) {
                std::cout << "unknot, " << sigI << "\n";
            } else if (result.genus && *result.genus != -1) {
                std::cout << "\n[!] WARNING! Recognized as a genus "
                          << *result.genus << " handlebody, ";
                std::cout << "\n[!] This is almost definitely a bug, "
                             "please "
                             "report it!\n";
            } else if (result.censusName) {
                std::cout << "      recognized as " << *result.censusName
                          << ", " << sigI << "\n";
            } else {
                std::cout << "NOT unknot, " << sigI << "\n";
            }
        }
    }
}

bool operator<(const Link &l1, const Link &l2) {
    if (l1.comps_.size() != l2.comps_.size())
        return l1.comps_.size() < l2.comps_.size();

    return static_cast<EdgeComplement>(l1) <
           static_cast<EdgeComplement>(l2);
}

std::ostream &operator<<(std::ostream &os, const Link &l) {
    os << "[";
    for (int i = 0; i < l.comps_.size(); ++i) {
        os << l.comps_[i];
        if (i != l.comps_.size() - 1) {
            os << ", ";
        }
    }
    return os << "]";
}
