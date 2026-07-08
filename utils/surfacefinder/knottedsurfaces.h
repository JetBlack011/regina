//
//  knottedsurfaces.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef KNOTTED_SURFACES_H

#define KNOTTED_SURFACES_H

#include <algorithm>
#include <memory>
#include <sstream>
#include <unordered_map>
#include <unordered_set>

#include <census/census.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>

#include "gluing.h"

enum class SurfaceCondition : std::uint8_t { all, boundary, closed };

/** Knot/Link implementation, specialized for taking complements */
class EdgeComplement {
  private:
    const regina::Triangulation<3> *tri_;
    std::vector<const regina::Edge<3> *> edges_;
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        tetEdges_;

  public:
    EdgeComplement(const regina::Triangulation<3> &tri,
                   const std::vector<const regina::Edge<3> *> &edges)
        : tri_(&tri), edges_(edges) {
        // Add the edge to the correct tetrahedra
        for (const regina::Edge<3> *edge : edges) {
            for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
                tetEdges_[emb.tetrahedron()].insert(emb.face());
            }
        }
    }

    EdgeComplement &operator=(const EdgeComplement &other) {
        if (this != &other) {
            tri_ = other.tri_;
            edges_ = other.edges_;
            tetEdges_ = other.tetEdges_;
        }
        return *this;
    }

    regina::Triangulation<3> buildComplement() const {
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
        complement.simplify();

        return complement;
    }

    bool recognizeComplement() const {
        auto complement = buildComplement();
        auto hits = regina::Census::lookup(complement);
        ssize_t genus = complement.recogniseHandlebody();

        if (!hits.empty()) {
            std::cout << "      recognized as " << hits.front().name() << ", "
                      << complement.isoSig() << "\n";
            return true;
        } else if (genus == 1) {
            std::cout << "      unknot, " << complement.isoSig() << "\n";
            return true;
        }
        return false;
    }

    friend bool operator<(const EdgeComplement &e1, const EdgeComplement &e2) {
        return e1.edges_.size() < e2.edges_.size();
    }

    friend std::ostream &operator<<(std::ostream &os, const EdgeComplement &e) {
        os << "{";
        for (int i = 0; i < e.edges_.size(); ++i) {
            os << e.edges_[i]->index();
            if (i != e.edges_.size() - 1) {
                os << ", ";
            }
        }
        return os << "}";
    }
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

    Link &operator=(const Link &other) {
        if (this != &other) {
            EdgeComplement::operator=(other);
            comps_ = other.comps_;
        }
        return *this;
    }

    regina::Triangulation<3> buildComplement(int component) const {
        return comps_[component].buildComplement();
    }

    regina::Triangulation<3> buildComplement() const {
        return EdgeComplement::buildComplement();
    }

    int countComponents() const { return comps_.size(); }

    void recognizeComplement() const {
        if (EdgeComplement::recognizeComplement()) {
            return;
        }

        auto complement = buildComplement();
        ssize_t genus = complement.recogniseHandlebody();
        int numComponents = countComponents();

        if (genus != -1 && numComponents == 1) {
            std::cout << "[!] WARNING! Recognized as a genus " << genus
                      << " handlebody, " << complement.isoSig() << "\n";
            std::cout << "[!] This is almost definitely a bug, please "
                         "report it!\n";
        } else if (numComponents == 1) {
            std::cout << "      NOT unknot, " << complement.isoSig() << "\n";
        } else if (numComponents > 1) {
            for (int i = 0; i < numComponents; ++i) {
                regina::Triangulation<3> complement = buildComplement(i);
                auto hits = regina::Census::lookup(complement);
                ssize_t genus = complement.recogniseHandlebody();

                std::cout << "    Component " << i + 1 << ": ";
                if (!hits.empty()) {
                    std::cout << "      recognized as " << hits.front().name()
                              << ", " << complement.isoSig() << "\n";
                } else if (genus == 1) {
                    std::cout << "unknot, " << complement.isoSig() << "\n";
                } else if (genus != -1) {
                    std::cout << "\n[!] WARNING! Recognized as a genus "
                              << genus << " handlebody, ";
                    std::cout << "\n[!] This is almost definitely a bug, "
                                 "please "
                                 "report it!\n";
                } else {
                    std::cout << "NOT unknot, " << complement.isoSig() << "\n";
                }
            }
        }
    }

    friend bool operator<(const Link &l1, const Link &l2) {
        if (l1.comps_.size() != l2.comps_.size())
            return l1.comps_.size() < l2.comps_.size();

        return static_cast<EdgeComplement>(l1) <
               static_cast<EdgeComplement>(l2);
    }

    friend std::ostream &operator<<(std::ostream &os, const Link &l) {
        os << "[";
        for (int i = 0; i < l.comps_.size(); ++i) {
            os << l.comps_[i];
            if (i != l.comps_.size() - 1) {
                os << ", ";
            }
        }
        return os << "]";
    }
};

/** KnottedSurface Implementation */

template <int dim>
class KnottedSurface {
  private:
    const regina::Triangulation<dim> *tri_;
    /**< The dim-manifold triangulation in which this subdim-manifold is
     * embedded */
    std::vector<std::shared_ptr<const regina::Triangulation<dim - 1>>>
        bdryComponents_;
    /**< bdryComponents_[c] is a standalone triangulation of
     * tri_->boundaryComponent(c), built once per component (in the very
     * first KnottedSurface constructed for a given tri_) and shared by
     * pointer across every subsequent copy -- these depend only on the
     * fixed ambient triangulation, never on anything the search mutates,
     * so there's no reason to rebuild them every time a surface gets
     * copied into the results set. A surface's boundary can touch more
     * than one of tri_'s boundary components (e.g. a cylinder built by
     * CobordismBuilder::thicken() naturally has a bottom and a top copy
     * of the base manifold), so this is a vector, not a single
     * triangulation -- see boundary(). */
    regina::Triangulation<2> surface_;
    /**< The triangulation of the sub-triangulation induced by the embedding
     */

    std::vector<const regina::Triangle<dim> *> emb_;
    /**< emb_[t->index()] is the tri_-triangle that surface_-triangle t maps
     * to. Grows/shrinks in lockstep with surface_ itself (push_back on
     * addTriangle, pop_back on removeTriangle) -- correct because
     * surface_'s own newSimplex()/removeSimplex() calls are always properly
     * nested (the DFS only ever removes the triangle it most recently
     * added), so surface_'s triangle indices behave exactly like a stack. */
    std::vector<regina::Triangle<2> *> inv_;
    /**< inv_[f->index()] is the surface_ triangle currently mapping to
     * tri_-triangle f, or nullptr if f isn't part of the surface right now.
     * Sized to tri_->countTriangles() once at construction, since tri_'s
     * triangle indices are fixed for the lifetime of the search -- this
     * replaces what used to be an unordered_map keyed on f, so every lookup
     * in the addTriangle/removeTriangle hot path is a plain array access
     * instead of a hash. */

    mutable std::array<int, 3> invariants_;
    mutable bool invariantsValid_ = false;
    /**< Regina's isOrientable()/countBoundaryComponents()/eulerChar() force
     * a full recompute of surface_'s skeleton (its cache is invalidated by
     * every addTriangle/removeTriangle mutation), so we no longer call them
     * eagerly on every DFS step -- only lazily, the first time detail() (or
     * any future invariant accessor) actually needs them. */

    std::set<int> indices_;

    // --- Incremental self-intersection tracking ---
    //
    // Rather than rescanning every vertex of the surface built so far on
    // every single addTriangle() attempt (an O(current surface size) cost
    // repeated at every DFS node), we maintain a rollback-capable union-find
    // over "raw corners" (a triangle together with one of its 3 local
    // vertices). Two raw corners end up in the same component exactly when
    // they're connected by a chain of edge-gluings among the triangles
    // currently in the surface -- i.e. one component per abstract vertex of
    // the surface being built, same as regina::Triangulation<2>'s own
    // vertex identification, just maintained incrementally by us instead of
    // recomputed from scratch after every mutation.
    //
    // imageRootCount_ tracks, for each ambient vertex, how many *distinct*
    // components currently map to it; selfIntersections_ counts how many
    // ambient vertices have more than one such component (i.e. two
    // genuinely different points of the abstract surface landing on the
    // same ambient point -- a self-intersection). A rollback union-find (no
    // path compression, union by size, undo log) is correct here because
    // DFS add/remove is properly nested (stack-like): whenever a triangle
    // is removed it's always the most recently added triangle that hasn't
    // already been removed, so undoing its union-find operations in LIFO
    // order always exactly reverses them.
    //
    // All three of ufParent_/ufSize_/imageRootCount_ are indexed by a fixed,
    // known-bounded key derived from tri_ (corner keys range over
    // [0, 3*tri_->countTriangles()), vertex indices over
    // [0, tri_->countVertices())), so -- like emb_/inv_ above -- plain
    // arrays sized once at construction replace what used to be
    // unordered_maps, eliminating hashing from the DFS hot path.
    std::vector<size_t> ufParent_;
    std::vector<int> ufSize_;
    std::vector<int> imageRootCount_;
    int selfIntersections_ = 0;

    struct UfUndoEntry {
        size_t childRoot;
        size_t parentRoot;
        int parentOldSize;
        const regina::Vertex<dim> *image;
    };
    std::vector<UfUndoEntry> ufUndoLog_;
    // Per-triangle (indexed by f->index()) snapshot of ufUndoLog_.size()
    // taken just before that triangle's unions were performed, so
    // removeTriangle() knows exactly how many log entries to undo.
    std::vector<size_t> ufCheckpoints_;

    // Number of (triangle, local facet) pairs in the surface not currently
    // glued to a neighbour -- lets isClosed() be O(1) instead of going
    // through regina::Triangulation<2>'s (skeleton-dependent) isClosed().
    int numBoundaryFacets_ = 0;

    // --- Unresolvable-conflict tracking (for early DFS pruning) ---
    //
    // hasSelfIntersection() (above) answers "is there a self-intersection
    // *right now*" -- true the instant two different abstract-surface
    // points land on the same ambient vertex, even if a still-undecided
    // triangle would go on to merge them into one point later (a genuine
    // branch point: 3+ triangles meeting at one ambient vertex without
    // being pairwise edge-adjacent yet). Rejecting on that signal alone is
    // what silently dropped valid surfaces with real branch points (see
    // the investigation this supports). hasUnresolvableConflict() answers
    // the *stronger, tighter* question actually needed for early pruning:
    // is there a conflict that can *never* be resolved by anything still
    // possible? It's still sound to reject on (never lets a truly-bad
    // surface through undetected -- checkWin_'s unconditional
    // hasSelfIntersection() call remains the final gate regardless), but
    // strictly weaker than hasSelfIntersection(), so it lets the DFS carry
    // a transient, still-resolvable conflict forward instead of discarding
    // it immediately.
    //
    // A conflict at ambient vertex v is judged permanent only once at
    // least two of v's currently-distinct corner-components (union-find
    // roots) are individually "closed" -- meaning every corner in that
    // component has exhausted every facet that could ever grow or merge
    // it further. Two closed components can never become one, so if two
    // exist simultaneously the conflict there truly can't be undone by
    // anything that happens later. This is deliberately scoped to just
    // the corners actually involved in a live conflict (checked only when
    // selfIntersections_ > 0, i.e. essentially never on the overwhelmingly
    // common conflict-free path) rather than the whole ambient vertex's
    // total degree in the *entire* triangulation -- the latter (tried and
    // reverted before this) is sound but can take arbitrarily long to
    // converge on a high-degree vertex, regardless of how small the actual
    // conflict is, which is what caused a 76,000x DFS-call blowup on one
    // of the existing test cases.
    //
    // "Exhausted every facet that could ever grow or merge it further"
    // needs two per-triangle-facet pieces of information:
    //  - excluded_[i]: has triangle i been permanently decided *against*
    //    for the rest of the current DFS path? (Toggled by
    //    decideExclude()/undecideExclude(), called from extend_'s
    //    exclude-branch -- see gluing.h's `queued` field for why this is
    //    now safe to treat as truly permanent within its scope.) A
    //    triangle nobody has decided about yet (never entered the
    //    frontier, or entered but not yet reached) is *not* excluded --
    //    it might still show up and resolve something, so treating
    //    "undecided" as "still a candidate" is the sound default.
    //  - candidatesByFacet_[i][j]: the ambient triangle indices that could
    //    ever be glued to triangle i's local facet j, i.e. every neighbour
    //    addTriangle() saw in `adj` with that srcFacet, recorded as plain
    //    indices (not pointers into the caller's AdjList -- some callers,
    //    e.g. direct KnottedSurface tests, pass a temporary AdjList that
    //    doesn't outlive the call, so a stored pointer would dangle).
    //    Recomputed fresh every addTriangle() call for that triangle,
    //    since it's cheap (same order as the union loop right below it)
    //    and this way it can never go stale between a triangle's removal
    //    and re-addition.
    std::vector<bool> excluded_;
    std::vector<std::array<std::vector<size_t>, 3>> candidatesByFacet_;

  public:
    std::unordered_set<const regina::Edge<dim> *> improperEdges_;
    /**< Notes whether an embedding is proper, in the sense that the image
     * of the boundary of the sub-triangulation is entirely contained within
     * the boundary of the dim-manifold triangulation */

    KnottedSurface(const regina::Triangulation<dim> *tri)
        : tri_(tri), inv_(tri_->countTriangles(), nullptr),
          ufParent_(3 * tri_->countTriangles()),
          ufSize_(3 * tri_->countTriangles()),
          imageRootCount_(tri_->countVertices(), 0),
          ufCheckpoints_(tri_->countTriangles()),
          excluded_(tri_->countTriangles(), false),
          candidatesByFacet_(tri_->countTriangles()) {
        if (!tri_->isClosed()) {
            bdryComponents_.reserve(tri_->countBoundaryComponents());
            for (size_t c = 0; c < tri_->countBoundaryComponents(); ++c) {
                bdryComponents_.push_back(
                    std::make_shared<const regina::Triangulation<dim - 1>>(
                        tri_->boundaryComponent(c)->build()));
            }
        }
    }

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_), bdryComponents_(other.bdryComponents_),
          surface_(other.surface_), emb_(other.emb_),
          inv_(other.inv_.size(), nullptr), invariants_(other.invariants_),
          invariantsValid_(other.invariantsValid_), indices_(other.indices_),
          ufParent_(other.ufParent_), ufSize_(other.ufSize_),
          imageRootCount_(other.imageRootCount_),
          selfIntersections_(other.selfIntersections_),
          ufUndoLog_(other.ufUndoLog_), ufCheckpoints_(other.ufCheckpoints_),
          numBoundaryFacets_(other.numBoundaryFacets_),
          excluded_(other.excluded_),
          candidatesByFacet_(other.candidatesByFacet_),
          improperEdges_(other.improperEdges_) {
        // emb_ was blitted straight from other (its values point into the
        // shared, immutable tri_, so they need no adjustment), but inv_'s
        // values point into surface_ itself -- a fresh, distinct copy of
        // other.surface_ -- so those pointers have to be rebuilt from our
        // own surface_ rather than copied from other.inv_.
        for (regina::Triangle<2> *t : surface_.triangles()) {
            inv_[emb_[t->index()]->index()] = t;
        }
    }

    KnottedSurface &operator=(const KnottedSurface &other) {
        if (this != &other) {
            tri_ = other.tri_;
            bdryComponents_ = other.bdryComponents_;
            surface_ = other.surface_;
            emb_ = other.emb_;
            improperEdges_ = other.improperEdges_;
            invariants_ = other.invariants_;
            invariantsValid_ = other.invariantsValid_;
            indices_ = other.indices_;
            ufParent_ = other.ufParent_;
            ufSize_ = other.ufSize_;
            imageRootCount_ = other.imageRootCount_;
            selfIntersections_ = other.selfIntersections_;
            ufUndoLog_ = other.ufUndoLog_;
            ufCheckpoints_ = other.ufCheckpoints_;
            numBoundaryFacets_ = other.numBoundaryFacets_;
            excluded_ = other.excluded_;
            candidatesByFacet_ = other.candidatesByFacet_;

            // See the copy constructor: inv_'s values are per-instance
            // surface_ pointers, so they must be rebuilt, not blitted.
            inv_.assign(other.inv_.size(), nullptr);
            for (regina::Triangle<2> *t : surface_.triangles()) {
                inv_[emb_[t->index()]->index()] = t;
            }
        }

        return *this;
    }

    const regina::Triangulation<2> &surface() const { return surface_; }

    bool isClosed() const { return numBoundaryFacets_ == 0; }

    bool isProper() const { return isClosed() || improperEdges_.empty(); }

    bool hasSelfIntersection() const { return selfIntersections_ > 0; }

    // Records that triangle w has been permanently decided *against* for
    // the rest of the current DFS path -- i.e. nothing will ever glue to
    // it again down this branch. Called from extend_'s exclude-branch,
    // matched by a later undecideExclude() when that branch backtracks.
    // Only safe to treat as truly permanent because extend_ now also
    // deduplicates the frontier (see gluing.h's `queued` field) -- without
    // that, the same candidate could still be re-pushed and included
    // later in the same subtree, which would make this bookkeeping wrong.
    void decideExclude(size_t triangleIndex) {
        excluded_[triangleIndex] = true;
    }

    void undecideExclude(size_t triangleIndex) {
        excluded_[triangleIndex] = false;
    }

    // See the class-level comment above excluded_/candidatesByFacet_ for
    // the full rationale. Cheap (O(1)) on the overwhelmingly common path
    // where there's no live conflict at all; the more expensive scan below
    // only runs while selfIntersections_ > 0.
    bool hasUnresolvableConflict() const {
        if (selfIntersections_ == 0)
            return false;

        // Group currently-mapped-to-v corners by (vertex, union-find
        // root), but only for vertices that actually have 2+ roots right
        // now -- everywhere else there's no conflict to resolve.
        std::unordered_map<
            size_t,
            std::unordered_map<
                size_t,
                std::vector<std::pair<const regina::Triangle<dim> *, int>>>>
            byVertex;

        for (int fi : indices_) {
            const regina::Triangle<dim> *f = tri_->triangle(fi);
            for (int c = 0; c < 3; ++c) {
                size_t vIdx = f->vertex(c)->index();
                if (imageRootCount_[vIdx] < 2)
                    continue;
                size_t root = ufFind_(cornerKey_(f, c));
                byVertex[vIdx][root].emplace_back(f, c);
            }
        }

        for (const auto &[vIdx, roots] : byVertex) {
            if (roots.size() < 2)
                continue;

            int closedCount = 0;
            for (const auto &[root, members] : roots) {
                if (isRootClosed_(members))
                    ++closedCount;
            }
            if (closedCount >= 2)
                return true;
        }
        return false;
    }

    // Returns one Link per tri_ boundary component that this surface's
    // boundary actually touches, paired with that component's index (see
    // bdryComponents_'s class-level comment for why this isn't just a
    // single Link). A connected boundary circle of the abstract surface
    // is a continuous image, so it always lands entirely within one
    // ambient boundary component -- it can never straddle two -- which is
    // what makes attributing edges to a single component well-defined
    // here.
    std::vector<std::pair<size_t, Link>> boundary() const {
        std::vector<std::pair<size_t, Link>> result;

        // TODO: Make this O(1) somehow
        for (size_t c = 0; c < bdryComponents_.size(); ++c) {
            const regina::BoundaryComponent<dim> *triBoundaryComp =
                tri_->boundaryComponent(c);
            std::unordered_set<const regina::Edge<dim - 1> *> edges;

            for (const regina::BoundaryComponent<2> *surfaceBoundaryComp :
                 surface_.boundaryComponents()) {
                for (const regina::Edge<2> *e : surfaceBoundaryComp->edges()) {
                    const regina::Edge<dim> *im = image(e);
                    for (int i = 0; i < triBoundaryComp->countEdges(); ++i) {
                        if (im == triBoundaryComp->edge(i)) {
                            edges.insert(bdryComponents_[c]->edge(i));
                            break;
                        }
                    }
                }
            }

            if (edges.empty())
                continue;

            std::vector<const regina::Edge<dim - 1> *> edgeList(edges.begin(),
                                                                edges.end());
            result.emplace_back(c, Link(*bdryComponents_[c], edgeList));
        }

        return result;
    }

    template <int facedim>
    regina::Face<dim, facedim> *image(const regina::Face<2, facedim> *f) const {
        static_assert(0 <= facedim && facedim <= 2,
                      "Must have 0 <= facedim <= 2");

        const regina::FaceEmbedding<2, facedim> &fEmb = f->front();
        return emb_[fEmb.simplex()->index()]->template face<facedim>(
            fEmb.face());
    }

    const regina::Triangle<dim> *image(regina::Triangle<2> *t) const {
        return emb_[t->index()];
    }

    regina::Triangle<2> *preimage(const regina::Triangle<dim> *f) const {
        return inv_[f->index()];
    }

    // Clears all per-search mutable state, as if freshly constructed --
    // but leaves tri_/bdryComponents_ untouched, since those depend only
    // on the fixed ambient triangulation. Lets a search reuse a single
    // KnottedSurface across every seed without re-deriving (and
    // re-copying) every boundary component's triangulation each time,
    // which `surface_ = {&tri_};`-style reassignment would otherwise do.
    void reset() {
        surface_ = regina::Triangulation<2>();
        emb_.clear();
        inv_.assign(tri_->countTriangles(), nullptr);
        invariantsValid_ = false;
        indices_.clear();
        std::fill(ufParent_.begin(), ufParent_.end(), size_t{0});
        std::fill(ufSize_.begin(), ufSize_.end(), 0);
        std::fill(imageRootCount_.begin(), imageRootCount_.end(), 0);
        selfIntersections_ = 0;
        ufUndoLog_.clear();
        std::fill(ufCheckpoints_.begin(), ufCheckpoints_.end(), size_t{0});
        numBoundaryFacets_ = 0;
        std::fill(excluded_.begin(), excluded_.end(), false);
        // candidatesByFacet_ is purely structural (derived from `adj`,
        // never from DFS state) and gets fully overwritten the next time
        // each triangle is added, so it doesn't need resetting here.
        improperEdges_.clear();
    }

    /** Mutation methods */

    // checkSelfIntersection=false defers the self-intersection check (still
    // performed, just not acted on here -- see hasSelfIntersection()) to
    // the caller. This matters for building a fixed multi-triangle base
    // (as opposed to the single-triangle-at-a-time DFS search, which wants
    // the immediate check for early pruning): a base with an internal
    // branch point can need triangle A and triangle B to both individually
    // touch some ambient vertex v before the *connecting* triangle between
    // them has been added, even though the full base is perfectly valid
    // once everything is in. Checking (and rejecting) after every single
    // insertion can therefore reject a valid base purely because of
    // insertion order, when what actually matters is whether a
    // self-intersection remains once the whole base is present.
    bool addTriangle(const regina::Triangle<dim> *f,
                     const typename GluingNode<dim>::AdjList &adj,
                     bool checkSelfIntersection = true) {
        // TODO: Figure out a nice way for this method not to know what an
        // AdjList is
        regina::Triangle<2> *src = surface_.newSimplex();
        std::array<bool, 3> isBoundaryEdge = {true, true, true};

        // src is always the newest surface_ triangle (index == emb_.size()
        // before this push), so emb_ can just grow in lockstep instead of
        // being keyed by src's pointer.
        emb_.push_back(f);
        inv_[f->index()] = src;

        // Structural (doesn't depend on DFS state): which ambient triangles
        // could ever be glued to each of f's 3 facets. See
        // candidatesByFacet_'s class-level comment for why this is
        // recomputed from `adj` here rather than just storing a pointer to
        // it.
        auto &byFacet = candidatesByFacet_[f->index()];
        for (auto &slot : byFacet)
            slot.clear();
        for (const auto &[adjNode, gluings] : adj) {
            for (const auto &g : gluings) {
                byFacet[g.srcFacet].push_back(adjNode->f->index());
            }
        }

        // Register f's 3 raw corners as fresh singleton union-find
        // components before attempting any joins, so ufUnite_() below always
        // has something to look up. ufCheckpoint marks where in the undo log
        // this triangle's own unions start, for removeTriangle() to roll
        // back to later.
        size_t ufCheckpoint = ufUndoLog_.size();
        for (int i = 0; i < 3; ++i) {
            size_t key = cornerKey_(f, i);
            ufParent_[key] = key;
            ufSize_[key] = 1;
            addCornerImage_(f->vertex(i));
        }

        int boundaryFacetsConsumed = 0;
        // Of boundaryFacetsConsumed, how many were glued to a genuinely
        // different, already-present neighbour triangle (as opposed to a
        // self-fold, where both "sides" are f's own brand-new facets, and
        // there's no pre-existing neighbour slot being flipped) -- see the
        // comment above numBoundaryFacets_'s update below for why this
        // needs to be tracked separately from boundaryFacetsConsumed.
        int externalFacetsConsumed = 0;

        // adjNode may be reachable via more than one gluing (two triangles
        // can share more than one edge -- see GluingNode::AdjList's class
        // comment), so every gluing in the vector needs to be joined and
        // unioned, not just one.
        for (const auto &[adjNode, gluings] : adj) {
            regina::Triangle<2> *dst = inv_[adjNode->f->index()];
            if (dst == nullptr)
                continue;
            bool selfFold = (dst == src);

            for (const auto &g : gluings) {
                int dstFacet = g.gluing[g.srcFacet];

                // A self-folded triangle (two of f's own facets identified
                // to the same ambient edge) gets *two* adjList entries for
                // this one physical gluing -- one per direction (see
                // buildGluingEdges_'s comment). Whichever is processed
                // first joins normally below; Regina's own join() then
                // sets *both* facets' adjacentSimplex reciprocally (dst
                // and src are the same triangle here), so the second
                // entry arrives seeing its own facet already glued -- that
                // is not a genuine double-glued-facet conflict, just this
                // gluing's other direction catching up. Its facet is
                // genuinely interior now, so isBoundaryEdge/
                // boundaryFacetsConsumed still need updating, but the join
                // and vertex-identification unions were already done from
                // the first direction (the gluing permutation is a full
                // bijection, so one direction's ufUnite_ calls already
                // cover every identification this gluing produces).
                if (selfFold && src->adjacentSimplex(g.srcFacet) != nullptr) {
                    isBoundaryEdge[g.srcFacet] = false;
                    ++boundaryFacetsConsumed;
                    continue;
                }

                // Saves us from catching an error
                if (dst->adjacentSimplex(dstFacet) != nullptr ||
                    src->adjacentSimplex(g.srcFacet) != nullptr) {
                    undoFailedAdd_(f, ufCheckpoint, src);
                    return false;
                }

                isBoundaryEdge[g.srcFacet] = false;
                ++boundaryFacetsConsumed;
                if (!selfFold)
                    ++externalFacetsConsumed;
                src->join(g.srcFacet, dst, g.gluing);

                // The gluing perm maps src's local vertex numbering to dst's;
                // the two vertices on the shared edge (everything but the
                // facet itself) get identified as the same abstract surface
                // vertex.
                for (int corner = 0; corner < 3; ++corner) {
                    if (corner == g.srcFacet)
                        continue;
                    ufUnite_(cornerKey_(f, corner),
                             cornerKey_(adjNode->f, g.gluing[corner]),
                             f->vertex(corner));
                }
            }
        }

        // f needs to be in indices_ *before* hasUnresolvableConflict() runs
        // (moved up from below, where it used to sit after this check):
        // that scan walks indices_ to find every corner currently mapping
        // to a conflicted vertex, and f's own corners are exactly what a
        // conflict just created -- leaving f out would make its brand-new
        // root at that vertex invisible to the closedness check.
        indices_.insert(f->index());

        // hasUnresolvableConflict(), not hasSelfIntersection(): see the
        // class-level comment above excluded_/candidatesByFacet_. The main
        // per-triangle DFS (extend_ in surfacefinder.h) is the only caller
        // that reaches here with checkSelfIntersection=true; checkWin_
        // still gates every recorded surface on the unconditional, real
        // hasSelfIntersection() regardless of what happens here.
        if (checkSelfIntersection && hasUnresolvableConflict()) {
            undoFailedAdd_(f, ufCheckpoint, src);
            return false;
        }

        for (int i = 0; i < 3; ++i) {
            if (!isBoundaryEdge[i]) {
                improperEdges_.erase(f->edge(i));
            } else if (isBoundaryEdge[i] && !f->edge(i)->isBoundary()) {
                improperEdges_.insert(f->edge(i));
            }
        }

        // Every consumed facet removes one boundary slot on our own side
        // (it was never counted as boundary to begin with, since f is
        // brand new); externally-consumed facets *additionally* flip one
        // slot on the neighbour's side from boundary to interior, since
        // that slot already existed and was already counted. A self-fold
        // has no such neighbour -- both "sides" are f's own new facets --
        // so it only accounts for the first "-1", not a second one (see
        // the comment above externalFacetsConsumed's declaration). The
        // remaining (3 - boundaryFacetsConsumed) facets are new boundary
        // slots of our own.
        numBoundaryFacets_ +=
            3 - boundaryFacetsConsumed - externalFacetsConsumed;

        ufCheckpoints_[f->index()] = ufCheckpoint;
        invariantsValid_ = false;

        return true;
    }

    void removeTriangle(const regina::Triangle<dim> *f) {
        regina::Triangle<2> *src = inv_[f->index()];

        for (int i = 0; i < 3; ++i) {
            regina::Triangle<2> *dst = src->adjacentSimplex(i);

            if (dst == src) {
                // A self-folded facet: both sides of this gluing are f's
                // own facets, and f is being removed entirely, so unlike a
                // genuinely external neighbour there's no *other*
                // triangle left standing to expose a newly-un-glued
                // facet on -- this pair should net to zero, matching
                // addTriangle's own numBoundaryFacets_ accounting for a
                // self-fold (see externalFacetsConsumed's comment there).
                // Never flagged improper while f was present (addTriangle
                // never inserts a self-fold's edge -- see its own
                // isBoundaryEdge loop), so nothing to erase there either,
                // but do it anyway for symmetry/defensiveness.
                regina::Edge<dim> *edge = f->edge(i);
                if (!edge->isBoundary()) {
                    improperEdges_.erase(edge);
                }
            } else if (dst != nullptr) {
                regina::Edge<dim> *edge = f->edge(i);
                if (!edge->isBoundary()) {
                    improperEdges_.insert(edge);
                }
                ++numBoundaryFacets_;
            } else {
                // f was the sole surface triangle exposing this facet as
                // un-glued (any other triangle sharing the same ambient edge
                // would have been joined to it in addTriangle instead), so
                // if it was flagged improper there, that flag must be
                // cleared now -- otherwise it lingers in improperEdges_
                // after f is gone, permanently poisoning isProper() for
                // whatever the surface goes on to become.
                regina::Edge<dim> *edge = f->edge(i);
                if (!edge->isBoundary()) {
                    improperEdges_.erase(edge);
                }
                --numBoundaryFacets_;
            }
        }

        ufUndoTo_(ufCheckpoints_[f->index()]);
        for (int i = 0; i < 3; ++i) {
            removeCornerImage_(f->vertex(i));
        }

        // src is always surface_'s current last triangle -- the DFS only
        // ever removes the triangle it most recently added -- so emb_ can
        // just shrink in lockstep instead of being keyed by src's pointer.
        emb_.pop_back();
        inv_[f->index()] = nullptr;
        indices_.erase(f->index());
        surface_.removeSimplex(src); // deletes src
        invariantsValid_ = false;
    }

    std::string detail() const {
        // if (!detail_.empty()) return detail_;

        ensureInvariants_();

        /* Generate surface details */
        bool isOrientable = invariants_[0];
        int punctures = invariants_[1];
        int genus = invariants_[2];
        std::ostringstream ans;

        if (!surface_.isConnected()) {
            throw regina::InvalidArgument("Surface must be connceted!");
        }

        if (isOrientable) {
            // Special names for surfaces with boundary:
            if (genus == 0 && punctures == 1)
                ans << "Disc";
            else if (genus == 0 && punctures == 2)
                ans << "Annulus";
            else {
                if (genus == 0)
                    ans << "Sphere";
                else if (genus == 1)
                    ans << "Torus";
                else
                    ans << "Orientable genus " << genus << " surface";

                if (punctures == 1)
                    ans << ", 1 puncture";
                else if (punctures > 1)
                    ans << ", " << punctures << " punctures";
            }
        } else {
            // Special names for surfaces with boundary:
            if (genus == 1 && punctures == 1)
                ans << "Möbius band";
            else {
                if (genus == 1)
                    ans << "Projective plane";
                else if (genus == 2)
                    ans << "Klein bottle";
                else
                    ans << "Non-orientable genus " << genus << " surface";

                if (punctures == 1)
                    ans << ", 1 puncture";
                else if (punctures > 1)
                    ans << ", " << punctures << " punctures";
            }
        }

        return ans.str();
    }

    // Independently recomputes every piece of this object's bookkeeping
    // from just tri_/surface_/emb_/inv_ (never trusting any of the other
    // incrementally-maintained state while doing so), and cross-checks it
    // against what's actually stored. Intended for tests, not the DFS hot
    // path -- most of this is O(current surface size) or worse. Returns
    // true iff everything agrees; on any mismatch, writes a description of
    // *every* mismatch found (not just the first) to `diag` and still
    // returns false, so a single failing test run shows the whole picture.
    bool checkInvariants(std::ostream &diag = std::cerr) const {
        bool ok = true;
        auto fail = [&](const std::string &msg) {
            diag << "[checkInvariants] " << msg << "\n";
            ok = false;
        };

        // --- emb_/inv_/indices_/surface_ mutual consistency ---
        if (emb_.size() != indices_.size())
            fail("emb_.size()=" + std::to_string(emb_.size()) +
                 " != indices_.size()=" + std::to_string(indices_.size()));
        if (surface_.countTriangles() != emb_.size())
            fail("surface_.countTriangles()=" +
                 std::to_string(surface_.countTriangles()) +
                 " != emb_.size()=" + std::to_string(emb_.size()));

        for (regina::Triangle<dim> *f : tri_->triangles()) {
            bool inIndices = indices_.count(f->index()) != 0;
            bool invSet = inv_[f->index()] != nullptr;
            if (inIndices != invSet) {
                fail("triangle " + std::to_string(f->index()) +
                     ": indices_ membership (" + std::to_string(inIndices) +
                     ") disagrees with inv_ being set (" +
                     std::to_string(invSet) + ")");
                continue;
            }
            if (invSet) {
                regina::Triangle<2> *t = inv_[f->index()];
                if (t->index() >= (int)emb_.size() || emb_[t->index()] != f) {
                    fail("triangle " + std::to_string(f->index()) +
                         ": inv_/emb_ round trip broken (inv_ points to "
                         "surface_ triangle " +
                         std::to_string(t->index()) +
                         ", which emb_ maps back to " +
                         (t->index() < (int)emb_.size()
                              ? std::to_string(emb_[t->index()]->index())
                              : "out-of-range"));
                }
            }
        }
        for (regina::Triangle<2> *t : surface_.triangles()) {
            if (t->index() >= (int)emb_.size()) {
                fail("surface_ triangle " + std::to_string(t->index()) +
                     " has no corresponding emb_ entry");
                continue;
            }
            const regina::Triangle<dim> *f = emb_[t->index()];
            if (inv_[f->index()] != t)
                fail("surface_ triangle " + std::to_string(t->index()) +
                     " maps via emb_ to ambient triangle " +
                     std::to_string(f->index()) +
                     ", but inv_ of that ambient triangle doesn't point "
                     "back to this surface_ triangle");
        }

        // --- numBoundaryFacets_ ---
        if (numBoundaryFacets_ != (int)surface_.countBoundaryFacets())
            fail("numBoundaryFacets_=" + std::to_string(numBoundaryFacets_) +
                 " != surface_.countBoundaryFacets()=" +
                 std::to_string(surface_.countBoundaryFacets()));

        // --- improperEdges_: an ambient edge belongs here iff it's not
        // itself on tri_'s boundary, but underlies at least one currently
        // un-glued facet of the surface. ---
        std::unordered_set<const regina::Edge<dim> *> groundTruthImproper;
        for (regina::Triangle<2> *t : surface_.triangles()) {
            const regina::Triangle<dim> *f = emb_[t->index()];
            for (int i = 0; i < 3; ++i) {
                if (t->adjacentSimplex(i) != nullptr)
                    continue;
                const regina::Edge<dim> *e = f->edge(i);
                if (!e->isBoundary())
                    groundTruthImproper.insert(e);
            }
        }
        if (groundTruthImproper.size() != improperEdges_.size()) {
            fail("improperEdges_ has " + std::to_string(improperEdges_.size()) +
                 " entries, ground truth scan found " +
                 std::to_string(groundTruthImproper.size()));
        } else {
            for (const regina::Edge<dim> *e : groundTruthImproper) {
                if (!improperEdges_.count(e))
                    fail("improperEdges_ is missing ambient edge " +
                         std::to_string(e->index()) +
                         " (currently un-glued in surface_, not on tri_'s "
                         "boundary)");
            }
        }

        // --- Self-intersection ground truth ---
        //
        // Deliberately NOT a second hand-rolled union-find: re-deriving the
        // "2+ roots at a vertex" idea with another from-scratch DSU would
        // only catch *drift* between the incremental and batch versions of
        // the *same* algorithm -- a conceptual bug in the algorithm itself
        // (e.g. if "same abstract vertex" were subtly the wrong notion)
        // would sail through both unchanged. Instead this cross-checks
        // against Regina's *own*, separately-implemented vertex-skeleton
        // computation on surface_ (regina::Triangulation<2>::vertices(),
        // the same machinery isOrientable()/eulerChar() rely on) -- a
        // genuinely independent source of "which corners are the same
        // abstract point", not a second copy of ufUnite_/ufFind_.
        std::unordered_map<size_t, size_t> maintainedToRegina,
            reginaToMaintained;
        std::unordered_map<size_t, const regina::Vertex<dim> *>
            maintainedVertex;
        std::unordered_map<size_t, std::unordered_set<size_t>> rootsAtVertex;
        for (regina::Vertex<2> *sv : surface_.vertices()) {
            // Well-definedness check: every embedding of sv (every
            // (surface_ triangle, local vertex) pair Regina considers the
            // same abstract point) should map, via emb_, to the *same*
            // ambient vertex -- this is what addTriangle's gluing is
            // supposed to guarantee, checked here directly against
            // Regina's own embeddings() rather than assumed.
            const regina::Vertex<dim> *ambient = nullptr;
            for (const auto &e : sv->embeddings()) {
                const regina::Triangle<dim> *f = emb_[e.simplex()->index()];
                const regina::Vertex<dim> *thisAmbient = f->vertex(e.face());
                if (ambient == nullptr) {
                    ambient = thisAmbient;
                } else if (ambient != thisAmbient) {
                    fail("surface_ vertex " + std::to_string(sv->index()) +
                         " (a single abstract point per Regina's own "
                         "skeleton) has embeddings mapping to two "
                         "different ambient vertices (" +
                         std::to_string(ambient->index()) + " and " +
                         std::to_string(thisAmbient->index()) + ")");
                }
            }
            if (ambient == nullptr)
                continue;
            rootsAtVertex[ambient->index()].insert(sv->index());

            // Cross-check against the maintained union-find on every
            // corner Regina attributes to sv: same two-way-consistency
            // shape as the emb_/inv_ round trip above.
            for (const auto &e : sv->embeddings()) {
                const regina::Triangle<dim> *f = emb_[e.simplex()->index()];
                size_t maintained = ufFind_(cornerKey_(f, e.face()));

                auto [mIt, mNew] =
                    maintainedToRegina.emplace(maintained, sv->index());
                if (!mNew && mIt->second != sv->index())
                    fail("maintained union-find group " +
                         std::to_string(maintained) +
                         " contains corners Regina's own skeleton "
                         "considers *different* abstract vertices (" +
                         std::to_string(mIt->second) + " and " +
                         std::to_string(sv->index()) +
                         ") -- ufUnite_ over-merged");
                auto [rIt, rNew] =
                    reginaToMaintained.emplace(sv->index(), maintained);
                if (!rNew && rIt->second != maintained)
                    fail("Regina vertex " + std::to_string(sv->index()) +
                         " (one abstract point) is split across two "
                         "different maintained union-find groups (" +
                         std::to_string(rIt->second) + " and " +
                         std::to_string(maintained) +
                         ") -- ufUnite_ under-merged");

                auto [vIt, vNew] =
                    maintainedVertex.emplace(maintained, ambient);
                if (!vNew && vIt->second != ambient)
                    fail("maintained union-find group " +
                         std::to_string(maintained) +
                         " contains corners mapping to two different "
                         "ambient vertices (" +
                         std::to_string(vIt->second->index()) + " and " +
                         std::to_string(ambient->index()) + ")");
            }
        }
        for (size_t vIdx = 0; vIdx < tri_->countVertices(); ++vIdx) {
            int groundTruthCount = 0;
            auto it = rootsAtVertex.find(vIdx);
            if (it != rootsAtVertex.end())
                groundTruthCount = (int)it->second.size();
            if (groundTruthCount != imageRootCount_[vIdx])
                fail("imageRootCount_[" + std::to_string(vIdx) +
                     "]=" + std::to_string(imageRootCount_[vIdx]) +
                     " != ground truth (Regina vertex count) " +
                     std::to_string(groundTruthCount));
        }
        int groundTruthSelfIntersections = 0;
        for (auto &[vIdx, roots] : rootsAtVertex)
            if (roots.size() > 1)
                ++groundTruthSelfIntersections;
        if (groundTruthSelfIntersections != selfIntersections_)
            fail("selfIntersections_=" + std::to_string(selfIntersections_) +
                 " != ground truth (Regina vertex count) " +
                 std::to_string(groundTruthSelfIntersections));

        // --- candidatesByFacet_: for every currently-included triangle's
        // every facet, an independent from-scratch scan of *every*
        // triangle in tri_ (using only the basic Triangle<dim>::edge()
        // primitive -- deliberately not Edge<dim>::embeddings(), which for
        // dim > 2 gives top-dimensional-simplex-level embeddings, not
        // triangle-level ones, and not any reuse of buildGluingEdges_'s own
        // O(n^2) double loop), compared as a multiset (duplicates matter --
        // see the class-level comment on candidatesByFacet_ re: an edge
        // shared by the same neighbour more than once). This is also,
        // incidentally, a direct regression check for the AdjList
        // multi-edge-sharing bug's failure mode, one layer down from where
        // that bug actually lived. ---
        for (int fi : indices_) {
            const regina::Triangle<dim> *f = tri_->triangle(fi);
            for (int j = 0; j < 3; ++j) {
                std::vector<size_t> groundTruth;
                const regina::Edge<dim> *e = f->edge(j);
                for (regina::Triangle<dim> *g : tri_->triangles()) {
                    for (int k = 0; k < 3; ++k) {
                        if (g == f && k == j)
                            continue;
                        if (g->edge(k) == e)
                            groundTruth.push_back(g->index());
                    }
                }
                std::vector<size_t> actual = candidatesByFacet_[fi][j];
                std::sort(groundTruth.begin(), groundTruth.end());
                std::sort(actual.begin(), actual.end());
                if (groundTruth != actual) {
                    std::ostringstream gt, ac;
                    for (size_t x : groundTruth)
                        gt << x << " ";
                    for (size_t x : actual)
                        ac << x << " ";
                    fail("candidatesByFacet_[" + std::to_string(fi) + "][" +
                         std::to_string(j) + "] = { " + ac.str() +
                         "} but ground truth (scanning tri_ directly) is "
                         "{ " +
                         gt.str() + "}");
                }
            }
        }

        return ok;
    }

    // The triangle index set already uniquely identifies a surface (the same
    // set of triangles is definitionally the same surface, and no two
    // different sets are "the same" surface), so it alone is sufficient for
    // std::set's ordering/dedup -- invariants_ is purely derived/display
    // data and doesn't need to participate here. That's what lets
    // computing invariants_ be deferred out of the DFS hot path entirely.
    friend bool operator==(const KnottedSurface &lhs,
                           const KnottedSurface &rhs) {
        return lhs.indices_ == rhs.indices_;
    }

    friend bool operator<(const KnottedSurface &lhs,
                          const KnottedSurface &rhs) {
        return lhs.indices_ < rhs.indices_;
    }

  private:
    static size_t cornerKey_(const regina::Triangle<dim> *f, int i) {
        return f->index() * 3 + static_cast<size_t>(i);
    }

    size_t ufFind_(size_t x) const {
        while (true) {
            size_t parent = ufParent_[x];
            if (parent == x)
                return x;
            x = parent;
        }
    }

    // Is facet `facet` of triangle f still able to ever be glued to
    // something? False (closed) if it's already glued (nothing left to
    // decide -- already resolved into whatever root it merged into) or if
    // every triangle that could ever occupy that facet has been
    // permanently excluded. True (open) if it's un-glued and at least one
    // not-yet-excluded candidate remains -- f might still end up glued
    // there later, potentially merging f's root with whatever that
    // neighbour's root turns out to be.
    bool isFacetOpen_(const regina::Triangle<dim> *f, int facet) const {
        regina::Triangle<2> *src = inv_[f->index()];
        if (src->adjacentSimplex(facet) != nullptr)
            return false;

        for (size_t neighborIdx : candidatesByFacet_[f->index()][facet]) {
            if (!excluded_[neighborIdx])
                return true;
        }
        return false;
    }

    // A union-find root (given as its member raw corners, all currently
    // mapping to the same ambient vertex) is "closed" -- can never grow or
    // merge with anything else again -- once every member corner's two
    // vertex-incident facets (the two facets other than the corner's own
    // local vertex index) are closed. See hasUnresolvableConflict()'s
    // class-level comment for why two simultaneously-closed roots at the
    // same vertex means a permanent, unresolvable conflict.
    bool isRootClosed_(
        const std::vector<std::pair<const regina::Triangle<dim> *, int>>
            &members) const {
        for (const auto &[f, c] : members) {
            for (int facet = 0; facet < 3; ++facet) {
                if (facet == c)
                    continue;
                if (isFacetOpen_(f, facet))
                    return false;
            }
        }
        return true;
    }

    // Registers a freshly-created raw corner mapping to ambient vertex v,
    // updating the self-intersection count if this creates (or grows) a
    // conflict.
    void addCornerImage_(const regina::Vertex<dim> *v) {
        int oldCount = imageRootCount_[v->index()]++;
        if (oldCount == 1)
            ++selfIntersections_;
    }

    // Inverse of addCornerImage_(), used both when unwinding a failed add
    // and when a corner's owning triangle is fully removed.
    void removeCornerImage_(const regina::Vertex<dim> *v) {
        int oldCount = imageRootCount_[v->index()]--;
        if (oldCount == 2)
            --selfIntersections_;
    }

    // Unions the union-find components of raw corners x and y, which must
    // both currently map to ambient vertex v (guaranteed by construction:
    // the gluing permutation used to reach here was derived directly from
    // the ambient triangulation's own edge identifications). No-op if
    // they're already in the same component (e.g. a triangle glued to
    // itself along another edge that already connected them).
    void ufUnite_(size_t x, size_t y, const regina::Vertex<dim> *v) {
        size_t rx = ufFind_(x), ry = ufFind_(y);
        if (rx == ry)
            return;
        if (ufSize_[rx] < ufSize_[ry])
            std::swap(rx, ry);

        int oldParentSize = ufSize_[rx];
        ufUndoLog_.push_back({ry, rx, oldParentSize, v});
        ufParent_[ry] = rx;
        ufSize_[rx] = oldParentSize + ufSize_[ry];

        removeCornerImage_(v);
    }

    // Rolls ufUndoLog_ back to the given checkpoint (a prior
    // ufUndoLog_.size()), in LIFO order -- correct because DFS add/remove is
    // properly nested (see the class-level comment on ufParent_).
    void ufUndoTo_(size_t checkpoint) {
        while (ufUndoLog_.size() > checkpoint) {
            const UfUndoEntry &e = ufUndoLog_.back();
            ufParent_[e.childRoot] = e.childRoot;
            ufSize_[e.parentRoot] = e.parentOldSize;
            addCornerImage_(e.image);
            ufUndoLog_.pop_back();
        }
    }

    // Common cleanup when addTriangle() must bail out (a double-glued
    // facet, or a self-intersection): undo this attempt's union-find
    // operations and singleton corners, then the surface_/emb_/inv_ state
    // that addTriangle() had already set up for f.
    void undoFailedAdd_(const regina::Triangle<dim> *f, size_t ufCheckpoint,
                        regina::Triangle<2> *src) {
        ufUndoTo_(ufCheckpoint);
        for (int i = 0; i < 3; ++i) {
            removeCornerImage_(f->vertex(i));
        }
        surface_.removeSimplex(src);
        emb_.pop_back();
        inv_[f->index()] = nullptr;
        // No-op if the double-glued-facet check (earlier in addTriangle,
        // before f ever enters indices_) is what called us; undoes the
        // hasUnresolvableConflict()-triggered insertion otherwise.
        indices_.erase(f->index());
    }

    void ensureInvariants_() const {
        if (invariantsValid_)
            return;

        bool isOrientable = surface_.isOrientable();
        int punctures = surface_.countBoundaryComponents();
        int genus = isOrientable ? (2 - surface_.eulerChar() - punctures) / 2
                                 : 2 - surface_.eulerChar() - punctures;

        invariants_[0] = isOrientable;
        invariants_[1] = punctures;
        invariants_[2] = genus;
        invariantsValid_ = true;
    }
};

#endif
