//
//  knottedsurfaces.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef KNOTTED_SURFACES_H

#define KNOTTED_SURFACES_H

#include <memory>
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
    std::shared_ptr<const regina::Triangulation<dim - 1>> bdry_;
    /**< tri_'s boundary, built once (in the very first KnottedSurface
     * constructed for a given tri_) and shared by pointer across every
     * subsequent copy -- it depends only on the fixed ambient triangulation,
     * never on anything the search mutates, so there's no reason to rebuild
     * it from scratch every time a surface gets copied into the results
     * set. */
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
          ufCheckpoints_(tri_->countTriangles()) {
        if (!tri_->isClosed()) {
            bdry_ = std::make_shared<const regina::Triangulation<dim - 1>>(
                tri_->boundaryComponent(0)->build());
        }
    }

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_), bdry_(other.bdry_), surface_(other.surface_),
          emb_(other.emb_), inv_(other.inv_.size(), nullptr),
          invariants_(other.invariants_),
          invariantsValid_(other.invariantsValid_), indices_(other.indices_),
          ufParent_(other.ufParent_), ufSize_(other.ufSize_),
          imageRootCount_(other.imageRootCount_),
          selfIntersections_(other.selfIntersections_),
          ufUndoLog_(other.ufUndoLog_), ufCheckpoints_(other.ufCheckpoints_),
          numBoundaryFacets_(other.numBoundaryFacets_),
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
            bdry_ = other.bdry_;
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

    Link boundary() const {
        std::unordered_set<const regina::Edge<dim - 1> *> edges;

        // TODO: Make this O(1) somehow
        for (const regina::BoundaryComponent<dim> *triBoundaryComp :
             tri_->boundaryComponents()) {

            for (const regina::BoundaryComponent<2> *surfaceBoundaryComp :
                 surface_.boundaryComponents()) {
                for (const regina::Edge<2> *e : surfaceBoundaryComp->edges()) {
                    const regina::Edge<dim> *im = image(e);
                    for (int i = 0; i < triBoundaryComp->countEdges(); ++i) {
                        if (im == triBoundaryComp->edge(i)) {
                            edges.insert(bdry_->edge(i));
                            break;
                        }
                    }
                }
            }
        }

        std::vector<const regina::Edge<dim - 1> *> edgeList(edges.begin(),
                                                            edges.end());
        return {*bdry_, edgeList};
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

        for (const auto &[adjNode, g] : adj) {
            regina::Triangle<2> *dst = inv_[adjNode->f->index()];
            if (dst == nullptr)
                continue;

            int dstFacet = g.gluing[g.srcFacet];

            // Saves us from catching an error
            if (dst->adjacentSimplex(dstFacet) != nullptr ||
                src->adjacentSimplex(g.srcFacet) != nullptr) {
                undoFailedAdd_(f, ufCheckpoint, src);
                return false;
            }

            isBoundaryEdge[g.srcFacet] = false;
            ++boundaryFacetsConsumed;
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

        if (checkSelfIntersection && hasSelfIntersection()) {
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

        // Each consumed facet removes one boundary slot on our own side (it
        // was never counted) and flips one on the neighbour's side from
        // boundary to interior; the remaining (3 - boundaryFacetsConsumed)
        // facets are new boundary slots of our own.
        numBoundaryFacets_ += 3 - 2 * boundaryFacetsConsumed;

        ufCheckpoints_[f->index()] = ufCheckpoint;
        indices_.insert(f->index());
        invariantsValid_ = false;

        return true;
    }

    void removeTriangle(const regina::Triangle<dim> *f) {
        regina::Triangle<2> *src = inv_[f->index()];

        for (int i = 0; i < 3; ++i) {
            regina::Triangle<2> *dst = src->adjacentSimplex(i);

            if (dst != nullptr) {
                regina::Edge<dim> *edge = f->edge(i);
                if (!edge->isBoundary()) {
                    improperEdges_.insert(edge);
                }
                ++numBoundaryFacets_;
            } else {
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
