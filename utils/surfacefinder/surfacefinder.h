//
//  surfacefinder.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef SURFACEFINDER_H

#define SURFACEFINDER_H

#include "knottedsurfaces.h"

template <int dim>
class SurfaceFinder {
  private:
    size_t calls_ = 0;

    SurfaceCondition cond_;

    const regina::Triangulation<dim> &tri_;
    /**< The given triangulation of a dim-manifold */

    std::vector<GluingNode<dim>> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    std::set<KnottedSurface<dim>> surfaces_;

    KnottedSurface<dim> surface_;

  public:
    SurfaceFinder(const regina::Triangulation<dim> &tri, SurfaceCondition cond)
        : tri_(tri), surface_(&tri), cond_(cond) {
        buildGluingNodes_();
        std::cout << "[+] Built gluing nodes\n";
        buildGluingEdges_();
        std::cout << "[+] Built gluing edges\n";
    }

    std::set<KnottedSurface<dim>> &surfaces() { return surfaces_; }

    std::set<KnottedSurface<dim>> &findSurfaces() {
        surfaces_.clear();

        std::cout << "\n";
        for (auto &seed : nodes_) {
            // Reset the surface and visited flags
            surface_.reset();
            for (auto &node : nodes_) {
                node.visited = false;
            }

            if (!surface_.addTriangle(seed.f, seed.adjList))
                continue;
            seed.visited = true;
            checkWin_();

            std::vector<GluingNode<dim> *> frontier;
            for (auto &[neighbor, g] : seed.adjList) {
                if (!neighbor->visited)
                    frontier.push_back(neighbor);
            }

            extend_(frontier, 0, seed.f->index());

            // Symmetric with extend_'s own add/remove discipline: undo the
            // seed itself before the next iteration's reset.
            seed.visited = false;
            surface_.removeTriangle(seed.f);
        }
        std::cout << "\n\n";
        std::cout << "[+] Total search calls = " << calls_ << "\n\n";

        return surfaces_;
    }

    // Runs the same search as findSurfaces(), but with a set of triangles
    // that must all already be part of the surface (e.g. an annulus whose
    // boundary is a known knot/link) -- the search only ever considers
    // *extending* this fixed base, never removing from it.
    std::set<KnottedSurface<dim>> &findSurfaces(
        std::unordered_set<regina::Triangle<dim> *> &startingTriangles) {
        surfaces_.clear();

        std::unordered_set<GluingNode<dim> *> startingNodes;
        // Find the nodes corresponding to the starting triangles (can be made
        // O(1) if we use a map, but this is fine for now)
        for (auto &node : nodes_) {
            if (startingTriangles.find(node.f) != startingTriangles.end()) {
                startingNodes.insert(&node);
            }
        }
        if (startingNodes.size() != startingTriangles.size()) {
            throw regina::InvalidArgument(
                "SurfaceFinder::findSurfaces: Starting "
                "triangle(s) not found in the graph");
        }

        surface_.reset();
        for (auto &node : nodes_) {
            node.visited = false;
        }

        std::cout << "\n";
        for (auto *base : startingNodes) {
            // Defer the self-intersection check (see addTriangle's
            // documentation): a base with an internal branch point can
            // legitimately have two of its own triangles both touching an
            // ambient vertex before the triangle connecting them has been
            // inserted, even though the full base is perfectly valid once
            // everything is in. Rejecting per-insertion here would make
            // the result depend on startingNodes' (hash-based, arbitrary)
            // iteration order instead of on whether the base is actually
            // valid. A double-glued facet is still a hard, immediate
            // failure -- that can never be "fixed" by a later insertion.
            if (!surface_.addTriangle(base->f, base->adjList,
                                      /*checkSelfIntersection=*/false))
                throw regina::InvalidArgument(
                    "SurfaceFinder::findSurfaces: Starting triangles do not "
                    "form a valid embedded partial surface (an edge is "
                    "glued twice)");
            base->visited = true;
        }

        // Now that every base triangle is in, check self-intersection
        // exactly once, over the whole base. This is the check that was
        // deferred above -- see addTriangle's documentation for why doing
        // it per-triangle during the loop is wrong for a base with
        // internal branch points.
        if (surface_.hasSelfIntersection()) {
            throw regina::InvalidArgument(
                "SurfaceFinder::findSurfaces: Starting triangles do not "
                "form a valid embedded partial surface (self-intersects)");
        }

        if (!surface_.surface().isConnected()) {
            throw regina::InvalidArgument(
                "SurfaceFinder::findSurfaces: Starting triangles do not form "
                "a connected partial surface");
        }

        checkWin_();

        std::vector<GluingNode<dim> *> frontier;
        for (auto *base : startingNodes) {
            for (auto &[neighbor, g] : base->adjList) {
                if (!neighbor->visited)
                    frontier.push_back(neighbor);
            }
        }

        extend_(frontier, 0);

        std::cout << "\n\n";
        std::cout << "[+] Total search calls = " << calls_ << "\n\n";

        return surfaces_;
    }

    int countNodes() const { return nodes_.size(); }

    int countEdges() const {
        int count = 0;
        for (const auto &node : nodes_) {
            count += node.adjList.size();
        }
        return count;
    }

    // Exposes a triangle's real gluing-graph adjacency, so tests can drive
    // KnottedSurface::addTriangle directly with genuine adjacency data
    // instead of re-deriving gluing permutations by hand.
    const typename GluingNode<dim>::AdjList &
    adjacencyOf(const regina::Triangle<dim> *f) const {
        for (const auto &node : nodes_) {
            if (node.f == f)
                return node.adjList;
        }
        throw regina::InvalidArgument(
            "SurfaceFinder::adjacencyOf: triangle not found in the graph");
    }

    friend std::ostream &operator<<(std::ostream &os,
                                    const SurfaceFinder<dim> &graph) {
        os << "SurfaceFinder<" << dim << ">(nodes = " << graph.countNodes()
           << ", edges = " << graph.countEdges() << ")\n";
        for (const auto &node : graph.nodes_) {
            os << "  Triangle: " << node.f->index() << "\n";
            for (const auto &[adjNode, gluing] : node.adjList) {
                os << "    Adjacent: " << adjNode->f->index() << ", " << gluing
                   << "\n";
            }
        }
        return os;
    }

  private:
    /**
     * Start building our gluing graph. Associate a node to each triangle in the
     * triangulation.
     */
    void buildGluingNodes_() {
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            nodes_.emplace_back(triangle);
        }
    }

    /**
     * Populate the node's adjacency lists.
     */
    void buildGluingEdges_() {
        for (auto &node : nodes_) {
            std::unordered_set<int> selfGluingEdges;
            for (int i = 0; i < 3; ++i) {
                regina::Triangle<dim> *triangle = node.f;
                regina::Edge<dim> *edge = triangle->edge(i);
                regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(i);

                // TODO: Would be nice for this to be O(1), not sure if possible
                for (auto &otherNode : nodes_) {
                    regina::Triangle<dim> *otherTriangle = otherNode.f;
                    for (int j = 0; j < 3; ++j) {
                        // Only glue identified edges and ignore identity
                        // mappings/opposite self gluings
                        if (otherTriangle->edge(j) != edge ||
                            (triangle == otherTriangle &&
                             (i == j || selfGluingEdges.find(i) !=
                                            selfGluingEdges.end())))
                            continue;

                        if (triangle == otherTriangle) {
                            selfGluingEdges.insert(j);
                        }

                        regina::Perm<dim + 1> edgeToOther =
                            otherTriangle->edgeMapping(j);
                        // regina::Perm is immutable, use std::array instead
                        std::array<int, 3> p;
                        for (int k = 0; k < 3; ++k) {
                            p[edgeToTriangle[k]] = edgeToOther[k];
                        }

                        node.adjList[&otherNode] = {triangle, i, otherTriangle,
                                                    p};
                    }
                }
            }
        }
    }

    // Win condition: we've found a surface that meets the given boundary
    // requirement and has no self-intersection. Called after every triangle
    // is added to surface_, regardless of whether the search goes on to add
    // more triangles afterward.
    //
    // The hasSelfIntersection() check here is a deliberate, unconditional
    // final gate, independent of whatever checks already happened upstream
    // (extend_'s per-triangle check, or findSurfaces(startingTriangles)'s
    // deferred whole-base check): this is the one place a surface actually
    // gets recorded, so it's the one place that must never let a
    // self-intersecting surface through, regardless of how surface_ got
    // into its current state. Embedded-ness is the entire point of this
    // tool -- this check should never be relied on to be redundant.
    void checkWin_() {
        if (surface_.hasSelfIntersection())
            return;

        if ((cond_ == SurfaceCondition::all) ||
            (cond_ == SurfaceCondition::closed && surface_.isClosed()) ||
            (cond_ == SurfaceCondition::boundary && surface_.isProper())) {
            surfaces_.insert(surface_);
        }
    }

    // Enumerates every way to extend the surface currently held in
    // surface_ using some subset of `frontier` (candidate triangles
    // adjacent to the surface, not yet decided one way or the other).
    //
    // Unlike a plain "add one triangle, recurse into its neighbours,
    // backtrack" DFS -- which can only ever discover connected triangle
    // sets that happen to admit a Hamiltonian path through the ambient
    // adjacency graph, since it only ever looks at the most-recently-added
    // triangle's neighbours -- this tries every remaining candidate against
    // the *entire* frontier accumulated so far, so a triangle with two or
    // more simultaneously-attached neighbours that aren't themselves
    // adjacent (a genuine branch point) is still found.
    //
    // frontier[0..idx) have already been decided (included or excluded) by
    // an ancestor call; frontier[idx..] are still undecided. Every call
    // restores `frontier` to the state it was given before returning, so
    // sibling branches (and the caller) always see a consistent frontier.
    //
    // seedIdx is the index of the triangle that seeded this whole search
    // (0 -- i.e. no restriction -- for the fixed-base findSurfaces()
    // overload, which isn't looped over seeds and so has no redundancy to
    // avoid). Every embedded surface has a unique minimum-index triangle, so
    // refusing to include any candidate with index < seedIdx means a given
    // surface is only ever assembled starting from that one triangle as
    // seed, instead of being rediscovered once per triangle it contains.
    void extend_(std::vector<GluingNode<dim> *> &frontier, size_t idx,
                 size_t seedIdx = 0) {
        ++calls_;
        if (calls_ % 1000000 == 0)
            std::cout << "\r" << std::flush << "[*] Call " << (calls_ / 1000000)
                      << "M, Surface size " << surface_.surface().size()
                      << ", Surfaces " << surfaces_.size() << "          ";

        if (idx >= frontier.size())
            return;

        GluingNode<dim> *w = frontier[idx];

        // Branch 1: permanently exclude w from the surface for the rest of
        // this search path, and move on to the next candidate.
        extend_(frontier, idx + 1, seedIdx);

        // Branch 2: include w (if it isn't already part of the surface via
        // some other branch, its index doesn't fall below the seed's, and
        // it can be validly added).
        if (!w->visited && w->f->index() >= seedIdx &&
            surface_.addTriangle(w->f, w->adjList)) {
            w->visited = true;
            checkWin_();

            size_t frontierSizeBefore = frontier.size();
            for (auto &[neighbor, g] : w->adjList) {
                if (!neighbor->visited)
                    frontier.push_back(neighbor);
            }

            extend_(frontier, idx + 1, seedIdx);

            frontier.resize(frontierSizeBefore);
            w->visited = false;
            surface_.removeTriangle(w->f);
        }
    }
};

#endif
