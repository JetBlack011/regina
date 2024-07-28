//
//  surfacefinder.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <array>
#include <cstring>
#include <iostream>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>

#include "maths/perm.h"
#include "triangulation/example2.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/faceembedding.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

/* Constants */
static const regina::Triangulation<2> GENUS_2_SURFACE =
    regina::Triangulation<2>::fromGluings(6, {{0, 2, 5, {2, 1, 0}},
                                              {0, 1, 1, {0, 2, 1}},
                                              {0, 0, 5, {1, 0, 2}},
                                              {1, 1, 2, {0, 2, 1}},
                                              {1, 0, 3, {0, 2, 1}},
                                              {2, 1, 3, {0, 2, 1}},
                                              {2, 0, 4, {0, 2, 1}},
                                              {3, 1, 4, {0, 2, 1}},
                                              {4, 1, 5, {0, 2, 1}}});

// template <int n>
// regina::Triangulation<3> knotComplement(const regina::Triangulation<3> t,
//                                         const Knot &k) {
//     return linkComplement(t, {k});
// }

enum SurfaceCondition { all, boundary, closed };

template <int n>
using Curve = std::vector<regina::Edge<n> *>;

template <int n>
using DualCurve = std::vector<regina::Triangle<n> *>;

/** Gluing Implementation */

template <int dim>
struct Gluing {
    regina::Face<dim, 2> *src;
    int srcFacet;
    regina::Face<dim, 2> *dst;
    regina::Perm<3> gluing;

    Gluing() = default;

    Gluing(regina::Face<dim, 2> *src, int srcFacet, regina::Face<dim, 2> *dst,
           regina::Perm<3> gluing)
        : src(src), srcFacet(srcFacet), dst(dst), gluing(gluing) {}

    friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
        return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
                  << g.dst->index() << ", " << g.gluing << ")";
    }
};

/** Knot implementation */

class Knot {
   private:
    regina::Triangulation<3> &tri_;
    Curve<3> edges_;
    std::unordered_set<const regina::Edge<3> *> edgeSet_;
    std::unordered_map<const regina::Tetrahedron<3> *, int> tetEdgeCount_;

    bool isUnknot_ = false;

    /**
     * Finds a push off of the given knot onto the dual 1-skeleton of the
     * triangulation it's sitting in. Note that this may change the underlying
     * triangulation.
     */
    DualCurve<3> toDualCurve_(const Knot &knot);

   public:
    Knot(regina::Triangulation<3> &tri, const Curve<3> &edges)
        : tri_(tri), edges_(edges) {
        // Precompute the set of edges and the count of edges in each
        // tetrahedron
        for (const regina::Edge<3> *edge : edges) {
            edgeSet_.insert(edge);
            for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
                auto search = tetEdgeCount_.find(emb.tetrahedron());
                if (search == tetEdgeCount_.end()) {
                    tetEdgeCount_.insert({emb.tetrahedron(), 1});
                } else {
                    ++search->second;
                }
            }
        }
    }

    /**
     * Simplifies a given knot so that each tetrahedron it runs through contains
     * at most one of the knot's edges.
     *
     * Algorithm:
     * 0) If we're the unknot (e.g. all edges run through a single tetrahedron)
     * or if every tetrahedron contains exactly one edge of the knot, then we're
     * done.
     * 1) For each tetrahedron containing an edge, if it contains more than one,
     * we have two cases:
     *    i. The edges are connected. In this case, replace these edges with the
     * single edge running from the two endpoints of the path formed in this
     * tetrahedron. Remember to update the tetrahedron edge counts.
     *     ii. There are exactly two edges using all of the vertices of the
     * tetrahedron. In this case, do a 1-4 move. Remember to update the
     * tetrahedron edge counts, including modifying which tetrahedra edges run
     * through.
     * 2) Repeat from 0). It is very possible that we've reduced to only one
     * terahedron, in which case we've arrived at the unknot.
     */
    void simplify() {
        bool isDone = false;
        while (!isUnknot_ && !isDone) {
            for (const auto &[tet, count] : tetEdgeCount_) {
                if (count == 2 && numVerticesUsed_(tet) == 4) {
                    // Separated case

                } else if (count >= 2) {
                    // Connected case
                }
            }
        }
    }

   private:
    int numVerticesUsed_(const regina::Tetrahedron<3> *tet) {
        std::unordered_set<regina::Vertex<3> *> verts;
        // For each edge in tet
        for (int i = 0; i < 6; ++i) {
            auto search = edgeSet_.find(tet->edge(i));
            if (search != edgeSet_.end()) {
                const regina::Edge<3> *e = *search;
                verts.insert(e->vertex(0));
                verts.insert(e->vertex(1));
            }
        }
        return verts.size();
    }

    std::array<regina::Tetrahedron<3> *, 4> fourOneTet() {
        auto tets = tri_.newTetrahedra<4>();

        for (int i = 0; i < 3; ++i) {
            for (int j = i + 1; j < 4; ++j) {
                tets[i]->join(j, tets[j], {i, j});
            }
        }

        // Same as
        // tets[0]->join(1, tets[1], {1, 0, 2, 3});
        // tets[0]->join(2, tets[2], {2, 1, 0, 3});
        // tets[0]->join(3, tets[3], {3, 1, 2, 0});
        // tets[1]->join(2, tets[2], {0, 2, 1, 3});
        // tets[1]->join(3, tets[3], {0, 3, 2, 1});
        // tets[2]->join(3, tets[3], {0, 1, 3, 2});

        return tets;
    }

    void fourOneMove(regina::Tetrahedron<3> *tet) {
        // std::array<std::optional<Gluing>, 4> gluings;
        // std::vector<regina::Edge<3> *> edgeMap;
        // std::unordered_map<regina::Edge<3> *, std::pair<int, int>> edges;

        // int numTets = tri_.countTetrahedra();

        //// For each triangle in tet
        // for (int i = 0; i < 4; ++i) {
        //     regina::Tetrahedron<3> *adjTet = tet->adjacentSimplex(i);
        //     if (adjTet != nullptr) {
        //         gluings[i] = {tet, i, adjTet, tet->adjacentGluing(i)};
        //     }
        // }

        // tri_.removeTetrahedron(tet);
        // auto tets = fourOneTet();

        // for (int i = 0; i < 4; ++i) {
        //     if (auto g = gluings[i]) {
        //         tets[i]->join(i, g->dst, g->gluing);
        //     }
        // }

        //// Update edge data
        // for (int i = 0; i < numTets - 1; ++i) {
        // }
    }
};

/** Link Implementation */

// std::pair<EdgeMap, regina::Triangulation<3>> drillDualCurve(
//     const regina::Triangulation<3> &tri, const DualCurve<3> &curve) {
//
// }

class Link {
   private:
    regina::Triangulation<3> &tri_;
    std::vector<Knot> components_;

   public:
    Link(regina::Triangulation<3> &tri, const std::vector<Curve<3>> &components)
        : tri_(tri) {
        for (const Curve<3> &edges : components) {
            components_.emplace_back(tri, edges);
        }
    }

    // regina::Triangulation<3> linkComplement(const
    // regina::Triangulation<3> &t,
    //                                         const Link &l) {
    //     regina::Triangulation<3> tri = t;
    //     Link link = l;
    //
    //     while (!link.empty()) {
    //         const Knot &knot = link.front();
    //         // 1) Push off the current link component
    //         DualCurve<3> curve = knotToDualCurve(knot);
    //
    //         /// 2) Drill, keep track of where the remaining link
    //         components end
    //
    //         const auto &[edgeMap, newTri] = drillDualCurve(tri, curve);
    //
    //         Link newLink;
    //         newLink.reserve(link.size() - 1);
    //         for (const Knot &knot2 : link) {
    //             if (knot2 == knot) continue;  // Ignore the now
    //             drilled-out
    //     ot
    //
    //             Knot newKnot;
    //             newKnot.reserve(knot.size());
    //             for (const regina::Edge<3> *edge : knot2) {
    //                 auto search = edgeMap.find(edge->index());
    //                 if (search == edgeMap.end()) {
    //                     throw regina::InvalidArgument(
    //                         "Could not find corresponding edge in
    //                         edgeMap!");
    //                 }
    //                 newKnot.push_back(newTri.edge(search->second));
    //             }
    //             newLink.push_back(newKnot);
    //         }
    //
    //         link = newLink;
    //         tri = newTri;
    //         // Repeat with newLink in the drilled triangulation
    //     }
    //
    //     // Now tri is a triangulation of the complement of the original
    //     link std::cout << tri.detail() << "\n";
    //
    //     return tri;
    // }
};

template <int dim>
class GluingNode {
   public:
    using AdjList = std::unordered_map<GluingNode *, Gluing<dim>>;

    regina::Triangle<dim> *f;
    AdjList adjList;
    AdjList invalids;
    bool valid = true;
    bool visited = false;

    GluingNode(regina::Face<dim, 2> *f) : f(f) {}

    friend std::ostream &operator<<(std::ostream &os, const GluingNode &node) {
        os << "(" << node.f->index() << " { ";
        for (const auto &[_, gluing] : node.adjList) {
            os << gluing << " ";
        }

        os << "})";

        return os;
    }
};

/** KnottedSurface Implementation */

template <int dim>
class KnottedSurface {
   private:
    const regina::Triangulation<dim> *tri_;
    /**< The dim-manifold triangulation in which this subdim-manifold is
     * embedded */
    regina::Triangulation<2> surface_;
    /**< The triangulation of the sub-triangulation induced by the embedding
     */

    std::unordered_map<regina::Triangle<2> *, const regina::Triangle<dim> *>
        emb_;
    /**< A mapping from the top-dimensional simplices of the
     * sub-triangulation into the subdim-simplices of the dim-triangulation
     */
    std::unordered_map<const regina::Triangle<dim> *, regina::Triangle<2> *>
        inv_;
    /**< A mapping from the subdim-simplices of the dim-triangulation to the
     * top-dimensional simplices of the sub-triangulation. Note that
     * inv_.at(emb_.at(face)) = face */

    std::array<int, 3> invariants_;

    std::set<int> indices_;

   public:
    std::unordered_set<const regina::Edge<dim> *> improperEdges_;
    /**< Notes whether an embedding is proper, in the sense that the image
     * of the boundary of the sub-triangulation is entirely contained within
     * the boundary of the dim-manifold triangulation */
    KnottedSurface(const regina::Triangulation<dim> *tri) : tri_(tri) {}

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_),
          surface_(other.surface_),
          improperEdges_(other.improperEdges_),
          invariants_(other.invariants_),
          indices_(other.indices_) {
        // Connect the wires...
        for (regina::Triangle<2> *t : other.surface_.triangles()) {
            const regina::Triangle<dim> *f = other.emb_.at(t);
            regina::Triangle<2> *newT = surface_.triangle(t->index());
            emb_[newT] = f;
            inv_[f] = newT;
        }
    }

    KnottedSurface &operator=(const KnottedSurface &other) {
        if (this != &other) {
            tri_ = other.tri_;
            surface_ = other.surface_;
            improperEdges_ = other.improperEdges_;
            invariants_ = other.invariants_;
            indices_ = other.indices_;

            // Connect the wires...
            for (regina::Triangle<2> *t : other.surface_.triangles()) {
                const regina::Triangle<dim> *f = other.emb_.at(t);
                regina::Triangle<2> *newT = surface_.triangle(t->index());
                emb_[newT] = f;
                inv_[f] = newT;
            }
        }

        return *this;
    }

    const regina::Triangulation<2> &surface() const { return surface_; }

    bool isProper() const { return improperEdges_.empty(); }

    template <int facedim>
    regina::Face<dim, facedim> *image(const regina::Face<2, facedim> *f) const {
        static_assert(0 <= dim <= 2, "Must have 0 <= facedim <= 2");

        // for (auto &[t, f] : emb_) {
        //     std::cout << "(" << t->index() << ", " << f->index() << ") ";
        // }
        // std::cout << "\n";

        const regina::FaceEmbedding<2, facedim> &fEmb = f->front();
        return emb_.at(fEmb.simplex())->template face<facedim>(fEmb.face());
    }

    const regina::Triangle<dim> *image(regina::Triangle<2> *t) const {
        return emb_.at(t);
    }

    regina::Triangle<2> *preimage(const regina::Triangle<dim> *f) const {
        auto search = inv_.find(f);
        if (search == inv_.end()) return nullptr;
        return search->second;
    }

    /** Mutation methods */

    bool addTriangle(const regina::Triangle<dim> *f,
                     const typename GluingNode<dim>::AdjList &adj) {
        regina::Triangle<2> *src = surface_.newSimplex();
        std::array<bool, 3> isBoundaryEdge = {true, true, true};

        emb_[src] = f;
        inv_[f] = src;

        for (const auto &[adjNode, g] : adj) {
            auto search = inv_.find(adjNode->f);
            if (search == inv_.end()) continue;

            regina::Triangle<2> *dst = search->second;
            int dstFacet = g.gluing[g.srcFacet];

            // Saves us from catching an error
            if (dst->adjacentSimplex(dstFacet) != nullptr ||
                src->adjacentSimplex(g.srcFacet) != nullptr) {
                surface_.removeSimplex(src);
                emb_.erase(src);
                inv_.erase(f);
                return false;
            }

            isBoundaryEdge[g.srcFacet] = false;
            src->join(g.srcFacet, dst, g.gluing);
        }

        for (int i = 0; i < 3; ++i) {
            if (!isBoundaryEdge[i]) {
                improperEdges_.erase(f->edge(i));
            } else if (isBoundaryEdge[i] && !f->edge(i)->isBoundary()) {
                improperEdges_.insert(f->edge(i));
            }
        }

        indices_.insert(f->index());
        updateInvariants_();

        return true;
    }

    void removeTriangle(const regina::Triangle<dim> *f) {
        regina::Triangle<2> *src = inv_.at(f);

        for (int i = 0; i < 3; ++i) {
            regina::Triangle<2> *dst = src->adjacentSimplex(i);

            if (dst != nullptr) {
                regina::Edge<dim> *edge = f->edge(i);
                if (!edge->isBoundary()) {
                    improperEdges_.insert(edge);
                }
            }
        }

        surface_.removeSimplex(src);
        emb_.erase(src);
        inv_.erase(f);
        indices_.erase(f->index());
        updateInvariants_();
    }

    std::string detail() const {
        // if (!detail_.empty()) return detail_;

        /* Generate surface details */
        bool isOrientable = invariants_[0];
        int punctures = invariants_[1];
        int genus = invariants_[2];
        std::ostringstream ans;

        if (!surface_.isConnected()) {
            throw regina::InvalidArgument("Surface must be connceted!");
        }

        if (isOrientable) {
            // Special names for surface_s with boundary:
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
                    ans << "genus " << genus << " torus";

                if (punctures == 1)
                    ans << ", 1 puncture";
                else if (punctures > 1)
                    ans << ", " << punctures << " punctures";
            }
        } else {
            // Special names for surface_s with boundary:
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

    friend bool operator==(const KnottedSurface &lhs,
                           const KnottedSurface &rhs) {
        return lhs.invariants_ == rhs.invariants_ &&
               lhs.indices_ == rhs.indices_;
    }

    friend bool operator<(const KnottedSurface &lhs,
                          const KnottedSurface &rhs) {
        return lhs.invariants_ < rhs.invariants_ ||
               (lhs.invariants_ == rhs.invariants_ &&
                lhs.indices_ < rhs.indices_);
    }

   private:
    void updateInvariants_() {
        bool isOrientable = surface_.isOrientable();
        int punctures = surface_.countBoundaryComponents();
        int genus = isOrientable ? (2 - surface_.eulerChar() - punctures) / 2
                                 : 2 - surface_.eulerChar() - punctures;

        invariants_[0] = isOrientable;
        invariants_[1] = punctures;
        invariants_[2] = genus;
    }
};

template <int dim>
class GluingGraph {
   private:
    size_t calls_ = 0;

    SurfaceCondition cond_;

    const regina::Triangulation<dim> &tri_;
    /**< The given triangulation of a dim-manifold */

    std::map<regina::Triangle<dim> *, GluingNode<dim>> nodes_;
    /**< Guaranteed not to be nodes with the same gluing data */

    std::set<KnottedSurface<dim>> surfaces_;

    KnottedSurface<dim> surface_;

   public:
    GluingGraph(const regina::Triangulation<dim> &tri, SurfaceCondition cond)
        : tri_(tri), surface_(&tri), cond_(cond) {
        buildGluingNodes_();
        buildGluingEdges_();
        // buildGluingInvalids_();
    }

    std::set<KnottedSurface<dim>> &findSurfaces() {
        if (!surfaces_.empty()) return surfaces_;

        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            reset_();
            dfs_(&nodes_.at(triangle), 0);
            // std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";
        }
        std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";

        return surfaces_;
    }

   private:
    void nspaces_(int m) {
        for (int i = 0; i < m; ++i) {
            std::cout << "  ";
        }
    }

    void buildGluingNodes_() {
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            nodes_.insert({triangle, triangle});
        }
    }

    void buildGluingEdges_() {
        for (auto &[triangle, node] : nodes_) {
            std::unordered_set<int> selfGluingEdges;
            for (int i = 0; i < 3; ++i) {
                regina::Edge<dim> *edge = triangle->edge(i);
                regina::Perm<dim + 1> edgeToTriangle = triangle->edgeMapping(i);

                for (auto &[other, otherNode] : nodes_) {
                    for (int j = 0; j < 3; ++j) {
                        // Only glue identified edges and ignore identity
                        // mappings/opposite self gluings
                        if (other->edge(j) != edge ||
                            (triangle == other &&
                             (i == j || selfGluingEdges.find(i) !=
                                            selfGluingEdges.end())))
                            continue;

                        if (triangle == other) {
                            selfGluingEdges.insert(j);
                        }

                        regina::Perm<dim + 1> edgeToOther =
                            other->edgeMapping(j);
                        // regina::Perm is immut able, use std::array instead
                        std::array<int, 3> p;
                        for (int k = 0; k < 3; ++k) {
                            p[edgeToTriangle[k]] = edgeToOther[k];
                        }

                        node.adjList[&otherNode] = {triangle, i, other, p};
                    }
                }
            }
        }
    }

    void dfs_(GluingNode<dim> *node, int layer) {
        ++calls_;
        if (calls_ % 100000 == 0)
            std::cout << "Call " << calls_ << ", Layer " << layer
                      << ", Surfaces " << surfaces_.size() << "\n";

        if (node->visited || !node->valid) {
            return;
        }

        // nspaces_(layer);
        // std::cout << layer << ": ";
        if (!surface_.addTriangle(node->f, node->adjList)) {
            // std::cout << "Failed to add " << *node << "\n";
            return;
        }

        // std::cout << "Added " << *node << "\n";

        node->visited = true;

        if (surface_.isProper()) {
            // nspaces_(layer);
            // std::cout << "IMPROPER EDGES: { ";
            // for (const regina::Edge<dim> *edge : surface_.improperEdges_) {
            //    std::cout << edge->index() << " ";
            //}
            // std::cout << "}\n";
            surfaces_.insert(surface_);
        }

        for (auto &[nextNode, g] : node->adjList) {
            ////nspaces_(layer);
            // std::cout << layer << ": Gluing " << node->f->index() << " to "
            //           << nextNode->f->index() << " along (" << g.srcFacet
            //           << ", " << g.gluing[g.srcFacet] << ")\n";
            dfs_(nextNode, layer + 1);
        }

        surface_.removeTriangle(node->f);
    }

    void reset_() {
        surface_ = {&tri_};

        for (auto &[_, node] : nodes_) {
            node.valid = true;
            node.visited = false;
        }
    }
};

void usage(const char *progName, const std::string &error = std::string()) {
    if (!error.empty()) std::cerr << error << "\n\n";

    std::cerr << "Usage:\n";
    std::cerr << "    " << progName
              << " { -a, --all | -b, --boundary | -c, --closed | -l, --links } "
                 "<isosig>\n"
                 "    "
              << progName << " [ -v, --version | -h, --help ]\n\n";
    std::cerr << "    -a, --all      : Find all surfaces, regardless of "
                 "boundary conditions\n";
    std::cerr
        << "    -b, --boundary : Find surfaces such that their boundary is "
           "contained entirely in\n"
           "                     the boundary of the given 4-manifold (Note, "
           "if the 4-manifold\n                     is closed, this is "
           "equivalent to --closed)\n";
    std::cerr << "    -c, --closed   : Find only closed surfaces in the given "
                 "4-manifold\n";
    std::cerr << "    -l, --links    : Same as --boundary, but also gives "
                 "a list of link types which\n"
                 "                     are the boundaries of the surfaces "
                 "we find\n\n";
    std::cerr << "    -v, --version  : Show which version of Regina is "
                 "being used\n";
    std::cerr << "    -h, --help     : Display this help\n";
    exit(1);
}

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        if (strcmp(argv[i], "-h") == 0 || strcmp(argv[i], "--help") == 0)
            usage(argv[0]);
        if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
            if (argc != 2)
                usage(argv[0],
                      "Option --version cannot be used with "
                      "any other arguments.");
            std::cout << PACKAGE_BUILD_STRING << "\n";
            exit(0);
        }
    }

    if (argc != 3) {
        usage(argv[0],
              "Please specify a surface condition and provide an "
              "isomorphism signature.");
    }

    std::string isoSig = argv[2];
    SurfaceCondition cond;

    if (strcmp(argv[1], "-a") == 0 || strcmp(argv[1], "--all") == 0) {
        cond = SurfaceCondition::all;
    } else if (strcmp(argv[1], "-b") == 0 ||
               strcmp(argv[1], "--boundary") == 0 ||
               strcmp(argv[1], "-l") == 0 || strcmp(argv[1], "--links") == 0) {
        cond = SurfaceCondition::boundary;
    } else if (strcmp(argv[1], "-c") == 0 || strcmp(argv[1], "--closed") == 0) {
        cond = SurfaceCondition::closed;
    } else {
        usage(argv[0], "Please specify a valid surface condition.");
    }

    // regina::Triangulation<4> tri(isoSig);
    regina::Triangulation<4> tri;
    //regina::Triangulation<2> tri = regina::Example<2>::orientable(5, 1);
    tri.newSimplex();
    tri.subdivide();
    //tri.subdivide();


    std::cout << tri.detail() << "\n";
    // Knot<3> k;

    // knotComplement(tri, k);

    std::cout << "--- "
              << (cond == SurfaceCondition::all
                      ? ""
                      : (cond == SurfaceCondition::boundary
                             ? "Proper "
                             : (cond == SurfaceCondition::closed ? "Closed "
                                                                 : "")))
              << "Surfaces ---\n";
    int surfaceCount = 0;
    int closedCount = 0;
    GluingGraph graph(tri, cond);

    std::cout << "Boundary = ";
    for (const auto comp : tri.boundaryComponents()) {
        for (const auto edge : comp->edges()) {
            std::cout << edge->index() << " ";
        }
    }
    std::cout << "\n\n";

    auto surfaces = graph.findSurfaces();
    std::map<std::string, int> countMap;

    for (const auto &surface : surfaces) {
        ++surfaceCount;
        if (surface.surface().isClosed()) {
            ++closedCount;
        }

        bool isProper = true;
        for (const regina::BoundaryComponent<2> *comp :
             surface.surface().boundaryComponents()) {
            for (const regina::Edge<2> *edge : comp->edges()) {
                if (edge->isBoundary() && !surface.image(edge)->isBoundary()) {
                    std::cout << "NOTE: " << surface.detail()
                              << " is not proper!!\n";
                    isProper = false;
                }
            }
        }

        if (isProper) {
            std::cout << surface.detail() << " is PROPER!"
                      << "\n";
        }

        auto search = countMap.find(surface.detail());
        if (search == countMap.end()) {
            countMap.insert({surface.detail(), 1});
        } else {
            ++search->second;
        }
    }

    std::cout << "\n";

    std::vector<std::string> descriptions;
    descriptions.reserve(countMap.size());
    for (const auto &[description, count] : countMap) {
        descriptions.push_back(description);
    }

    std::sort(descriptions.begin(), descriptions.end(),
              [](const std::string &first, const std::string &second) {
                  return first.size() < second.size() ||
                         (first.size() == second.size() && first < second);
              });

    for (std::string &description : descriptions) {
        int count = countMap.find(description)->second;
        std::cout << "- " << count << (count == 1 ? " copy of " : " copies of ")
                  << description << "\n";
    }

    std::cout << "\n";
    if (cond != SurfaceCondition::closed) {
        std::cout << "Total admissible surfaces = " << surfaceCount << "\n";
    }
    std::cout << "Total closed surfaces = " << closedCount << "\n\n";

    return 0;
}
