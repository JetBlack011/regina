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

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <iterator>
#include <map>
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "maths/perm.h"
#include "maths/vector.h"
#include "triangulation/example2.h"
#include "triangulation/facenumbering.h"
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

template <int dim, int subdim>
struct Gluing {
    regina::Face<dim, subdim> *src;
    int srcFacet;
    regina::Face<dim, subdim> *dst;
    regina::Perm<subdim + 1> gluing;

    Gluing() = default;

    Gluing(regina::Face<dim, subdim> *src, int srcFacet,
           regina::Face<dim, subdim> *dst, regina::Perm<subdim + 1> gluing)
        : src(src), srcFacet(srcFacet), dst(dst), gluing(gluing) {}

    friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
        return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
                  << g.dst->index() << ", " << g.gluing << ")";
    }
};

/** Knot implementation */

class Knot {
   public:
    regina::Triangulation<3> tri_;
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>>
        tetEdges_;
    size_t numEdges_ = 0;

    bool isUnknot_ = false;

    enum JoinType { face, edge, vertex };

   public:
    Knot(regina::Triangulation<3> &tri, const Curve<3> &edges) : tri_(tri) {
        for (const regina::Edge<3> *edge : edges) {
            addEdge_(tri_.edge(edge->index()));
        }

        if (edges.size() <= 4) isUnknot_ = true;
    }

    bool isUnknot() const { return isUnknot_; }

    /**
     * Simplifies a given knot so that each tetrahedron it runs through contains
     * at most one of the knot's edges.
     *
     * Algorithm:
     * 0) If we're the unknot (e.g. all edges run through a single tetrahedron)
     * or if every tetrahedron contains exactly one edge of the knot, then we're
     * done.
     *
     * 1) Perform shrinkage. That is, suppose a tetrahedron contains more than
     * one edge. If those edges are connected end to end, we can isotope that
     * edge to go from its start point to its end point.
     *
     * 2) Separate the remaining edges. After step 1), a tetrahedron containing
     * an edge of the knot contains either a single edge or two disconnected
     * edges. In the second case, perform a 1-4 move to separate the edges.
     */
    void simplify() {
        if (isUnknot_) return;
        shrink_();
        separate_();
    }

    //void map(const std::function<void(regina::Tetrahedron<3> *,
    //                                  regina::Edge<3> *)> &f) {
    //    map([f](regina::Tetrahedron<3> *currTet, regina::Edge<3> *currEdge,
    //            regina::Tetrahedron<3> *, regina::Edge<3> *,
    //            JoinType) { f(currTet, currEdge); });
    //}

    //void map(const std::function<
    //         void(regina::Tetrahedron<3> *, regina::Edge<3> *,
    //              regina::Tetrahedron<3> *, regina::Edge<3> *, JoinType)> &f) {
    //    auto it = tetEdges_.begin();
    //    regina::Tetrahedron<3> *firstTet = it->first;
    //    regina::Tetrahedron<3> *prevTet = nullptr;
    //    regina::Tetrahedron<3> *currTet = firstTet;
    //    regina::Edge<3> *currEdge = currTet->edge(*it->second.begin());

    //    do {
    //        // TODO: Optimization, fix later
    //        // bool isFaceJoin = false;

    //        // for (int i = 0; i < 4; ++i) {
    //        //     // If there is no adjacent tetrahedron along this face or if
    //        //     // it's the previous tetrahedron, continue
    //        //     regina::Tetrahedron<3> *nextTet =
    //        //     currTet->adjacentSimplex(i); if (nextTet == nullptr ||
    //        //     nextTet == prevTet) continue;

    //        //    // If the adjacent tetrahedron has no edges, continue
    //        //    auto search = tetEdges_.find(nextTet);
    //        //    if (search == tetEdges_.end()) continue;

    //        //    // Otherwise, check if the edge inside the adjacent
    //        //    tetrahedron
    //        //    // is a continuation of the edge in the current tetrahedron
    //        //    ---
    //        //    // if so, this is our new current tetrahedron
    //        //    regina::Edge<3> *nextEdge =
    //        //        nextTet->edge(*search->second.begin());

    //        //    if (hasSharedVertex_(nextEdge, currEdge)) {
    //        //        f(currTet, currEdge, nextTet, nextEdge, JoinType::face);
    //        //        isFaceJoin = true;
    //        //        prevTet = currTet;
    //        //        currTet = nextTet;
    //        //        currEdge = nextEdge;
    //        //        break;
    //        //    }
    //        //}

    //        // if (isFaceJoin) continue;

    //        for (auto [nextTet, edge] : tetEdges_) {
    //            regina::Edge<3> *nextEdge = nextTet->edge(*edge.begin());

    //            if (nextTet == prevTet || !hasSharedVertex_(nextEdge, currEdge))
    //                continue;

    //            JoinType joinType;

    //            if (sharedFace_(nextTet, currTet) != nullptr)
    //                joinType = JoinType::face;
    //            else if (hasSharedEdge_(nextTet, currTet))
    //                joinType = JoinType::edge;
    //            else
    //                joinType = JoinType::vertex;

    //            f(currTet, currEdge, nextTet, nextEdge, joinType);
    //        }
    //    } while (currTet != firstTet);
    //}

    //void resolveSingularities() {
    //    if (isUnknot_) return;

    //    map([this](regina::Tetrahedron<3> *currTet, regina::Edge<3> *currEdge,
    //           regina::Tetrahedron<3> *nextTet, regina::Edge<3> *nextEdge,
    //           JoinType joinType) {
    //            if (joinType == JoinType::edge) {
    //                this->resolveEdgeSingularity_(currTet, nextTet);
    //            } else if (joinType == JoinType::vertex) {
    //            this->resolveVertexSingularity_(currTet, nextTet);
    //            }
    //    });

    //    return ans;
    //}

    void resolveEdgeSingularity_(regina::Tetrahedron<3> *tet1, regina::Tetrahedron<3> *tet2) {
        std::vector<regina::Tetrahedron<3> *> sharedEdgeSeq;

        
    }

    void resolveVertexSingularity_(regina::Vertex<3> *tet1, regina::Vertex<3> *tet2) {

    }

    /**
     * Pushes off the knot into the dual 1-skeleton of the given triangulation.
     *
     * Algorithm:
     *
     * 1) First, we need to isotope the knot so that any contiguous edges e1, e2
     * belong to tetrahedra that are joined along a common face. There are 2
     * cases where this is not true a priori:
     *
     *   i. e1 and e2 only belong to tetrahedra joined along a common edge. In
     * this case, there exists a sequence of tetrahedra {T_i} sharing this edge
     * such that T_0 contains e1, T_1 contains e2, and T_i shares a face with
     * T_i+1. Subdividing if necessary, we can isotope the knot to around this
     * shared edge.
     *
     *   ii. e1 and e2 only belong to tetrahedra joined along a common vertex.
     * Once again, there exists a sequence of tetrahedra {T_i} as above, where
     * now each tetrahedron shares this common vertex. Subdividing as necessary,
     * we can repeat the shared edge algorithm above and build a "bridge"
     * between shared edges.
     *
     * 2) Once every knot satisfies the above, iterate over these shared faces
     * and add them to our dual curve
     */
    //DualCurve<3> dual() {
    //    simplify();
    //    resolveSingularities();
    //    return getDualCurve_();
    //}

    /**
     * Computes a triangulation of the complement of this knot in the given
     * triangulation.
     *
     * Algorithm:
     * 1) If we're the unknot, simply replace a tetrahedron in the triangulation
     * with a known drilling of the unknot. 2)
     */
    // regina::Triangulation<3> complement() {}

    friend std::ostream &operator<<(std::ostream &os, const Knot &k) {
        if (k.isUnknot_) {
            return os << "Unknot";
        }

        std::unordered_set<const regina::Edge<3> *> edges;
        for (const auto &[tet, tetEdges] : k.tetEdges_) {
            for (size_t e : tetEdges) {
                edges.insert(tet->edge(e));
            }
        }

        os << "{ ";
        for (const regina::Edge<3> *e : edges) {
            os << "(" << e->vertex(0)->index() << ", " << e->vertex(1)->index()
               << ") ";
        }
        return os << "}";

        //// wtf...
        // if (k.isUnknot_) {
        //     return os << "Unknot";
        // }
        // os << "{";
        // auto it = k.edges_.begin();
        // int next;
        // for (auto it = k.edges_.begin(); it != k.edges_.end() - 1; ++it) {
        //     int v0 = (*it)->vertex(0)->index();
        //     int v1 = (*it)->vertex(1)->index();
        //     int u0 = (*(it + 1))->vertex(0)->index();
        //     int u1 = (*(it + 1))->vertex(1)->index();
        //     if (v0 != u0 && v0 != u1) {
        //         os << v0 << ", ";
        //         next = v1 == u0 ? u0 : u1;
        //     } else {
        //         os << v1 << ", ";
        //         next = v0 == u0 ? u0 : u1;
        //     }
        // }

        // return os << next << "}";
    }

   private:
    bool hasSharedVertex_(const regina::Edge<3> *e1,
                          const regina::Edge<3> *e2) {
        return e1->vertex(0) == e2->vertex(0) ||
               e1->vertex(0) == e2->vertex(1) ||
               e1->vertex(1) == e2->vertex(0) || e1->vertex(1) == e2->vertex(1);
    }

    regina::Triangle<3> *sharedFace_(const regina::Tetrahedron<3> *tet1,
                                     const regina::Tetrahedron<3> *tet2) {
        for (int i = 0; i < 4; ++i) {
            if (tet1->adjacentSimplex(i) == tet2) return tet1->triangle(i);
        }

        return nullptr;
    }

    bool hasSharedEdge_(const regina::Tetrahedron<3> *tet1,
                        const regina::Tetrahedron<3> *tet2) {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                if (tet1->edge(i) == tet2->edge(j)) return true;
            }
        }

        return false;
    }

    bool hasSharedVertex_(const regina::Tetrahedron<3> *tet1,
                          const regina::Tetrahedron<3> *tet2) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (tet1->vertex(i) == tet2->vertex(j)) return true;
            }
        }

        return false;
    }

    //DualCurve<3> getDualCurve_() {
    //    DualCurve<3> ans;

    //    map([&ans, this](regina::Tetrahedron<3> *currTet, regina::Edge<3> *,
    //                     regina::Tetrahedron<3> *nextTet, regina::Edge<3> *,
    //                     JoinType joinType) {
    //        if (joinType != JoinType::face)
    //            throw regina::InvalidArgument(
    //                "All adjacent edges must be joined along faces!");
    //        ans.push_back(this->sharedFace_(currTet, nextTet));
    //    });

    //    return ans;
    //}

    //Curve<3> getCurve_() {
    //    Curve<3> ans;

    //    map([&ans](regina::Tetrahedron<3> *, regina::Edge<3> *edge) {
    //        ans.push_back(edge);
    //    });

    //    return ans;
    //}

    bool updateIsUnknot_() {
        if (tetEdges_.size() <= 1) {
            return isUnknot_ = true;
        }

        bool isAllInOneTet = false;

        for (const auto &[_, edges] : tetEdges_) {
            if (edges.size() == numEdges_) {
                return isUnknot_ = true;
            }
        }

        return false;
    }

    void shrink_() {
        bool isDone = false;
        while (!updateIsUnknot_() && !isDone) {
            isDone = true;
            for (const auto &[tet, _] : tetEdges_) {
                // If we find a tetrahedron we can shrink, we need to keep going
                if (shrink_(tet)) {
                    isDone = false;
                    // Since at least one tetrahedron's edge information has
                    // been updated, we must start over
                    break;
                }
            }
        }
    }

    void separate_() {
        bool isDone = false;
        while (!updateIsUnknot_() && !isDone) {
            // TODO: Obvious optimizations
            isDone = true;
            for (auto &[tet, edges] : tetEdges_) {
                assert(edges.size() <= 2);

                if (edges.size() == 2) {
                    fourOneMove_(tet);
                    isDone = false;
                    break;
                }
            }
        }
    }

    bool shrink_(regina::Tetrahedron<3> *tet) {
        std::unordered_set<size_t> edges = tetEdges_.at(tet);

        std::unordered_map<size_t, size_t> vertCounts;
        for (size_t e : edges) {
            regina::Perm<4> p = tet->edgeMapping(e);
            ++vertCounts[p[0]];
            ++vertCounts[p[1]];
        }

        if (edges.size() < 2 || (edges.size() == 2 && vertCounts.size() == 4))
            return false;

        for (size_t e : edges) {
            removeEdge_(tet->edge(e));
        }

        /* Add the new edge */
        // In this case, there are exactly two endpoints in this
        // tet
        std::array<size_t, 2> endVerts;
        size_t i = 0;
        for (auto [vert, count] : vertCounts) {
            if (count == 1) endVerts[i++] = vert;
            if (i >= 2) break;
        }

        addEdge_(tet->edge(endVerts[0], endVerts[1]));

        return true;
    }

    void addEdge_(const regina::Edge<3> *edge) {
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            tetEdges_[emb.tetrahedron()].insert(emb.face());
        }

        ++numEdges_;
    }

    void removeEdge_(const regina::Edge<3> *edge) {
        for (const regina::EdgeEmbedding<3> &emb : edge->embeddings()) {
            regina::Tetrahedron<3> *tet = emb.tetrahedron();
            tetEdges_[tet].erase(emb.face());
            if (tetEdges_[tet].empty()) {
                tetEdges_.erase(tet);
            }
        }

        --numEdges_;
    }

    std::array<regina::Tetrahedron<3> *, 4> fourOneTets_() {
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

    void fourOneMove_(regina::Tetrahedron<3> *tet) {
        std::array<std::optional<Gluing<3, 3>>, 4> gluings;

        // For each face of tet, track gluings
        for (int i = 0; i < 4; ++i) {
            regina::Tetrahedron<3> *adjTet = tet->adjacentSimplex(i);
            if (adjTet != nullptr) {
                gluings[i] = {tet, i, adjTet, tet->adjacentGluing(i)};
            }
        }

        std::array<regina::Tetrahedron<3> *, 4> newTets = fourOneTets_();
        // Add the edges contained in tet to their appropriate position in
        // newTets
        for (size_t e : tetEdges_.at(tet)) {
            regina::Perm<4> p = tet->edgeMapping(e);
            tetEdges_[newTets[p[2]]].insert(e);
            tetEdges_[newTets[p[3]]].insert(e);
        }

        tetEdges_.erase(tet);
        tri_.removeTetrahedron(tet);

        // Glue in the new tetrahedra
        for (int i = 0; i < 4; ++i) {
            if (auto g = gluings[i]) {
                newTets[i]->join(i, g->dst, g->gluing);
            }
        }
    }

    std::array<std::unordered_map<int, regina::Tetrahedron<3> *>, 2>
    twoSixTets_(Gluing<3, 3> g) {
        std::array<std::unordered_map<int, regina::Tetrahedron<3> *>, 2> tets;

        for (int i = 0; i < 4; ++i) {
            if (i != g.srcFacet) {
                tets[0][i] = tri_.newTetrahedron();
            }
        }
        for (int i = 0; i < 4; ++i) {
            if (i != g.gluing[g.srcFacet]) {
                tets[1][i] = tri_.newTetrahedron();
            }
        }

        // Inner gluings
        for (int i = 0; i < 2; ++i) {
            int gluingFace = i == 0 ? g.srcFacet : g.gluing[g.srcFacet];
            for (int j = 0; j < 4; ++j) {
                if (j == gluingFace) continue;
                for (int k = j + 1; k < 4; ++k) {
                    if (k == gluingFace) continue;
                    tets[i][j]->join(k, tets[i][k], {j, k});
                }
            }
        }

        // Gluings between tet0 and tet1
        for (int i = 0; i < 4; ++i) {
            if (i != g.srcFacet) {
                tets[0][i]->join(g.srcFacet, tets[1][g.gluing[i]], g.gluing);
            }
        }

        return tets;
    }

    std::array<std::unordered_map<int, regina::Tetrahedron<3> *>, 2>
    twoSixMove_(regina::Tetrahedron<3> *tet0, regina::Tetrahedron<3> *tet1) {
        std::array<regina::Tetrahedron<3> *, 2> tets = {tet0, tet1};
        auto tetGluing = getGluing_(tet0, tet1);
        std::array<std::array<std::optional<Gluing<3, 3>>, 4>, 2> gluings;

        // For each face of tet, track gluings
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (j == (i == 0 ? tetGluing.srcFacet
                                 : tetGluing.gluing[tetGluing.srcFacet]))
                    continue;
                regina::Tetrahedron<3> *adjTet = tets[i]->adjacentSimplex(j);
                if (adjTet != nullptr)
                    gluings[i][j] = {tets[i], j, adjTet,
                                     tets[i]->adjacentGluing(j)};
            }
        }

        auto newTets = twoSixTets_(tetGluing);
        std::cout << "New tets: " << tri_.isoSig() << "\n";
        // Add the edges contained in tet to their appropriate position in tets
        for (int i = 0; i < 2; ++i) {
            if (tetEdges_.find(tets[i]) != tetEdges_.end()) {
                for (size_t e : tetEdges_.at(tets[i])) {
                    regina::Perm<4> p = tets[i]->edgeMapping(e);
                    if (p[2] != tetGluing.srcFacet)
                        tetEdges_[newTets[i][p[2]]].insert(e);
                    if (p[3] != tetGluing.srcFacet)
                        tetEdges_[newTets[i][p[3]]].insert(e);
                }
                tetEdges_.erase(tets[i]);
            }
        }
        std::cout << "Before gluing: " << tri_.isoSig() << "\n";

        tri_.removeTetrahedron(tet0);
        tri_.removeTetrahedron(tet1);

        // Glue in the new tetrahedra
        for (int i = 0; i < 2; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (auto g = gluings[i][j]) {
                    newTets[i][j]->join(j, g->dst, g->gluing);
                }
            }
        }
        std::cout << "After gluing: " << tri_.isoSig() << "\n";

        std::cout << "\n";
        return newTets;
    }

   public:
    Gluing<3, 3> getGluing_(regina::Tetrahedron<3> *tet0,
                            regina::Tetrahedron<3> *tet1) const {
        Gluing<3, 3> g;
        for (int i = 0; i < 4; ++i) {
            if (tet0->adjacentSimplex(i) == tet1)
                return {tet0, i, tet1, tet0->adjacentGluing(i)};
        }

        throw regina::InvalidArgument(
            "Given tetrahedra are not glued together!");
    }

    void subdivideSharedEdgeSequence_(
        std::vector<regina::Tetrahedron<3> *> &tets) {
        regina::Perm<4> p;
        for (const regina::EdgeEmbedding<3> &emb :
             sharedEdge_(tets)->embeddings()) {
            if (emb.tetrahedron() == tets[0]) {
                p = emb.vertices();
                break;
            }
        }

        for (int i = 0; i < tets.size() - 1; ++i) {
            auto g = getGluing_(tets[i], tets[i + 1]);
            int sharedFace = p[2] != g.srcFacet ? p[2] : p[3];

            auto newTets = twoSixMove_(tets[i], tets[i + 1]);

            tets[i] = newTets[0][sharedFace];
            tets[i + 1] = newTets[0][sharedFace]->adjacentSimplex(g.srcFacet);

            p = tets[i + 1]->edgeMapping(
                regina::FaceNumbering<3, 1>::edgeNumber[g.gluing[p[0]]]
                                                       [g.gluing[p[1]]]);
        }
    }

    bool isGlued_(const regina::Tetrahedron<3> *tet0,
                  const regina::Tetrahedron<3> *tet1) const {
        for (int i = 0; i < 4; ++i) {
            if (tet0->adjacentSimplex(i) == tet1) return true;
        }
        return false;
    }

    regina::Edge<3> *sharedEdge_(
        const std::vector<regina::Tetrahedron<3> *> &tets) {
        auto ts = lastThree_(tets);
        regina::Triangle<3> *t0 =
            ts[0]->triangle(getGluing_(ts[0], ts[1]).srcFacet);
        regina::Triangle<3> *t1 =
            ts[1]->triangle(getGluing_(ts[1], ts[2]).srcFacet);
        std::unordered_set<regina::Edge<3> *> edges = {t0->edge(0), t0->edge(1),
                                                       t0->edge(2)};

        for (int i = 0; i < 3; ++i) {
            if (edges.find(t1->edge(i)) != edges.end()) return t1->edge(i);
        }

        throw regina::InvalidArgument(
            "Given tetrahedra don't have a shared edge!");
    }

    std::array<regina::Tetrahedron<3> *, 3> lastThree_(
        const std::vector<regina::Tetrahedron<3> *> &tets) {
        std::array<regina::Tetrahedron<3> *, 3> ans;
        int j = 0;
        for (int i = tets.size() - 3; i < tets.size(); ++i) {
            ans[j++] = tets[i];
        }
        return ans;
    }

    regina::Vertex<3> *sharedVertex_(
        const std::vector<regina::Tetrahedron<3> *> &tets) const {
        std::set<regina::Vertex<3> *> verts;
        for (int j = 0; j < 4; ++j) {
            verts.insert(tets[0]->vertex(j));
        }

        for (int i = 1; i < tets.size(); ++i) {
            std::set<regina::Vertex<3> *> newVerts;
            for (int j = 0; j < 4; ++j) {
                newVerts.insert(tets[i]->vertex(j));
            }

            std::set<regina::Vertex<3> *> intersection;
            std::set_intersection(
                verts.begin(), verts.end(), newVerts.begin(), newVerts.end(),
                std::inserter(intersection, intersection.begin()));
            verts = intersection;

            if (verts.size() == 1) return *(verts.begin());
        }

        if (verts.size() >= 1) return *(verts.begin());

        throw regina::InvalidArgument(
            "Given tetrahedra do not share a common vertex!");
    }

    void subdivideSharedVertexSequence_(
        std::vector<regina::Tetrahedron<3> *> &tets) {
        const regina::Vertex<3> *sharedVertex = sharedVertex_(tets);
        int v;
        for (const regina::VertexEmbedding<3> &emb :
             sharedVertex->embeddings()) {
            if (emb.simplex() == tets[0]) {
                v = emb.vertex();
                break;
            }
        }

        for (int i = 0; i < tets.size() - 1;) {
            int start = i;
            std::vector<regina::Tetrahedron<3> *> sharedEdgeSubseq = {
                tets[i], tets[++i], tets[++i]};
            regina::Edge<3> *sharedEdge = sharedEdge_(sharedEdgeSubseq);

            while (i < tets.size() - 1 &&
                   sharedEdge ==
                       sharedEdge_({tets[i - 1], tets[i], tets[i + 1]})) {
                sharedEdgeSubseq.push_back(tets[++i]);
            }

            // std::cout << "ISO = " << tri_.isoSig() << "\n";
            subdivideSharedEdgeSequence_(sharedEdgeSubseq);
            // std::cout << "SUBDIVIDED = " << tri_.isoSig() << "\n";
            // for (regina::Tetrahedron<3> *tet : sharedEdgeSubseq) {
            //     std::cout << "Tet " << tet->index() << ": ";
            //     for (int i = 0; i < 4; ++i) {
            //         std::cout << tet->vertex(i)->index() << " ";
            //     }
            //     std::cout << "\n";
            // }

            for (int j = start; j <= i; ++j) {
                tets[j] = sharedEdgeSubseq[j - start];
            }

            if (i >= tets.size() - 1) return;

            sharedEdge = sharedEdge_(sharedEdgeSubseq);
            sharedVertex = tets[0]->vertex(v);
            regina::Vertex<3> *otherVertex =
                sharedEdge->vertex(0) != sharedVertex ? sharedEdge->vertex(0)
                                                      : sharedEdge->vertex(1);

            for (const regina::VertexEmbedding<3> &emb :
                 otherVertex->embeddings()) {
                if (emb.simplex() == tets[i]) {
                    tets[i] = tets[i]->adjacentSimplex(emb.vertex());
                    break;
                }
            }
            for (const regina::VertexEmbedding<3> &emb :
                 otherVertex->embeddings()) {
                if (emb.simplex() == tets[i - 1]) {
                    tets.insert(tets.begin() + i,
                                tets[i - 1]->adjacentSimplex(emb.vertex()));
                    break;
                }
            }
        }
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
    std::vector<Knot> comps_;

   public:
    Link(regina::Triangulation<3> &tri, const std::vector<Curve<3>> &components)
        : tri_(tri) {
        for (const Curve<3> &edges : components) {
            comps_.emplace_back(tri, edges);
        }
    }

    bool isUnlink() const {
        for (const Knot &k : comps_) {
            if (!k.isUnknot()) return false;
        }
        return true;
    }

    void simplify() {
        for (Knot &k : comps_) {
            k.simplify();
            break;
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

    friend std::ostream &operator<<(std::ostream &os, const Link &l) {
        for (const Knot &k : l.comps_) {
            os << k << "\n";
        }
        return os;
    }
};

template <int dim>
class GluingNode {
   public:
    using AdjList = std::unordered_map<GluingNode *, Gluing<dim, 2>>;

    regina::Triangle<dim> *f;
    AdjList adjList;
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
    regina::Triangulation<dim - 1> bdry_;
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
    KnottedSurface(const regina::Triangulation<dim> *tri) : tri_(tri) {
        if (!tri_->isClosed()) bdry_ = tri_->boundaryComponent(0)->build();
    }

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_),
          surface_(other.surface_),
          improperEdges_(other.improperEdges_),
          invariants_(other.invariants_),
          indices_(other.indices_) {
        if (!tri_->isClosed()) bdry_ = tri_->boundaryComponent(0)->build();
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

    bool isProper() const {
        return surface_.isClosed() || improperEdges_.empty();
    }

    bool hasSelfIntersection() const {
        std::unordered_set<const regina::Vertex<dim> *> images;
        for (const regina::Vertex<2> *v : surface_.vertices()) {
            const regina::Vertex<dim> *im = image(v);
            if (images.find(im) != images.end()) return true;
            images.insert(im);
        }
        return false;
    }

    Link boundary() {
        const regina::BoundaryComponent<dim> *b = tri_->boundaryComponent(0);
        std::unordered_map<const regina::Edge<dim> *,
                           const regina::Edge<dim - 1> *>
            bdryMap;
        std::unordered_map<const regina::Edge<2> *,
                           const regina::Edge<dim - 1> *>
            edgeMap;

        for (const regina::Edge<dim - 1> *e : bdry_.edges()) {
            bdryMap[b->edge(e->index())] = e;
        }

        for (const regina::BoundaryComponent<2> *comp :
             surface_.boundaryComponents()) {
            for (const regina::Edge<2> *e : comp->edges()) {
                edgeMap[e] = bdryMap.at(image(e));
            }
        }

        std::vector<Curve<dim - 1>> comps;

        // TODO: Likely lots of optimizations can be made here
        // for (const regina::BoundaryComponent<2> *comp :
        //     surface_.boundaryComponents()) {
        //    const regina::Edge<2> *first = comp->edge(0);
        //    const regina::Edge<2> *next;
        //    std::unordered_set<const regina::Vertex<2> *> verts = {
        //        first->vertex(0), first->vertex(1)};

        //    Curve<dim - 1> edges = {edgeMap.at(first)};

        //    while (edges.size() < comp->countEdges()) {
        //        for (const regina::Edge<2> *e : comp->edges()) {
        //            if (std::find(edges.begin(), edges.end(), edgeMap.at(e))
        //            ==
        //                    edges.end() &&
        //                (verts.find(e->vertex(0)) != verts.end() ||
        //                 verts.find(e->vertex(1)) != verts.end())) {
        //                next = e;
        //                break;
        //            }
        //        }

        //        verts = {next->vertex(0), next->vertex(1)};
        //        edges.push_back(edgeMap.at(next));
        //    }

        //    comps.push_back(edges);
        //}

        for (const regina::BoundaryComponent<2> *comp :
             surface_.boundaryComponents()) {
            Curve<dim - 1> edges;
            //for (const regina::Edge<2> *e : comp->edges()) {
            //    edges.push_back(edgeMap.at(e));
            //}
            comps.push_back(edges);
        }

        return {bdry_, comps};
    }

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
        std::cout << "[+] Built gluing nodes\n";
        buildGluingEdges_();
        std::cout << "[+] Built gluing edges\n";
        // buildGluingInvalids_();
        // std::cout << "Built gluing invalids\n";
    }

    void pruneSurfaces() {
        // TODO: Do this during DFS
        for (auto it = surfaces_.begin(); it != surfaces_.end();) {
            if ((*it).hasSelfIntersection()) {
                it = surfaces_.erase(it);
                continue;
            }
            ++it;
        }
    }

    std::set<KnottedSurface<dim>> &findSurfaces() {
        if (!surfaces_.empty()) return surfaces_;

        std::cout << "\n";
        for (regina::Triangle<dim> *triangle : tri_.triangles()) {
            reset_();
            dfs_(&nodes_.at(triangle), 0);
            // std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";
        }
        std::cout << "\n\n";
        std::cout << "[+] Total search calls = " << calls_ << "\n\n";

        pruneSurfaces();

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

                // for (const regina::TriangleEmbedding<dim> &emb :
                //      triangle->embeddings()) {
                //     // TODO: GAH
                //     const regina::Simplex<dim> *s = emb.simplex();
                //     for (int j = 0; j < 10; ++j) {
                //         regina::Triangle<dim> *other = s->triangle(j);
                //         GluingNode<dim> &otherNode = nodes_.at(other);

                //        for (int k = 0; k < 3; ++k) {
                //            // Only glue identified edges and ignore identity
                //            // mappings/opposite self gluings
                //            if (other->edge(k) != edge ||
                //                (triangle == other &&
                //                 (i == k || selfGluingEdges.find(i) !=
                //                                selfGluingEdges.end())))
                //                continue;

                //            if (triangle == other) {
                //                selfGluingEdges.insert(k);
                //            }

                //            regina::Perm<dim + 1> edgeToOther =
                //                other->edgeMapping(k);
                //            // regina::Perm is immutable, use std::array
                //            instead std::array<int, 3> p; for (int k = 0; k <
                //            3; ++k) {
                //                p[edgeToTriangle[k]] = edgeToOther[k];
                //            }

                //            node.adjList[&otherNode] = {triangle, i, other,
                //            p};
                //        }
                //    }
                //}

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
                        // regina::Perm is immutable, use std::array instead
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

    // void buildGluingInvalids_() {
    //     for (auto &[f, node] : nodes_) {
    //         std::unordered_set<const regina::Triangle<dim> *> adjTriangles;
    //         std::unordered_set<const regina::Vertex<dim> *> vertices;

    //        for (const auto [adj, _] : node.adjList) {
    //            adjTriangles.insert(adj->f);
    //        }

    //        for (int i = 0; i < 3; ++i) {
    //            vertices.insert(f->vertex(i));
    //        }

    //        // TODO: Optimize
    //        for (auto &[adjTriangle, adjNode] : nodes_) {

    //            if (adjTriangles.find(adjTriangle) != adjTriangles.end())
    //                continue;

    //            bool containsVertex = false;
    //            for (int i = 0; i < 3; ++i) {
    //                if (vertices.find(adjTriangle->vertex(i)) !=
    //                vertices.end())
    //                    containsVertex = true;
    //            }

    //            if (containsVertex) node.invalids.push_back(&adjNode);
    //        }
    //    }
    //}

    void dfs_(GluingNode<dim> *node, int layer) {
        ++calls_;
        if (calls_ % 100000 == 0)
            std::cout << "\r" << std::flush << "[*] Call " << calls_
                      << ", Layer " << layer << ", Surfaces "
                      << surfaces_.size() << "          ";

        if (node->visited) {
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
            //  nspaces_(layer);
            //  std::cout << "IMPROPER EDGES: { ";
            //  for (const regina::Edge<dim> *edge : surface_.improperEdges_) {
            //     std::cout << edge->index() << " ";
            // }
            //  std::cout << "}\n";
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

template <int dim>
void surfacesDetail(std::set<KnottedSurface<dim>> &surfaces,
                    SurfaceCondition cond) {
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

    std::map<std::string, int> countMap;

    for (const auto &surface : surfaces) {
        ++surfaceCount;
        if (surface.surface().isClosed()) {
            ++closedCount;
        }

        // For checking if properness works
        // bool isProper = true;
        // for (const regina::BoundaryComponent<2> *comp :
        //     surface.surface().boundaryComponents()) {
        //    for (const regina::Edge<2> *edge : comp->edges()) {
        //        if (edge->isBoundary() && !surface.image(edge)->isBoundary())
        //        {
        //            std::cout << "NOTE: " << surface.detail()
        //                      << " is not proper!!\n";
        //            isProper = false;
        //        }
        //    }
        //}

        // if (isProper) {
        //     std::cout << surface.detail() << " is PROPER!"
        //               << "\n";
        // }

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
}

int main(int argc, char *argv[]) {
    // Check for standard arguments:
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "-h" || arg == "--help") usage(argv[0]);
        if (arg == "-v" || arg == "--version") {
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

    std::string arg = argv[1];
    std::string isoSig = argv[2];
    SurfaceCondition cond;

    if (arg == "-a" || arg == "--all") {
        cond = SurfaceCondition::all;
    } else if (arg == "-b" || arg == "--boundary" || arg == "-l" ||
               arg == "--links") {
        cond = SurfaceCondition::boundary;
    } else if (arg == "-c" || arg == "--closed") {
        cond = SurfaceCondition::closed;
    } else {
        usage(argv[0], "Please specify a valid surface condition.");
    }

    //    regina::Triangulation<4> tri;
    //    tri.newSimplex();
    //    tri.pachner(tri.pentachoron(0));
    //
    // regina::Triangulation<3> tri("caba");
    // for (const regina::Tetrahedron<3> *tet : tri.tetrahedra()) {
    //     std::cout << "Tetrahedron " << tet->index() << "\n";
    //     for (int i = 0; i < 4; ++i) {
    //         std::cout << "Vertex " << i << " = " << tet->vertex(i)->index()
    //                   << "\n";
    //     }
    //     std::cout << "\n";
    // }
    // std::vector<int> edgeIndices = {2, 0, 7, 8, 5};
    // Curve<3> c;
    // for (int i : edgeIndices) {
    //     c.push_back(tri.edge(i));
    // }
    // Knot k(tri, c);
    // std::cout << k << "\n";
    // k.simplify();
    // std::cout << k;

    regina::Triangulation<4> tri(isoSig);
    // regina::Triangulation<3> tri("eabbba");
    // regina::Triangulation<3> tri("gabbbbba");
    // regina::Triangulation<3> tri("fabbfaa");
    //regina::Triangulation<3> tri;
    //tri.newSimplex();
    //tri.newSimplex();
    //tri.newSimplex();
    //tri.newSimplex();
    //tri.newSimplex();
    //tri.newSimplex();
    //tri.newSimplex();
    //tri.tetrahedron(0)->join(0, tri.tetrahedron(1), {});
    //tri.tetrahedron(1)->join(1, tri.tetrahedron(2), {});
    //tri.tetrahedron(2)->join(0, tri.tetrahedron(3), {});
    //tri.tetrahedron(3)->join(1, tri.tetrahedron(4), {});
    //tri.tetrahedron(4)->join(2, tri.tetrahedron(5), {});
    //tri.tetrahedron(5)->join(1, tri.tetrahedron(6), {});

    //std::cout << "Original = " << tri.isoSig() << "\n";
    //Knot k = {tri, {}};
    //std::vector<regina::Tetrahedron<3> *> tets;
    //for (regina::Tetrahedron<3> *tet : k.tri_.tetrahedra()) {
    //    tets.push_back(tet);
    //}

    //k.subdivideSharedVertexSequence_(tets);
    //std::cout << "Subdivided = " << k.tri_.isoSig() << "\n";
    // tri.newSimplex();

    // regina::Triangulation<2> tri = regina::Example<2>::orientable(6, 0);
    // tri.newSimplex();
    // tri.subdivide();
    // tri.subdivide();
    //    tri.subdivide();

    std::cout << tri.detail() << "\n";

    GluingGraph graph(tri, cond);
    auto &surfaces = graph.findSurfaces();

    surfacesDetail(surfaces, cond);

    for (const auto *comp : tri.boundaryComponents()) {
        for (const auto *v : comp->vertices()) {
            std::cout << v->index() << " ";
        }
        std::cout << "\n";
    }
    std::cout << "\n";

    int numNonUnlinks = 0;
    for (auto surface : surfaces) {
        if (surface.surface().isClosed()) continue;
        Link l = surface.boundary();
        // std::cout << l << "\n";
        // l.simplify();
        if (!l.isUnlink()) {
            ++numNonUnlinks;
            std::cout << "DETAIL: " << surface.detail() << "\n";
            std::cout << l << "\n";
        }
    }

    std::cout << "NUMBER OF NONTRIVIAL LINKS = " << numNonUnlinks;

    return 0;
}
