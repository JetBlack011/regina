//
//  knottedsurfaces.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef KNOTTED_SURFACES_H

#define KNOTTED_SURFACES_H

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
#include <ostream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "maths/perm.h"
#include "triangulation/facenumbering.h"
#include "triangulation/forward.h"
#include "triangulation/generic/boundarycomponent.h"
#include "triangulation/generic/face.h"
#include "triangulation/generic/faceembedding.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

#include "gluing.h"

enum class SurfaceCondition { all, boundary, closed };

template <int n>
using Curve = std::vector<regina::Edge<n> *>;

template <int n>
using DualCurve = std::vector<regina::Triangle<n> *>;

/** Knot implementation */

class Knot {
  public:
    regina::Triangulation<3> tri_;
    std::unordered_map<regina::Tetrahedron<3> *, std::unordered_set<size_t>> tetEdges_;
    size_t numEdges_ = 0;

    bool isUnknot_ = false;

    enum class JoinType { face, edge, vertex };

  public:
    Knot(regina::Triangulation<3> &tri, const Curve<3> &edges) : tri_(tri) {
        for (const regina::Edge<3> *edge : edges) {
            addEdge_(tri_.edge(edge->index()));
        }

        if (edges.size() <= 4)
            isUnknot_ = true;
    }

    bool isUnknot() const { return isUnknot_; }

    /**
     * Simplifies a given knot so that each tetrahedron it runs through contains at most one of the
     * knot's edges.
     *
     * Algorithm:
     * 0) If we're the unknot (e.g. all edges run through a single tetrahedron) or if every
     * tetrahedron contains exactly one edge of the knot, then we're done.
     *
     * 1) Perform shrinkage. That is, suppose a tetrahedron contains more than one edge. If those
     * edges are connected end to end, we can isotope that edge to go from its start point to its
     * end point.
     *
     * 2) Separate the remaining edges. After step 1), a tetrahedron containing an edge of the knot
     * contains either a single edge or two disconnected edges. In the second case, perform a 1-4
     * move to separate the edges.
     */
    void simplify() {
        if (isUnknot_)
            return;
        shrink_();
        separate_();
    }

    // void map(const std::function<void(regina::Tetrahedron<3> *,
    //                                   regina::Edge<3> *)> &f) {
    //     map([f](regina::Tetrahedron<3> *currTet, regina::Edge<3> *currEdge,
    //             regina::Tetrahedron<3> *, regina::Edge<3> *,
    //             JoinType) { f(currTet, currEdge); });
    // }

    // void map(const std::function<
    //          void(regina::Tetrahedron<3> *, regina::Edge<3> *,
    //               regina::Tetrahedron<3> *, regina::Edge<3> *, JoinType)> &f) {
    //     auto it = tetEdges_.begin();
    //     regina::Tetrahedron<3> *firstTet = it->first;
    //     regina::Tetrahedron<3> *prevTet = nullptr;
    //     regina::Tetrahedron<3> *currTet = firstTet;
    //     regina::Edge<3> *currEdge = currTet->edge(*it->second.begin());

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

    // void resolveSingularities() {
    //     if (isUnknot_) return;

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

    void resolveVertexSingularity_(regina::Vertex<3> *tet1, regina::Vertex<3> *tet2) {}

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
    // DualCurve<3> dual() {
    //     simplify();
    //     resolveSingularities();
    //     return getDualCurve_();
    // }

    /**
     * Computes a triangulation of the complement of this knot in the given
     * triangulation.
     *
     * Algorithm:
     * 1) If we're the unknot, simply replace a tetrahedron in the triangulation
     * with a known drilling of the unknot. 2)
     */
    // regina::Triangulation<3> complement() {}

    /**
     * Prints the knot in the form { v1, v2, ..., vn } where v1, v2, ..., vn are the indices of the
     * vertices appearing in the knot
     */
    friend std::ostream &operator<<(std::ostream &os, const Knot &k) {
        //if (k.isUnknot_) {
        //    return os << "Unknot";
        //}

        std::unordered_set<const regina::Edge<3> *> edgeSet;
        for (const auto &[tet, tetEdges] : k.tetEdges_) {
            for (size_t e : tetEdges) {
                edgeSet.insert(tet->edge(e));
            }
        }

        std::vector<const regina::Vertex<3> *> vertices;
        const auto e = edgeSet.begin();
        vertices.push_back((*e)->vertex(0));
        auto currVert = (*e)->vertex(1);
        edgeSet.erase(e);

        while (!edgeSet.empty()) {
            for (const regina::Edge<3> *edge : edgeSet) {
                if (edge->vertex(0) == currVert) {
                    vertices.push_back(currVert);
                    currVert = edge->vertex(1);
                    edgeSet.erase(edge);
                    break;
                } else if (edge->vertex(1) == currVert) {
                    vertices.push_back(currVert);
                    currVert = edge->vertex(0);
                    edgeSet.erase(edge);
                    break;
                }
            }
        }

        os << "{ ";
        for (const regina::Vertex<3> *v : vertices) {
            os << v->index() << ", ";
        }
        return os << vertices[0]->index() << " }";
    }

  public:
    bool hasSharedVertex_(const regina::Edge<3> *e1, const regina::Edge<3> *e2) const {
        return e1->vertex(0) == e2->vertex(0) || e1->vertex(0) == e2->vertex(1) ||
               e1->vertex(1) == e2->vertex(0) || e1->vertex(1) == e2->vertex(1);
    }

    regina::Triangle<3> *sharedFace_(const regina::Tetrahedron<3> *tet1,
                                     const regina::Tetrahedron<3> *tet2) const {
        for (int i = 0; i < 4; ++i) {
            if (tet1->adjacentSimplex(i) == tet2)
                return tet1->triangle(i);
        }

        return nullptr;
    }

    bool hasSharedEdge_(const regina::Tetrahedron<3> *tet1, const regina::Tetrahedron<3> *tet2) {
        for (int i = 0; i < 6; ++i) {
            for (int j = 0; j < 6; ++j) {
                if (tet1->edge(i) == tet2->edge(j))
                    return true;
            }
        }

        return false;
    }

    bool hasSharedVertex_(const regina::Tetrahedron<3> *tet1, const regina::Tetrahedron<3> *tet2) {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                if (tet1->vertex(i) == tet2->vertex(j))
                    return true;
            }
        }

        return false;
    }

    // DualCurve<3> getDualCurve_() {
    //     DualCurve<3> ans;

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

    // Curve<3> getCurve_() {
    //     Curve<3> ans;

    //    map([&ans](regina::Tetrahedron<3> *, regina::Edge<3> *edge) {
    //        ans.push_back(edge);
    //    });

    //    return ans;
    //}

    bool updateIsUnknot_() {
        if (tetEdges_.size() <= 1) {
            return isUnknot_ = true;
        }

        for (const auto &[_, edges] : tetEdges_) {
            if (edges.size() == numEdges_) {
                return isUnknot_ = true;
            }
        }

        return false;
    }

    /**
     * Locally "cleans up" the knot by isotoping away redundant edges in every tetrahedron
     * containing an edge of the knot.
     */
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

    /**
     * In the case where there are two edges in a given tetrahedron which are disjoint, we separate
     * them by performing a 1-4 move. This should only be run after shrink_() has been called.
     */
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
            // TODO: Don't be stupid (impossible!)
            removeEdge_(tet->edge(e));
        }

        /* Add the new edge */
        // In this case, there are exactly two endpoints in this tet
        std::array<size_t, 2> endVerts;
        size_t i = 0;
        for (auto [vert, count] : vertCounts) {
            if (count == 1)
                endVerts[i++] = vert;
            if (i >= 2)
                break;
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

    std::array<std::unordered_map<int, regina::Tetrahedron<3> *>, 2> twoSixTets_(Gluing<3, 3> g) {
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
                if (j == gluingFace)
                    continue;
                for (int k = j + 1; k < 4; ++k) {
                    if (k == gluingFace)
                        continue;
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
                if (j == (i == 0 ? tetGluing.srcFacet : tetGluing.gluing[tetGluing.srcFacet]))
                    continue;
                regina::Tetrahedron<3> *adjTet = tets[i]->adjacentSimplex(j);
                if (adjTet != nullptr)
                    gluings[i][j] = {tets[i], j, adjTet, tets[i]->adjacentGluing(j)};
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
    Gluing<3, 3> getGluing_(regina::Tetrahedron<3> *tet0, regina::Tetrahedron<3> *tet1) const {
        Gluing<3, 3> g;
        for (int i = 0; i < 4; ++i) {
            if (tet0->adjacentSimplex(i) == tet1)
                return {tet0, i, tet1, tet0->adjacentGluing(i)};
        }

        throw regina::InvalidArgument("Given tetrahedra are not glued together!");
    }

    void subdivideSharedEdgeSequence_(std::vector<regina::Tetrahedron<3> *> &tets) {
        regina::Perm<4> p;
        for (const regina::EdgeEmbedding<3> &emb : sharedEdge_(tets)->embeddings()) {
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
                regina::FaceNumbering<3, 1>::edgeNumber[g.gluing[p[0]]][g.gluing[p[1]]]);
        }
    }

    bool isGlued_(const regina::Tetrahedron<3> *tet0, const regina::Tetrahedron<3> *tet1) const {
        for (int i = 0; i < 4; ++i) {
            if (tet0->adjacentSimplex(i) == tet1)
                return true;
        }
        return false;
    }

    regina::Edge<3> *sharedEdge_(const std::vector<regina::Tetrahedron<3> *> &tets) {
        auto ts = lastThree_(tets);
        regina::Triangle<3> *t0 = ts[0]->triangle(getGluing_(ts[0], ts[1]).srcFacet);
        regina::Triangle<3> *t1 = ts[1]->triangle(getGluing_(ts[1], ts[2]).srcFacet);
        std::unordered_set<regina::Edge<3> *> edges = {t0->edge(0), t0->edge(1), t0->edge(2)};

        for (int i = 0; i < 3; ++i) {
            if (edges.find(t1->edge(i)) != edges.end())
                return t1->edge(i);
        }

        throw regina::InvalidArgument("Given tetrahedra don't have a shared edge!");
    }

    std::array<regina::Tetrahedron<3> *, 3>
    lastThree_(const std::vector<regina::Tetrahedron<3> *> &tets) {
        std::array<regina::Tetrahedron<3> *, 3> ans;
        int j = 0;
        for (int i = tets.size() - 3; i < tets.size(); ++i) {
            ans[j++] = tets[i];
        }
        return ans;
    }

    regina::Vertex<3> *sharedVertex_(const std::vector<regina::Tetrahedron<3> *> &tets) const {
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
            std::set_intersection(verts.begin(), verts.end(), newVerts.begin(), newVerts.end(),
                                  std::inserter(intersection, intersection.begin()));
            verts = intersection;

            if (verts.size() == 1)
                return *(verts.begin());
        }

        if (verts.size() >= 1)
            return *(verts.begin());

        throw regina::InvalidArgument("Given tetrahedra do not share a common vertex!");
    }

    void subdivideSharedVertexSequence_(std::vector<regina::Tetrahedron<3> *> &tets) {
        const regina::Vertex<3> *sharedVertex = sharedVertex_(tets);
        int v;
        for (const regina::VertexEmbedding<3> &emb : sharedVertex->embeddings()) {
            if (emb.simplex() == tets[0]) {
                v = emb.vertex();
                break;
            }
        }

        for (int i = 0; i < tets.size() - 1;) {
            int start = i;
            std::vector<regina::Tetrahedron<3> *> sharedEdgeSubseq = {tets[i], tets[++i],
                                                                      tets[++i]};
            regina::Edge<3> *sharedEdge = sharedEdge_(sharedEdgeSubseq);

            while (i < tets.size() - 1 &&
                   sharedEdge == sharedEdge_({tets[i - 1], tets[i], tets[i + 1]})) {
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

            if (i >= tets.size() - 1)
                return;

            sharedEdge = sharedEdge_(sharedEdgeSubseq);
            sharedVertex = tets[0]->vertex(v);
            regina::Vertex<3> *otherVertex = sharedEdge->vertex(0) != sharedVertex
                                                 ? sharedEdge->vertex(0)
                                                 : sharedEdge->vertex(1);

            for (const regina::VertexEmbedding<3> &emb : otherVertex->embeddings()) {
                if (emb.simplex() == tets[i]) {
                    tets[i] = tets[i]->adjacentSimplex(emb.vertex());
                    break;
                }
            }
            for (const regina::VertexEmbedding<3> &emb : otherVertex->embeddings()) {
                if (emb.simplex() == tets[i - 1]) {
                    tets.insert(tets.begin() + i, tets[i - 1]->adjacentSimplex(emb.vertex()));
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
    Link(regina::Triangulation<3> &tri, const std::vector<Curve<3>> &components) : tri_(tri) {
        for (const Curve<3> &edges : components) {
            comps_.emplace_back(tri, edges);
        }
    }

    bool isUnlink() const {
        for (const Knot &k : comps_) {
            if (!k.isUnknot())
                return false;
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

    std::unordered_map<regina::Triangle<2> *, const regina::Triangle<dim> *> emb_;
    /**< A mapping from the top-dimensional simplices of the
     * sub-triangulation into the subdim-simplices of the dim-triangulation
     */
    std::unordered_map<const regina::Triangle<dim> *, regina::Triangle<2> *> inv_;
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
        if (!tri_->isClosed())
            bdry_ = tri_->boundaryComponent(0)->build();
    }

    KnottedSurface(const KnottedSurface &other)
        : tri_(other.tri_), surface_(other.surface_), improperEdges_(other.improperEdges_),
          invariants_(other.invariants_), indices_(other.indices_) {
        if (!tri_->isClosed())
            bdry_ = tri_->boundaryComponent(0)->build();
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

    bool isProper() const { return surface_.isClosed() || improperEdges_.empty(); }

    bool hasSelfIntersection() const {
        std::unordered_set<const regina::Vertex<dim> *> images;
        for (const regina::Vertex<2> *v : surface_.vertices()) {
            const regina::Vertex<dim> *im = image(v);
            if (images.find(im) != images.end())
                return true;
            images.insert(im);
        }
        return false;
    }

    Link boundary() {
        const regina::BoundaryComponent<dim> *b = tri_->boundaryComponent(0);
        std::unordered_map<const regina::Edge<dim> *, const regina::Edge<dim - 1> *> bdryMap;
        std::unordered_map<const regina::Edge<2> *, const regina::Edge<dim - 1> *> edgeMap;

        for (const regina::Edge<dim - 1> *e : bdry_.edges()) {
            bdryMap[b->edge(e->index())] = e;
        }

        for (const regina::BoundaryComponent<2> *comp : surface_.boundaryComponents()) {
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

        for (const regina::BoundaryComponent<2> *comp : surface_.boundaryComponents()) {
            Curve<dim - 1> edges;
            // for (const regina::Edge<2> *e : comp->edges()) {
            //     edges.push_back(edgeMap.at(e));
            // }
            comps.push_back(edges);
        }

        return {bdry_, comps};
    }

    template <int facedim>
    regina::Face<dim, facedim> *image(const regina::Face<2, facedim> *f) const {
        static_assert(0 <= facedim && facedim <= 2, "Must have 0 <= facedim <= 2");

        // for (auto &[t, f] : emb_) {
        //     std::cout << "(" << t->index() << ", " << f->index() << ") ";
        // }
        // std::cout << "\n";

        const regina::FaceEmbedding<2, facedim> &fEmb = f->front();
        return emb_.at(fEmb.simplex())->template face<facedim>(fEmb.face());
    }

    const regina::Triangle<dim> *image(regina::Triangle<2> *t) const { return emb_.at(t); }

    regina::Triangle<2> *preimage(const regina::Triangle<dim> *f) const {
        auto search = inv_.find(f);
        if (search == inv_.end())
            return nullptr;
        return search->second;
    }

    /** Mutation methods */

    bool addTriangle(const regina::Triangle<dim> *f, const typename GluingNode<dim>::AdjList &adj) {
        // TODO: Figure out a nice way for this method not to know what an
        // AdjList is
        regina::Triangle<2> *src = surface_.newSimplex();
        std::array<bool, 3> isBoundaryEdge = {true, true, true};

        emb_[src] = f;
        inv_[f] = src;

        for (const auto &[adjNode, g] : adj) {
            auto search = inv_.find(adjNode->f);
            if (search == inv_.end())
                continue;

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

    friend bool operator==(const KnottedSurface &lhs, const KnottedSurface &rhs) {
        return lhs.invariants_ == rhs.invariants_ && lhs.indices_ == rhs.indices_;
    }

    friend bool operator<(const KnottedSurface &lhs, const KnottedSurface &rhs) {
        return lhs.invariants_ < rhs.invariants_ ||
               (lhs.invariants_ == rhs.invariants_ && lhs.indices_ < rhs.indices_);
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

// template <int n>
// regina::Triangulation<3> knotComplement(const regina::Triangulation<3> t,
//                                         const Knot &k) {
//     return linkComplement(t, {k});
// }

#endif
