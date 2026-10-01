//
//  linkingnumber.cpp
//

#include "surfer/submanifold/linkingnumber.h"

#include <algorithm>
#include <deque>
#include <map>
#include <unordered_map>
#include <utility>

namespace linkingnumber {

std::atomic<bool> auditLinkingNumbers{false};

namespace {

// Arithmetic in GF(p), p = 2^61 - 1.
constexpr uint64_t P = (uint64_t{1} << 61) - 1;

uint64_t reduce(int64_t v) {
    int64_t r = v % static_cast<int64_t>(P);
    return static_cast<uint64_t>(r < 0 ? r + static_cast<int64_t>(P) : r);
}
uint64_t add(uint64_t a, uint64_t b) {
    uint64_t s = a + b;
    return s >= P ? s - P : s;
}
uint64_t sub(uint64_t a, uint64_t b) { return a >= b ? a - b : a + P - b; }
uint64_t mul(uint64_t a, uint64_t b) {
    __uint128_t z = static_cast<__uint128_t>(a) * b;
    uint64_t s = static_cast<uint64_t>(z & P) + static_cast<uint64_t>(z >> 61);
    return s >= P ? s - P : s;
}
uint64_t inverse(uint64_t a) { // a != 0; Fermat
    uint64_t result = 1, base = a, e = P - 2;
    while (e) {
        if (e & 1)
            result = mul(result, base);
        base = mul(base, base);
        e >>= 1;
    }
    return result;
}
// The representative in (-p/2, p/2).
int64_t lift(uint64_t v) {
    return v > P / 2 ? static_cast<int64_t>(v) - static_cast<int64_t>(P)
                     : static_cast<int64_t>(v);
}

// One edge of an oriented cycle: `edge` run from `from` to `to`, `dir` = +1
// along the edge's own vertex(0) -> vertex(1), -1 against it.
struct Step {
    int edge, dir, from, to;
};

// Walks `curve` into an oriented cycle, or nullopt if it is not a single
// simple closed curve of `tri`.
std::optional<std::vector<Step>>
orient(const regina::Triangulation<3> &tri,
       const std::vector<const regina::Edge<3> *> &curve) {
    if (curve.empty())
        return std::nullopt;
    std::unordered_map<int, std::vector<int>> at; // vertex -> slots in curve
    for (size_t i = 0; i < curve.size(); ++i) {
        if (&curve[i]->triangulation() != &tri)
            return std::nullopt;
        at[static_cast<int>(curve[i]->vertex(0)->index())].push_back(
            static_cast<int>(i));
        at[static_cast<int>(curve[i]->vertex(1)->index())].push_back(
            static_cast<int>(i));
    }
    for (const auto &[v, slots] : at)
        if (slots.size() != 2)
            return std::nullopt; // not a simple closed curve
    std::vector<Step> steps;
    std::vector<char> used(curve.size(), 0);
    int slot = 0;
    int from = static_cast<int>(curve[0]->vertex(0)->index());
    const int start = from;
    for (size_t n = 0; n < curve.size(); ++n) {
        used[slot] = 1;
        const regina::Edge<3> *e = curve[slot];
        const int v0 = static_cast<int>(e->vertex(0)->index());
        const int v1 = static_cast<int>(e->vertex(1)->index());
        const int dir = (v0 == from) ? 1 : -1;
        const int to = dir == 1 ? v1 : v0;
        steps.push_back({static_cast<int>(e->index()), dir, from, to});
        from = to;
        if (n + 1 == curve.size())
            break;
        slot = -1;
        for (int s : at[from])
            if (!used[s]) {
                slot = s;
                break;
            }
        if (slot < 0)
            return std::nullopt;
    }
    if (from != start)
        return std::nullopt;
    return steps;
}

} // namespace

Complex::Complex(const regina::Triangulation<3> &tri) : tri_(&tri) {
    nV_ = static_cast<int>(tri.countVertices());
    nE_ = static_cast<int>(tri.countEdges());
    nF_ = static_cast<int>(tri.countTriangles());
    nT_ = static_cast<int>(tri.size());
    if (nT_ == 0 || !tri.isOrientable() || tri.hasBoundaryFacets() ||
        !tri.isConnected()) {
        valid_ = false;
        return;
    }

    edgeEnds_.resize(nE_);
    for (int e = 0; e < nE_; ++e)
        edgeEnds_[e] = {static_cast<int>(tri.edge(e)->vertex(0)->index()),
                        static_cast<int>(tri.edge(e)->vertex(1)->index())};

    // (delta x)(f) = sum_k (-1)^k x([f without its vertex k]), each edge taken
    // in f's induced order and compared with the edge's own.
    faceEdge_.resize(nF_);
    faceEdgeSign_.resize(nF_);
    for (int f = 0; f < nF_; ++f) {
        const regina::Triangle<3> *t = tri.triangle(f);
        for (int k = 0; k < 3; ++k) {
            faceEdge_[f][k] = static_cast<int>(t->edge(k)->index());
            // q[0], q[1]: the triangle's vertices at the edge's vertex(0),
            // vertex(1).
            const regina::Perm<4> q = t->edgeMapping(k);
            const int along = q[0] < q[1] ? 1 : -1;
            faceEdgeSign_[f][k] = static_cast<int8_t>((k % 2 ? -1 : 1) * along);
        }
    }

    // The coefficient of each tetrahedron's faces in the boundary of the
    // oriented fundamental cycle sum_t s(t) t.
    tetFace_.resize(nT_);
    tetSign_.resize(nT_);
    tetAdj_.resize(nT_);
    tetGluing_.resize(nT_);
    tetVertex_.resize(nT_);
    for (int i = 0; i < nT_; ++i) {
        const regina::Tetrahedron<3> *tet = tri.tetrahedron(i);
        const int s = tet->orientation();
        for (int c = 0; c < 4; ++c)
            tetVertex_[i][c] = static_cast<int>(tet->vertex(c)->index());
        for (int j = 0; j < 4; ++j) {
            tetFace_[i][j] = static_cast<int>(tet->triangle(j)->index());
            const regina::Perm<4> p = tet->triangleMapping(j);
            const int inversions =
                (p[0] > p[1]) + (p[0] > p[2]) + (p[1] > p[2]);
            const int o = inversions % 2 ? -1 : 1;
            tetSign_[i][j] = static_cast<int8_t>(s * (j % 2 ? -1 : 1) * o);
            const regina::Tetrahedron<3> *adj = tet->adjacentTetrahedron(j);
            if (!adj) {
                valid_ = false;
                return;
            }
            tetAdj_[i][j] = static_cast<int>(adj->index());
            tetGluing_[i][j] = tet->adjacentGluing(j);
        }
    }
    // A closed, consistently oriented triangulation has sum_t s(t) t a
    // cycle: the two sides of every triangle cancel. Checked, since
    // everything below rests on it.
    for (int i = 0; i < nT_; ++i)
        for (int j = 0; j < 4; ++j) {
            const int other = tetAdj_[i][j];
            const int otherFace = tetGluing_[i][j][j];
            if (tetSign_[i][j] != -tetSign_[other][otherFace]) {
                valid_ = false;
                return;
            }
        }
}

std::optional<long>
linkingNumber(const Complex &cx, const std::vector<const regina::Edge<3> *> &A,
              const std::vector<const regina::Edge<3> *> &B) {
    if (!cx.valid_)
        return std::nullopt;
    const auto stepsA = orient(*cx.tri_, A);
    const auto stepsB = orient(*cx.tri_, B);
    if (!stepsA || !stepsB)
        return std::nullopt;
    {
        std::vector<char> onA(cx.nV_, 0);
        for (const Step &s : *stepsA)
            onA[s.from] = 1;
        for (const Step &s : *stepsB)
            if (onA[s.from])
                return std::nullopt; // not vertex-disjoint
    }

    // 1. beta = PD(B*): walk B's push-off through the dual cells.
    //
    // Each vertex b_k of B needs a path of tetrahedra from t_{k-1} (on the
    // edge before it) to t_k (on the edge after it), through triangles
    // containing b_k. The whole star of b_k is searched, not just until t_k
    // turns up, because the homotopy argument needs EVERY tetrahedron there
    // to meet b_k at one corner (see linkingnumber.h).
    const std::vector<Step> &b = *stepsB;
    const int m = static_cast<int>(b.size());
    std::vector<int> onEdge(m); // t_k
    for (int k = 0; k < m; ++k)
        onEdge[k] = static_cast<int>(
            cx.tri_->edge(b[k].edge)->front().simplex()->index());

    auto cornerOf = [&](int tet, int vertex) {
        int corner = -1, count = 0;
        for (int c = 0; c < 4; ++c)
            if (cx.tetVertex_[tet][c] == vertex) {
                corner = c;
                ++count;
            }
        return count == 1 ? corner : -1;
    };

    std::vector<int64_t> beta(cx.nF_, 0);
    std::vector<int> seenStamp(cx.nT_, -1);
    std::vector<std::pair<int, int>> parent(cx.nT_); // (tetrahedron, face) we came through
    std::deque<int> queue;
    for (int k = 0; k < m; ++k) {
        const int vertex = b[k].from;
        const int source = onEdge[(k + m - 1) % m];
        const int target = onEdge[k];
        if (cornerOf(source, vertex) < 0)
            return std::nullopt;
        queue.clear();
        queue.push_back(source);
        seenStamp[source] = k;
        while (!queue.empty()) {
            const int t = queue.front();
            queue.pop_front();
            const int c = cornerOf(t, vertex);
            for (int i = 0; i < 4; ++i) {
                if (i == c)
                    continue; // the one face missing the vertex
                const int next = cx.tetAdj_[t][i];
                const int nextCorner = cx.tetGluing_[t][i][c];
                if (cornerOf(next, vertex) != nextCorner)
                    return std::nullopt; // meets the vertex twice
                if (seenStamp[next] == k)
                    continue;
                seenStamp[next] = k;
                parent[next] = {t, i};
                queue.push_back(next);
            }
        }
        if (seenStamp[target] != k)
            return std::nullopt;
        // Crossing from t through its face i into next (entering by face j)
        // adds the entered side's sign; see Complex::tetSign_.
        for (int t = target; t != source;) {
            const auto [prev, face] = parent[t];
            const int entered = cx.tetGluing_[prev][face][face];
            beta[cx.tetFace_[prev][face]] += cx.tetSign_[t][entered];
            t = prev;
        }
    }
    // Self-check: B* is a cycle, i.e. beta a cocycle.
    for (int t = 0; t < cx.nT_; ++t) {
        int64_t sum = 0;
        for (int i = 0; i < 4; ++i)
            sum += cx.tetSign_[t][i] * beta[cx.tetFace_[t][i]];
        if (sum != 0)
            return std::nullopt;
    }

    // 2. Solve delta x = beta over GF(p). Gauge: x = 0 on a spanning tree of
    // the 1-skeleton, which leaves exactly one solution. Then peel triangles
    // with one unknown edge; eliminate whatever is left.
    std::vector<uint64_t> x(cx.nE_, 0);
    std::vector<char> known(cx.nE_, 0);
    std::vector<std::vector<int>> vertexEdges(cx.nV_);
    std::vector<std::vector<int>> edgeFaces(cx.nE_);
    for (int e = 0; e < cx.nE_; ++e) {
        vertexEdges[cx.edgeEnds_[e][0]].push_back(e);
        if (cx.edgeEnds_[e][1] != cx.edgeEnds_[e][0])
            vertexEdges[cx.edgeEnds_[e][1]].push_back(e);
    }
    for (int f = 0; f < cx.nF_; ++f)
        for (int k = 0; k < 3; ++k) {
            auto &list = edgeFaces[cx.faceEdge_[f][k]];
            if (list.empty() || list.back() != f)
                list.push_back(f);
        }
    {
        std::vector<char> reached(cx.nV_, 0);
        std::deque<int> q{0};
        reached[0] = 1;
        while (!q.empty()) {
            const int u = q.front();
            q.pop_front();
            for (int e : vertexEdges[u]) {
                const int w = cx.edgeEnds_[e][0] == u ? cx.edgeEnds_[e][1]
                                                      : cx.edgeEnds_[e][0];
                if (reached[w])
                    continue;
                reached[w] = 1;
                known[e] = 1; // x(e) = 0
                q.push_back(w);
            }
        }
    }

    auto unknownEdgesOf = [&](int f) {
        // Distinct unknown edges of f, with their summed coefficients.
        std::array<std::pair<int, int>, 3> out{};
        int n = 0;
        for (int k = 0; k < 3; ++k) {
            const int e = cx.faceEdge_[f][k];
            if (known[e])
                continue;
            int slot = 0;
            while (slot < n && out[slot].first != e)
                ++slot;
            if (slot == n)
                out[n++] = {e, 0};
            out[slot].second += cx.faceEdgeSign_[f][k];
        }
        return std::make_pair(out, n);
    };
    auto knownPart = [&](int f) { // beta(f) - sum over known edges
        uint64_t rhs = reduce(beta[f]);
        for (int k = 0; k < 3; ++k) {
            const int e = cx.faceEdge_[f][k];
            if (known[e])
                rhs = sub(rhs, mul(reduce(cx.faceEdgeSign_[f][k]), x[e]));
        }
        return rhs;
    };

    std::deque<int> peel;
    for (int f = 0; f < cx.nF_; ++f)
        peel.push_back(f);
    while (!peel.empty()) {
        const int f = peel.front();
        peel.pop_front();
        const auto [unknowns, n] = unknownEdgesOf(f);
        if (n != 1 || unknowns[0].second == 0)
            continue;
        const int e = unknowns[0].first;
        x[e] = mul(knownPart(f), inverse(reduce(unknowns[0].second)));
        known[e] = 1;
        for (int g : edgeFaces[e])
            peel.push_back(g);
    }

    // Whatever peeling left: sparse Gaussian elimination, pivoting on each
    // row's smallest column, so every stored row has its pivot first.
    using Row = std::vector<std::pair<int, uint64_t>>; // (edge, coefficient), sorted
    std::map<int, std::pair<Row, uint64_t>> pivots;    // pivot edge -> (row, rhs)
    size_t work = 0;
    constexpr size_t WORK_LIMIT = 200'000'000;
    for (int f = 0; f < cx.nF_; ++f) {
        const auto [unknowns, n] = unknownEdgesOf(f);
        if (n == 0)
            continue;
        Row row;
        for (int s = 0; s < n; ++s)
            if (unknowns[s].second != 0)
                row.emplace_back(unknowns[s].first,
                                 reduce(unknowns[s].second));
        std::sort(row.begin(), row.end());
        uint64_t rhs = knownPart(f);
        while (!row.empty()) {
            const auto it = pivots.find(row.front().first);
            if (it == pivots.end()) {
                const uint64_t scale = inverse(row.front().second);
                for (auto &[col, coef] : row)
                    coef = mul(coef, scale);
                rhs = mul(rhs, scale);
                pivots.emplace(row.front().first,
                               std::make_pair(std::move(row), rhs));
                row.clear();
                break;
            }
            // row -= row[0] * pivotRow (whose leading coefficient is 1).
            const uint64_t factor = row.front().second;
            const Row &pr = it->second.first;
            Row merged;
            merged.reserve(row.size() + pr.size());
            size_t a = 0, z = 0;
            while (a < row.size() || z < pr.size()) {
                if (z == pr.size() ||
                    (a < row.size() && row[a].first < pr[z].first)) {
                    merged.push_back(row[a++]);
                } else if (a == row.size() || pr[z].first < row[a].first) {
                    merged.emplace_back(pr[z].first,
                                        sub(0, mul(factor, pr[z].second)));
                    ++z;
                } else {
                    const uint64_t v =
                        sub(row[a].second, mul(factor, pr[z].second));
                    if (v != 0)
                        merged.emplace_back(row[a].first, v);
                    ++a;
                    ++z;
                }
            }
            rhs = sub(rhs, mul(factor, it->second.second));
            work += merged.size();
            if (work > WORK_LIMIT)
                return std::nullopt; // fill-in out of hand; let the caller fall back
            row = std::move(merged);
        }
        // A row reduced to nothing is a dependent equation. If its rhs is
        // not 0 the system is inconsistent -- impossible for a cocycle beta
        // -- and the final check below refuses the answer.
    }
    // Back-substitution, highest pivot first; an unknown edge with no pivot
    // is free (none should be) and stays 0.
    for (auto it = pivots.rbegin(); it != pivots.rend(); ++it) {
        const auto &[row, rhs] = it->second;
        uint64_t v = rhs;
        for (size_t s = 1; s < row.size(); ++s)
            v = sub(v, mul(row[s].second, x[row[s].first]));
        x[it->first] = v;
        known[it->first] = 1;
    }

    // Self-check: x really solves delta x = beta, on every triangle. This is
    // what makes any failure above harmless -- no x that fails it is used.
    for (int f = 0; f < cx.nF_; ++f) {
        uint64_t lhs = 0;
        for (int k = 0; k < 3; ++k)
            lhs = add(lhs, mul(reduce(cx.faceEdgeSign_[f][k]),
                               x[cx.faceEdge_[f][k]]));
        if (lhs != reduce(beta[f]))
            return std::nullopt;
    }

    // 3. lk(A, B) = +-x(A).
    uint64_t value = 0;
    for (const Step &s : *stepsA)
        value = s.dir > 0 ? add(value, x[s.edge]) : sub(value, x[s.edge]);
    const int64_t lk = lift(value);
    return static_cast<long>(lk < 0 ? -lk : lk);
}

} // namespace linkingnumber
