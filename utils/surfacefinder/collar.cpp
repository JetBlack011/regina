//
//  collar.cpp
//

#include "collar.h"

void CollarBuilder::addLayer(const CobordismBuilder<3> &cob) {
    const regina::Triangulation<3> &baseTri = cob.baseTriangulation();
    for (int edgeIdx : edgeIndices_) {
        const regina::Edge<3> *e = baseTri.edge(edgeIdx);

        // Any tetrahedron incident to e gives the same ambient
        // triangles once gluing is accounted for (verified empirically
        // against a shared-face example), so the first embedding is as
        // good as any other.
        const regina::EdgeEmbedding<3> &emb = e->front();
        regina::Tetrahedron<3> *tet = emb.tetrahedron();
        regina::Perm<4> m = emb.vertices();
        int a = m[0], b = m[1];
        if (a > b)
            std::swap(a, b);

        addEdgeSweep_(cob, tet, a, b);
    }
}

std::unordered_set<regina::Triangle<4> *> CollarBuilder::resolve() const {
    std::unordered_set<regina::Triangle<4> *> result;
    for (const FaceDesc &f : faces_) {
        result.insert(resolveFace_(f));
    }
    return result;
}

int CollarBuilder::encodeVertex_(int v, bool isTop) {
    constexpr int D = 4;
    if (!isTop)
        return v;
    return v == 0 ? D : v - 1;
}

void CollarBuilder::addEdgeSweep_(const CobordismBuilder<3> &cob,
                                  regina::Tetrahedron<3> *tet, int a, int b) {
    regina::Simplex<4> *simpA = cob.currentTopSimplex(tet, a);
    regina::Simplex<4> *simpB = cob.currentTopSimplex(tet, b);

    faces_.push_back({simpA,
                      {encodeVertex_(a, false), encodeVertex_(b, false),
                       encodeVertex_(a, true)}});
    faces_.push_back({simpB,
                      {encodeVertex_(b, false), encodeVertex_(a, true),
                       encodeVertex_(b, true)}});
}

regina::Triangle<4> *CollarBuilder::resolveFace_(const FaceDesc &f) {
    std::array<int, 5> full;
    std::unordered_set<int> used(f.verts.begin(), f.verts.end());
    full[0] = f.verts[0];
    full[1] = f.verts[1];
    full[2] = f.verts[2];
    int idx = 3;
    for (int i = 0; i <= 4; ++i) {
        if (!used.count(i))
            full[idx++] = i;
    }
    return f.simplex->triangle(
        regina::Face<4, 2>::faceNumber(regina::Perm<5>(full)));
}
