//
//  incoming.cpp
//

#include "cobound/search/incoming.h"

#include "linknaming/complement/edgecycles.h"

#include <algorithm>
#include <optional>
#include <tuple>
#include <unordered_map>

namespace cobordismgraph {

namespace {
// Maps ambient vertex `v` to its corresponding vertex in `dest`, via `iso`
// (which must map `v`'s own triangulation to `dest`).
size_t mapVertexIndex(const regina::Vertex<3> *v,
                      const regina::Triangulation<3> &dest,
                      const regina::Isomorphism<3> &iso) {
    auto emb = v->front();
    size_t destTet = iso.simpImage(emb.tetrahedron()->index());
    int destLocal = iso.facetPerm(emb.tetrahedron()->index())[emb.vertex()];
    return dest.tetrahedron(destTet)->vertex(destLocal)->index();
}
} // namespace

namespace {
// Maps edge `e` to its image's index in `dest` under `iso`.
size_t mapEdgeIndex(const regina::Edge<3> *e,
                    const regina::Triangulation<3> &dest,
                    const regina::Isomorphism<3> &iso) {
    auto emb = e->front();
    size_t tet = emb.tetrahedron()->index();
    regina::Perm<4> p = iso.facetPerm(tet);
    regina::Perm<4> v = emb.vertices();
    return dest.tetrahedron(iso.simpImage(tet))->edge(p[v[0]], p[v[1]])->index();
}

// The row's directed link under `iso`: edge index -> (tail, head) vertex
// indices in `dest`, and the edge's position in rowEdges.
using RowImage = std::unordered_map<size_t, std::tuple<size_t, size_t, size_t>>;
RowImage
directedImage(const std::vector<const regina::Edge<3> *> &rowEdges,
              const std::vector<bool> &rowReversed,
              const regina::Triangulation<3> &dest,
              const regina::Isomorphism<3> &iso) {
    RowImage image;
    for (size_t i = 0; i < rowEdges.size(); ++i) {
        const regina::Vertex<3> *tail =
            rowReversed[i] ? rowEdges[i]->vertex(1) : rowEdges[i]->vertex(0);
        const regina::Vertex<3> *head =
            rowReversed[i] ? rowEdges[i]->vertex(0) : rowEdges[i]->vertex(1);
        image[mapEdgeIndex(rowEdges[i], dest, iso)] = {
            mapVertexIndex(tail, dest, iso), mapVertexIndex(head, dest, iso), i};
    }
    return image;
}

std::vector<size_t> sortedKeys(const RowImage &image) {
    std::vector<size_t> keys;
    keys.reserve(image.size());
    for (const auto &[e, ends] : image)
        keys.push_back(e);
    std::ranges::sort(keys);
    return keys;
}
} // namespace

RowOrientation
buildRowOrientation(const std::vector<const regina::Edge<3> *> &rowEdges,
                    const std::vector<bool> &rowReversed,
                    const regina::Triangulation<3> &searchSideTri,
                    const std::vector<size_t> *requiredEdges) {
    if (rowEdges.empty())
        throw regina::InvalidArgument(
            "buildRowOrientation(): rowEdges must not be empty");
    if (rowReversed.size() != rowEdges.size())
        throw regina::InvalidArgument(
            "buildRowOrientation(): rowReversed and rowEdges differ in size");

    const regina::Triangulation<3> &rowTri = rowEdges.front()->triangulation();

    // What the old code used: whichever isomorphism isIsomorphicTo() returns.
    std::optional<regina::Isomorphism<3>> legacy =
        rowTri.isIsomorphicTo(searchSideTri);
    if (!legacy)
        throw regina::InvalidArgument(
            "buildRowOrientation(): the row's own triangulation is not "
            "isomorphic to searchSideTri");
    auto legacyImage =
        directedImage(rowEdges, rowReversed, searchSideTri, *legacy);

    std::optional<RowImage> chosen;
    if (!requiredEdges || sortedKeys(legacyImage) == *requiredEdges) {
        chosen = legacyImage;
    } else {
        rowTri.findAllIsomorphisms(
            searchSideTri, [&](const regina::Isomorphism<3> &iso) {
                auto image =
                    directedImage(rowEdges, rowReversed, searchSideTri, iso);
                if (sortedKeys(image) != *requiredEdges)
                    return false; // keep looking
                chosen = std::move(image);
                return true;
            });
        if (!chosen)
            throw regina::InvalidArgument(
                "buildRowOrientation(): no isomorphism takes the row's link "
                "onto the seed's edges in the search-side boundary");
    }

    RowOrientation result;
    std::vector<edgecycles::EdgeEnds> directed;
    directed.reserve(chosen->size());
    for (const auto &[e, ends] : *chosen) {
        const auto &[tail, head, rowIndex] = ends;
        result.tailOf[e] = tail;
        result.rowIndexOf[e] = rowIndex;
        directed.push_back({e, tail, head});
    }
    const std::optional<size_t> cycles = edgecycles::countDirectedCycles(directed);
    if (!cycles)
        throw regina::InvalidArgument(
            "buildRowOrientation(): the link's edges do not chain into "
            "closed directed curves (one edge leaving and one arriving at "
            "every vertex)");
    result.components = *cycles;
    result.edges = sortedKeys(*chosen);
    result.divergedFromDefaultIsomorphism = (*chosen != legacyImage);
    return result;
}

} // namespace cobordismgraph

namespace farside {

std::vector<size_t> boundaryEdgesOf(const regina::Triangulation<4> &tri,
                                    const std::vector<int> &faces,
                                    size_t bcIndex) {
    const regina::BoundaryComponent<4> *bc = tri.boundaryComponent(bcIndex);
    std::unordered_map<const regina::Edge<4> *, size_t> local;
    for (size_t k = 0; k < bc->countEdges(); ++k) local.emplace(bc->edge(k), k);
    std::vector<size_t> edges;
    for (int f : faces) {
        const regina::Triangle<4> *t = tri.triangle(f);
        for (int i = 0; i < 3; ++i)
            if (auto it = local.find(t->edge(i)); it != local.end())
                edges.push_back(it->second);
    }
    std::ranges::sort(edges);
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

} // namespace farside

namespace rowsearch {

void orientRow(RowBuild &row) {
    const auto &edges2 = row.link.edges;
    const auto &reversed2 = row.link.reversed;
    // The seed's own edges on the search side: exactly L x {0}.
    if (!row.seedFaces.empty())
        row.searchEdges = farside::boundaryEdgesOf(row.tri, row.seedFaces,
                                                   row.searchSideBC);
    row.orientation = cobordismgraph::buildRowOrientation(
        edges2, reversed2, row.tri.boundaryComponent(row.searchSideBC)->build(),
        row.seedFaces.empty() ? nullptr : &row.searchEdges);
    if (row.seedFaces.empty())
        row.searchEdges = row.orientation->edges;

    // Setup-time checks on the row's own link, in place of any per-surface
    // ones: the search side is fixed from here on.
    if (row.searchEdges.size() != edges2.size())
        throw regina::InvalidArgument(
            "the search side holds " + std::to_string(row.searchEdges.size()) +
            " link edges, the diagram " + std::to_string(edges2.size()));
    if (row.orientation->components !=
        static_cast<size_t>(row.componentCount))
        throw regina::InvalidArgument(
            "the search-side link has " +
            std::to_string(row.orientation->components) +
            " components, the diagram " + std::to_string(row.componentCount));
}

void buildRow(const std::string &pdNotation, int thickenLayers,
              int collarLayers, RowBuild &row) {
    buildAmbient(pdNotation, thickenLayers, collarLayers, row);
    orientRow(row);
}
} // namespace rowsearch
