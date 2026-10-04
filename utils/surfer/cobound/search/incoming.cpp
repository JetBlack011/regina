//
//  incoming.cpp
//

#include "cobound/search/incoming.h"

#include "linknaming/complement/edgecycles.h"
#include "linknaming/diagrams/diagramiso.h"
#include "linknaming/tables.h"

#include <algorithm>
#include <numeric>
#include <optional>
#include <tuple>
#include <unordered_map>

namespace search {

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

// The incoming link's directed edges under `iso`: edge index -> (tail, head)
// vertex indices in `dest`, and the edge's position in diagramEdges.
using IncomingImage = std::unordered_map<size_t, std::tuple<size_t, size_t, size_t>>;
IncomingImage
directedImage(const std::vector<const regina::Edge<3> *> &diagramEdges,
              const std::vector<bool> &diagramReversed,
              const regina::Triangulation<3> &dest,
              const regina::Isomorphism<3> &iso) {
    IncomingImage image;
    for (size_t i = 0; i < diagramEdges.size(); ++i) {
        const regina::Vertex<3> *tail =
            diagramReversed[i] ? diagramEdges[i]->vertex(1) : diagramEdges[i]->vertex(0);
        const regina::Vertex<3> *head =
            diagramReversed[i] ? diagramEdges[i]->vertex(0) : diagramEdges[i]->vertex(1);
        image[mapEdgeIndex(diagramEdges[i], dest, iso)] = {
            mapVertexIndex(tail, dest, iso), mapVertexIndex(head, dest, iso), i};
    }
    return image;
}

std::vector<size_t> sortedKeys(const IncomingImage &image) {
    std::vector<size_t> keys;
    keys.reserve(image.size());
    for (const auto &[e, ends] : image)
        keys.push_back(e);
    std::ranges::sort(keys);
    return keys;
}
} // namespace

IncomingOrientation
buildIncomingOrientation(const std::vector<const regina::Edge<3> *> &diagramEdges,
                    const std::vector<bool> &diagramReversed,
                    const regina::Triangulation<3> &incomingTri,
                    const std::vector<size_t> *requiredEdges) {
    if (diagramEdges.empty())
        throw regina::InvalidArgument(
            "buildIncomingOrientation(): diagramEdges must not be empty");
    if (diagramReversed.size() != diagramEdges.size())
        throw regina::InvalidArgument(
            "buildIncomingOrientation(): diagramReversed and diagramEdges differ in size");

    const regina::Triangulation<3> &diagramTri = diagramEdges.front()->triangulation();

    // What the old code used: whichever isomorphism isIsomorphicTo() returns.
    std::optional<regina::Isomorphism<3>> legacy =
        diagramTri.isIsomorphicTo(incomingTri);
    if (!legacy)
        throw regina::InvalidArgument(
            "buildIncomingOrientation(): the diagram's own triangulation is not "
            "isomorphic to incomingTri");
    auto legacyImage =
        directedImage(diagramEdges, diagramReversed, incomingTri, *legacy);

    std::optional<IncomingImage> chosen;
    if (!requiredEdges || sortedKeys(legacyImage) == *requiredEdges) {
        chosen = legacyImage;
    } else {
        diagramTri.findAllIsomorphisms(
            incomingTri, [&](const regina::Isomorphism<3> &iso) {
                auto image =
                    directedImage(diagramEdges, diagramReversed, incomingTri, iso);
                if (sortedKeys(image) != *requiredEdges)
                    return false; // keep looking
                chosen = std::move(image);
                return true;
            });
        if (!chosen)
            throw regina::InvalidArgument(
                "buildIncomingOrientation(): no isomorphism takes the diagram's link "
                "onto the seed's edges in the incoming boundary");
    }

    IncomingOrientation result;
    std::vector<edgecycles::EdgeEnds> directed;
    directed.reserve(chosen->size());
    for (const auto &[e, ends] : *chosen) {
        const auto &[tail, head, diagramIndex] = ends;
        result.tailOf[e] = tail;
        result.incomingIndexOf[e] = diagramIndex;
        directed.push_back({e, tail, head});
    }
    const std::optional<size_t> cycles = edgecycles::countDirectedCycles(directed);
    if (!cycles)
        throw regina::InvalidArgument(
            "buildIncomingOrientation(): the link's edges do not chain into "
            "closed directed curves (one edge leaving and one arriving at "
            "every vertex)");
    result.components = *cycles;
    result.edges = sortedKeys(*chosen);
    result.divergedFromDefaultIsomorphism = (*chosen != legacyImage);
    return result;
}

} // namespace search

namespace search {

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

} // namespace search

namespace search {

void orientIncoming(IncomingThickening &thickened) {
    const auto &edges2 = thickened.link.edges;
    const auto &reversed2 = thickened.link.reversed;
    // The seed's own edges on the incoming side: exactly L x {0}.
    if (!thickened.seedFaces.empty())
        thickened.incomingEdges = search::boundaryEdgesOf(thickened.tri, thickened.seedFaces,
                                                   thickened.incomingBC);
    thickened.orientation = search::buildIncomingOrientation(
        edges2, reversed2, thickened.tri.boundaryComponent(thickened.incomingBC)->build(),
        thickened.seedFaces.empty() ? nullptr : &thickened.incomingEdges);
    if (thickened.seedFaces.empty())
        thickened.incomingEdges = thickened.orientation->edges;

    // Setup-time checks on the incoming link, in place of any per-surface
    // ones: the incoming side is fixed from here on.
    if (thickened.incomingEdges.size() != edges2.size())
        throw regina::InvalidArgument(
            "the search side holds " + std::to_string(thickened.incomingEdges.size()) +
            " link edges, the diagram " + std::to_string(edges2.size()));
    if (thickened.orientation->components !=
        static_cast<size_t>(thickened.componentCount))
        throw regina::InvalidArgument(
            "the search-side link has " +
            std::to_string(thickened.orientation->components) +
            " components, the diagram " + std::to_string(thickened.componentCount));
}

void buildIncoming(const std::string &pdNotation, int thickenLayers,
              int collarLayers, IncomingThickening &thickened) {
    buildAmbient(pdNotation, thickenLayers, collarLayers, thickened);
    orientIncoming(thickened);
}

linknaming::GaussDiagram gaussOf(const diagramtriangulation::Diagram &d) {
    linknaming::GaussDiagram g;
    for (const auto &c : d.crossings) g.signs.push_back(c.sign);
    g.comps = d.gauss;
    g.origin.resize(d.components);
    for (size_t i = 0; i < d.components; ++i) g.origin[i] = i;
    if (g.comps.size() != d.components)
        throw std::logic_error("gaussOf: gauss codes do not cover every component");
    return g;
}

std::vector<int> certifyIncoming(const diagramtriangulation::DiagramDrawer &drawer,
                                 const std::vector<diagramtriangulation::EdgeCycle> &cycles,
                                 const linknaming::GaussDiagram &given) {
    const linknaming::GaussDiagram drawn = gaussOf(drawer.draw(cycles));
    auto iso = linknaming::findDiagramIsomorphism(drawn, given, /*allowMirror=*/false,
                                                  /*allowReverse=*/false);
    if (!iso)
        throw IncomingNotCertified(
            "the triangulated incoming link does not redraw as its diagram (no "
            "orientation-preserving isomorphism): not certified, refused");
    return iso->componentMap;
}

linknaming::GaussDiagram diagramOfPD(const std::string &pd) {
    const regina::Link link = linknaming::linkFromTablePD(pd);
    std::vector<size_t> origin(link.countComponents());
    std::iota(origin.begin(), origin.end(), size_t(0));
    return linknaming::GaussDiagram::of(link, origin);
}
} // namespace search
