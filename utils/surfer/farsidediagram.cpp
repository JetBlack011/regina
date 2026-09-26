//
//  farsidediagram.cpp
//
//  Oriented diagrams of witnesses' outgoing links, from their pair
//  signatures.
//
//  Usage:
//    farsidediagram [--layers N] '<row PD code>' < pairsigs
//
//  --layers is the witnesses' thicken_layers (cobordisms.csv): 2, the
//  default, for everything since early September; 1 for the earliest runs.
//
//  stdin: one "<id> <pair signature>" per line, every witness from the same
//  row (the row whose PD code is the argument). stdout, one line per
//  witness, plus a ROW line first:
//
//    ROW components=<n> lk=<matrix>
//    W <id> ok components=<m> crossingless=<list> pd=<[[a,b,c,d],...]>
//        lk=<m x m matrix> surface=<surface component of each outgoing curve>
//        incoming=<for each surface component: the row components it meets>
//    W <id> FAILED <reason>
//
//  The row's link is knotbuilder's, oriented as its PD code says; `lk` of
//  the ROW line is its linking matrix in DiagramDrawer::cyclesOf()'s
//  component order, which is also the order `incoming` refers to.
//
//  How: the row is thickened exactly as verifyslicegenus thickens it
//  (knotbuilder, CobordismBuilder x2, CollarBuilder), the witness's decoded
//  pair is carried onto that thickening by an isomorphism sending its
//  incoming curve onto the row's L x {0}, and from there on everything is
//  what the search itself would have had: farside::orientedOutgoingLink()
//  and knotbuilder::DiagramDrawer. The isomorphism is the only search in
//  the pipeline, and it is pinned by L x {0}, so the outgoing side is read
//  through the thickening's own product structure -- never mirrored by an
//  automorphism of T.
//

#include <algorithm>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "cobordismgraph.h"
#include "collar.h"
#include "embeddedsubmanifold.h"
#include "farsidecurves.h"
#include "knotbuilder/diagramdrawer.h"
#include "knotbuilder/knotbuilder.h"
#include "pairsig.h"
#include "skeleton.h"

namespace {

std::string matrix(const knotbuilder::Diagram &d) {
    std::ostringstream o;
    o << '[';
    for (size_t i = 0; i < d.components; ++i) {
        o << (i ? ",[" : "[");
        for (size_t j = 0; j < d.components; ++j)
            o << (j ? "," : "") << (i == j ? 0 : d.linkingNumber(i, j));
        o << ']';
    }
    return o.str() + ']';
}

template <typename T>
std::string list(const std::vector<T> &v) {
    std::ostringstream o;
    o << '[';
    for (size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
    return o.str() + ']';
}

// Face index of `t` (a triangle of `from`) under `iso`, in `to`.
int carryTriangle(const regina::Triangle<4> *t, const regina::Isomorphism<4> &iso,
                  const regina::Triangulation<4> &to) {
    auto emb = t->front();
    size_t p = emb.pentachoron()->index();
    regina::Perm<5> v = iso.facetPerm(p) * emb.vertices();
    return static_cast<int>(
        to.pentachoron(iso.simpImage(p))->triangle(regina::Face<4, 2>::faceNumber(v))->index());
}

} // namespace

int main(int argc, char **argv) {
    int layers = 2;
    int arg = 1;
    if (argc == 4 && std::string(argv[1]) == "--layers") {
        layers = std::stoi(argv[2]);
        arg = 3;
    }
    if (argc != arg + 1 || layers < 1) {
        std::cerr << "usage: farsidediagram [--layers N] '<row PD code>' < pairsigs\n";
        return 2;
    }
    const knotbuilder::PDCode pd = knotbuilder::parsePDCode(argv[arg]);
    auto [T, edges, reversed] = knotbuilder::buildLink(pd);

    // The row's thickening, exactly as verifyslicegenus builds it.
    std::vector<int> edgeIndices;
    for (const regina::Edge<3> *e : edges) edgeIndices.push_back(static_cast<int>(e->index()));
    CobordismBuilder<3> cob(T);
    CollarBuilder collar(edgeIndices);
    for (int i = 0; i < layers; ++i) {
        cob.thicken();
        collar.addLayer(cob);
    }
    const regina::Triangulation<4> &W = cob.getCobordism();
    std::vector<int> seed;
    for (regina::Triangle<4> *t : collar.resolve()) seed.push_back(static_cast<int>(t->index()));
    const size_t incomingBC = cob.baseBoundaryComponent()->index();
    const std::vector<size_t> rowEdges = farside::boundaryEdgesOf(W, seed, incomingBC);
    const cobordismgraph::RowOrientation row = cobordismgraph::buildRowOrientation(
        edges, reversed, W.boundaryComponent(incomingBC)->build(), &rowEdges);
    const farside::OutgoingMap outgoing(T, cob);
    const knotbuilder::DiagramDrawer drawer(T, pd.size());

    // The row's own components, in cyclesOf() order, and each row edge's component.
    const auto rowCycles = knotbuilder::DiagramDrawer::cyclesOf(edges, reversed);
    std::vector<size_t> componentOfRowEdge(edges.size());
    {
        std::unordered_map<size_t, size_t> compOfT;
        for (size_t c = 0; c < rowCycles.size(); ++c)
            for (const auto &de : rowCycles[c]) compOfT[de.edge] = c;
        for (size_t i = 0; i < edges.size(); ++i)
            componentOfRowEdge[i] = compOfT.at(edges[i]->index());
    }
    const knotbuilder::Diagram rowDiagram = drawer.draw(rowCycles);
    std::cout << "ROW components=" << rowCycles.size() << " lk=" << matrix(rowDiagram) << "\n";

    Skeleton<4, 2> skeleton(W);
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string id, sig;
        if (!(in >> id >> sig)) continue;
        try {
            DecodedKnottedSurfaceSig dec = fromKnottedSurfaceSig(sig);
            const std::vector<int> faces = dec.surface->markedFaces();

            // The isomorphism onto W sending the witness's incoming curve onto L x {0}.
            std::vector<int> carried;
            bool found = false;
            dec.ambient->findAllIsomorphisms(W, [&](const regina::Isomorphism<4> &iso) {
                std::vector<int> image;
                image.reserve(faces.size());
                for (int f : faces) image.push_back(carryTriangle(dec.ambient->triangle(f), iso, W));
                if (farside::boundaryEdgesOf(W, image, incomingBC) != rowEdges) return false;
                carried = std::move(image);
                found = true;
                return true;
            });
            if (!found) {
                // Say why: no isomorphism at all (another thickening), or
                // isomorphisms whose incoming curve is not L x {0}.
                size_t isos = 0;
                std::string sizes;
                dec.ambient->findAllIsomorphisms(W, [&](const regina::Isomorphism<4> &iso) {
                    std::vector<int> image;
                    for (int f : faces) image.push_back(carryTriangle(dec.ambient->triangle(f), iso, W));
                    if (isos < 4) {
                        auto in = farside::boundaryEdgesOf(W, image, incomingBC);
                        auto out = farside::boundaryEdgesOf(W, image, outgoing.boundaryComponent());
                        std::vector<size_t> common;
                        std::ranges::set_intersection(in, rowEdges, std::back_inserter(common));
                        sizes += " [in " + std::to_string(in.size()) + " edges, " +
                                 std::to_string(common.size()) + " on L; out " +
                                 std::to_string(out.size()) + "]";
                    }
                    ++isos;
                    return false;
                });
                std::cout << "W " << id << " FAILED no isomorphism carries its incoming curve onto L ("
                          << isos << " isomorphisms onto the thickening; L has " << rowEdges.size()
                          << " edges;" << sizes << ")\n";
                continue;
            }
            KnottedSurface surface(skeleton, carried);
            auto link = farside::orientedOutgoingLink(surface, outgoing, row, incomingBC);
            if (!link) {
                std::cout << "W " << id << " FAILED incoming orientation is inconsistent\n";
                continue;
            }
            knotbuilder::Diagram d = drawer.draw(link->curves);

            // Which row components each surface component meets.
            std::map<size_t, std::vector<size_t>> incoming;
            const auto surfaceOf = surface.boundaryEdgeSurfaceComponent();
            for (const auto &[bc, curves] : surface.orientedBoundaryLinks()) {
                if (bc != incomingBC) continue;
                for (const OrientedCurve &curve : curves) {
                    if (curve.empty()) continue;
                    size_t rowEdge = row.rowIndexOf.at(curve.front().edge->index());
                    incoming[surfaceOf.at(curve.front().edge)].push_back(componentOfRowEdge[rowEdge]);
                }
            }
            std::ostringstream inc;
            inc << '{';
            bool first = true;
            for (const auto &[sc, comps] : incoming) {
                inc << (first ? "" : ",") << '"' << sc << "\":" << list(comps);
                first = false;
            }
            inc << '}';
            std::ostringstream pdOut;
            pdOut << '[';
            for (size_t i = 0; i < d.pd.size(); ++i)
                pdOut << (i ? "," : "") << '[' << d.pd[i][0] << ',' << d.pd[i][1] << ','
                      << d.pd[i][2] << ',' << d.pd[i][3] << ']';
            pdOut << ']';
            std::cout << "W " << id << " ok components=" << d.components
                      << " crossingless=" << list(d.crossingless) << " pd=" << pdOut.str()
                      << " lk=" << matrix(d) << " surface=" << list(link->surfaceComponent)
                      << " incoming=" << inc.str() << "\n";
        } catch (const std::exception &e) {
            std::cout << "W " << id << " FAILED " << e.what() << "\n";
        }
    }
    return 0;
}
