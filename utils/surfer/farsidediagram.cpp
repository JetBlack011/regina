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
#include <map>
#include <optional>
#include <unordered_map>
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
#include "farsideredraw.h"
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
    const farside::WitnessRedrawer redraw(argv[arg], layers);
    const knotbuilder::TriangulationWithLink &built = redraw.built();

    // The row's own components, in cyclesOf() order, and each row edge's component.
    const auto &rowCycles = redraw.rowCycles();
    std::vector<size_t> componentOfRowEdge(built.edges.size());
    {
        std::unordered_map<size_t, size_t> compOfT;
        for (size_t c = 0; c < rowCycles.size(); ++c)
            for (const auto &de : rowCycles[c]) compOfT[de.edge] = c;
        for (size_t i = 0; i < built.edges.size(); ++i)
            componentOfRowEdge[i] = compOfT.at(built.edges[i]->index());
    }
    const knotbuilder::Diagram rowDiagram = redraw.drawer().draw(rowCycles);
    std::cout << "ROW components=" << rowCycles.size() << " lk=" << matrix(rowDiagram) << "\n";

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string id, sig;
        if (!(in >> id >> sig)) continue;
        try {
            std::string why;
            std::optional<std::vector<int>> carried = redraw.carry(sig, why);
            if (!carried) {
                std::cout << "W " << id << " FAILED " << why << "\n";
                continue;
            }
            KnottedSurface surface(redraw.skeleton(), *carried);
            auto link = farside::orientedOutgoingLink(surface, redraw.outgoing(), redraw.row(),
                                                      redraw.incomingBC());
            if (!link) {
                std::cout << "W " << id << " FAILED incoming orientation is inconsistent\n";
                continue;
            }
            knotbuilder::Diagram d = redraw.drawer().draw(link->curves);

            // Which row components each surface component meets.
            std::map<size_t, std::vector<size_t>> incoming;
            const auto surfaceOf = surface.boundaryEdgeSurfaceComponent();
            for (const auto &[bc, curves] : surface.orientedBoundaryLinks()) {
                if (bc != redraw.incomingBC()) continue;
                for (const OrientedCurve &curve : curves) {
                    if (curve.empty()) continue;
                    size_t rowEdge = redraw.row().rowIndexOf.at(curve.front().edge->index());
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
