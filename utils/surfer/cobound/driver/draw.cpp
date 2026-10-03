//
//  draw.cpp
//
//  cobound draw (was farsidediagram): oriented diagrams of stored
//  cobordisms' outgoing links, from their pair signatures.
//
//  Usage:
//    cobound draw [--layers N] [--gauss] [--faces] '<row PD code>' < pairsigs
//
//  The flags are farsidediagram's, which the cascade's checker and recorder
//  pass; each is a config key (layers, draw_gauss, draw_faces, draw_pairsig,
//  pair_sig_cache), so --config and --set work too.
//
//  --gauss appends each diagram's signed Gauss data (signs=, gauss=) to the
//  ROW line and to every W line: per component, in curve order, the
//  crossings it passes. cascade_check.py reads it; without the flag the
//  output is unchanged. The ROW line then also carries build=, the
//  thickening's digest (WitnessRedrawer::buildChecksum()).
//
//  --faces reads "<id> <f1,f2,...>" lines instead: a surface given by its
//  triangles in this row's thickening, as a cascade certificate records an
//  in-process witness. It is rebuilt face by face with the search's own
//  checks (WitnessRedrawer::rebuild()) before anything is read from it.
//
//  --pairsig (with --faces) appends each rebuilt surface's pair signature,
//  as verifyslicegenus would record it (pairsig=, over the row's own
//  thickening), whether it is connected (connected=) and its resolved
//  vertices (resolved=): what turns a certificate's in-process witness into
//  an ordinary cobordisms.csv row. --sig-cache DIR (with --pairsig) keeps the
//  thickening's own part of the signature -- 99.8% of its cost -- in DIR by
//  the thickening's digest, checked on every use (pairsig.h, Detail), so a
//  row met again costs a fraction of a second instead of tens.
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
//  what the search itself would have had: outgoing::orientedOutgoingLink()
//  and knotbuilder::DiagramDrawer. The isomorphism is the only search in
//  the pipeline, and it is pinned by L x {0}, so the outgoing side is read
//  through the thickening's own product structure -- never mirrored by an
//  automorphism of T.
//

#include <algorithm>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <iterator>
#include <map>
#include <optional>
#include <unordered_map>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "diagramtriangulation/thickening/collar.h"
#include "surfer/submanifold/submanifold.h"
#include "cobound/outgoing/outgoinglink.h"
#include "cobound/outgoing/fromdatabase.h"
#include "diagramtriangulation/todiagram.h"
#include "diagramtriangulation/fromdiagram.h"
#include "diagramtriangulation/pdcode.h"
#include "cobound/json.h"
#include "cobound/cobordisms/pairsigner.h"
#include "surfer/pairsig/pairsig.h"
#include "surfer/submanifold/skeleton.h"
#include "surfer/submanifold/vertexlinks.h"
#include "cobound/driver/commands.h"
#include "cobound/driver/config.h"
#include "cobound/frozen.h"

namespace {

std::string matrix(const knotbuilder::Diagram &d) {
    std::vector<std::vector<long>> lk(d.components, std::vector<long>(d.components, 0));
    for (size_t i = 0; i < d.components; ++i)
        for (size_t j = 0; j < d.components; ++j)
            if (i != j) lk[i][j] = d.linkingNumber(i, j);
    return json::matrix(lk);
}

template <typename T>
std::string list(const std::vector<T> &v) {
    return json::array(v);
}

// Signed Gauss data (--gauss): the crossing signs, then per component the
// crossings it passes in order (+(k+1) over crossing k, -(k+1) under), which
// is the whole oriented diagram and names every curve by its index.
std::string gaussFields(const knotbuilder::Diagram &d) {
    std::ostringstream o;
    o << " signs=[";
    for (size_t k = 0; k < d.crossings.size(); ++k) o << (k ? "," : "") << d.crossings[k].sign;
    o << "] gauss=[";
    for (size_t c = 0; c < d.gauss.size(); ++c) o << (c ? "," : "") << list(d.gauss[c]);
    return o.str() + "]";
}

// Each far-side curve's edges of T, sorted (--gauss): what identifies a curve
// across two reads that list the curves in different orders.
std::string curveEdges(const std::vector<knotbuilder::EdgeCycle> &curves) {
    std::ostringstream o;
    o << " edges=[";
    for (size_t c = 0; c < curves.size(); ++c) {
        std::vector<size_t> es;
        for (const auto &de : curves[c]) es.push_back(de.edge);
        std::sort(es.begin(), es.end());
        o << (c ? "," : "") << list(es);
    }
    return o.str() + "]";
}

} // namespace

int commands::draw(const std::vector<std::string> &args) {
    const char *usage = "usage: cobound draw [--layers N] [--gauss] "
                        "[--faces [--pairsig [--sig-cache DIR]]] '<row PD code>' < pairsigs (or faces)\n";
    int layers = 2;
    bool gauss = false, facesInput = false, pairsigOut = false;
    std::string sigCache;
    std::vector<std::string> positional;
    try {
        const config::Config cfg = config::forCommand(
            "draw", config::Context::draw, args,
            {{"--layers", "layers"},
             {"--gauss", "draw_gauss", false, "1"},
             {"--faces", "draw_faces", false, "1"},
             {"--pairsig", "draw_pairsig", false, "1"},
             {"--sig-cache", "pair_sig_cache"}},
            &positional);
        layers = static_cast<int>(cfg.integer("layers"));
        gauss = cfg.flag("draw_gauss");
        facesInput = cfg.flag("draw_faces");
        pairsigOut = cfg.flag("draw_pairsig");
        sigCache = cfg.text("pair_sig_cache");
    } catch (const config::Error &e) {
        std::cerr << e.what() << "\n" << usage;
        return 2;
    }
    if (positional.size() != 1 || layers < 1 || (pairsigOut && !facesInput) ||
        (!sigCache.empty() && !pairsigOut)) {
        std::cerr << usage;
        return 2;
    }
    const outgoing::WitnessRedrawer redraw(positional.front(), layers);
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
    std::cout << kFrozenDrawRowLine << rowCycles.size() << " lk=" << matrix(rowDiagram)
              << (gauss ? gaussFields(rowDiagram) + " build=" + redraw.buildChecksum()
                        : std::string())
              << "\n";

    // Pair signatures (--pairsig): the thickening's own part is computed once,
    // at the first surface, as pairSigsOf() does.
    std::unique_ptr<PairSigContext<4, 2>> sigContext;
    // With --sig-cache, the ambient part is kept per thickening in the one
    // context cache (pairSigContextFor(): PairSigContext::cached(), a
    // .pairsigctx file per ambient, checked on every use). Before phase 3
    // this tool kept its own <build>.detail files instead; those are no
    // longer read, which costs a rebuild once per thickening.
    auto makeContext = [&] {
        bool loaded = false;
        sigContext = cobordisms::pairSigContextFor(redraw.thickening(), sigCache, 1, &loaded);
        if (!sigCache.empty())
            std::cerr << (loaded ? "sig-cache hit " : "sig-cache miss ")
                      << redraw.buildChecksum() << "\n";
    };
    std::vector<int> currentFaces;

    // One witness's W line, from its surface in the thickening.
    auto describe = [&](const std::string &id, KnottedSurface &surface) {
            auto link = outgoing::orientedOutgoingLink(surface, redraw.outgoing(), redraw.row(),
                                                      redraw.incomingBC());
            if (!link) {
                std::cout << "W " << id << " FAILED incoming orientation is inconsistent\n";
                return;
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
            std::cout << "W " << id << " ok components=" << d.components
                      << " crossingless=" << list(d.crossingless) << " pd="
                      << knotbuilder::formatPDCode(d.pd, knotbuilder::PDSpelling::commas)
                      << " lk=" << matrix(d) << " surface=" << list(link->surfaceComponent)
                      << " incoming=" << inc.str()
                      << (gauss ? gaussFields(d) + curveEdges(link->curves) + " genus=" +
                                      std::to_string(KnottedSurface::tubedSurfaceType(
                                                         surface.triangulation())
                                                         .genus)
                                : std::string());
            if (pairsigOut) {
                if (!sigContext) makeContext();
                std::cout << " resolved=" << surface.singularVertexCount()
                          << " connected=" << (surface.triangulation().isConnected() ? 1 : 0)
                          << " pairsig=" << sigContext->sig(currentFaces);
            }
            std::cout << "\n";
    };

    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string id, data;
        if (!(in >> id >> data)) continue;
        try {
            std::string why;
            if (facesInput) {
                // Faces of thickening(), rebuilt with the search's own checks.
                std::vector<int> faces;
                std::istringstream fs(data);
                for (std::string f; std::getline(fs, f, ',');) faces.push_back(std::stoi(f));
                PetalCache cache;
                KnottedSurface::SelfIntersectionOptions options;
                options.resolveUnlinked = true;
                KnottedSurface surface(options, redraw.skeleton(), cache);
                if (!redraw.rebuild(faces, surface, why)) {
                    std::cout << "W " << id << " FAILED " << why << "\n";
                    continue;
                }
                currentFaces = faces;
                describe(id, surface);
            } else {
                std::optional<std::vector<int>> carried = redraw.carry(data, why);
                if (!carried) {
                    std::cout << "W " << id << " FAILED " << why << "\n";
                    continue;
                }
                KnottedSurface surface(redraw.skeleton(), *carried);
                describe(id, surface);
            }
        } catch (const std::exception &e) {
            std::cout << "W " << id << " FAILED " << e.what() << "\n";
        }
    }
    return 0;
}
