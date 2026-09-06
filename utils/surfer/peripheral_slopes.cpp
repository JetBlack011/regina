//
//  peripheral_slopes.cpp
//
//  Created by John Teague on 09/03/2026.
//
//  The C++ half of the far-side identification pass (see peripheral.h and
//  cobordism-atlas/tools/identify_far_sides.py). Everything topological
//  happens here; the SnapPy half only installs a basis and asks whether two
//  manifolds are isometric.
//
//  Two subcommands, both batch (one process for a whole run, not one per
//  witness), reading and writing the record format described below:
//
//    dump    stdin:  "<id> <pairsig>" per line
//            stdout: one RECORD per boundary component of each surface,
//                    carrying the drilled complement as a SnapPea data file
//                    with all curves zero, plus our SIGNED meridians
//
//    dump-link
//            stdin:  "<id> <pd code>" per line
//            stdout: the same records, built from a diagram instead of from a
//                    witness surface -- so the reference table and the far
//                    sides are drilled by one code path with one sign
//                    convention, and no comparison rests on the two agreeing
//                    by luck
//
//    slope   stdin:  those records, with SnapPea's basis now installed in the
//                    TRI block (the Python half round-trips them through
//                    snappy.Triangulation(..., remove_finite_vertices=False))
//            stdout: "SLOPE <id> <bc> <component> <cusp> <a> <b> <c> <d>",
//                    where (a, b) is our meridian in SnapPea's basis and
//                    [(a,b), (c,d)] is the change-of-basis matrix to install
//
//  Record format (blank-line free, so it survives being piped):
//
//    RECORD <id> <boundary component> <components>
//    ORIENT <component> <surface component index, or -1 if unknown>
//    MERIDIAN <component> <count> [<tet> <vertex> <face> <sign>]*
//    TRI <line count>
//    <that many lines>
//    ENDRECORD
//
//  The ORIENT lines are what make the meridian signs readable. A meridian is
//  oriented by lk(mu, K) = +1, so its sign is meaningful only once the
//  component is directed; dump takes that direction from the surface's own
//  orientation (KnottedSurface::orientedBoundaryLinks()). But that orientation
//  is chosen per CONNECTED COMPONENT of the surface, so two meridians' signs
//  are mutually meaningful exactly when their curves bound the same surface
//  component. ORIENT reports that component per link component, so the reader
//  can tell "these signs are related" from "these are independently flippable"
//  instead of assuming the first. For dump-link, where the diagram directs
//  every component and nothing is independently flippable, every ORIENT line
//  reads 0.
//
//  A boundary component whose drilling or meridian construction fails is
//  reported as
//
//    FAILED <id> <boundary component> <reason>
//
//  rather than skipped, so the Python half can account for every input.

#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include <map>

#include "embeddedsubmanifold.h"
#include "knotbuilder.h"
#include "linkcomplement.h"
#include "pairsig.h"
#include "peripheral.h"

namespace {

void usage(const char *progName) {
    std::cerr << "Usage:\n"
              << "    " << progName << " dump      < ids-and-pairsigs\n"
              << "    " << progName << " dump-link < ids-and-pd-codes\n"
              << "    " << progName << " slope     < records\n"
              << "    " << progName << " sig       < ids-and-pairsigs\n";
    exit(1);
}

/** Writes one record's TRI block, which is the bulk of the output. */
void writeTri(std::ostream &out, const std::string &text) {
    std::vector<std::string> lines;
    std::istringstream in(text);
    std::string line;
    while (std::getline(in, line))
        lines.push_back(line);
    out << "TRI " << lines.size() << "\n";
    for (const std::string &l : lines)
        out << l << "\n";
}

// Emits the isomorphism signature of each boundary component's complement,
// as "SIG <id> <boundary-component> <components> <isosig>".
//
// A far side reaches cobordisms.csv as a NAME -- often a census name that no
// literature table knows -- and a name is not a signature, so it cannot be
// looked up in a table of complements we built ourselves. This recovers the
// signature from the pair signature, which every witness carries, so a far
// side named only "L109021" or "m129 : #3" can still be matched against
// cobordism-atlas/results/reference_complements.json.
//
// Unlike dump, this simplifies: the caller wants the same canonical form
// Census::lookup() and the reference table use, not the drilled
// triangulation with its peripheral curves intact.
int runSig() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty())
            continue;
        std::istringstream fields(line);
        std::string id, sig;
        if (!(fields >> id >> sig)) {
            std::cerr << "skipping malformed input line\n";
            continue;
        }

        DecodedKnottedSurfaceSig decoded;
        try {
            decoded = fromKnottedSurfaceSig(sig);
        } catch (const std::exception &e) {
            std::cout << "FAILED " << id << " - decode:" << e.what() << "\n";
            continue;
        }

        for (const auto &[bc, link] : decoded.surface->boundaryLinks()) {
            try {
                regina::Triangulation<3> complement = link.buildComplement();
                complement.simplify();
                std::cout << "SIG " << id << ' ' << bc << ' '
                          << link.countComponents() << ' '
                          << complement.isoSig() << "\n";
            } catch (const std::exception &e) {
                std::cout << "FAILED " << id << ' ' << bc << " drill:"
                          << e.what() << "\n";
            }
        }
    }
    return 0;
}

/** Writes one record: the drilled complement, its meridians, and provenance. */
void writeRecord(const std::string &id, size_t bc,
                 const peripheral::DrilledWithMeridians &drilled,
                 const std::vector<long> &surfaceComponent) {
    const std::string tri = peripheral::snapPeaOriented(drilled.tri);
    std::cout << "RECORD " << id << ' ' << bc << ' '
              << drilled.meridians.size() << "\n";
    for (size_t c = 0; c < drilled.meridians.size(); ++c)
        std::cout << "ORIENT " << c << ' '
                  << (c < surfaceComponent.size() ? surfaceComponent[c] : -1)
                  << "\n";
    for (size_t c = 0; c < drilled.meridians.size(); ++c) {
        const peripheral::Curve &mu = drilled.meridians[c];
        std::cout << "MERIDIAN " << c << ' ' << mu.size();
        for (const peripheral::Crossing &x : mu)
            std::cout << ' ' << x.tet << ' ' << x.vertex << ' ' << x.face
                      << ' ' << x.sign;
        std::cout << "\n";
    }
    writeTri(std::cout, tri);
    std::cout << "ENDRECORD\n";
}

/**
 * Puts each of `link`'s components into the direction `curves` gives it.
 *
 * The surface's own orientation directs its boundary curves, but it does so in
 * the boundary component's numbering and in whatever order the curves came
 * out; Link groups the same edges into its own components in its own order.
 * This matches the two up by edge identity -- every edge of a curve belongs to
 * exactly one Link component -- so the meridians come back in Link's component
 * order with the surface's directions on them.
 *
 * Returns false if any component is left undirected, which means the two views
 * disagree about the edge set and no signed answer should be reported at all.
 */
bool directComponents(
    const Link &link, const std::vector<OrientedCurve> &curves,
    const std::map<const regina::Edge<3> *, size_t> &surfaceOf,
    std::vector<std::vector<peripheral::DirectedEdge>> &directions,
    std::vector<long> &surfaceComponent) {
    std::map<const regina::Edge<3> *, bool> reversedOf;
    for (const OrientedCurve &curve : curves)
        for (const OrientedEdge &oe : curve)
            reversedOf[oe.edge] = oe.reversed;

    directions.assign(link.countComponents(), {});
    surfaceComponent.assign(link.countComponents(), -1);
    for (int c = 0; c < link.countComponents(); ++c) {
        for (const regina::Edge<3> *e : link.comps_[c].edges()) {
            auto found = reversedOf.find(e);
            if (found == reversedOf.end())
                return false;
            directions[c].push_back({e, found->second});
            auto owner = surfaceOf.find(e);
            if (owner != surfaceOf.end())
                surfaceComponent[c] = static_cast<long>(owner->second);
        }
    }
    return true;
}

int runDump() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty())
            continue;
        std::istringstream fields(line);
        std::string id, sig;
        if (!(fields >> id >> sig)) {
            std::cerr << "skipping malformed input line\n";
            continue;
        }

        DecodedKnottedSurfaceSig decoded;
        try {
            decoded = fromKnottedSurfaceSig(sig);
        } catch (const std::exception &e) {
            std::cout << "FAILED " << id << " - decode:" << e.what() << "\n";
            continue;
        }

        // The surface's own orientation is what signs the meridians. It is
        // only defined up to a flip per connected component of the surface,
        // which is why the component index is reported alongside rather than
        // quietly dropped.
        std::vector<std::pair<size_t, std::vector<OrientedCurve>>> oriented;
        std::map<const regina::Edge<3> *, size_t> surfaceOf;
        bool haveDirections = false;
        try {
            oriented = decoded.surface->orientedBoundaryLinks();
            surfaceOf = decoded.surface->boundaryEdgeSurfaceComponent();
            haveDirections = true;
        } catch (const std::exception &e) {
            std::cout << "FAILED " << id << " - orient:" << e.what() << "\n";
        }

        for (const auto &[bc, link] : decoded.surface->boundaryLinks()) {
            try {
                std::vector<std::vector<peripheral::DirectedEdge>> directions;
                std::vector<long> surfaceComponent;
                bool directed = false;
                if (haveDirections) {
                    for (const auto &[obc, curves] : oriented) {
                        if (obc != bc)
                            continue;
                        directed = directComponents(link, curves, surfaceOf,
                                                    directions,
                                                    surfaceComponent);
                        break;
                    }
                }

                // Undirected is a strictly weaker answer, not a wrong one:
                // the meridians come back signed only up to an independent
                // +- per component, and the -1 in ORIENT says so. Falling
                // back silently to a *signed-looking* result would be the
                // unsound thing.
                peripheral::DrilledWithMeridians drilled =
                    directed ? link.buildComplementWithPeripheral(directions)
                             : link.buildComplementWithPeripheral();
                if (!directed)
                    surfaceComponent.assign(link.countComponents(), -1);

                writeRecord(id, bc, drilled, surfaceComponent);
            } catch (const std::exception &e) {
                std::cout << "FAILED " << id << ' ' << bc << " drill:"
                          << e.what() << "\n";
            }
        }
    }
    return 0;
}

/**
 * Records built from a diagram rather than from a witness surface.
 *
 * The reference table has to be drilled by the same code with the same sign
 * convention as the far sides it will be compared against, or the comparison
 * measures the difference between two conventions rather than between two
 * links. A PD code directs every component outright (knotbuilder returns the
 * traversal direction per edge), so nothing here is independently flippable
 * and every ORIENT line reads 0.
 */
int runDumpLink() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.empty())
            continue;
        const size_t space = line.find(' ');
        if (space == std::string::npos) {
            std::cerr << "skipping malformed input line\n";
            continue;
        }
        const std::string id = line.substr(0, space);
        const std::string pd = line.substr(space + 1);

        try {
            knotbuilder::TriangulationWithLink built =
                knotbuilder::buildLink(knotbuilder::parsePDCode(pd));
            Link link(built.tri, built.edges);

            // knotbuilder hands back edges and their traversal directions as
            // parallel arrays over the whole diagram; regroup them into Link's
            // own component order by edge identity.
            std::map<const regina::Edge<3> *, bool> reversedOf;
            for (size_t i = 0; i < built.edges.size(); ++i)
                reversedOf[built.edges[i]] =
                    i < built.reversed.size() ? built.reversed[i] : false;

            std::vector<std::vector<peripheral::DirectedEdge>> directions(
                link.countComponents());
            for (int c = 0; c < link.countComponents(); ++c)
                for (const regina::Edge<3> *e : link.comps_[c].edges())
                    directions[c].push_back({e, reversedOf[e]});

            peripheral::DrilledWithMeridians drilled =
                link.buildComplementWithPeripheral(directions);
            std::vector<long> surfaceComponent(link.countComponents(), 0);
            writeRecord(id, 0, drilled, surfaceComponent);
        } catch (const std::exception &e) {
            std::cout << "FAILED " << id << " 0 build:" << e.what() << "\n";
        }
    }
    return 0;
}

int runSlope() {
    std::string line;
    while (std::getline(std::cin, line)) {
        if (line.rfind("RECORD ", 0) != 0)
            continue;

        std::istringstream header(line.substr(7));
        std::string id;
        size_t bc = 0, numComponents = 0;
        header >> id >> bc >> numComponents;

        std::vector<peripheral::Curve> meridians(numComponents);
        std::string tri;
        bool ok = true;

        while (std::getline(std::cin, line) && line != "ENDRECORD") {
            if (line.rfind("MERIDIAN ", 0) == 0) {
                std::istringstream in(line.substr(9));
                size_t comp = 0, count = 0;
                in >> comp >> count;
                if (comp >= meridians.size()) {
                    ok = false;
                    continue;
                }
                for (size_t k = 0; k < count; ++k) {
                    peripheral::Crossing x{};
                    in >> x.tet >> x.vertex >> x.face >> x.sign;
                    meridians[comp].push_back(x);
                }
            } else if (line.rfind("TRI ", 0) == 0) {
                size_t count = std::stoul(line.substr(4));
                std::ostringstream body;
                for (size_t k = 0; k < count; ++k) {
                    std::string l;
                    if (!std::getline(std::cin, l)) {
                        ok = false;
                        break;
                    }
                    body << l << "\n";
                }
                tri = body.str();
            }
        }

        if (!ok || tri.empty()) {
            std::cout << "FAILED " << id << ' ' << bc << " record:truncated\n";
            continue;
        }

        try {
            peripheral::SnapPeaFile file = peripheral::parseSnapPea(tri);
            for (size_t c = 0; c < meridians.size(); ++c) {
                // Which cusp this component became is SnapPea's choice, not
                // ours: we hand it a file with no cusp indices and it works
                // them out. Read the answer back off the meridian's own
                // crossings, which must all sit on one cusp.
                int cusp = -1;
                bool consistent = true;
                for (const peripheral::Crossing &x : meridians[c]) {
                    if (x.tet >= file.tets.size()) {
                        consistent = false;
                        break;
                    }
                    const int here = file.tets[x.tet].cusp[x.vertex];
                    if (cusp == -1)
                        cusp = here;
                    else if (cusp != here)
                        consistent = false;
                }
                if (!consistent || cusp < 0) {
                    std::cout << "FAILED " << id << ' ' << bc
                              << " cusp:meridian " << c
                              << " does not lie on a single cusp\n";
                    continue;
                }

                auto [a, b] = peripheral::slope(
                    file, cusp,
                    peripheral::toField(meridians[c], file.tets.size()));
                auto [cc, d] = peripheral::completeBasis(a, b);
                std::cout << "SLOPE " << id << ' ' << bc << ' ' << c << ' '
                          << cusp << ' ' << a << ' ' << b << ' ' << cc << ' '
                          << d << "\n";
            }
        } catch (const std::exception &e) {
            std::cout << "FAILED " << id << ' ' << bc << " slope:" << e.what()
                      << "\n";
        }
    }
    return 0;
}

} // namespace

int main(int argc, char *argv[]) {
    std::ios::sync_with_stdio(false);
    if (argc != 2)
        usage(argv[0]);

    const std::string mode = argv[1];
    if (mode == "dump")
        return runDump();
    if (mode == "dump-link")
        return runDumpLink();
    if (mode == "slope")
        return runSlope();
    if (mode == "sig")
        return runSig();
    usage(argv[0]);
    return 1;
}
