//
//  peripheral.cpp
//
//  Created by John Teague on 09/03/2026.
//

#include "peripheral.h"

#include <algorithm>
#include <cstdlib>
#include <sstream>

namespace peripheral {

namespace {

/**
 * remaining_face[i][j] from engine/snappea/kernel/tables.cpp: the face such
 * that i, j and it are arranged counterclockwise around the missing vertex.
 * intersectionNumber() needs SnapPea's own table here, not a rederivation --
 * the sign of every interior crossing depends on it.
 */
constexpr int remainingFace[4][4] = {{9, 3, 1, 2},
                                     {2, 9, 3, 0},
                                     {3, 0, 9, 1},
                                     {1, 2, 0, 9}};

/**
 * FLOW(A, B) from engine/snappea/kernel/kernel_typedefs.h: given that a curve
 * meets one side of a triangle A times and a different side B times, the
 * number of strands running from the first side to the second.
 */
int flow(int a, int b) {
    if ((a < 0) != (b < 0))
        return ((a < 0) != (a + b < 0)) ? a : -b;
    return 0;
}

/** The next line of `in`, skipping blank ones. Throws at end of input. */
std::string nextNonEmpty(std::istringstream &in) {
    std::string line;
    while (std::getline(in, line)) {
        if (line.find_first_not_of(" \t\r") != std::string::npos)
            return line;
    }
    throw regina::InvalidArgument(
        "peripheral::parseSnapPea(): file ended early");
}

} // namespace

DrilledWithMeridians drillWithMeridians(
    const regina::Triangulation<3> &tri,
    const std::vector<std::vector<const regina::Edge<3> *>> &components) {
    // Every edge in its own edgeVertex direction: the undirected behaviour
    // this function had before meridians carried a meaningful sign.
    std::vector<std::vector<DirectedEdge>> directed;
    directed.reserve(components.size());
    for (const auto &comp : components) {
        std::vector<DirectedEdge> edges;
        edges.reserve(comp.size());
        for (const regina::Edge<3> *e : comp)
            edges.push_back(DirectedEdge{e, false});
        directed.push_back(std::move(edges));
    }
    return drillWithMeridians(tri, directed);
}

DrilledWithMeridians drillWithMeridians(
    const regina::Triangulation<3> &tri,
    const std::vector<std::vector<DirectedEdge>> &components) {
    if (!tri.isOrientable())
        throw regina::InvalidArgument(
            "peripheral::drillWithMeridians(): triangulation is not "
            "orientable, so there is no right-handed sheet to put the "
            "meridians on");

    regina::Triangulation<3> complement(tri);

    // Follow the caller's edges across orient() exactly rather than
    // re-looking them up by index afterwards. orient() flips vertices 2 and 3
    // of precisely those tetrahedra whose orientation() is -1 (see
    // TriangulationBase<dim>::orient(), which swaps vertices dim-1 and dim)
    // and never renumbers tetrahedra, so recording the pre-orient
    // orientations is enough to map every (tetrahedron, local edge)
    // descriptor to its image.
    std::vector<int> preOrientation(complement.size());
    for (size_t i = 0; i < complement.size(); ++i)
        preOrientation[i] = complement.tetrahedron(i)->orientation();

    // Record each edge as (tetrahedron, tail local vertex, head local vertex),
    // with tail/head in the direction the caller's traversal runs, rather than
    // as a bare local edge number. Carrying the two endpoints separately is
    // what lets the direction survive orient() below: a local edge number
    // alone cannot say which way round the edge is being used.
    struct Desc {
        size_t tet;
        int tail;
        int head;
    };
    std::vector<std::vector<Desc>> tracked;
    tracked.reserve(components.size());
    for (const auto &comp : components) {
        if (comp.empty())
            throw regina::InvalidArgument(
                "peripheral::drillWithMeridians(): empty component");
        std::vector<Desc> descs;
        descs.reserve(comp.size());
        for (const DirectedEdge &de : comp) {
            auto emb = de.edge->front();
            // emb.vertices() maps 0, 1 to the local vertices realizing the
            // edge's own vertex(0), vertex(1).
            const int from = emb.vertices()[de.reversed ? 1 : 0];
            const int to = emb.vertices()[de.reversed ? 0 : 1];
            descs.push_back(Desc{emb.tetrahedron()->index(), from, to});
        }
        tracked.push_back(std::move(descs));
    }

    complement.orient();

    for (auto &descs : tracked)
        for (Desc &d : descs) {
            if (preOrientation[d.tet] != -1)
                continue;
            // This tetrahedron had vertices 2 and 3 swapped. Mapping both
            // endpoints (rather than just re-deriving an edge number) keeps
            // the direction right: the swap reverses edge {2,3}, and leaves
            // the direction of every other edge alone.
            auto swap23 = [](int v) { return v == 2 ? 3 : (v == 3 ? 2 : v); };
            d.tail = swap23(d.tail);
            d.head = swap23(d.head);
        }

    // Drill, recording each component's meridian from the two tetrahedra its
    // *first* pinch inserts. pinchEdge() appends exactly two tetrahedra and
    // never removes or renumbers any (the same invariant
    // EdgeComplement::drillTrackingEdges_() relies on), so those two indices
    // stay valid through every later pinch, and so do the descriptors above.
    std::vector<Curve> meridians;
    meridians.reserve(tracked.size());
    for (const auto &descs : tracked) {
        Curve meridian;
        bool first = true;
        for (const Desc &d : descs) {
            const int localEdge = regina::Edge<3>::edgeNumber[d.tail][d.head];
            size_t before = complement.size();
            complement.pinchEdge(complement.tetrahedron(d.tet)->edge(localEdge));
            if (complement.size() != before + 2)
                throw regina::InvalidArgument(
                    "peripheral::drillWithMeridians(): pinchEdge() did not "
                    "insert exactly two tetrahedra");
            if (first) {
                // pinchEdge() lays its annulus out relative to the edge's own
                // edgeVertex direction, so the tabulated core runs the wrong
                // way round exactly when the traversal does.
                const int sign =
                    (regina::Edge<3>::edgeVertex[localEdge][0] == d.tail) ? 1
                                                                          : -1;
                const size_t inserted[2] = {before, before + 1};
                for (const auto &entry : pinchMeridian)
                    meridian.push_back(Crossing{inserted[entry[0]], entry[1],
                                                entry[2], sign * entry[3]});
                first = false;
            }
        }
        meridians.push_back(std::move(meridian));
    }

    return {std::move(complement), std::move(meridians)};
}

std::string snapPeaOriented(const regina::Triangulation<3> &tri) {
    if (!tri.isOriented())
        throw regina::InvalidArgument(
            "peripheral::snapPeaOriented(): triangulation is not oriented");

    // Regina's writer hard-codes "unknown_orientability", which makes
    // SnapPea run its own orient() on load. That flips vertices 2 and 3 of
    // some tetrahedra and would silently invalidate every Crossing we
    // computed, so declare the orientation we already have.
    std::string text = tri.snapPea();
    static const std::string unknown = "unknown_orientability";
    const size_t at = text.find(unknown);
    if (at == std::string::npos)
        throw regina::InvalidArgument(
            "peripheral::snapPeaOriented(): Regina's SnapPea export no longer "
            "declares unknown_orientability; the header patch needs updating");
    text.replace(at, unknown.size(), "oriented_manifold");
    return text;
}

SnapPeaFile parseSnapPea(const std::string &text) {
    std::istringstream in(text);
    SnapPeaFile file;
    std::string line;

    if (!std::getline(in, line) || line.rfind("%", 0) != 0)
        throw regina::InvalidArgument(
            "peripheral::parseSnapPea(): missing '% Triangulation' header");
    std::getline(in, file.name);
    std::getline(in, file.solution);
    std::getline(in, file.orientability);
    std::getline(in, file.chernSimons);
    // Trim the orientability line, which callers compare against.
    file.orientability.erase(file.orientability.find_last_not_of(" \t\r") + 1);

    {
        std::istringstream counts(nextNonEmpty(in));
        if (!(counts >> file.numOrCusps >> file.numNonOrCusps))
            throw regina::InvalidArgument(
                "peripheral::parseSnapPea(): bad cusp counts");
    }
    const int numCusps = file.numOrCusps + file.numNonOrCusps;
    for (int i = 0; i < numCusps; ++i)
        file.cuspLines.push_back(nextNonEmpty(in));

    size_t numTet = 0;
    {
        std::istringstream count(nextNonEmpty(in));
        if (!(count >> numTet))
            throw regina::InvalidArgument(
                "peripheral::parseSnapPea(): bad tetrahedron count");
    }

    file.meridian = CurveField(numTet);
    file.longitude = CurveField(numTet);
    file.tets.resize(numTet);

    for (size_t t = 0; t < numTet; ++t) {
        SnapPeaFile::Tet &tet = file.tets[t];
        {
            std::istringstream nbr(nextNonEmpty(in));
            for (int f = 0; f < 4; ++f)
                if (!(nbr >> tet.neighbour[f]))
                    throw regina::InvalidArgument(
                        "peripheral::parseSnapPea(): bad neighbour row");
        }
        {
            std::istringstream glu(nextNonEmpty(in));
            for (int f = 0; f < 4; ++f)
                if (!(glu >> tet.gluing[f]))
                    throw regina::InvalidArgument(
                        "peripheral::parseSnapPea(): bad gluing row");
        }
        {
            std::istringstream cusp(nextNonEmpty(in));
            for (int v = 0; v < 4; ++v)
                if (!(cusp >> tet.cusp[v]))
                    throw regina::InvalidArgument(
                        "peripheral::parseSnapPea(): bad cusp-index row");
        }
        // Four rows of sixteen: (meridian, longitude) x (right, left) sheets,
        // each row running over vertex then face.
        for (int which = 0; which < 2; ++which)
            for (int sheet = 0; sheet < 2; ++sheet) {
                std::istringstream row(nextNonEmpty(in));
                CurveField &field = which == 0 ? file.meridian : file.longitude;
                for (int v = 0; v < 4; ++v)
                    for (int f = 0; f < 4; ++f)
                        if (!(row >> field(t, sheet, v, f)))
                            throw regina::InvalidArgument(
                                "peripheral::parseSnapPea(): bad curve row");
            }
        tet.shape = nextNonEmpty(in);
    }

    return file;
}

CurveField toField(const Curve &curve, size_t numTet) {
    CurveField field(numTet);
    for (const Crossing &c : curve) {
        if (c.tet >= numTet)
            throw regina::InvalidArgument(
                "peripheral::toField(): crossing names a tetrahedron outside "
                "the triangulation");
        field(c.tet, 0, c.vertex, c.face) += c.sign;
    }
    return field;
}

long intersectionNumber(const SnapPeaFile &file, int cusp, const CurveField &a,
                        const CurveField &b) {
    const size_t numTet = file.tets.size();

    // copy_curves_to_scratch()'s double copy: on a torus cusp, curve[0] is
    // summed across the two sheets and written to both, so that curve[0] and
    // curve[1] can never end up on opposite sheets of the orientation double
    // cover (which would make the intersection number wrong rather than
    // merely negated).
    CurveField scratch(numTet);
    for (size_t t = 0; t < numTet; ++t)
        for (int v = 0; v < 4; ++v)
            for (int f = 0; f < 4; ++f) {
                const int total = a(t, 0, v, f) + a(t, 1, v, f);
                scratch(t, 0, v, f) = total;
                scratch(t, 1, v, f) = total;
            }

    long total = 0;

    // Crossings on the edges. Counted only where curve[0] is entering, so
    // that each edge crossing is counted once rather than once per incident
    // triangle.
    for (size_t t = 0; t < numTet; ++t)
        for (int v = 0; v < 4; ++v) {
            if (file.tets[t].cusp[v] != cusp)
                continue;
            for (int f = 0; f < 4; ++f) {
                if (f == v)
                    continue;
                for (int sheet = 0; sheet < 2; ++sheet)
                    if (scratch(t, sheet, v, f) > 0)
                        total += static_cast<long>(scratch(t, sheet, v, f)) *
                                 b(t, sheet, v, f);
            }
        }

    // Crossings in the interiors of the triangles, counting only the strand
    // of curve[0] running from the current side towards the right; the other
    // strands are picked up by the other values of f.
    for (size_t t = 0; t < numTet; ++t)
        for (int v = 0; v < 4; ++v) {
            if (file.tets[t].cusp[v] != cusp)
                continue;
            for (int f = 0; f < 4; ++f) {
                if (f == v)
                    continue;
                const int onLeft = remainingFace[v][f];
                const int onRight = remainingFace[f][v];
                total += static_cast<long>(
                             flow(scratch(t, 0, v, f), scratch(t, 0, v, onRight))) *
                         b(t, 0, v, onRight);
                total += static_cast<long>(
                             flow(scratch(t, 1, v, f), scratch(t, 1, v, onLeft))) *
                         b(t, 1, v, onLeft);
            }
        }

    return total;
}

std::pair<long, long> slope(const SnapPeaFile &file, int cusp,
                            const CurveField &mu) {
    const long s = intersectionNumber(file, cusp, file.meridian, file.longitude);
    if (s != 1 && s != -1)
        throw regina::InvalidArgument(
            "peripheral::slope(): the file's stored meridian and longitude do "
            "not intersect once, so they are not a peripheral basis");
    const long a = intersectionNumber(file, cusp, mu, file.longitude);
    const long b = intersectionNumber(file, cusp, mu, file.meridian);
    return {a / s, -(b / s)};
}

std::pair<long, long> completeBasis(long a, long b) {
    // Extended Euclid on (a, b): x, y with a*x + b*y = gcd. Setting
    // (c, d) = (-y, x) gives a*d - b*c = a*x + b*y = 1.
    long oldR = a, r = b;
    long oldX = 1, x = 0;
    long oldY = 0, y = 1;
    while (r != 0) {
        const long q = oldR / r;
        std::tie(oldR, r) = std::make_pair(r, oldR - q * r);
        std::tie(oldX, x) = std::make_pair(x, oldX - q * x);
        std::tie(oldY, y) = std::make_pair(y, oldY - q * y);
    }
    if (oldR != 1 && oldR != -1)
        throw regina::InvalidArgument(
            "peripheral::completeBasis(): slope is not primitive, so it is "
            "not a meridian");
    if (oldR == -1) {
        oldX = -oldX;
        oldY = -oldY;
    }
    return {-oldY, oldX};
}

} // namespace peripheral
