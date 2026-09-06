//
//  peripheral.h
//
//  Created by John Teague on 09/03/2026.
//

#ifndef PERIPHERAL_H

#define PERIPHERAL_H

#include <array>
#include <string>
#include <utility>
#include <vector>

#include <triangulation/dim3.h>

/*! \file utils/surfer/peripheral.h
 *  \brief Drills a link out of a triangulated 3-manifold while retaining its
 *  meridians, and expresses them as slopes in SnapPea's peripheral basis.
 *
 *  \section per_why Why this exists
 *
 *  identify::identify() names a link by its *complement*, and a complement
 *  does not determine a link: Rolfsen twisting along an unknotted component
 *  changes the link while preserving the exterior. Knots are exempt
 *  (Gordon-Luecke); links are not, and the ambiguity is real in our own data
 *  (L7n2, L5a1, L8n2, L9n3 and L10n9 all share one exterior, with slice
 *  genus 0, 1, 0, 1 and 1 respectively).
 *
 *  A link *is* determined by its complement together with the oriented
 *  peripheral system, since filling along the meridians recovers (S^3, L).
 *  EdgeComplement::buildComplement() throws the meridians away; this file
 *  keeps them.
 *
 *  \section per_frame The common frame
 *
 *  Naming the far side means comparing against SnapPy, which needs our
 *  meridian as a pair of coefficients against *some* basis. Which basis is
 *  irrelevant -- SnapPea's is arbitrary (the meridian is (0,+-1) for the
 *  trefoil and (+-1,0) for the figure-8) -- but our meridian and the
 *  reference basis must live in a *common combinatorial frame*.
 *
 *  \warning Do **not** obtain that frame by wrapping the drilled complement
 *  in regina::SnapPeaTriangulation. That constructor calls SnapPea's
 *  remove_finite_vertices() whenever
 *  `countVertices() > countBoundaryComponents()`
 *  (engine/snappea/snappeatriangulation.cpp), which a drilled complement
 *  always satisfies -- the ambient triangulation's non-link vertices survive
 *  the drilling as finite vertices. The triangulation is then retriangulated
 *  (measured: 30 tetrahedra in, 3 out) and the peripheral curves it reports
 *  describe a manifold we can no longer point at. The trap is that the
 *  condition is *false* for regina::Example3::whiteheadLink(), so the obvious
 *  spot check passes.
 *
 *  Instead we go through Triangulation<3>::snapPea(), which writes our
 *  triangulation verbatim with all curves zero, and let SnapPea install a
 *  basis into *that* file (peripheral_curves_as_needed() only touches cusps
 *  whose curves are all zero). See tools/identify_far_sides.py for the
 *  round trip.
 */

namespace peripheral {

/**
 * One signed crossing of a curve with one side of one triangle of a cusp
 * cross-section, in the coordinates SnapPea's data files use.
 *
 * The cusp cross-section is triangulated by the corners of the ambient
 * tetrahedra: tetrahedron \a tet contributes one triangle at each ideal
 * vertex \a vertex, whose three sides are numbered by the tetrahedron faces
 * other than \a vertex. A positive \a sign means the curve is *entering* the
 * triangle across side \a face, negative that it is leaving.
 */
struct Crossing {
    size_t tet;  /**< Ambient tetrahedron index. */
    int vertex;  /**< 0-3: which ideal vertex's corner triangle. */
    int face;    /**< 0-3, != vertex: which side of that corner. */
    int sign;    /**< +1 entering, -1 leaving. */
};

/** A closed curve on a cusp cross-section, as its crossings. */
using Curve = std::vector<Crossing>;

/**
 * The core of the annulus that Triangulation<3>::pinchEdge() leaves behind,
 * as (which inserted tetrahedron, vertex, face, sign).
 *
 * pinchEdge() replaces a triangular face containing \a e with a
 * two-tetrahedron "pinched ball" -- Regina's own words, "a 3-ball in which
 * some internal curve joining two distinct boundary points is collapsed to a
 * point, whose link then becomes an annulus". The collapsed curve is parallel
 * to \a e, so the core of that annulus encircles the drilled component
 * exactly once: it is a meridian.
 *
 * The two inserted tetrahedra are glued to each other by a *constant*
 * pattern (`t0 f0 <-> t1 f0` by (1 2), `t0 f3 <-> t1 f3` by (0 1), and
 * `t1 f1 <-> t1 f2` by (1 2)), so the annulus and its core are the same
 * every time and can be tabulated once. Working the identifications through,
 * the pinch point is `t0` vertices 0, 1, 2 together with `t1` vertices 0, 1,
 * 2, and the annulus is those six corners; of the two independent cycles in
 * its dual graph, the self-loop at `(t1, 0)` is null-homologous and the
 * 3-cycle below is the core.
 *
 * Entries are `{which inserted tetrahedron (0 = t0, 1 = t1), vertex, face,
 * sign}`. The overall sign is a convention: only +-meridian is needed to
 * decide whether two links are the same, and pinning the sign down (which the
 * *oriented* comparison additionally needs) is the job of the linking-number
 * regression, not of this table.
 */
inline constexpr std::array<std::array<int, 4>, 6> pinchMeridian = {{
    {0, 2, 3, +1}, {0, 2, 0, -1},
    {1, 1, 0, +1}, {1, 1, 2, -1},
    {1, 2, 1, +1}, {1, 2, 3, -1},
}};

/**
 * A drilled complement together with one meridian per drilled component.
 *
 * \a tri is oriented, ideal at each drilled component, and deliberately
 * **not** simplified: simplify() would destroy the tetrahedra the meridians
 * are expressed in. It is large, and meant to be transient -- SnapPy carries
 * peripheral data through its own simplification once the curves are
 * installed, so nothing downstream needs to keep it.
 */
struct DrilledWithMeridians {
    regina::Triangulation<3> tri;  /**< The drilled complement. */
    std::vector<Curve> meridians;  /**< One per component, input order. */
};

/**
 * One edge of a link component, together with the direction a traversal of
 * that component runs along it.
 *
 * \a reversed means the traversal runs `edge->vertex(1)` -> `edge->vertex(0)`,
 * the same convention as knotbuilder::TriangulationWithLink::reversed and
 * OrientedEdge in embeddedsubmanifold.h, so a PD-tagged diagram edge or a
 * surface's own induced boundary direction can be handed here unchanged.
 */
struct DirectedEdge {
    const regina::Edge<3> *edge;
    bool reversed;
};

/**
 * Drills every edge of every component, recording each component's meridian
 * from the two tetrahedra its *first* pinch inserts.
 *
 * The ambient triangulation is copied and oriented before anything else, so
 * that the curves land on the right-handed sheet SnapPea expects. orient()
 * flips vertices 2 and 3 of exactly those tetrahedra whose orientation() is
 * -1 and never renumbers tetrahedra, so the caller's edges are followed
 * across it exactly rather than re-looked-up by index.
 *
 * \section per_sign What the sign means
 *
 * A meridian's orientation is not free: with the ambient orientation fixed,
 * orienting a component \a K orients its meridian by `lk(mu, K) = +1`, and
 * reversing \a K negates \a mu. So the *directions* in \a components are what
 * make the returned Curves signed, and a link is only pinned down as an
 * ORIENTED link once they are.
 *
 * pinchEdge() builds its annulus relative to the drilled edge's own
 * `Edge<3>::edgeVertex` direction, so this negates a component's Curve exactly
 * when the traversal runs against that direction on the edge the meridian is
 * read from. Two corrections are folded in: the caller's \a reversed flag, and
 * the fact that orient()'s vertex 2/3 swap reverses edge {2,3} (and only that
 * edge) in a negatively oriented tetrahedron.
 *
 * One residual convention is deliberately NOT pinned down: whether the
 * tabulated core runs with or against the edge. Getting it backwards negates
 * every component's meridian at once, which is exactly the ambiguity between a
 * link and its reverse (or its mirror), and both preserve the slice genus and
 * land on the same LinkInfo orientation variant -- that notation already
 * records orientations relative to the first component. What matters, and what
 * this signs correctly, is the meridians' signs *relative to each other*.
 *
 * \pre Every component is a cycle of distinct edges of \a tri, traversed
 * head-to-tail in the given order, and distinct components are disjoint.
 *
 * \exception regina::InvalidArgument \a tri is not orientable, or a component
 * is empty.
 */
DrilledWithMeridians drillWithMeridians(
    const regina::Triangulation<3> &tri,
    const std::vector<std::vector<DirectedEdge>> &components);

/**
 * As above, with every edge taken in its own `Edge<3>::edgeVertex` direction.
 *
 * The components carry no direction, so neither do the meridians: each one is
 * `+-mu` with the sign fixed by how Regina happens to number the first edge
 * drilled, independently per component. That is enough to decide whether two
 * links are the same UNORIENTED link (test `mu` against `+-mu`, which is what
 * SnapPy's Isometry::extends_to_link() does) and is not enough to decide
 * orientation. Prefer the directed overload wherever the direction is known.
 */
DrilledWithMeridians drillWithMeridians(
    const regina::Triangulation<3> &tri,
    const std::vector<std::vector<const regina::Edge<3> *>> &components);

/**
 * Triangulation<3>::snapPea() with `oriented_manifold` declared in place of
 * `unknown_orientability`.
 *
 * Regina always writes the latter, which makes SnapPea run its own orient()
 * on load; that flips vertices 2 and 3 of some tetrahedra and would silently
 * invalidate every Crossing we computed. Declaring the orientation we already
 * established keeps SnapPea's hands off the labelling.
 *
 * \pre \a tri is oriented (Triangulation<3>::isOriented()).
 */
std::string snapPeaOriented(const regina::Triangulation<3> &tri);

/**
 * A curve in SnapPea's crossing-number coordinates, indexed
 * `(tetrahedron, sheet, vertex, face)`.
 *
 * Sheet 0 is right_handed, sheet 1 left_handed; for an oriented manifold
 * everything lives on sheet 0 and sheet 1 stays zero.
 */
class CurveField {
  private:
    std::vector<std::array<int, 32>> data_; /**< [sheet*16 + vertex*4 + face] */

  public:
    /** All-zero over `numTet` tetrahedra. */
    explicit CurveField(size_t numTet) : data_(numTet, std::array<int, 32>{}) {}

    /** The number of tetrahedra this field spans. */
    size_t size() const { return data_.size(); }

    /** The crossing number at `(tet, sheet, vertex, face)`. */
    int &operator()(size_t tet, int sheet, int v, int f) {
        return data_[tet][sheet * 16 + v * 4 + f];
    }

    /** The crossing number at `(tet, sheet, vertex, face)`. */
    int operator()(size_t tet, int sheet, int v, int f) const {
        return data_[tet][sheet * 16 + v * 4 + f];
    }
};

/** Everything peripheral we read back out of a SnapPea data file. */
struct SnapPeaFile {
    /** One tetrahedron's block. */
    struct Tet {
        std::array<int, 4> neighbour;      /**< Adjacent tetrahedron indices. */
        std::array<std::string, 4> gluing; /**< Gluings, as 4-digit strings. */
        std::array<int, 4> cusp;           /**< Cusp index per vertex, -1 if finite. */
        std::string shape;                 /**< The shape line, passed through. */
    };

    std::string name;                   /**< The file's name line. */
    std::string solution;               /**< The solution-type + volume line. */
    std::string orientability;          /**< e.g. `oriented_manifold`. */
    std::string chernSimons;            /**< The CS line. */
    int numOrCusps = 0;                 /**< Orientable cusps. */
    int numNonOrCusps = 0;              /**< Non-orientable cusps. */
    std::vector<std::string> cuspLines; /**< One topology/filling line per cusp. */
    std::vector<Tet> tets;              /**< Per-tetrahedron blocks. */
    CurveField meridian{0};             /**< curve[M], as loaded. */
    CurveField longitude{0};            /**< curve[L], as loaded. */
};

/**
 * Parses a SnapPea data file, keeping the peripheral-curve block.
 *
 * \exception regina::InvalidArgument \a text is not a SnapPea data file we
 * recognise.
 */
SnapPeaFile parseSnapPea(const std::string &text);

/** A Curve as a CurveField over `numTet` tetrahedra, on the right-handed sheet. */
CurveField toField(const Curve &curve, size_t numTet);

/**
 * The algebraic intersection number `<A, B>` on cusp \a cusp.
 *
 * Mirrors compute_intersection_numbers() in
 * engine/snappea/kernel/intersection_numbers.cpp, including its convention
 * that edge crossings are counted only where \a A is entering (not leaving)
 * a triangle, so that each edge crossing is counted once rather than once per
 * incident triangle. \a A is double-copied onto both sheets exactly as
 * copy_curves_to_scratch() does for torus cusps.
 */
long intersectionNumber(const SnapPeaFile &file, int cusp, const CurveField &a,
                        const CurveField &b);

/**
 * Expresses \a mu as `(a, b)` against the file's stored `(meridian,
 * longitude)` basis on cusp \a cusp.
 *
 * \exception regina::InvalidArgument the stored basis does not intersect
 * once, i.e. the file does not actually carry a peripheral basis.
 */
std::pair<long, long> slope(const SnapPeaFile &file, int cusp,
                            const CurveField &mu);

/**
 * Completes `(a, b)` to a change-of-basis matrix `[(a,b), (c,d)]` with
 * determinant 1, as SnapPy's set_peripheral_curves() wants.
 *
 * \pre gcd(a, b) == 1, which holds for any meridian (it is primitive).
 */
std::pair<long, long> completeBasis(long a, long b);

} // namespace peripheral

#endif // PERIPHERAL_H
