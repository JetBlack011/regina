//
//  incoming.h
//
//  The incoming link of a search, in the search's own terms.
//

#ifndef SURFER_COBOUND_INCOMING_H
#define SURFER_COBOUND_INCOMING_H

#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "diagramtriangulation/thickening/thickening.h"
#include "diagramtriangulation/todiagram.h"
#include "linknaming/diagrams/gaussdiagram.h"

/*! \file utils/surfer/cobound/search/incoming.h
 *  \brief The incoming link of a search -- the link L searched from, on the
 *  incoming boundary T x {0} -- in the search's own terms: its edges and its
 *  PD orientation there (the incoming map), and the certification that T
 *  carries the diagram the search was given (certifyIncoming()).
 */

namespace search {

/**
 * The incoming diagram's PD-tagged edges (diagramtriangulation::TriangulationWithLink's
 * `edges`/`reversed`), translated into directed pairs of *vertex indices*
 * within some triangulation combinatorially isomorphic to the diagram's own
 * diagram triangulation -- normally the ambient incoming boundary
 * component, rebuilt via `BoundaryComponent<4>::build()` (see
 * buildIncomingOrientation()).
 */
struct IncomingOrientation {
    std::unordered_map<size_t, size_t> tailOf;
    /**< Edge index of L in the incoming triangulation -> index of that
         edge's tail vertex under the diagram's PD orientation. Keyed by edge,
         not by vertex pair, so two edges joining the same pair of vertices
         can never be confused. */
    std::vector<size_t> edges; /**< Sorted keys of tailOf: L's edge set. */
    std::unordered_map<size_t, size_t> incomingIndexOf;
    /**< Edge index of L in the incoming triangulation -> that edge's
         position in the diagramEdges given to buildIncomingOrientation(), i.e. which
         edge of the diagram's own link it is. Lets a caller tell which
         component of L an incoming curve is. */
    size_t components = 0;
    /**< How many closed curves `edges` forms, each checked to chain head to
         tail under the PD orientation (buildIncomingOrientation() throws
         otherwise). */
    bool divergedFromDefaultIsomorphism = false;
    /**< Whether the isomorphism isIsomorphicTo() would have returned maps L
         differently -- onto other edges, or with other directions. That is
         the map this code used to trust; see buildIncomingOrientation(). */
};

/**
 * Builds `diagramEdges`/`diagramReversed`'s IncomingOrientation against `incomingTri`.
 *
 * The diagram's own triangulation and `incomingTri` are isomorphic, but PD
 * triangulations have automorphisms, so "an isomorphism" does not determine
 * where L goes. When `requiredEdges` is given (a seeded search: the seed's
 * own edges in that boundary component, i.e. L x {0} exactly), the
 * isomorphism used is one taking L's edges onto exactly that set. Any such
 * choice is sound: two of them differ by an automorphism of the diagram's
 * triangulation preserving L, a PL homeomorphism of S^3 taking L to itself,
 * and the slice genus is invariant under homeomorphism and mirroring.
 *
 * \throws regina::InvalidArgument if `diagramEdges` is empty, if the two
 * triangulations are not isomorphic, if no isomorphism takes L onto
 * `requiredEdges`, or if the image of L fails to chain into closed directed
 * curves.
 */
IncomingOrientation
buildIncomingOrientation(const std::vector<const regina::Edge<3> *> &diagramEdges,
                    const std::vector<bool> &diagramReversed,
                    const regina::Triangulation<3> &incomingTri,
                    const std::vector<size_t> *requiredEdges = nullptr);
} // namespace search

namespace search {

/**
 * The edges of `faces` (triangle indices of `tri`) lying in boundary
 * component `bcIndex`, as sorted local edge indices of that component --
 * the numbering its built triangulation uses.
 */
std::vector<size_t> boundaryEdgesOf(const regina::Triangulation<4> &tri,
                                    const std::vector<int> &faces,
                                    size_t bcIndex);
} // namespace search


namespace search {

/**
 * The incoming link's ambient S^3 x I and its seed (ThickenedLink), with the incoming map.
 * Filled in place by buildIncoming() and never moved: an OutgoingNamer and a
 * SurfaceSearch built from it hold pointers into `link.tri` and `cob`.
 */
struct IncomingThickening : ThickenedLink {
    std::optional<search::IncomingOrientation> orientation;
    /**< The incoming map: L's edges and PD orientation on the incoming boundary. */
    std::vector<size_t> incomingEdges;
    /**< The incoming link, as sorted edge indices of the incoming boundary
         component's built triangulation. Seeded, the seed's own edges there
         (L x {0}); unseeded, the image of L under the incoming map. */
};

/**
 * Builds `thickened` for PD code `pdNotation`: buildAmbient(), then orientIncoming().
 *
 * \throws regina::InvalidArgument as either does.
 */
void buildIncoming(const std::string &pdNotation, int thickenLayers,
              int collarLayers, IncomingThickening &thickened);

/**
 * The incoming map for an ambient built by buildAmbient(). Checks, once, that the
 * incoming side holds exactly L's edges in L's number of components.
 *
 * \throws regina::InvalidArgument for an incoming map that cannot be built or
 * fails those checks.
 */
void orientIncoming(IncomingThickening &thickened);

/** The linknaming::GaussDiagram of a drawn diagram (todiagram.h): its
 *  crossings' signs and Gauss codes, `origin` the drawn component index. */
linknaming::GaussDiagram gaussOf(const diagramtriangulation::Diagram &d);

/** Thrown by certifyIncoming(): T does not carry the diagram the search was
 *  given. */
struct IncomingNotCertified : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/**
 * The certification that a search's triangulation carries the diagram it was
 * given (plan, phase 7.2): the incoming link as T holds it -- `cycles`, the
 * incoming link's edge cycles in T's component order -- drawn back
 * from T by `drawer` (diagramtriangulation::DiagramDrawer on T) is
 * isomorphic to `given` as a diagram, orientation kept and no mirror
 * (linknaming::findDiagramIsomorphism()). Returns that isomorphism's
 * component map: incoming component i is `given`'s component map[i].
 *
 * Every search certifies its incoming link against the PD it was given
 * (Searcher::run()), at depth 0 as with a goal; a goal run's assembler maps
 * the incoming curves onto its graph link's diagram through it too
 * (bounds::CobordismAssembler). Never a fallback: a table row that does not
 * certify is refused, since searching any other diagram would silently
 * change T. (Only an untabulated goal target may fall back to its
 * simplified diagram, and says so: scheduler.h.)
 *
 * \throws IncomingNotCertified when there is no such isomorphism.
 */
std::vector<int> certifyIncoming(const diagramtriangulation::DiagramDrawer &drawer,
                                 const std::vector<diagramtriangulation::EdgeCycle> &cycles,
                                 const linknaming::GaussDiagram &given);

/** `pd`'s own diagram, as certifyIncoming() checks a search against it:
 *  linknaming::linkFromTablePD(), components in the PD's order. */
linknaming::GaussDiagram diagramOfPD(const std::string &pd);

} // namespace search

#endif // SURFER_COBOUND_INCOMING_H
