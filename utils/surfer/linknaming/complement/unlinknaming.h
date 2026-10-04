//
//  unlinknaming.h
//
//  The unknot and unlinks, named from their complements.
//

#ifndef SURFER_LINKNAMING_UNLINKNAMING_H
#define SURFER_LINKNAMING_UNLINKNAMING_H

#include <optional>
#include <string>
#include <vector>

#include <triangulation/dim3.h>

#include "linknaming/complement/linkcomplement.h"

/*! \file utils/surfer/linknaming/complement/unlinknaming.h
 *  \brief Names the unknot and unlinks from their complements, with
 *  proof and without any census: a handlebody-genus check (cached in
 *  complementcache.h), a free-group certificate, and the capping that turns
 *  an arc in a ball into a closed curve. This is what KnottedSurface's local
 *  flatness and resolution checks rest on; naming anything else is
 *  censusnaming.h's.
 */

namespace complement {

/**
 * The handlebody genus of `complement` (whose isoSig is `sig`), cached in
 * the complement cache: 1 for a solid torus (the unknot's complement), -1
 * for "not any handlebody". Tries two fast, sound one-sided checks first
 * (a fundamental group that simplifies to Z proves 1; a hyperbolic
 * structure proves -1) before recogniseHandlebody(). Never touches the
 * census, so it is safe on isUnknot()'s hot, highly parallel path.
 */
ssize_t cachedGenus(const regina::Triangulation<3> &complement,
                    const std::string &sig);

/**
 * The rank of `t`'s fundamental group when the presentation group() reaches
 * has no relations -- a free group, of rank countGenerators(), however the
 * simplification got there -- else nullopt, which is only ever inconclusive.
 * The one group test behind groupProvesUnlink(), the unknot's fast path in
 * cachedGenus() (rank 1) and certifiesUnlink() (rank m).
 */
std::optional<size_t> freeGroupRank(const regina::Triangulation<3> &t);

/**
 * Fast, sound, one-sided proof that `t` (a link complement, possibly of
 * several components) is the complement of a split unlink: its fundamental
 * group simplifies to a presentation with no relations, i.e. it is free
 * (freeGroupRank()). A \c false is only ever inconclusive.
 */
bool groupProvesUnlink(const regina::Triangulation<3> &t);

/**
 * The name of the `components`-component unlink: `"Unknot"` for one
 * component, else `"<n>-component unlink"`. The one spelling of both; every
 * namer writes it through here.
 */
std::string unlinkName(size_t components);

/** Whether `name` is `"<n>-component unlink"` (an unlink of two or more
 *  components; not `"Unknot"`). */
bool isMultiComponentUnlinkName(const std::string &name);

/**
 * Whether `name` is an unlink's name, `"Unknot"` or `"<n>-component unlink"`
 * -- and so (as returned by census::nameComplement(const Link&)) safe to use for a genus
 * deduction on a MULTI-component outgoing link; false for everything else.
 *
 * Those two are the only multi-curve names nameComplement() produces from a
 * structural proof of the link itself (free pi_1 => split unlink), rather
 * than from a lookup of the complement's homeomorphism type. Every other
 * multi-component name -- a Thistlethwaite name, a census name, a bare
 * isoSig -- records the complement, and a link complement belongs to
 * infinitely many non-isotopic links (Rolfsen twisting along an unknotted
 * component), with different slice genera. So such a name does not tell you
 * which link was found, and NOT merely which orientation of it: the
 * candidate set is not enumerable at all. Single-curve names need no such
 * check, since Gordon-Luecke makes a knot's complement determine it.
 *
 * The solver applies this through solver::outgoingBearsBound(). Kept
 * here, next to where these strings are actually produced, rather than
 * pattern-matched elsewhere, so the two stay in sync if the format ever
 * changes.
 */
bool isUnlinkName(const std::string &name);

/**
 * Cheap test for whether `e`'s complement is a genus-1 handlebody (a solid
 * torus) -- unlike census::nameComplement(), this never falls back to the slower
 * Census::lookup().
 *
 * \return \c true if and only if the complement is a genus-1 handlebody.
 */
bool isUnknot(const EdgeComplement &e);

/**
 * `e`'s name without any census: `"Unknot"` if its complement is a solid
 * torus, else the complement's isoSig -- what census::nameComplement() names it when the
 * census has no name for it.
 */
std::string unlinkNameOrIsoSig(const EdgeComplement &e);

/**
 * `l`'s name without any census: `"<n>-component unlink"` for n > 1
 * components whose complement's group is free (groupProvesUnlink()),
 * `"Unknot"` for a solid-torus complement, else the complement's isoSig --
 * what census::nameComplement(const Link&) names it when the census has no name for it.
 */
std::string unlinkNameOrIsoSig(const Link &l);

/**
 * Sound, one-sided certificate that `edges` is the `m`-component unlink in
 * `tri`: returns \c true only if the edges form exactly `m` pairwise
 * disjoint cycles and the fundamental group of their complement simplifies
 * to a presentation with `m` generators and no relations. A link in S^3 is
 * the unlink if and only if its group is free (pl_enumeration_draft,
 * Lemma "unlink-free"), so \c true is conclusive. \c false is conclusive
 * only about malformed input; for well-formed input it means "not proven"
 * (the presentation did not simplify all the way), never "proven linked".
 *
 * Unlike groupProvesUnlink(), this checks the rank against `m` explicitly
 * and validates the edge set itself, rather than relying on Link's
 * component walk and on the ambient being a homology sphere.
 *
 * \pre `tri` is a triangulated 3-sphere. The certificate is about links in
 * S^3; a free group of the right rank proves nothing elsewhere.
 */
bool certifiesUnlink(const regina::Triangulation<3> &tri,
                     const std::vector<const regina::Edge<3> *> &edges,
                     size_t m);

/**
 * A triangulated 3-ball with its boundary sphere coned off to a single
 * apex (so a 3-sphere), together with a set of curves carried across into
 * it; see capInCone().
 *
 * \warning `edges` points into `tri`, so an instance must stay where
 * capInCone() filled it: copying or moving it would leave `edges` pointing
 * into the old triangulation.
 */
struct CappedCurves {
    regina::Triangulation<3> tri; /**< The coned-off ball: a 3-sphere. */
    std::vector<const regina::Edge<3> *> edges;
    /**< The input edges, re-resolved in `tri`, plus -- if the input had an
         open arc -- the two cone edges closing it up through the apex. */
    size_t components = 0; /**< Number of closed curves `edges` forms. */

    CappedCurves() = default;
    CappedCurves(const CappedCurves &) = delete;
    CappedCurves &operator=(const CappedCurves &) = delete;
};

/**
 * Cones off the boundary of `ball` (a triangulated 3-ball, e.g. the link of
 * a boundary vertex of a 4-manifold) into `out.tri` and carries `edges`
 * across into `out.edges`. `edges` must be a disjoint union of cycles plus
 * at most one arc whose two endpoints lie on the boundary sphere; the arc
 * is closed up through the apex. This is Definition "petal-knotted"'s
 * capping construction: the two cone edges through the apex form an arc
 * isotopic, rel endpoints, to any arc in the boundary sphere pushed slightly
 * outward, so the closed-up curve is unknotted exactly when the arc is.
 *
 * \return \c false (leaving `out` unspecified) if `ball` has no boundary, or
 * `edges` is not of the stated form: a repeated edge, a vertex of degree
 * greater than 2, more than one arc, or an arc endpoint off the boundary.
 */
bool capInCone(const regina::Triangulation<3> &ball,
               const std::vector<const regina::Edge<3> *> &edges,
               CappedCurves &out);

} // namespace complement

#endif // SURFER_LINKNAMING_UNLINKNAMING_H
