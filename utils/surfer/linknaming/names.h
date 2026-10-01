//
//  names.h
//
//  The grammar of far-side names.
//

#ifndef SURFER_LINKNAMING_NAMES_H
#define SURFER_LINKNAMING_NAMES_H

#include <optional>
#include <string>
#include <utility>
#include <vector>

/*! \file utils/surfer/linknaming/names.h
 *  \brief The grammar of the names a far side is recorded under, read back:
 *  component counts, bases, splits (`A u B`), alternatives (`A|B`),
 *  composites (`K #_c L`, `A#B`), sums (`#{...}`) and knot marks.
 *
 *  Pure string functions, shared by both solvers' C++ side and the search.
 */

namespace cobordismgraph {

/**
 * The number of components of whatever `name` names, worked out from the
 * name itself.
 *
 * \note This is a fallback for names we have no better information about.
 * When a component count is *observed* (the number of boundary curves a
 * found surface actually put on an ambient boundary component), prefer the
 * observation: it is a fact about this surface, whereas this is an
 * inference from a string.
 */
int componentsFromName(const std::string &name);

/** `name` with any trailing orientation tag removed: `L6a3{0}` -> `L6a3`. */
std::string baseName(const std::string &name);

/**
 * A COMPOSITE far-side name, "K #_c L": the knot K connect-summed into
 * component c of the link L (L a base name, c in L's PD component order).
 *
 * The component is part of the name because K # L is only well defined once
 * it is given. It does not enter the bound: g_4(K #_c L) <= g_4(K) + g_4(L)
 * for every c, by boundary-connect-summing the two minimal surfaces.
 * \a knot keeps its chirality ("m3_1"); g_4 is mirror-invariant, so the
 * solver looks it up without the "m".
 */
struct CompositeName {
    std::string knot;
    int component = 0;  /**< -1 for "?": the namer did not compute it */
    std::string link;   /**< a base, or a tagged variant ("L7n1{1}") from an exact name */
};

/** The parts of a composite name, or nullopt if `name` is not one. */
std::optional<CompositeName> compositeParts(const std::string &name);

/**
 * The summands of a composite KNOT name, "A#B#..." (each "K", "mK" or
 * "Unknot"), or an empty vector if `name` is not one. A census name
 * ("m129 : #2"), a knot-into-link sum ("3_1 #_0 L2a1") and a split name
 * ("A u B") are all refused.
 */
std::vector<std::string> knotSummands(const std::string &name);

/**
 * A (possibly marked) prime knot name without its marks: "mr8_17" -> "8_17".
 * An exact name marks a summand "m" (mirrored) and/or "r" (reversed)
 * relative to the table's diagram, where its symmetry type makes the mark
 * matter; g_4 sees neither. Any other name is returned unchanged.
 */
std::string stripKnotMarks(const std::string &name);

/**
 * The pieces of a sum along components, each with its component count:
 * "K1 #_? K2 #_? L" (knots summed into components of the link L), or
 * "#{A[?] # B[?] ; ...}" (sum sites, one per shared component). nullopt for
 * anything else, or when a piece's name does not state its component count.
 *
 * In "#{...}" a piece summed at two sites is written at both, so every
 * occurrence counts as a piece of its own. That can only over-count, which
 * weakens both bounds built from the pieces and never makes one unsound.
 */
std::optional<std::vector<std::pair<std::string, int>>> sumPieces(const std::string &name);

/**
 * What a far side named EXACTLY (--far-side-exact) may be: the name itself,
 * or -- where the namer proved it one of a few orientation variants it could
 * not tell apart -- those alternatives "A|B". Never widened to a base's
 * variants. A split or a sum is one candidate: its alternatives live inside
 * its factors and pieces.
 */
std::vector<std::string> exactCandidates(const std::string &name);

/**
 * An identify() result reduced to the name the graph should key on.
 *
 * Strips the trailing parenthetical if there is one; leaves everything else
 * (`"Unknot"`, `"3-component unlink"`, a bare isoSig, an undecorated census
 * name) untouched.
 */
std::string normalizeIdentifiedName(const std::string &name);

/**
 * The factors of a split name, "A u B u ...", or an empty vector if `name`
 * is not one.
 *
 * A far side is split exactly when its exterior is reducible, which is the
 * commonest thing an unnameable multi-component far side turns out to be:
 * measured over 70 sampled unidentified link far sides, 43 of them. The
 * factors are named separately, each by its own exterior, and the link is
 * their disjoint union.
 */
std::vector<std::string> splitFactors(const std::string &name);

/** The alternatives of a split factor, "A|B|...": `{factor}` when there is
 *  only one. */
std::vector<std::string> factorAlternatives(const std::string &factor);

/** A tagged table link, "L7n1{1}": its component count is in the tag. */
bool isTaggedLinkName(const std::string &s);

} // namespace cobordismgraph

#endif // SURFER_LINKNAMING_NAMES_H
