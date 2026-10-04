//
//  cobordism.h
//
//  One cobordism found by a search.
//

#ifndef SURFER_COBOUND_COBORDISM_H
#define SURFER_COBOUND_COBORDISM_H

#include <string>
#include <vector>

/*! \file utils/surfer/cobound/cobordisms/cobordism.h
 *  \brief One surface a search found, as the record it is stored under, and
 *  the identity a search deduplicates by.
 */

namespace cobordisms {

/** Whether a cobordism bounds its subject on its own, or only relative to
 * another name. */
enum class CobordismKind { direct, cobordism };

/**
 * One surface we actually found.
 */
struct Cobordism {
    CobordismKind kind = CobordismKind::cobordism;

    std::string
        subject; /**< The name whose own boundary this surface realizes. */
    int subjectComponents = 1; /**< Observed curve count on `subject`'s side. */

    std::string other; /**< The outgoing link's name (or description); empty iff `kind ==
                          direct`. */
    std::vector<std::string> otherCandidates;
    /**< The oriented variants `other` could be (see NameTable::candidates).
         Empty iff `kind == direct`. */
    int otherComponents = 1; /**< Observed curve count on the outgoing link. */

    int genus = 0;
    /**< The genus of the connected surface this cobordism provides. For a
         disconnected find, the tubed genus (SurfaceFoundInfo::tubedGenus),
         not the meaningless whole-complex figure. */
    bool tubed = false;
    /**< Whether the surface as found was disconnected and had to be tubed
         to give `genus`. Recorded so a reader can tell that the pairSig
         names a disconnected complex, with the tubing as the step from it
         to the surface the bound is about. */

    std::string pairSig;
    /**< The found surface's pair signature, if captured -- held in memory
         only until the cobordism is on disk (see fileOffset). */
    long long fileOffset = -1;
    /**< Byte offset of this cobordism's line in the database it was loaded
         from or appended to, or -1 while it exists only in memory. A loaded
         cobordism never holds its ~11 KB pair signature in memory; anything
         that needs it reads it back through this offset. */
    std::string pairSigKey;
    /**< cobordisms::cobordismKey(pairSig), sha1(pairsig)[:12] -- the key the
         per-cobordism outgoing resolution file is keyed on. Filled at load
         only when resolutions are in use, and for every new cobordism. */
    /**
     * Whether this cobordism's outgoing link has been PROVED, per cobordism, to be
     * the link named in `other` -- by the peripheral system (an isometry or
     * curve-respecting isomorphism carrying our meridians), or by a split
     * decomposition whose every factor is a knot (Gordon-Luecke). Set only
     * on the solver-side copy by applyOutgoingResolutions(); never recorded.
     *
     * This is what lets outgoingBearsBound() accept a multi-component outgoing
     * link: the complement alone never determines a link, but complement
     * plus meridians does, and then the oriented variants of `other` ARE a
     * complete candidate set, so the max/min over them is sound.
     */
    bool outgoingProved = false;
    /**
     * Whether the outgoing link was named (outgoing_names_file, --far-side-exact): redrawn
     * from this cobordism's own pair signature, oriented by its surface, and
     * named with a proof that is an identity (one link up to mirror and
     * global reversal). Such an outgoing link may also RECEIVE a bound from the
     * subject, whatever its component count. Solver-side only; never
     * recorded.
     */
    bool outgoingNamed = false;

    // Provenance: which search produced this, and under what budget.
    std::string sourceSearch;
    int thickenLayers = 1;
    long long maxFaces = 0; /**< 0 means the search was unbounded. */

    int resolvedVertices = 0;
    /**< How many ambient vertices the found surface meets itself at. 0 for
         an embedded surface. Positive only under --resolve-unlinked, where
         each such vertex is an unlinked self-intersection, and a
         perturbation near those vertices turns the surface into an embedded
         one of the same topology and boundary (paper §4.5), so the cobordism
         bounds exactly as an embedded one would. The solver ignores it, but
         whether it is zero is part of cobordismIdentity(): an embedded
         cobordism must never be discarded as a duplicate of one that rests on
         the resolution theorem. */
};

/** The dedup identity of a cobordism: kind, subject and its component count,
 * outgoing name, outgoing component count, genus, tubed, and whether it is
 * resolved. Only the first cobordism of each identity is recorded.
 *
 * For an outgoing link that bears a bound (outgoingBearsBound()) equal identity
 * means the cobordisms bound exactly the same thing. For any other outgoing link
 * the name comes from its complement, which several links can share, so two
 * surfaces reaching DIFFERENT links can collapse here. That is deliberate:
 * keying such outgoing links by their own edges instead (tried 2026-09-26) made
 * nearly every surface its own cobordism -- 19,690 cobordisms from 23,672
 * surfaces for L9a43{1;0} at cap 3, against 109 by name -- which no database
 * can hold at a million surfaces per search, while the per-cobordism pipeline
 * has found such a collapse in 21 of 228,580 cobordisms, and those outgoing links
 * bear nothing in the solvers. A cheap invariant of the outgoing LINK
 * (per-curve knot types, linking numbers) would separate them boundedly. */
std::string cobordismIdentity(const Cobordism &w);

/** Whether `cobordisms` already contains a cobordism with `w`'s
 * cobordismIdentity(). Linear; the search keeps a hashed set of identities
 * instead, and this remains for tests and small callers. */
bool haveCobordism(const std::vector<Cobordism> &cobordisms, const Cobordism &w);
} // namespace cobordisms

#endif // SURFER_COBOUND_COBORDISM_H
