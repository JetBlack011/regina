//
//  solverinputs.h
//
//  What the atlas solver reads besides the database and the tables.
//

#ifndef SURFER_COBOUND_SOLVERINPUTS_H
#define SURFER_COBOUND_SOLVERINPUTS_H

#include <cstddef>
#include <filesystem>
#include <functional>
#include <string>
#include <unordered_map>
#include <vector>

#include "cobound/cobordisms/cobordism.h"
#include "cobound/solver/literature.h"
#include "cobound/solver/solver.h"

/*! \file utils/surfer/cobound/solver/solverinputs.h
 *  \brief The solver's inputs beyond the database and the tables (`solve`):
 *  name aliases, per-cobordism resolutions and outgoing names, the
 *  table's link classes and certified bounds. Each is applied to a SEPARATE
 *  copy of the cobordisms the solver reads (the database keeps what the
 *  search observed). The resolutions and --sum-rules stay until the atlas
 *  task cuts them (plan hand-off item 6).
 */

namespace solverinputs {

/**
 * Loads the observed-name -> classical-name table (name_aliases).
 *
 * An outgoing link is named by census::nameComplement(), from its complement alone, and that
 * often lands on something no literature table knows: a bare isomorphism
 * signature, a Christy census name ("L108019"), or a SnapPy census manifold
 * name ("m129 : #3"). Such an outgoing link bounds nothing. Where we have since
 * PROVED what one of those is -- by a Pachner match against a complement
 * built from a PD code, or by the peripheral test for a link -- this table
 * records it.
 *
 * Deliberately a separate file rather than a rewrite of cobordisms.csv:
 * that file records what the search observed, and must keep doing so. See
 * applyNameAliases() for why the distinction has to survive into memory too.
 */
std::unordered_map<std::string, std::string> loadNameAliases(const std::filesystem::path &path);

/** `name` as loadNameAliases() keys it: any " : #N" census-hit suffix
 *  stripped, matching linknames::name()'s own convention. */
std::string aliasKey(const std::string &name);

/**
 * Resolves every cobordism's outgoing link through the alias table, returning a
 * SEPARATE vector for the solver to consume.
 *
 * Returning a copy rather than mutating in place is the whole point: the
 * vector the search holds is what the database records, so aliasing it
 * would quietly bake resolved names into the observation record -- exactly
 * what keeping a separate alias table was meant to avoid.
 *
 * `otherCandidates` is re-derived rather than carried across: propagate()
 * consumes the stored candidate list, so leaving it keyed to the old name
 * would let a cobordism claim an outgoing link of one name and the variants of
 * another.
 */
std::vector<cobordisms::Cobordism>
applyNameAliases(const std::vector<cobordisms::Cobordism> &cobordisms,
                 const std::unordered_map<std::string, std::string> &aliases,
                 const solver::NameTable &names, size_t &appliedOut);

/** One proved outgoing identity, keyed on the cobordism rather than the name. */
struct OutgoingResolution {
  std::string boundaryComponent; // "0" or "1", as peripheral_slopes reports it
  std::string name;              // the ORIENTED name we have proved it to be
};

/** One row of far_side_exact: the outgoing link redrawn from the cobordism's own
 *  pair signature, oriented by its surface and named with a proof (name). */
struct OutgoingName {
  std::string name;
  bool isName = false; // an identity (may receive a bound), not a description
  int components = 0; // curves drawn: must equal the cobordism's observed count
};

/** Loads far_side_exact: witness,name,exact,pinned,components,proof. A key
 *  seen with two different names is a contradiction and is dropped. */
std::unordered_map<std::string, OutgoingName>
loadOutgoingNames(const std::filesystem::path &path, size_t &clashes);

/**
 * Applies far_side_exact to the solver's copy of the cobordisms, last, so it
 * outranks aliases and resolutions: it is the outgoing link drawn from this
 * cobordism's own surface. Refused (and counted) when the drawing's curve count
 * is not the count the search observed. The candidates are the name alone,
 * or its proved alternatives -- never a base's variants.
 */
std::vector<cobordisms::Cobordism>
applyOutgoingNames(std::vector<cobordisms::Cobordism> cobordisms,
                  const std::unordered_map<std::string, OutgoingName> &outgoingNames, size_t &applied,
                  size_t &refused);

/**
 * Loads the per-cobordism outgoing resolution table (far_side_resolutions).
 *
 * WHY THIS EXISTS SEPARATELY FROM name_aliases. An alias is keyed on the
 * observed NAME, which is sound only where a name determines the object.
 * For a knot it does: Gordon-Luecke makes the complement determine the knot
 * up to mirroring, and g_4 is mirror-invariant. For a LINK it does not --
 * one complement belongs to infinitely many links (Rolfsen twisting), and
 * in our own data one observed census name is a dozen different links
 * across different cobordisms. A name-keyed row for such an outgoing link would be
 * wrong on most of the cobordisms it matched.
 *
 * The pair signature does determine the outgoing link, so outgoing links
 * of several components are keyed on it (via cobordisms::cobordismKey) plus
 * which boundary component of that cobordism is meant.
 */
std::unordered_map<std::string, std::vector<OutgoingResolution>>
loadOutgoingResolutions(const std::filesystem::path &path);

/**
 * Resolves outgoing links cobordism by cobordism, returning a SEPARATE vector for the
 * same reason applyNameAliases() does.
 *
 * Applied AFTER applyNameAliases(), and strictly more specific than it: a
 * resolution names one cobordism's outgoing link, where an alias can only speak
 * about a name. For a KNOT outgoing link the two must agree -- a knot is
 * determined by its complement (Gordon-Luecke), so an alias is an identity
 * and a disagreement is a bug, and the run stops (std::runtime_error). For a
 * LINK outgoing link an alias can only say which COMPLEMENT was observed, and one
 * complement is many links: on 2026-09-24, 21 cobordisms whose outgoing link was
 * aliased from a census name (e.g. 9^2_55 -> L9n6) were proved per cobordism,
 * with their own meridians, to be another link with the same complement
 * (L9n8). There the resolution wins, and the count of such overrides is
 * reported. This mirrors cobordism-atlas/tools/frontier.py's
 * load_witnesses(); the two implementations are deliberately independent,
 * and `frontier.py --check` is only a check while they stay that way.
 */
std::vector<cobordisms::Cobordism> applyOutgoingResolutions(
    std::vector<cobordisms::Cobordism> resolved,
    const std::vector<cobordisms::Cobordism> &observed,
    const std::unordered_map<std::string, std::vector<OutgoingResolution>> &resolutions,
    const solver::NameTable &names, size_t &appliedOut);

/**
 * link_classes: table names that are one oriented link up to mirror and
 * global reversal (tableclasses), name -> its class's canonical name.
 * \throws std::runtime_error "could not open link classes <path>".
 */
std::unordered_map<std::string, std::string> loadLinkClasses(const std::filesystem::path &path);

/**
 * cascade_proofs: target,goal,bound,basis,support,verdict,source,... Only
 * CERTIFIED proofs of the connected goal bound g4; names (the target and
 * every literature leaf) read through `classOf`, as cobordisms'.
 * \throws std::runtime_error the file cannot be opened, or lacks a column.
 */
struct CertifiedBounds {
  std::vector<solver::ExternalProof> proofs;
  size_t skipped = 0; ///< rows not CERTIFIED proofs of the connected goal
};
CertifiedBounds loadCertifiedBounds(const std::filesystem::path &path,
                                  const std::function<const std::string &(const std::string &)> &classOf);

} // namespace solverinputs

#endif // SURFER_COBOUND_SOLVERINPUTS_H
