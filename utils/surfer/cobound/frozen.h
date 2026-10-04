//
//  frozen.h
//
//  Output tokens frozen in retired terms.
//

#ifndef SURFER_COBOUND_FROZEN_H
#define SURFER_COBOUND_FROZEN_H

/*! \file utils/surfer/cobound/frozen.h
 *  \brief The output tokens that still spell a retired term (witness, hop,
 *  node, row, far side, search side, cascade, store, exact name, profile,
 *  identification), named for what they are. The code says search, link,
 *  cobordism, outgoing, incoming, database, ...; the bytes stay the old ones
 *  (plan, Hard constraint 2), until the atlas task's format change renames
 *  them here.
 *
 *  Here are the tokens something reads back or keys on: file and directory
 *  names, prefixes parsed back, enumerated values, CSV headers, and the log
 *  tokens the atlas's tools parse. The rest of the frozen formats -- JSON keys
 *  and the free text of log lines -- are literals at the one place each is
 *  written, and the phase-6 token map lists every one with its future name.
 */

// ---- the run directory (driver/run.cpp, driver/scheduler.cpp, driver/runrecords.cpp)

/// A search's directory is <work>/hop_<k>_n<link>/: k the search's number in
/// the run, <link> its link's id in the cobordism graph (0 in a run without a
/// goal). Found again by its prefix and numbered by the digits after it
/// (cobordisms/pending.cpp, driver/run.cpp).
inline constexpr char kFrozenHopDirPrefix[] = "hop_";
inline constexpr char kFrozenHopDirNodeMark[] = "_n";
/// The goal run's records: one JSON line per search, its links' partition
/// genera, their lower bounds, and the untabulated links it met.
inline constexpr char kFrozenCascadeJsonl[] = "cascade.jsonl";
inline constexpr char kFrozenProfilesJsonl[] = "profiles.jsonl";
inline constexpr char kFrozenNodeBoundsJsonl[] = "node_bounds.jsonl";
inline constexpr char kFrozenNodesCsv[] = "nodes.csv";
/// nodes.csv's `label` value: "node <id>".
inline constexpr char kFrozenNodesCsvLabel[] = "node ";
/// An untabulated link's subject name: cascade:<run>/<target>/n<id>, in
/// cobordisms.csv (`subject`, `source_row`) and nodes.csv; nodes.csv lists
/// exactly the subjects with this prefix.
inline constexpr char kFrozenCascadeSubjectPrefix[] = "cascade:";
inline constexpr char kFrozenCascadeSubjectNodeMark[] = "/n";

// ---- the database (cobordisms/database.cpp, cobordisms/pending.cpp)

/// The sidecar beside a database: which diagram each cobordism was searched
/// on (`witness` its cobordism key, `row_pd` the incoming diagram's PD).
inline constexpr char kFrozenRowsSidecarSuffix[] = ".rows.csv";
inline constexpr char kFrozenRowsSidecarHeader[] = "witness,layers,row_pd\n";

// ---- the cobordism graph's values (labels, sources, keys and kinds)

/// A derivation leaf from a cobordism with no outgoing link: its source is
/// "direct witness <key>", and the certificate writer reads <key> back from
/// after the prefix.
inline constexpr char kFrozenDirectWitnessSource[] = "direct witness ";
/// A goal run's in-process cobordism key: hop<k>#<i>, the i-th kept surface
/// of search k (certificate.json `records[].witness`).
inline constexpr char kFrozenHopKeyPrefix[] = "hop";
/// Link labels: an outgoing link's, a split outgoing link's whole, a stored
/// cobordism's incoming link ("row <name>"), a summand ("summand <k> of node
/// <id>").
inline constexpr char kFrozenFarSideLabel[] = "far side of ";
inline constexpr char kFrozenSplitFarSideLabel[] = "split far side of ";
inline constexpr char kFrozenRowLabel[] = "row ";
inline constexpr char kFrozenSummandOfNodeLabel[] = " of node ";
/// Derivation kinds of a cobordism used forwards and in reverse
/// (certificate.json `records[].kind`).
inline constexpr char kFrozenKindWitnessForward[] = "witness-forward";
inline constexpr char kFrozenKindWitnessReverse[] = "witness-reverse";
/// A lower bound transported across a cobordism: node_bounds.jsonl's and
/// lower_certificate.json's `kind`.
inline constexpr char kFrozenLowerKindWitness[] = "witness";

// ---- a search's records (search/preconditions.cpp, search/searchreport.cpp, driver/run.cpp)

/// Gate reasons (rejection_sample_log `reason`): the incoming curves not the
/// incoming link's component count; more than one outgoing link.
inline constexpr char kFrozenReasonSearchSideBroken[] = "search-side-broken";
inline constexpr char kFrozenReasonMultiFarSide[] = "multi-far-side";
/// The accounting line keeps the bucket of the retired unseeded search, whose
/// count is always 0 (dispatch.py's RE_ACCOUNTING, the canaries' expected lines).
inline constexpr char kFrozenSearchSideElsewhere[] = ", search-side-elsewhere 0";
inline constexpr char kFrozenSurfaceStatsHeader[] =
    "row,max_faces,triangles,orientable,genus,punctures,tubed_genus,"
    "closed_components,connected,count\n";
inline constexpr char kFrozenSelfIntersectionCensusHeader[] =
    "row,max_faces,resolve_unlinked,satisfying,embedded,resolved,"
    "singular,interior_unlinked,interior_uncertified,"
    "boundary_unlinked,boundary_uncertified,multi_open,"
    "multi_open_search_side,multi_open_far,multi_open_far_clean,"
    "multi_open_far_simple,far_configs,far_clean_configs,"
    "configs_saturated,audited,audit_knotted,knotted_pairsigs\n";
inline constexpr char kFrozenRejectionSampleHeader[] =
    "row,reason,tubed_genus,connected,boundary,pairsig\n";

// ---- log tokens the atlas's tools parse

/// `[+] <name>: N new witnesses, outcome X` (dispatch.py's RE_OUTCOME; N the
/// search's or run's new cobordisms).
inline constexpr char kFrozenNewWitnessesOutcome[] = " new witnesses, outcome ";
/// `[+] witness store: K kept, F new, A appended to P` (signing into the database).
inline constexpr char kFrozenWitnessStoreLine[] = "[+] witness store: ";
/// `[+] exact far-side names: N table entries (T ms)`: the tables the outgoing
/// links are named against, loaded; every run without a goal prints it
/// (dispatch.py fails a search without it). A run that cannot load them exits
/// 1 before searching.
inline constexpr char kFrozenExactFarSideNamesLine[] = "[+] exact far-side names: ";
/// A goal run's per-search lines: `[+] hop <k> <name>: accounting: ...`,
/// `: diagram naming:`, `: breadth:`, `[!!] hop <k>: surface accounting
/// failed -- ...`, `[!] hop <k>: frontier not written: ...`, `[+] hop <k>:
/// node <id> (...)`.
inline constexpr char kFrozenHopLine[] = "hop ";
/// `[+] <name>: identification: ...` (the complement route's counters).
inline constexpr char kFrozenIdentificationLine[] = ": identification: ";
/// `[+] <name>: CONSTRUCTIVE witness found -- reaches genus G ...`.
inline constexpr char kFrozenConstructiveWitnessFound[] =
    ": CONSTRUCTIVE witness found -- reaches genus ";
/// `cobound draw`'s incoming-link line (the atlas's checker, cascade_check.py,
/// reads it).
inline constexpr char kFrozenDrawRowLine[] = "ROW components=";

#endif // SURFER_COBOUND_FROZEN_H
