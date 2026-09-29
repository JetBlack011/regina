// proofgraph.h
//
// The cascade's knowledge: links (nodes), the witnesses and split
// decompositions that relate them (edges), and immutable proof records for
// every (partition, genus) bound derived from them. See README.md, "The
// proof graph".
//
// The graph has cycles. Every witness is used in BOTH directions (a
// cobordism read backwards is a cobordism), a far side can be a link met
// before, and a search from a far side can find an edge back to an ancestor.
// Bounds are therefore a fixed point, recomputed incrementally whenever an
// edge or a leaf arrives (propagate()), never a tree evaluation.
//
// Proofs stay well-founded because records are immutable and created only
// on a strict improvement, from records that already exist: every record's
// children have smaller ids, so proof(r) is a DAG, and a node can never be
// justified by itself through a cycle (glue() never returns a genus below its
// input's, so a cycle cannot produce a strict improvement of its own start).

#pragma once

#include <map>
#include <optional>
#include <string>
#include <vector>

#include "profile.h"

namespace cascade {

using NodeId = int;
using EdgeId = int;
using RecordId = long;

enum class RecordKind {
  leaf,            ///< A fact from outside: a table, an anchor, a certificate.
  witnessForward,  ///< A witness's incoming link, from its outgoing link.
  witnessReverse,  ///< A witness's outgoing link, from its incoming link.
  splitCombine,    ///< A split link, from surfaces for its pieces.
  splitRestrict,   ///< A piece of a split link, from a surface for the whole.
};

const char *kindName(RecordKind k);

struct Record {
  RecordId id = -1;
  NodeId node = -1;
  Partition partition;
  int genus = 0;
  RecordKind kind = RecordKind::leaf;
  EdgeId edge = -1;                 ///< The witness or split edge used.
  std::vector<RecordId> children;   ///< Records used; all have smaller ids.
  std::string source;               ///< For leaves: where the fact comes from.
};

/**
 * A witness: a cobordism from node `in`'s link to node `out`'s link.
 * `inMap[i]` is the component of node `in` that the shape's incoming curve i
 * is, and `outMap[j]` likewise; both are bijections. Establishing them
 * (identity with orientation up to global reversal and mirror, README.md) is
 * the caller's job and is what soundness rests on.
 */
struct WitnessEdge {
  EdgeId id = -1;
  NodeId in = -1, out = -1;
  CobordismShape shape;
  std::vector<int> inMap, outMap;
  std::string key; ///< Provenance, e.g. the witness key sha1(pairsig)[:12].
};

/**
 * A split link and its pieces: node `whole` is the split union of the
 * `pieces` (separated by 2-spheres). `pieceMap[k][c]` is the component of
 * `whole` that component c of pieces[k] is; together they are a bijection
 * onto the components of `whole`.
 */
struct SplitEdge {
  EdgeId id = -1;
  NodeId whole = -1;
  std::vector<NodeId> pieces;
  std::vector<std::vector<int>> pieceMap;
};

struct Node {
  NodeId id = -1;
  int components = 0;
  std::string label;
  /// Linking matrix, when known: used only as a contradiction gate.
  std::optional<std::vector<std::vector<int>>> linking;
  /// A proved lower bound on the connected slice genus, when known: every
  /// entry's genus must be at least this (a surface with any partition tubes
  /// to a connected one of the same genus). Used only as a contradiction gate.
  std::optional<int> genusLowerBound;
  std::string lowerBoundSource;
  Profile profile;
  std::vector<EdgeId> witnessEdges; ///< Edges with this node at either end.
  std::vector<EdgeId> splitEdges;   ///< As whole or as a piece.
};

class ProofGraph {
public:
  NodeId addNode(int components, std::string label,
                 std::optional<std::vector<std::vector<int>>> linking = {});
  void setGenusLowerBound(NodeId n, int lo, std::string source);

  /// Records an outside fact; returns its record, or -1 if already implied.
  /// Call propagate() afterwards.
  RecordId addLeaf(NodeId n, const Partition &p, int genus,
                   std::string source);
  /// Adds a witness edge (validated); call propagate() afterwards.
  EdgeId addWitness(NodeId in, NodeId out, CobordismShape shape,
                    std::vector<int> inMap, std::vector<int> outMap,
                    std::string key);
  /// Adds a split edge (validated); call propagate() afterwards.
  EdgeId addSplit(NodeId whole, std::vector<NodeId> pieces,
                  std::vector<std::vector<int>> pieceMap);

  /// Derives every bound that follows, to a fixed point. Returns the number
  /// of records created.
  long propagate();
  /// Re-derives from every record, not just new ones, and propagates. At a
  /// fixed point this creates nothing; tests use it to check that.
  long saturate();

  const Node &node(NodeId n) const { return nodes_.at(n); }
  const Record &record(RecordId r) const { return records_.at(r); }
  const WitnessEdge &witness(EdgeId e) const { return witnesses_.at(e); }
  const SplitEdge &split(EdgeId e) const { return splits_.at(e); }
  size_t nodeCount() const { return nodes_.size(); }
  size_t recordCount() const { return records_.size(); }

  /// The least genus known for node n with a partition refining `target`.
  std::optional<ProfileEntry> best(NodeId n, const Partition &target) const;
  /// The connected slice genus bound for node n.
  std::optional<ProfileEntry> bestConnected(NodeId n) const;

  /// Every record `r` rests on, children before parents, `r` last.
  std::vector<RecordId> proof(RecordId r) const;

  /// Re-derives record r from its children with glue() and the edge's maps,
  /// returning an empty string if it reproduces exactly, else what differs.
  /// A bookkeeping check, not an independent one (cascade_check.py is).
  std::string recheck(RecordId r) const;

  /// Contradictions met so far (a derived bound below a proved lower bound,
  /// or a partition the linking numbers forbid). Each means a bug or a wrong
  /// identification somewhere; the driver halts on any.
  const std::vector<std::string> &contradictions() const {
    return contradictions_;
  }

  // ------------------------------------------------------------ lower bounds
  //
  // lower(n, Q): a lower bound on the total genus of EVERY surface for node
  // n whose partition refines Q (README.md, "Lower bounds"). Seeded by the
  // literature bound (setGenusLowerBound(), for every Q) and the linking
  // condition (a forbidden Q, and so each of its refinements, gets
  // kNoSurface); transported across witnesses in both directions and across
  // split edges, to a fixed point. Tracked only for nodes of at most
  // kMaxLowerComponents components (Bell(5) = 52 partitions).

  static constexpr int kNoSurface = 1 << 29;
  static constexpr int kMaxLowerComponents = 5;

  /// Relaxes every lower bound to a fixed point, then checks each against
  /// the upper bounds (a contradiction if some proved surface is below).
  /// Returns the number of improvements.
  long propagateLower();
  /// The lower bound for surfaces refining q (0 when nothing is known).
  int lower(NodeId n, const Partition &q) const;
  /// Forgets every lower bound: the literature seeds and everything
  /// propagated from them (the linking condition, read from the nodes, stays).
  /// For what-ifs that seed one fact and read where it reaches.
  void clearLowerBounds();

private:
  RecordId insert(NodeId n, const Partition &p, int genus, RecordKind kind,
                  EdgeId edge, std::vector<RecordId> children,
                  std::string source);
  void deriveFrom(RecordId r);

  // What glue() and the maps give for a witness edge from one entry; used by
  // both deriveFrom() and recheck().
  std::optional<std::pair<Partition, int>>
  throughWitness(const WitnessEdge &e, bool forward, const Partition &p,
                 int genus) const;
  std::optional<std::pair<Partition, int>>
  combineSplit(const SplitEdge &s,
               const std::vector<const Record *> &perPiece) const;
  std::optional<Partition> restrictSplit(const SplitEdge &s, size_t piece,
                                         const Partition &whole) const;

  // Lower bounds per node: partition labels -> bound (absent: 0).
  std::vector<std::map<std::vector<int>, int>> lower_;
  bool raiseLower(NodeId n, const Partition &q, int value);
  // The bound transported to `to` (partition q) across witness e from the
  // other end; minimum over refinements of q.
  int lowerAcross(const WitnessEdge &e, bool toIsIn, const Partition &q) const;

  std::vector<Node> nodes_;
  std::vector<Record> records_;
  std::vector<WitnessEdge> witnesses_;
  std::vector<SplitEdge> splits_;
  std::vector<RecordId> pending_;
  std::vector<std::string> contradictions_;
};

} // namespace cascade
