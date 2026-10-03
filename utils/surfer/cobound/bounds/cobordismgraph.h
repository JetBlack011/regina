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

#include "cobound/bounds/partitiongenera.h"

namespace bounds {

using NodeId = int;
using EdgeId = int;
using RecordId = long;

enum class RecordKind {
  leaf,            ///< A fact from outside: a table, an anchor, a certificate.
  witnessForward,  ///< A witness's incoming link, from its outgoing link.
  witnessReverse,  ///< A witness's outgoing link, from its incoming link.
  splitCombine,    ///< A split link, from surfaces for its pieces.
  splitRestrict,   ///< A piece of a split link, from a surface for the whole.
  sumCombine,      ///< A sum along components, from surfaces for its summands.
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

/**
 * A link that is a connected sum of its `pieces` along components (paper
 * lem:sum-along-components; README.md "Sums along components"): node
 * `whole` is obtained from the pieces by successive sums, and
 * `pieceMap[k][c]` is the component of `whole` that component c of
 * pieces[k] becomes part of. Every component of `whole` is the image of at
 * least one piece component; two piece components with one image were
 * summed together. Unlike a split, the map is not injective.
 */
struct SumEdge {
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
  std::vector<EdgeId> sumEdges;     ///< As whole or as a summand.
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
  /// Adds a sum-along-components edge (validated: every whole component
  /// is hit, at least two pieces); call propagate() afterwards.
  EdgeId addSum(NodeId whole, std::vector<NodeId> pieces,
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
  const SumEdge &sum(EdgeId e) const { return sums_.at(e); }
  size_t nodeCount() const { return nodes_.size(); }
  size_t recordCount() const { return records_.size(); }
  size_t witnessCount() const { return witnesses_.size(); }

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

  /// Why a lower bound holds: the last fact that raised it. Every raise is a
  /// strict increase and no rule increases what it reads (transport
  /// subtracts, a piece's bound subtracts proved genera, a whole's is a sum
  /// of its pieces'), so following reasons from any bound reaches literature
  /// or linking leaves and never cycles: a lower bound's proof is a tree.
  struct LowerReason {
    enum class Kind { none, literature, linking, witness, splitWhole, splitPiece, seed, sumPiece };
    Kind kind = Kind::none;
    /// witness: the edge and which end this is; the other end's partition
    /// (labels) whose bound `from` was read, and what the cap added.
    EdgeId edge = -1;
    bool toIsIn = false;
    std::vector<int> fromPartition;
    int from = 0, addition = 0;
    /// splitWhole: each piece's partition (labels) whose bounds were summed;
    /// splitPiece: the whole's partition (fromPartition) and the proved
    /// surfaces (records) of the other pieces whose genera were subtracted;
    /// sumPiece: the sum edge (edge), the summand (from = its connected
    /// bound; `piece` its index) and the records of the other summands'
    /// connected surfaces, whose genera plus components minus one were
    /// subtracted (addition).
    std::vector<std::vector<int>> pieces;
    std::vector<RecordId> records;
    int piece = -1;
  };
  /// The bound lower(n, q) reads, the partition it is stored for (q or a
  /// coarser one) and its reason. Kind none when nothing is known (value 0).
  struct LowerFact {
    int value = 0;
    Partition storedFor;
    LowerReason reason;
  };
  LowerFact lowerWhy(NodeId n, const Partition &q) const;

  /// What a transport read: the other end's partition and what the cap added.
  struct Transport {
    Partition otherPartition;
    int addition = 0;
    int other = 0; ///< lower() at the other end
  };
  /// What witness e transports to its `to` end (toIsIn: node e.in, else
  /// e.out) for surfaces with partition EXACTLY p: cap such a surface onto
  /// e, read the other end's partition and bound, and subtract what the cap
  /// added (glue() with a genus-0 cap). kNoSurface if p is forbidden or the
  /// other end has no surface of the capped partition; 0 if the other end
  /// has no curves. Monotone under refinement of p (README.md, "Lower
  /// bounds": lem:transport-monotone), which is why lowerAcross() need not
  /// minimise over refinements; proofgraph_test checks that on random graphs.
  int transportedLower(const WitnessEdge &e, bool toIsIn, const Partition &p,
                       Transport *detail = nullptr) const;
  /// Changes whenever a lower bound is raised or cleared: a cheap key for
  /// caching what-ifs across calls.
  long lowerVersion() const { return lowerVersion_; }

  /// A what-if: `seeds` raise the given nodes' bounds (each for surfaces
  /// refining its partition) in a copy that KEEPS every bound this graph
  /// already has, relax, and read lower(target, goal). nullopt if the copy
  /// meets a contradiction (a seed above a proved surface), in which case
  /// nothing it says is used. This graph is untouched.
  struct LowerSeed {
    NodeId node = -1;
    Partition partition;
    int value = 0;
  };
  std::optional<int> lowerIf(const std::vector<LowerSeed> &seeds, NodeId target,
                             const Partition &goal) const;
  /// Forgets every lower bound: the literature seeds and everything
  /// propagated from them (the linking condition, read from the nodes, stays).
  /// For what-ifs that seed one fact and read where it reaches.
  void clearLowerBounds();

  /// Node n's profile as JSON object fields, without the enclosing braces,
  /// for the driver's profiles.jsonl (README.md, "Profiles"):
  ///   "components", "linking" (when known), "genus_lower" (when seeded),
  ///   "entries": its Pareto set, [{"p": partition, "g": genus, "r": record}]
  ///   sorted by partition then genus, and, for at most kMaxLowerComponents
  ///   components, "lower": [{"p": Q, "lo": lower(n, Q)}] for every
  ///   partition Q with a positive bound, or {"p": Q, "forbidden": true}
  ///   where no surface can have partition Q (kNoSurface).
  std::string profileFields(NodeId n) const;

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
  // The whole's (partition, genus) from one surface per summand: genera add
  // and blocks merge where components were summed (paper lem:sum-partitions).
  std::pair<Partition, int> combineSum(const SumEdge &s,
                                       const std::vector<const Record *> &perPiece) const;

  // Lower bounds per node: partition labels -> bound (absent: 0), and why.
  std::vector<std::map<std::vector<int>, int>> lower_;
  std::vector<std::map<std::vector<int>, LowerReason>> lowerReason_;
  long lowerVersion_ = 0;
  bool raiseLower(NodeId n, const Partition &q, int value,
                  const LowerReason &why);
  // The bound transported to `to` for surfaces refining q across witness e:
  // transportedLower() at q itself, since that is monotone under refinement.
  int lowerAcross(const WitnessEdge &e, bool toIsIn, const Partition &q) const;

  std::vector<Node> nodes_;
  std::vector<Record> records_;
  std::vector<WitnessEdge> witnesses_;
  std::vector<SplitEdge> splits_;
  std::vector<SumEdge> sums_;
  std::vector<RecordId> pending_;
  std::vector<std::string> contradictions_;
};

} // namespace bounds
