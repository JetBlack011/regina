// nodes.h
//
// The cascade's links, each a node of the proof graph, identified exactly
// and with the component map that soundness rests on. See README.md,
// "Component maps".

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "exactnaming/gaussdiagram.h"
#include "exactnaming/snappeaisometry.h"
#include "proofgraph.h"

namespace cascade {

/// The linking matrix of a diagram, components in its own order.
std::vector<std::vector<int>> linkingMatrix(const exactnaming::GaussDiagram &d);

/**
 * Simplifies a diagram with regina::Link::simplify() (Reidemeister moves
 * only: never reflects or reverses, and keeps component indices, a zero-
 * crossing component in its slot), keeping `origin`, then removes every
 * nugatory crossing left (removeNugatoryCrossings()): a hop's row must be a
 * reduced diagram for knotbuilder's drawer to certify it. Throws std::logic_error
 * if the component count or any pairwise linking number changed: component
 * identity is what every profile is indexed by, so a violation must stop
 * everything rather than mislabel components.
 */
exactnaming::GaussDiagram simplifyKeepingComponents(const exactnaming::GaussDiagram &d);

/// A nugatory crossing of `d` (one whose removal disconnects its diagram:
/// always a self-crossing), or nullopt if `d` is reduced.
std::optional<size_t> nugatoryCrossing(const exactnaming::GaussDiagram &d);

/// `d` with every nugatory crossing removed: each by turning over the side it
/// cuts off, which deletes it, swaps over and under on that side's crossings
/// and keeps every sign. The same oriented link, components in the same
/// order, with every linking number kept.
exactnaming::GaussDiagram removeNugatoryCrossings(exactnaming::GaussDiagram d);

/// How a diagram was found to be a node.
struct NodeMatch {
  NodeId node = -1;
  /// The diagram's component i is the node's component componentMap[i].
  std::vector<int> componentMap;
  bool mirrored = false; ///< the diagram is the node's mirror image
  bool reversed = false; ///< every component of the diagram runs against the node's
  bool created = false;
  std::string method; ///< "new", "diagram", "isometry" or "unknot"
};

struct NodeInfo {
  exactnaming::GaussDiagram diagram; ///< the representative diagram, simplified
  std::vector<std::vector<int>> linking;
  bool hyperbolic = false;
  double volume = 0;
};

struct NodeRegistryStats {
  long lookups = 0, diagramHits = 0, isometryHits = 0, isometryRetries = 0,
       created = 0, unknots = 0;
};

/**
 * Interns connected (one split piece), simplified diagrams as proof-graph
 * nodes. A match is only ever claimed on an exact test:
 *   - a diagram isomorphism (cascade/diagramiso.h), up to mirror and global
 *     reversal, whose component map is the node's; or
 *   - for hyperbolic diagrams, an isometry of complements carrying meridians
 *     to meridians with ONE orientation sign on every component
 *     (exactnaming::KernelLink, exact when found), whose component map is
 *     the isometry's.
 * A miss creates a new node: two nodes for one link cost duplicated search,
 * never soundness. One node for two links is what must not happen, and
 * cannot, since every match above is a proof.
 */
class NodeRegistry {
public:
  explicit NodeRegistry(ProofGraph &graph);
  NodeRegistry(const NodeRegistry &) = delete;
  NodeRegistry &operator=(const NodeRegistry &) = delete;

  /// \pre `piece` is one split piece (exactnaming::splitPieces()), simplified.
  NodeMatch intern(const exactnaming::GaussDiagram &piece, const std::string &label);

  /// The unknot's node (created on first use, with its disc as a leaf).
  NodeId unknot();

  const NodeInfo &info(NodeId n) const { return info_.at(n); }
  bool known(NodeId n) const { return info_.count(n) > 0; }
  const NodeRegistryStats &stats() const { return stats_; }

private:
  static std::string diagramKey(const exactnaming::GaussDiagram &d);

  ProofGraph &g_;
  std::map<NodeId, NodeInfo> info_;
  std::map<NodeId, std::unique_ptr<exactnaming::KernelLink>> kernel_;
  std::multimap<std::string, NodeId> byDiagramKey_;
  // Hyperbolic nodes by (components, volume rounded to 1e-6); lookups scan
  // neighbouring buckets, so rounding never hides a match.
  std::multimap<std::pair<size_t, long long>, NodeId> byVolume_;
  NodeId unknot_ = -1;
  NodeRegistryStats stats_;
};

} // namespace cascade
