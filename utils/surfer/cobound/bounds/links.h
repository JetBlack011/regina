// links.h
//
// A goal run's links, each a link of the cobordism graph, matched exactly
// and with the component map that soundness rests on. See README.md,
// "Component maps".

#pragma once

#include <map>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "linknaming/diagrams/gaussdiagram.h"
#include "linknaming/diagrams/simplification.h"
#include "linknaming/isometry/isometry.h"
#include "cobound/bounds/cobordismgraph.h"

namespace bounds {

/// How a diagram was found to be a link.
struct LinkMatch {
  LinkId link = -1;
  /// The diagram's component i is the link's component componentMap[i].
  std::vector<int> componentMap;
  bool mirrored = false; ///< the diagram is the link's mirror image
  bool reversed = false; ///< every component of the diagram runs against the link's
  bool created = false;
  std::string method; ///< "new", "diagram", "isometry" or "unknot"
};

struct LinkInfo {
  linknaming::GaussDiagram diagram; ///< the representative diagram, simplified
  std::vector<std::vector<int>> linking;
  bool hyperbolic = false;
  double volume = 0;
};

struct LinkRegistryStats {
  long lookups = 0, diagramHits = 0, isometryHits = 0, isometryRetries = 0,
       created = 0, unknots = 0;
};

/**
 * Interns connected (one split piece), simplified diagrams as cobordism-graph
 * links. A match is only ever claimed on an exact test:
 *   - a diagram isomorphism (linknaming/diagrams/diagramiso.h), up to mirror and global
 *     reversal, whose component map is the link's; or
 *   - for hyperbolic diagrams, an isometry of complements carrying meridians
 *     to meridians with ONE orientation sign on every component
 *     (linknaming::KernelLink, exact when found), whose component map is
 *     the isometry's.
 * A miss creates a new graph link: two graph links for one link cost duplicated
 * search, never soundness. One graph link for two links is what must not happen, and
 * cannot, since every match above is a proof.
 */
class LinkRegistry {
public:
  explicit LinkRegistry(CobordismGraph &graph);
  LinkRegistry(const LinkRegistry &) = delete;
  LinkRegistry &operator=(const LinkRegistry &) = delete;

  /// \pre `piece` is one split piece (linknaming::splitPieces()), simplified.
  LinkMatch intern(const linknaming::GaussDiagram &piece, const std::string &label);

  /// The unknot's link (created on first use, with its disc as a leaf).
  LinkId unknot();

  const LinkInfo &info(LinkId n) const { return info_.at(n); }
  bool known(LinkId n) const { return info_.count(n) > 0; }
  const LinkRegistryStats &stats() const { return stats_; }

private:
  static std::string diagramKey(const linknaming::GaussDiagram &d);

  CobordismGraph &g_;
  std::map<LinkId, LinkInfo> info_;
  std::map<LinkId, std::unique_ptr<linknaming::KernelLink>> kernel_;
  std::multimap<std::string, LinkId> byDiagramKey_;
  // Hyperbolic links by (components, volume rounded to 1e-6); lookups scan
  // neighbouring buckets, so rounding never hides a match.
  std::multimap<std::pair<size_t, long long>, LinkId> byVolume_;
  LinkId unknot_ = -1;
  LinkRegistryStats stats_;
};

} // namespace bounds
