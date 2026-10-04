//
//  certificate.h
//
//  Upper and lower certificates: a goal's proof, as a checker replays it.
//

#ifndef SURFER_COBOUND_CERTIFICATE_H
#define SURFER_COBOUND_CERTIFICATE_H

#include <map>
#include <ostream>
#include <set>
#include <string>
#include <vector>

#include "cobound/bounds/cobordismgraph.h"
#include "cobound/bounds/links.h"
#include "cobound/bounds/partition.h"
#include "cobound/bounds/searchcobordisms.h"

/*! \file utils/surfer/cobound/bounds/certificate.h
 *  \brief certificate.json and lower_certificate.json (formats frozen:
 *  `records[].witness` as `hop<k>#<i>`, `records[].source`, `nodes[]`, every
 *  key), read by the atlas's cascade_check.py and cascade_record.py.
 *
 *  An upper certificate is the proof of the target's best genus on the goal
 *  partition: its derivations (`records[]`), children before parents, each cobordism
 *  derivation with what a checker needs to replay its surface (the searched
 *  diagram, the faces and
 *  the thickening's digest, or the pair signature), and every link it
 *  names. A lower certificate is the proof of the target's lower bound as a
 *  tree of facts.
 */

namespace bounds {

/// What a checker needs to replay one cobordism of the graph.
struct CobordismSource {
  std::string searchDir, incomingPD, key;
  AddedCobordism added;
  int layers = 2;
  std::string pairsig; ///< inline for master cobordisms (no search directory)
  /// The searched diagram's component i is the link's incomingLinkMap[i]:
  /// identity for a search on the link's own diagram, the registry's map for a
  /// master cobordism (searched on the table's diagram).
  std::vector<int> incomingLinkMap;
  /// An in-process search's surface, as triangles of its thickening,
  /// and that thickening's digest (OutgoingReader::buildChecksum()). A
  /// certificate carries both; the checker rebuilds the thickening, refuses a
  /// different digest, and rebuilds the surface from its faces. No pair
  /// signature is ever computed for it (that is only for a cobordism bound
  /// for the atlas; pairSigsOf()).
  std::vector<int> faces;
  std::string build;
};

/// Every cobordism's CobordismSource: graph cobordisms by relation, direct ones by key.
struct CobordismSources {
  std::map<RelationId, CobordismSource> byCobordism;
  std::map<std::string, CobordismSource> direct;
};

/// The goal a certificate proves.
struct CertificateGoal {
  std::string targetName, targetPD;
  int goalGenus = 0;
  int goalLower = -1;
  bool disjoint = false;
  LinkId target = -1;
  Partition partition; ///< the goal partition of the target
};

/// Writes a goal's certificates from the graph, its links, their table
/// names and how each cobordism's surface is found. All references must outlive
/// the writer.
class CertificateWriter {
public:
  CertificateWriter(const CobordismGraph &g, const LinkRegistry &reg,
                    const std::map<LinkId, std::string> &tableName, const CobordismSources &sources)
      : g_(g), reg_(reg), tableName_(tableName), sources_(sources) {}

  /// certificate.json at `path`: the proof of the target's best genus on
  /// the goal partition. Nothing when there is none. Throws
  /// std::runtime_error when the file cannot be opened or written (the run
  /// then halts and claims no goal: scheduler.cpp).
  void writeUpper(const std::string &path, const CertificateGoal &goal) const;
  /// lower_certificate.json at `path`: the proof of lower(target, goal)
  /// as a tree of facts, children before parents. Throws as writeUpper().
  void writeLower(const std::string &path, const CertificateGoal &goal) const;
  /// The proof of lower(n, q), readably: its reason, then the reasons of
  /// what it read, indented, down to literature and linking leaves.
  void describeLower(std::ostream &o, LinkId n, const Partition &q, int indent) const;

private:
  static void writeSurface(std::ostream &c, const CobordismSource &info);
  void writeCobordism(std::ostream &c, RelationId e, std::set<LinkId> &links) const;
  void writeDerivations(std::ostream &c, const std::vector<DerivationId> &ids,
                    std::set<LinkId> &links) const;
  void writeLinks(std::ostream &c, const std::set<LinkId> &links) const;

  const CobordismGraph &g_;
  const LinkRegistry &reg_;
  const std::map<LinkId, std::string> &tableName_;
  const CobordismSources &sources_;
};

} // namespace bounds

#endif // SURFER_COBOUND_CERTIFICATE_H
