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
 *  partition: its records, children before parents, each witness record
 *  with what a checker needs to replay its surface (the row, the faces and
 *  the thickening's digest, or the pair signature), and every node it
 *  names. A lower certificate is the proof of the target's lower bound as a
 *  tree of facts.
 */

namespace cascade {

/// What a checker needs to replay one cobordism edge of the graph.
struct EdgeInfo {
  std::string hopDir, rowPD, key;
  HopEdge he;
  int layers = 2;
  std::string pairsig; ///< inline for master witnesses (no hop directory)
  /// The row diagram's component i is the node's rowNodeMap[i]: identity
  /// for a hop on the node's own diagram, the registry's map for a master
  /// row (the table's diagram).
  std::vector<int> rowNodeMap;
  /// An in-process hop's surface, as triangles of the row's thickening,
  /// and that thickening's digest (WitnessRedrawer::buildChecksum()). A
  /// certificate carries both; the checker rebuilds the row, refuses a
  /// different digest, and rebuilds the surface from its faces. No pair
  /// signature is ever computed for it (that is only for a witness bound
  /// for the atlas; pairSigsOf()).
  std::vector<int> faces;
  std::string build;
};

/// Every edge's EdgeInfo: witness edges by edge, direct witnesses by key.
struct EdgeInfos {
  std::map<EdgeId, EdgeInfo> byEdge;
  std::map<std::string, EdgeInfo> direct;
};

/// The goal a certificate proves.
struct CertificateGoal {
  std::string targetName, targetPD;
  int goalGenus = 0;
  int goalLower = -1;
  bool disjoint = false;
  NodeId target = -1;
  Partition partition; ///< the goal partition of the target
};

/// Writes a goal's certificates from the graph, its links, their table
/// names and how each edge's surface is found. All references must outlive
/// the writer.
class CertificateWriter {
public:
  CertificateWriter(const ProofGraph &g, const NodeRegistry &reg,
                    const std::map<NodeId, std::string> &tableName, const EdgeInfos &edges)
      : g_(g), reg_(reg), tableName_(tableName), edges_(edges) {}

  /// certificate.json at `path`: the proof of the target's best genus on
  /// the goal partition. Nothing when there is none.
  void writeUpper(const std::string &path, const CertificateGoal &goal) const;
  /// lower_certificate.json at `path`: the proof of lower(target, goal)
  /// as a tree of facts, children before parents.
  void writeLower(const std::string &path, const CertificateGoal &goal) const;
  /// The proof of lower(n, q), readably: its reason, then the reasons of
  /// what it read, indented, down to literature and linking leaves.
  void describeLower(std::ostream &o, NodeId n, const Partition &q, int indent) const;

private:
  static void writeSurface(std::ostream &c, const EdgeInfo &info);
  void writeWitnessEdge(std::ostream &c, EdgeId e, std::set<NodeId> &nodes) const;
  void writeRecords(std::ostream &c, const std::vector<RecordId> &ids,
                    std::set<NodeId> &nodes) const;
  void writeNodes(std::ostream &c, const std::set<NodeId> &nodes) const;

  const ProofGraph &g_;
  const NodeRegistry &reg_;
  const std::map<NodeId, std::string> &tableName_;
  const EdgeInfos &edges_;
};

} // namespace cascade

#endif // SURFER_COBOUND_CERTIFICATE_H
