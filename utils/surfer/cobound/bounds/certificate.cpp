//
//  certificate.cpp
//

#include "cobound/bounds/certificate.h"

#include <algorithm>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <utility>

#include "cobound/json.h"
#include "cobound/frozen.h"

namespace bounds {

// How a certificate finds a witness's surface: a master witness's pair
// signature inline; an in-process one's faces in its row's thickening, with
// that thickening's digest.
void CertificateWriter::writeSurface(std::ostream &c, const EdgeInfo &info) {
  if (!info.pairsig.empty()) c << ",\"pairsig\":\"" << json::escape(info.pairsig) << "\"";
  if (!info.faces.empty()) {
    c << ",\"faces\":[";
    for (size_t i = 0; i < info.faces.size(); ++i) c << (i ? "," : "") << info.faces[i];
    c << "],\"build\":\"" << info.build << "\"";
  }
}

void CertificateWriter::writeWitnessEdge(std::ostream &c, EdgeId eid,
                                         std::set<NodeId> &nodes) const {
  // A witness edge as a checker replays it: its key, ends, shape and maps,
  // and (for an edge with a hop or master row) the row, the surface (faces
  // and build digest, or pair signature) and each far-side piece's match.
  const WitnessEdge &we = g_.witness(eid);
  nodes.insert(we.in);
  nodes.insert(we.out);
  const auto it = edges_.byEdge.find(eid);
  c << ",\"witness\":\"" << json::escape(we.key) << "\"";
  c << ",\"in\":" << we.in << ",\"out\":" << we.out
    << ",\"shape\":{\"components\":" << we.shape.components << ",\"genus\":" << we.shape.genus
    << ",\"inComponent\":" << json::array(we.shape.inComponent)
    << ",\"outComponent\":" << json::array(we.shape.outComponent) << "},\"inMap\":" << json::array(we.inMap)
    << ",\"outMap\":" << json::array(we.outMap);
  if (it == edges_.byEdge.end()) return;
  const HopEdge &he = it->second.he;
  c << ",\"hop_dir\":\"" << json::escape(it->second.hopDir) << "\",\"row_pd\":\""
    << json::escape(it->second.rowPD) << "\",\"layers\":" << it->second.layers
    << ",\"row_node_map\":" << json::array(it->second.rowNodeMap);
  writeSurface(c, it->second);
  c << ",\"split_edge\":" << he.splitEdge << ",\"farCurveEdges\":[";
  for (size_t j = 0; j < he.farCurveEdges.size(); ++j)
    c << (j ? "," : "") << json::array(he.farCurveEdges[j]);
  c << "],\"pieces\":[";
  for (size_t k = 0; k < he.pieces.size(); ++k) {
    const NodeMatch &m = he.pieces[k];
    c << (k ? "," : "") << "{\"node\":" << m.node << ",\"method\":\"" << m.method
      << "\",\"componentMap\":" << json::array(m.componentMap) << ",\"mirrored\":"
      << (m.mirrored ? "true" : "false") << ",\"reversed\":" << (m.reversed ? "true" : "false")
      << ",\"origins\":" << json::array(he.pieceOrigins[k]) << "}";
    nodes.insert(m.node);
  }
  c << "]";
}

void CertificateWriter::writeRecords(std::ostream &c, const std::vector<RecordId> &ids,
                                     std::set<NodeId> &nodes) const {
  bool first = true;
  for (RecordId r : ids) {
    const Record &rec = g_.record(r);
    nodes.insert(rec.node);
    c << (first ? "" : ",\n") << "{\"id\":" << r << ",\"node\":" << rec.node
      << ",\"partition\":\"" << rec.partition.str() << "\",\"genus\":" << rec.genus
      << ",\"kind\":\"" << kindName(rec.kind) << "\",\"edge\":" << rec.edge
      << ",\"source\":\"" << json::escape(rec.source) << "\",\"children\":[";
    for (size_t i = 0; i < rec.children.size(); ++i)
      c << (i ? "," : "") << rec.children[i];
    c << "]";
    if (rec.kind == RecordKind::witnessForward || rec.kind == RecordKind::witnessReverse)
      writeWitnessEdge(c, rec.edge, nodes);
    if (rec.kind == RecordKind::splitCombine || rec.kind == RecordKind::splitRestrict) {
      const SplitEdge &se = g_.split(rec.edge);
      nodes.insert(se.whole);
      nodes.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == RecordKind::sumCombine) {
      const SumEdge &se = g_.sum(rec.edge);
      nodes.insert(se.whole);
      nodes.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == RecordKind::leaf && rec.source.rfind(kFrozenDirectWitnessSource, 0) == 0) {
      const std::string key = rec.source.substr(sizeof kFrozenDirectWitnessSource - 1);
      if (auto it = edges_.direct.find(key); it != edges_.direct.end()) {
        c << ",\"witness\":\"" << json::escape(it->second.key)
          << "\",\"hop_dir\":\"" << json::escape(it->second.hopDir) << "\",\"row_pd\":\""
          << json::escape(it->second.rowPD) << "\",\"layers\":" << it->second.layers;
        writeSurface(c, it->second);
      }
    }
    c << "}";
    first = false;
  }
}

void CertificateWriter::writeLower(const std::string &path, const CertificateGoal &goal) const {
  // The proof of lower(target, goal) as a tree of facts, children before
  // parents (README.md, "Lower-bound mode"): each fact is a node, a
  // partition, the value the bound holds there, and its reason. A witness
  // fact carries the edge exactly as an upper record does, so the checker
  // replays the surface the same way, then recomputes the cap's addition
  // and the other end's partition itself. Split facts carry the upper
  // records they subtract, with those records' own proofs.
  using Kind = ProofGraph::LowerReason::Kind;
  struct Fact {
    NodeId node;
    Partition q;
    ProofGraph::LowerFact f;
    long from = -1;          ///< witness / split-piece: the fact read
    std::vector<long> pieces; ///< split-whole: the pieces' facts
  };
  std::vector<Fact> facts;
  std::map<std::pair<NodeId, std::vector<int>>, long> ids;
  std::set<std::pair<NodeId, std::vector<int>>> onStack;
  std::set<RecordId> records;
  std::function<long(NodeId, const Partition &)> visit = [&](NodeId n, const Partition &q) -> long {
    const ProofGraph::LowerFact f = g_.lowerWhy(n, q);
    const auto key = std::make_pair(n, f.storedFor.labels());
    if (auto it = ids.find(key); it != ids.end()) return it->second;
    if (!onStack.insert(key).second)
      throw std::logic_error("a lower bound's reasons cycle: node " + std::to_string(n));
    Fact fact{n, f.storedFor, f};
    switch (f.reason.kind) {
    case Kind::witness: {
      const WitnessEdge &e = g_.witness(f.reason.edge);
      fact.from = visit(f.reason.toIsIn ? e.out : e.in,
                        Partition::fromLabels(f.reason.fromPartition));
      break;
    }
    case Kind::splitWhole:
      for (EdgeId sid : g_.node(n).splitEdges) {
        const SplitEdge &s = g_.split(sid);
        if (s.whole != n) continue;
        for (size_t k = 0; k < s.pieces.size() && k < f.reason.pieces.size(); ++k)
          fact.pieces.push_back(visit(s.pieces[k], Partition::fromLabels(f.reason.pieces[k])));
        break;
      }
      break;
    case Kind::splitPiece:
      for (EdgeId sid : g_.node(n).splitEdges) {
        const SplitEdge &s = g_.split(sid);
        if (s.whole == n || std::find(s.pieces.begin(), s.pieces.end(), n) == s.pieces.end())
          continue;
        fact.from = visit(s.whole, Partition::fromLabels(f.reason.fromPartition));
        break;
      }
      for (RecordId r : f.reason.records)
        for (RecordId p : g_.proof(r)) records.insert(p);
      break;
    case Kind::sumPiece: {
      const SumEdge &s = g_.sum(f.reason.edge);
      const NodeId piece = s.pieces[static_cast<size_t>(f.reason.piece)];
      fact.from = visit(piece, Partition::coarsest(g_.node(piece).components));
      for (RecordId r : f.reason.records)
        for (RecordId p : g_.proof(r)) records.insert(p);
      break;
    }
    default:
      break;
    }
    onStack.erase(key);
    facts.push_back(std::move(fact));
    const long id = static_cast<long>(facts.size()) - 1;
    ids[key] = id;
    return id;
  };
  const long top = visit(goal.target, goal.partition);
  std::ofstream c(path);
  c << "{\"target\":\"" << json::escape(goal.targetName) << "\",\"target_pd\":\""
    << json::escape(goal.targetPD) << "\",\"goal_lower\":" << goal.goalLower << ",\"goal\":\""
    << (goal.disjoint ? "disjoint" : "connected") << "\",\"lower\":"
    << g_.lower(goal.target, goal.partition) << ",\"top\":" << top << ",\"facts\":[\n";
  std::set<NodeId> nodes;
  auto value = [](int v) {
    return v >= ProofGraph::kNoSurface ? std::string("\"inf\"") : std::to_string(v);
  };
  for (size_t i = 0; i < facts.size(); ++i) {
    const Fact &fact = facts[i];
    nodes.insert(fact.node);
    c << (i ? ",\n" : "") << "{\"id\":" << i << ",\"node\":" << fact.node << ",\"partition\":\""
      << fact.q.str() << "\",\"value\":" << value(fact.f.value);
    switch (fact.f.reason.kind) {
    case Kind::literature:
      c << ",\"kind\":\"literature\",\"source\":\"" << json::escape(g_.node(fact.node).lowerBoundSource)
        << "\"";
      break;
    case Kind::linking:
      c << ",\"kind\":\"linking\"";
      break;
    case Kind::witness: {
      const WitnessEdge &e = g_.witness(fact.f.reason.edge);
      c << ",\"kind\":\"" << kFrozenLowerKindWitness
        << "\",\"to_is_in\":" << (fact.f.reason.toIsIn ? "true" : "false")
        << ",\"from\":" << fact.from << ",\"from_partition\":\""
        << Partition::fromLabels(fact.f.reason.fromPartition).str() << "\",\"from_value\":"
        << value(fact.f.reason.from) << ",\"addition\":" << fact.f.reason.addition;
      writeWitnessEdge(c, e.id, nodes);
      break;
    }
    case Kind::splitWhole: {
      c << ",\"kind\":\"split-whole\",\"pieces\":[";
      for (size_t k = 0; k < fact.pieces.size(); ++k) c << (k ? "," : "") << fact.pieces[k];
      c << "]";
      for (EdgeId sid : g_.node(fact.node).splitEdges) {
        const SplitEdge &se = g_.split(sid);
        if (se.whole != fact.node) continue;
        c << ",\"whole\":" << se.whole << ",\"piece_nodes\":" << json::array(se.pieces)
          << ",\"pieceMap\":[";
        for (size_t k = 0; k < se.pieceMap.size(); ++k) c << (k ? "," : "") << json::array(se.pieceMap[k]);
        c << "]";
        break;
      }
      break;
    }
    case Kind::splitPiece: {
      c << ",\"kind\":\"split-piece\",\"from\":" << fact.from << ",\"from_partition\":\""
        << Partition::fromLabels(fact.f.reason.fromPartition).str() << "\",\"from_value\":"
        << value(fact.f.reason.from) << ",\"subtracted\":" << fact.f.reason.addition
        << ",\"records\":[";
      for (size_t k = 0; k < fact.f.reason.records.size(); ++k)
        c << (k ? "," : "") << fact.f.reason.records[k];
      c << "]";
      break;
    }
    case Kind::sumPiece: {
      const SumEdge &se = g_.sum(fact.f.reason.edge);
      c << ",\"kind\":\"sum-piece\",\"from\":" << fact.from << ",\"piece\":" << fact.f.reason.piece
        << ",\"from_value\":" << value(fact.f.reason.from) << ",\"subtracted\":"
        << fact.f.reason.addition << ",\"records\":[";
      for (size_t k = 0; k < fact.f.reason.records.size(); ++k)
        c << (k ? "," : "") << fact.f.reason.records[k];
      c << "],\"whole\":" << se.whole << ",\"piece_nodes\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k) c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
      for (NodeId p : se.pieces) nodes.insert(p);
      break;
    }
    default:
      c << ",\"kind\":\"none\"";
      break;
    }
    c << "}";
  }
  c << "\n],\"records\":[\n";
  writeRecords(c, std::vector<RecordId>(records.begin(), records.end()), nodes);
  c << "\n],\"nodes\":[\n";
  writeNodes(c, nodes);
  c << "\n]}\n";
}

void CertificateWriter::writeUpper(const std::string &path, const CertificateGoal &goal) const {
  auto best = g_.best(goal.target, goal.partition);
  if (!best) return;
  std::ofstream c(path);
  c << "{\"target\":\"" << json::escape(goal.targetName) << "\",\"target_pd\":\""
    << json::escape(goal.targetPD) << "\",\"goal_genus\":" << goal.goalGenus
    << ",\"goal\":\"" << (goal.disjoint ? "disjoint" : "connected") << "\",\"genus\":"
    << best->genus << ",\"records\":[\n";
  std::set<NodeId> nodes;
  writeRecords(c, g_.proof(best->record), nodes);
  c << "\n],\"nodes\":[\n";
  writeNodes(c, nodes);
  c << "\n]}\n";
}

void CertificateWriter::writeNodes(std::ostream &c, const std::set<NodeId> &nodes) const {
  bool first = true;
  for (NodeId n : nodes) {
    c << (first ? "" : ",\n") << "{\"id\":" << n << ",\"label\":\""
      << json::escape(g_.node(n).label) << "\",\"components\":" << g_.node(n).components;
    if (reg_.known(n) && reg_.info(n).diagram.crossings() > 0) {
      // The node's own diagram, as signed Gauss data: component maps refer
      // to ITS component order, which a PD round trip need not keep.
      const linknaming::GaussDiagram &d = reg_.info(n).diagram;
      c << ",\"pd\":\"" << json::escape(rowPD(d)) << "\",\"signs\":[";
      for (size_t k = 0; k < d.signs.size(); ++k) c << (k ? "," : "") << d.signs[k];
      c << "],\"gauss\":[";
      for (size_t i = 0; i < d.comps.size(); ++i) {
        c << (i ? ",[" : "[");
        for (size_t j = 0; j < d.comps[i].size(); ++j) c << (j ? "," : "") << d.comps[i][j];
        c << "]";
      }
      c << "]";
    }
    if (auto it = tableName_.find(n); it != tableName_.end())
      c << ",\"table\":\"" << json::escape(it->second) << "\"";
    c << "}";
    first = false;
  }
}

void CertificateWriter::describeLower(std::ostream &o, NodeId n, const Partition &q,
                                      int indent) const {
  using Kind = ProofGraph::LowerReason::Kind;
  const auto fact = g_.lowerWhy(n, q);
  const std::string pad(static_cast<size_t>(2 * indent + 4), ' ');
  auto name = [&](NodeId m) {
    auto t = tableName_.find(m);
    return "node " + std::to_string(m) +
           (t == tableName_.end() ? std::string(" (untabulated)") : " (" + t->second + ")");
  };
  o << pad << name(n) << " " << q.str() << " >= "
    << (fact.value >= ProofGraph::kNoSurface ? std::string("(no such surface)")
                                             : std::to_string(fact.value));
  if (!(fact.storedFor == q)) o << ", stored for " << fact.storedFor.str();
  if (indent > 40) {
    o << " ...\n";
    return;
  }
  switch (fact.reason.kind) {
  case Kind::none:
    o << ": nothing known\n";
    return;
  case Kind::literature:
    o << ": " << g_.node(n).lowerBoundSource << "\n";
    return;
  case Kind::linking:
    o << ": the linking numbers forbid this partition\n";
    return;
  case Kind::seed:
    o << ": what-if seed\n";
    return;
  case Kind::sumPiece: {
    const SumEdge &s = g_.sum(fact.reason.edge);
    const NodeId piece = s.pieces[static_cast<size_t>(fact.reason.piece)];
    o << ": a sum along components: summand " << name(piece) << " >= " << fact.reason.from
      << ", minus the other summands' connected genera plus components minus one ("
      << fact.reason.addition << "; records";
    for (RecordId r : fact.reason.records) o << " " << r;
    o << ")\n";
    describeLower(o, piece, Partition::coarsest(g_.node(piece).components), indent + 1);
    return;
  }
  case Kind::witness: {
    const WitnessEdge &e = g_.witness(fact.reason.edge);
    const NodeId other = fact.reason.toIsIn ? e.out : e.in;
    const Partition op = Partition::fromLabels(fact.reason.fromPartition);
    o << ": across witness " << e.key << " (genus " << e.shape.genus << ", "
      << e.shape.components << " pieces) from " << name(other) << " " << op.str() << " >= "
      << fact.reason.from << ", the cap adding " << fact.reason.addition << "\n";
    describeLower(o, other, op, indent + 1);
    return;
  }
  case Kind::splitWhole: {
    o << ": the split link's pieces, summed\n";
    for (EdgeId sid : g_.node(n).splitEdges) {
      const SplitEdge &s = g_.split(sid);
      if (s.whole != n) continue;
      for (size_t k = 0; k < s.pieces.size() && k < fact.reason.pieces.size(); ++k)
        describeLower(o, s.pieces[k], Partition::fromLabels(fact.reason.pieces[k]),
                      indent + 1);
      return;
    }
    return;
  }
  case Kind::splitPiece: {
    o << ": a piece of a split link whose whole is bounded, minus the proved genera ("
      << fact.reason.addition << ") of the other pieces (records";
    for (RecordId r : fact.reason.records) o << " " << r;
    o << ")\n";
    for (EdgeId sid : g_.node(n).splitEdges) {
      const SplitEdge &s = g_.split(sid);
      if (s.whole == n) continue;
      if (std::find(s.pieces.begin(), s.pieces.end(), n) == s.pieces.end()) continue;
      describeLower(o, s.whole, Partition::fromLabels(fact.reason.fromPartition), indent + 1);
      return;
    }
    return;
  }
  }
}

} // namespace bounds
