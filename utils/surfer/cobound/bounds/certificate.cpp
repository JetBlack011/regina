//
//  certificate.cpp
//

#include "cobound/bounds/certificate.h"

#include <algorithm>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <functional>
#include <stdexcept>
#include <utility>

#include "cobound/json.h"
#include "cobound/frozen.h"

namespace bounds {

namespace {
// A certificate's file, opened for writing; one that cannot be opened throws.
std::ofstream openCertificate(const std::string &path) {
  std::ofstream c(path);
  if (!c) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
  return c;
}
// Whether its writes reached the file: flushed first, since a failed write is
// only seen once the buffer reaches the file. A failure throws: a goal met
// with no certificate on disk is never reported as met.
void finishCertificate(std::ofstream &c, const std::string &path) {
  c.flush();
  if (!c) throw std::runtime_error("writing " + path + " failed: " + std::strerror(errno));
}
} // namespace

// How a certificate finds a cobordism's surface: a master cobordism's pair
// signature inline; an in-process one's faces in its searched link's thickening, with
// that thickening's digest.
void CertificateWriter::writeSurface(std::ostream &c, const CobordismSource &info) {
  if (!info.pairsig.empty()) c << ",\"pairsig\":\"" << json::escape(info.pairsig) << "\"";
  if (!info.faces.empty()) {
    c << ",\"faces\":[";
    for (size_t i = 0; i < info.faces.size(); ++i) c << (i ? "," : "") << info.faces[i];
    c << "],\"build\":\"" << info.build << "\"";
  }
}

void CertificateWriter::writeCobordism(std::ostream &c, RelationId eid,
                                         std::set<LinkId> &links) const {
  // A cobordism as a checker replays it: its key, ends, shape and maps,
  // and (for one found by a search or read from the master) the searched
  // diagram, the surface (faces and build digest, or pair signature) and
  // each outgoing piece's match.
  const LinkCobordism &cob = g_.cobordism(eid);
  links.insert(cob.in);
  links.insert(cob.out);
  const auto it = sources_.byCobordism.find(eid);
  c << ",\"witness\":\"" << json::escape(cob.key) << "\"";
  c << ",\"in\":" << cob.in << ",\"out\":" << cob.out
    << ",\"shape\":{\"components\":" << cob.shape.components << ",\"genus\":" << cob.shape.genus
    << ",\"inComponent\":" << json::array(cob.shape.inComponent)
    << ",\"outComponent\":" << json::array(cob.shape.outComponent) << "},\"inMap\":" << json::array(cob.inMap)
    << ",\"outMap\":" << json::array(cob.outMap);
  if (it == sources_.byCobordism.end()) return;
  const AddedCobordism &added = it->second.added;
  c << ",\"hop_dir\":\"" << json::escape(it->second.searchDir) << "\",\"row_pd\":\""
    << json::escape(it->second.incomingPD) << "\",\"layers\":" << it->second.layers
    << ",\"row_node_map\":" << json::array(it->second.incomingLinkMap);
  writeSurface(c, it->second);
  c << ",\"split_edge\":" << added.split << ",\"farCurveEdges\":[";
  for (size_t j = 0; j < added.outgoingCurveEdges.size(); ++j)
    c << (j ? "," : "") << json::array(added.outgoingCurveEdges[j]);
  c << "],\"pieces\":[";
  for (size_t k = 0; k < added.pieces.size(); ++k) {
    const LinkMatch &m = added.pieces[k];
    c << (k ? "," : "") << "{\"node\":" << m.link << ",\"method\":\"" << m.method
      << "\",\"componentMap\":" << json::array(m.componentMap) << ",\"mirrored\":"
      << (m.mirrored ? "true" : "false") << ",\"reversed\":" << (m.reversed ? "true" : "false")
      << ",\"origins\":" << json::array(added.pieceOrigins[k]) << "}";
    links.insert(m.link);
  }
  c << "]";
}

void CertificateWriter::writeDerivations(std::ostream &c, const std::vector<DerivationId> &ids,
                                     std::set<LinkId> &links) const {
  bool first = true;
  for (DerivationId r : ids) {
    const Derivation &rec = g_.derivation(r);
    links.insert(rec.link);
    c << (first ? "" : ",\n") << "{\"id\":" << r << ",\"node\":" << rec.link
      << ",\"partition\":\"" << rec.partition.str() << "\",\"genus\":" << rec.genus
      << ",\"kind\":\"" << kindName(rec.kind) << "\",\"edge\":" << rec.relation
      << ",\"source\":\"" << json::escape(rec.source) << "\",\"children\":[";
    for (size_t i = 0; i < rec.children.size(); ++i)
      c << (i ? "," : "") << rec.children[i];
    c << "]";
    if (rec.kind == DerivationKind::cobordismForward || rec.kind == DerivationKind::cobordismReverse)
      writeCobordism(c, rec.relation, links);
    if (rec.kind == DerivationKind::splitCombine || rec.kind == DerivationKind::splitRestrict) {
      const Split &se = g_.split(rec.relation);
      links.insert(se.whole);
      links.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == DerivationKind::sumCombine) {
      const Sum &se = g_.sum(rec.relation);
      links.insert(se.whole);
      links.insert(se.pieces.begin(), se.pieces.end());
      c << ",\"whole\":" << se.whole << ",\"pieces\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k)
        c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
    }
    if (rec.kind == DerivationKind::leaf && rec.source.rfind(kFrozenDirectWitnessSource, 0) == 0) {
      const std::string key = rec.source.substr(sizeof kFrozenDirectWitnessSource - 1);
      if (auto it = sources_.direct.find(key); it != sources_.direct.end()) {
        c << ",\"witness\":\"" << json::escape(it->second.key)
          << "\",\"hop_dir\":\"" << json::escape(it->second.searchDir) << "\",\"row_pd\":\""
          << json::escape(it->second.incomingPD) << "\",\"layers\":" << it->second.layers;
        writeSurface(c, it->second);
      }
    }
    c << "}";
    first = false;
  }
}

void CertificateWriter::writeLower(const std::string &path, const CertificateGoal &goal) const {
  // The proof of lower(target, goal) as a tree of facts, children before
  // parents (README.md, "Lower goals"): each fact is a link, a
  // partition, the value the bound holds there, and its reason. A cobordism
  // fact carries the cobordism exactly as an upper derivation does, so the checker
  // replays the surface the same way, then recomputes the cap's addition
  // and the other end's partition itself. Split facts carry the upper
  // derivations they subtract, with those derivations' own proofs.
  using Kind = CobordismGraph::LowerReason::Kind;
  struct Fact {
    LinkId link;
    Partition q;
    CobordismGraph::LowerFact f;
    long from = -1;          ///< cobordism / split-piece: the fact read
    std::vector<long> pieces; ///< split-whole: the pieces' facts
  };
  std::vector<Fact> facts;
  std::map<std::pair<LinkId, std::vector<int>>, long> ids;
  std::set<std::pair<LinkId, std::vector<int>>> onStack;
  std::set<DerivationId> derivations;
  std::function<long(LinkId, const Partition &)> visit = [&](LinkId n, const Partition &q) -> long {
    const CobordismGraph::LowerFact f = g_.lowerWhy(n, q);
    const auto key = std::make_pair(n, f.storedFor.labels());
    if (auto it = ids.find(key); it != ids.end()) return it->second;
    if (!onStack.insert(key).second)
      throw std::logic_error("a lower bound's reasons cycle: node " + std::to_string(n));
    Fact fact{n, f.storedFor, f};
    switch (f.reason.kind) {
    case Kind::cobordism: {
      const LinkCobordism &e = g_.cobordism(f.reason.relation);
      fact.from = visit(f.reason.toIsIn ? e.out : e.in,
                        Partition::fromLabels(f.reason.fromPartition));
      break;
    }
    case Kind::splitWhole:
      for (RelationId sid : g_.link(n).splits) {
        const Split &s = g_.split(sid);
        if (s.whole != n) continue;
        for (size_t k = 0; k < s.pieces.size() && k < f.reason.pieces.size(); ++k)
          fact.pieces.push_back(visit(s.pieces[k], Partition::fromLabels(f.reason.pieces[k])));
        break;
      }
      break;
    case Kind::splitPiece:
      for (RelationId sid : g_.link(n).splits) {
        const Split &s = g_.split(sid);
        if (s.whole == n || std::find(s.pieces.begin(), s.pieces.end(), n) == s.pieces.end())
          continue;
        fact.from = visit(s.whole, Partition::fromLabels(f.reason.fromPartition));
        break;
      }
      for (DerivationId r : f.reason.derivations)
        for (DerivationId p : g_.proof(r)) derivations.insert(p);
      break;
    case Kind::sumPiece: {
      const Sum &s = g_.sum(f.reason.relation);
      const LinkId piece = s.pieces[static_cast<size_t>(f.reason.piece)];
      fact.from = visit(piece, Partition::coarsest(g_.link(piece).components));
      for (DerivationId r : f.reason.derivations)
        for (DerivationId p : g_.proof(r)) derivations.insert(p);
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
  std::ofstream c = openCertificate(path);
  c << "{\"target\":\"" << json::escape(goal.targetName) << "\",\"target_pd\":\""
    << json::escape(goal.targetPD) << "\",\"goal_lower\":" << goal.goalLower << ",\"goal\":\""
    << (goal.disjoint ? "disjoint" : "connected") << "\",\"lower\":"
    << g_.lower(goal.target, goal.partition) << ",\"top\":" << top << ",\"facts\":[\n";
  std::set<LinkId> links;
  auto value = [](int v) {
    return v >= CobordismGraph::kNoSurface ? std::string("\"inf\"") : std::to_string(v);
  };
  for (size_t i = 0; i < facts.size(); ++i) {
    const Fact &fact = facts[i];
    links.insert(fact.link);
    c << (i ? ",\n" : "") << "{\"id\":" << i << ",\"node\":" << fact.link << ",\"partition\":\""
      << fact.q.str() << "\",\"value\":" << value(fact.f.value);
    switch (fact.f.reason.kind) {
    case Kind::literature:
      c << ",\"kind\":\"literature\",\"source\":\"" << json::escape(g_.link(fact.link).lowerBoundSource)
        << "\"";
      break;
    case Kind::linking:
      c << ",\"kind\":\"linking\"";
      break;
    case Kind::cobordism: {
      const LinkCobordism &e = g_.cobordism(fact.f.reason.relation);
      c << ",\"kind\":\"" << kFrozenLowerKindWitness
        << "\",\"to_is_in\":" << (fact.f.reason.toIsIn ? "true" : "false")
        << ",\"from\":" << fact.from << ",\"from_partition\":\""
        << Partition::fromLabels(fact.f.reason.fromPartition).str() << "\",\"from_value\":"
        << value(fact.f.reason.from) << ",\"addition\":" << fact.f.reason.addition;
      writeCobordism(c, e.id, links);
      break;
    }
    case Kind::splitWhole: {
      c << ",\"kind\":\"split-whole\",\"pieces\":[";
      for (size_t k = 0; k < fact.pieces.size(); ++k) c << (k ? "," : "") << fact.pieces[k];
      c << "]";
      for (RelationId sid : g_.link(fact.link).splits) {
        const Split &se = g_.split(sid);
        if (se.whole != fact.link) continue;
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
      for (size_t k = 0; k < fact.f.reason.derivations.size(); ++k)
        c << (k ? "," : "") << fact.f.reason.derivations[k];
      c << "]";
      break;
    }
    case Kind::sumPiece: {
      const Sum &se = g_.sum(fact.f.reason.relation);
      c << ",\"kind\":\"sum-piece\",\"from\":" << fact.from << ",\"piece\":" << fact.f.reason.piece
        << ",\"from_value\":" << value(fact.f.reason.from) << ",\"subtracted\":"
        << fact.f.reason.addition << ",\"records\":[";
      for (size_t k = 0; k < fact.f.reason.derivations.size(); ++k)
        c << (k ? "," : "") << fact.f.reason.derivations[k];
      c << "],\"whole\":" << se.whole << ",\"piece_nodes\":" << json::array(se.pieces) << ",\"pieceMap\":[";
      for (size_t k = 0; k < se.pieceMap.size(); ++k) c << (k ? "," : "") << json::array(se.pieceMap[k]);
      c << "]";
      for (LinkId p : se.pieces) links.insert(p);
      break;
    }
    default:
      c << ",\"kind\":\"none\"";
      break;
    }
    c << "}";
  }
  c << "\n],\"records\":[\n";
  writeDerivations(c, std::vector<DerivationId>(derivations.begin(), derivations.end()), links);
  c << "\n],\"nodes\":[\n";
  writeLinks(c, links);
  c << "\n]}\n";
  finishCertificate(c, path);
}

void CertificateWriter::writeUpper(const std::string &path, const CertificateGoal &goal) const {
  auto best = g_.best(goal.target, goal.partition);
  if (!best) return;
  std::ofstream c = openCertificate(path);
  c << "{\"target\":\"" << json::escape(goal.targetName) << "\",\"target_pd\":\""
    << json::escape(goal.targetPD) << "\",\"goal_genus\":" << goal.goalGenus
    << ",\"goal\":\"" << (goal.disjoint ? "disjoint" : "connected") << "\",\"genus\":"
    << best->genus << ",\"records\":[\n";
  std::set<LinkId> links;
  writeDerivations(c, g_.proof(best->derivation), links);
  c << "\n],\"nodes\":[\n";
  writeLinks(c, links);
  c << "\n]}\n";
  finishCertificate(c, path);
}

void CertificateWriter::writeLinks(std::ostream &c, const std::set<LinkId> &links) const {
  bool first = true;
  for (LinkId n : links) {
    c << (first ? "" : ",\n") << "{\"id\":" << n << ",\"label\":\""
      << json::escape(g_.link(n).label) << "\",\"components\":" << g_.link(n).components;
    if (reg_.known(n) && reg_.info(n).diagram.crossings() > 0) {
      // The link's own diagram, as signed Gauss data: component maps refer
      // to ITS component order, which a PD round trip need not keep.
      const linknaming::GaussDiagram &d = reg_.info(n).diagram;
      c << ",\"pd\":\"" << json::escape(diagramPD(d)) << "\",\"signs\":[";
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

void CertificateWriter::describeLower(std::ostream &o, LinkId n, const Partition &q,
                                      int indent) const {
  using Kind = CobordismGraph::LowerReason::Kind;
  const auto fact = g_.lowerWhy(n, q);
  const std::string pad(static_cast<size_t>(2 * indent + 4), ' ');
  auto name = [&](LinkId m) {
    auto t = tableName_.find(m);
    return "node " + std::to_string(m) +
           (t == tableName_.end() ? std::string(" (untabulated)") : " (" + t->second + ")");
  };
  o << pad << name(n) << " " << q.str() << " >= "
    << (fact.value >= CobordismGraph::kNoSurface ? std::string("(no such surface)")
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
    o << ": " << g_.link(n).lowerBoundSource << "\n";
    return;
  case Kind::linking:
    o << ": the linking numbers forbid this partition\n";
    return;
  case Kind::seed:
    o << ": what-if seed\n";
    return;
  case Kind::sumPiece: {
    const Sum &s = g_.sum(fact.reason.relation);
    const LinkId piece = s.pieces[static_cast<size_t>(fact.reason.piece)];
    o << ": a sum along components: summand " << name(piece) << " >= " << fact.reason.from
      << ", minus the other summands' connected genera plus components minus one ("
      << fact.reason.addition << "; records";
    for (DerivationId r : fact.reason.derivations) o << " " << r;
    o << ")\n";
    describeLower(o, piece, Partition::coarsest(g_.link(piece).components), indent + 1);
    return;
  }
  case Kind::cobordism: {
    const LinkCobordism &e = g_.cobordism(fact.reason.relation);
    const LinkId other = fact.reason.toIsIn ? e.out : e.in;
    const Partition op = Partition::fromLabels(fact.reason.fromPartition);
    o << ": across witness " << e.key << " (genus " << e.shape.genus << ", "
      << e.shape.components << " pieces) from " << name(other) << " " << op.str() << " >= "
      << fact.reason.from << ", the cap adding " << fact.reason.addition << "\n";
    describeLower(o, other, op, indent + 1);
    return;
  }
  case Kind::splitWhole: {
    o << ": the split link's pieces, summed\n";
    for (RelationId sid : g_.link(n).splits) {
      const Split &s = g_.split(sid);
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
    for (DerivationId r : fact.reason.derivations) o << " " << r;
    o << ")\n";
    for (RelationId sid : g_.link(n).splits) {
      const Split &s = g_.split(sid);
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
