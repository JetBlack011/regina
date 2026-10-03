// proofgraph.cpp

#include "cobound/bounds/cobordismgraph.h"

#include "cobound/json.h"
#include "cobound/frozen.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <set>
#include <sstream>
#include <stdexcept>

namespace bounds {

const char *kindName(RecordKind k) {
  switch (k) {
  case RecordKind::leaf: return "leaf";
  case RecordKind::witnessForward: return kFrozenKindWitnessForward;
  case RecordKind::witnessReverse: return kFrozenKindWitnessReverse;
  case RecordKind::splitCombine: return "split-combine";
  case RecordKind::splitRestrict: return "split-restrict";
  case RecordKind::sumCombine: return "sum-combine";
  }
  return "?";
}

namespace {

void requireBijection(const std::vector<int> &map, int n, const char *what) {
  if (static_cast<int>(map.size()) != n)
    throw std::invalid_argument(std::string(what) + ": wrong size");
  std::vector<char> seen(n, 0);
  for (int v : map) {
    if (v < 0 || v >= n || seen[v])
      throw std::invalid_argument(std::string(what) + ": not a bijection");
    seen[v] = 1;
  }
}

} // namespace

NodeId ProofGraph::addNode(int components, std::string label,
                           std::optional<std::vector<std::vector<int>>> lk) {
  if (components < 1)
    throw std::invalid_argument("addNode: a link has at least one component");
  if (lk && static_cast<int>(lk->size()) != components)
    throw std::invalid_argument("addNode: linking matrix size");
  Node n;
  n.id = static_cast<NodeId>(nodes_.size());
  n.components = components;
  n.label = std::move(label);
  n.linking = std::move(lk);
  n.profile = Profile(components);
  nodes_.push_back(std::move(n));
  lower_.emplace_back();
  lowerReason_.emplace_back();
  return nodes_.back().id;
}

void ProofGraph::setGenusLowerBound(NodeId n, int lo, std::string source) {
  Node &node = nodes_.at(n);
  node.genusLowerBound = lo;
  node.lowerBoundSource = std::move(source);
  for (const ProfileEntry &e : node.profile.entries())
    if (e.genus < lo)
      contradictions_.push_back(node.label + ": record " +
                                std::to_string(e.record) + " has genus " +
                                std::to_string(e.genus) +
                                " below the lower bound " + std::to_string(lo) +
                                " (" + node.lowerBoundSource + ")");
}

RecordId ProofGraph::addLeaf(NodeId n, const Partition &p, int genus,
                             std::string source) {
  if (p.size() != nodes_.at(n).components)
    throw std::invalid_argument("addLeaf: partition size");
  return insert(n, p, genus, RecordKind::leaf, -1, {}, std::move(source));
}

EdgeId ProofGraph::addWitness(NodeId in, NodeId out, CobordismShape shape,
                              std::vector<int> inMap, std::vector<int> outMap,
                              std::string key) {
  shape.validate();
  requireBijection(inMap, nodes_.at(in).components, "addWitness inMap");
  requireBijection(outMap, nodes_.at(out).components, "addWitness outMap");
  if (inMap.size() != shape.inComponent.size() ||
      outMap.size() != shape.outComponent.size())
    throw std::invalid_argument("addWitness: map size != curve count");
  WitnessEdge e;
  e.id = static_cast<EdgeId>(witnesses_.size());
  e.in = in;
  e.out = out;
  e.shape = std::move(shape);
  e.inMap = std::move(inMap);
  e.outMap = std::move(outMap);
  e.key = std::move(key);
  witnesses_.push_back(e);
  nodes_[in].witnessEdges.push_back(e.id);
  if (out != in)
    nodes_[out].witnessEdges.push_back(e.id);
  // Existing records at either end can now be pushed across it.
  for (NodeId x : {in, out})
    for (const ProfileEntry &pe : nodes_[x].profile.entries())
      pending_.push_back(pe.record);
  return e.id;
}

EdgeId ProofGraph::addSplit(NodeId whole, std::vector<NodeId> pieces,
                            std::vector<std::vector<int>> pieceMap) {
  if (pieces.size() != pieceMap.size() || pieces.empty())
    throw std::invalid_argument("addSplit: pieces/maps");
  std::vector<int> all;
  for (size_t k = 0; k < pieces.size(); ++k) {
    if (static_cast<int>(pieceMap[k].size()) != nodes_.at(pieces[k]).components)
      throw std::invalid_argument("addSplit: piece map size");
    if (pieces[k] == whole)
      throw std::invalid_argument("addSplit: a piece is the whole");
    all.insert(all.end(), pieceMap[k].begin(), pieceMap[k].end());
  }
  requireBijection(all, nodes_.at(whole).components, "addSplit pieceMap");
  SplitEdge s;
  s.id = static_cast<EdgeId>(splits_.size());
  s.whole = whole;
  s.pieces = std::move(pieces);
  s.pieceMap = std::move(pieceMap);
  splits_.push_back(s);
  nodes_[whole].splitEdges.push_back(s.id);
  std::set<NodeId> seen;
  for (NodeId p : s.pieces)
    if (seen.insert(p).second)
      nodes_[p].splitEdges.push_back(s.id);
  for (NodeId x : seen)
    for (const ProfileEntry &pe : nodes_[x].profile.entries())
      pending_.push_back(pe.record);
  for (const ProfileEntry &pe : nodes_[whole].profile.entries())
    pending_.push_back(pe.record);
  return s.id;
}

EdgeId ProofGraph::addSum(NodeId whole, std::vector<NodeId> pieces,
                          std::vector<std::vector<int>> pieceMap) {
  if (pieces.size() != pieceMap.size() || pieces.size() < 2)
    throw std::invalid_argument("addSum: at least two pieces, one map each");
  const int n = nodes_.at(whole).components;
  std::vector<int> hits(n, 0);
  for (size_t k = 0; k < pieces.size(); ++k) {
    if (static_cast<int>(pieceMap[k].size()) != nodes_.at(pieces[k]).components)
      throw std::invalid_argument("addSum: piece map size");
    if (pieces[k] == whole)
      throw std::invalid_argument("addSum: a piece is the whole");
    for (int w : pieceMap[k]) {
      if (w < 0 || w >= n)
        throw std::invalid_argument("addSum: piece map out of range");
      ++hits[w];
    }
  }
  for (int h : hits)
    if (h == 0)
      throw std::invalid_argument("addSum: a component of the whole is no piece's");
  // The sum sites must form a tree (a sphere decomposition always does):
  // the graph with a vertex per piece and per whole component, and an edge
  // per piece component, is connected with pieces + n - 1 edges. Otherwise
  // two pieces would be summed twice and the genus formula would be wrong.
  size_t edges = 0;
  std::vector<int> parent(pieces.size() + static_cast<size_t>(n));
  std::iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> find = [&](int x) { return parent[x] == x ? x : parent[x] = find(parent[x]); };
  for (size_t k = 0; k < pieces.size(); ++k)
    for (int w : pieceMap[k]) {
      ++edges;
      parent[find(static_cast<int>(k))] = find(static_cast<int>(pieces.size()) + w);
    }
  int roots = 0;
  for (size_t v = 0; v < parent.size(); ++v) roots += find(static_cast<int>(v)) == static_cast<int>(v);
  if (roots != 1 || edges != pieces.size() + static_cast<size_t>(n) - 1)
    throw std::invalid_argument("addSum: the sum sites do not form a tree");
  SumEdge s;
  s.id = static_cast<EdgeId>(sums_.size());
  s.whole = whole;
  s.pieces = std::move(pieces);
  s.pieceMap = std::move(pieceMap);
  sums_.push_back(s);
  nodes_[whole].sumEdges.push_back(s.id);
  std::set<NodeId> seen;
  for (NodeId p : s.pieces)
    if (seen.insert(p).second)
      nodes_[p].sumEdges.push_back(s.id);
  for (NodeId x : seen)
    for (const ProfileEntry &pe : nodes_[x].profile.entries())
      pending_.push_back(pe.record);
  return s.id;
}

RecordId ProofGraph::insert(NodeId n, const Partition &p, int genus,
                            RecordKind kind, EdgeId edge,
                            std::vector<RecordId> children,
                            std::string source) {
  Node &node = nodes_.at(n);
  if (node.profile.implies(p, genus))
    return -1;
  const RecordId id = static_cast<RecordId>(records_.size());
  for (RecordId c : children)
    if (c < 0 || c >= id)
      throw std::logic_error("insert: a child record is not older");
  Record r;
  r.id = id;
  r.node = n;
  r.partition = p;
  r.genus = genus;
  r.kind = kind;
  r.edge = edge;
  r.children = std::move(children);
  r.source = std::move(source);
  records_.push_back(std::move(r));
  node.profile.insert(p, genus, id);
  pending_.push_back(id);

  if (node.genusLowerBound && genus < *node.genusLowerBound)
    contradictions_.push_back(
        node.label + ": record " + std::to_string(id) + " (" +
        kindName(kind) + ") gives genus " + std::to_string(genus) +
        " below the proved lower bound " +
        std::to_string(*node.genusLowerBound) + " (" +
        node.lowerBoundSource + ")");
  if (node.linking && !linkingAllows(p, *node.linking))
    contradictions_.push_back(
        node.label + ": record " + std::to_string(id) + " (" +
        kindName(kind) + ") claims partition " + p.str() +
        ", which its linking numbers forbid: a component map is wrong");
  return id;
}

std::optional<std::pair<Partition, int>>
ProofGraph::throughWitness(const WitnessEdge &e, bool forward,
                           const Partition &p, int genus) const {
  // forward: p partitions node `out`; glue on the outgoing side and read the
  // result on the incoming curves. reverse: the other way round.
  const std::vector<int> &gluedMap = forward ? e.outMap : e.inMap;
  const std::vector<int> &freeMap = forward ? e.inMap : e.outMap;
  std::vector<int> curveLabels(gluedMap.size());
  for (size_t j = 0; j < gluedMap.size(); ++j)
    curveLabels[j] = p.blockOf(gluedMap[j]);
  auto g = glue(e.shape, forward ? Side::outgoing : Side::incoming,
                Partition::fromLabels(curveLabels), genus);
  if (!g)
    return std::nullopt;
  std::vector<int> nodeLabels(freeMap.size());
  for (size_t i = 0; i < freeMap.size(); ++i)
    nodeLabels[freeMap[i]] = g->partition.blockOf(i);
  return std::make_pair(Partition::fromLabels(nodeLabels), g->genus);
}

std::optional<std::pair<Partition, int>>
ProofGraph::combineSplit(const SplitEdge &s,
                         const std::vector<const Record *> &perPiece) const {
  // Surfaces for the pieces, placed in disjoint balls: blocks never merge
  // across pieces, and genera add.
  std::vector<int> labels(nodes_[s.whole].components, -1);
  int genus = 0;
  int offset = 0;
  for (size_t k = 0; k < s.pieces.size(); ++k) {
    const Record &r = *perPiece[k];
    for (size_t c = 0; c < s.pieceMap[k].size(); ++c)
      labels[s.pieceMap[k][c]] = offset + r.partition.blockOf(c);
    offset += r.partition.blocks();
    genus += r.genus;
  }
  return std::make_pair(Partition::fromLabels(labels), genus);
}

std::pair<Partition, int> ProofGraph::combineSum(
    const SumEdge &s, const std::vector<const Record *> &perPiece) const {
  // Boundary-connected-sum the surfaces at each sum site: the two pieces
  // containing the summed components become one, their genera add and no
  // genus is created (paper lem:sum-partitions). So the whole's blocks are
  // the unions, over the summands' blocks, of their images, merged wherever
  // two piece components share an image.
  const int n = nodes_[s.whole].components;
  std::vector<int> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> find = [&](int x) {
    return parent[x] == x ? x : parent[x] = find(parent[x]);
  };
  int genus = 0;
  for (size_t k = 0; k < s.pieces.size(); ++k) {
    const Record &r = *perPiece[k];
    genus += r.genus;
    for (int b = 0; b < r.partition.blocks(); ++b) {
      int first = -1;
      for (size_t c = 0; c < s.pieceMap[k].size(); ++c) {
        if (r.partition.blockOf(static_cast<int>(c)) != b) continue;
        const int w = s.pieceMap[k][c];
        if (first < 0) first = w;
        else parent[find(w)] = find(first);
      }
    }
  }
  std::vector<int> labels(n);
  for (int w = 0; w < n; ++w) labels[w] = find(w);
  return {Partition::fromLabels(labels), genus};
}

std::optional<Partition> ProofGraph::restrictSplit(const SplitEdge &s,
                                                   size_t piece,
                                                   const Partition &whole) const {
  // Only a partition whose blocks each lie within one piece restricts: then
  // the surface's pieces bounding this piece's components form a surface for
  // it, of total genus at most the whole's.
  std::vector<int> pieceOf(nodes_[s.whole].components, -1);
  for (size_t k = 0; k < s.pieces.size(); ++k)
    for (int w : s.pieceMap[k])
      pieceOf[w] = static_cast<int>(k);
  std::vector<int> blockPiece(whole.blocks(), -1);
  for (int w = 0; w < whole.size(); ++w) {
    int &bp = blockPiece[whole.blockOf(w)];
    if (bp < 0)
      bp = pieceOf[w];
    else if (bp != pieceOf[w])
      return std::nullopt;
  }
  std::vector<int> labels(s.pieceMap[piece].size());
  for (size_t c = 0; c < labels.size(); ++c)
    labels[c] = whole.blockOf(s.pieceMap[piece][c]);
  return Partition::fromLabels(labels);
}

void ProofGraph::deriveFrom(RecordId rid) {
  // Copy: insert() may reallocate records_ and nodes_ entries' vectors.
  const Record r = records_[rid];
  const Node &n = nodes_[r.node];
  const std::vector<EdgeId> wEdges = n.witnessEdges;
  const std::vector<EdgeId> sEdges = n.splitEdges;
  const std::vector<EdgeId> mEdges = n.sumEdges;

  // As a summand: combine with every current entry of the other summands.
  for (EdgeId mid : mEdges) {
    const SumEdge s = sums_[mid];
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      if (s.pieces[k] != r.node) continue;
      std::vector<std::vector<RecordId>> options(s.pieces.size());
      bool anyEmpty = false;
      for (size_t j = 0; j < s.pieces.size(); ++j) {
        if (j == k) {
          options[j] = {rid};
          continue;
        }
        for (const ProfileEntry &pe : nodes_[s.pieces[j]].profile.entries())
          options[j].push_back(pe.record);
        if (options[j].empty()) anyEmpty = true;
      }
      if (anyEmpty) continue;
      std::vector<std::vector<RecordId>> combos;
      std::vector<RecordId> pick(s.pieces.size());
      std::function<void(size_t)> rec = [&](size_t j) {
        if (j == s.pieces.size()) {
          combos.push_back(pick);
          return;
        }
        for (RecordId o : options[j]) {
          pick[j] = o;
          rec(j + 1);
        }
      };
      rec(0);
      for (const auto &ids : combos) {
        std::vector<const Record *> recs;
        for (RecordId id : ids) recs.push_back(&records_[id]);
        auto d = combineSum(s, recs);
        insert(s.whole, d.first, d.second, RecordKind::sumCombine, mid, ids, "");
      }
    }
  }

  for (EdgeId eid : wEdges) {
    const WitnessEdge e = witnesses_[eid];
    if (e.out == r.node)
      if (auto d = throughWitness(e, /*forward=*/true, r.partition, r.genus))
        insert(e.in, d->first, d->second, RecordKind::witnessForward, eid,
               {rid}, "");
    if (e.in == r.node)
      if (auto d = throughWitness(e, /*forward=*/false, r.partition, r.genus))
        insert(e.out, d->first, d->second, RecordKind::witnessReverse, eid,
               {rid}, "");
  }

  for (EdgeId sid : sEdges) {
    const SplitEdge s = splits_[sid];
    if (s.whole == r.node) {
      for (size_t k = 0; k < s.pieces.size(); ++k)
        if (auto p = restrictSplit(s, k, r.partition))
          insert(s.pieces[k], *p, r.genus, RecordKind::splitRestrict, sid,
                 {rid}, "");
    }
    // As a piece (possibly several times, if the same node appears twice):
    // combine with every current entry of the other pieces.
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      if (s.pieces[k] != r.node)
        continue;
      std::vector<std::vector<const Record *>> options(s.pieces.size());
      bool anyEmpty = false;
      for (size_t j = 0; j < s.pieces.size(); ++j) {
        if (j == k) {
          options[j] = {&records_[rid]};
          continue;
        }
        for (const ProfileEntry &pe : nodes_[s.pieces[j]].profile.entries())
          options[j].push_back(&records_[pe.record]);
        if (options[j].empty())
          anyEmpty = true;
      }
      if (anyEmpty)
        continue;
      // Collect every combination first: insert() can reallocate records_.
      std::vector<std::vector<RecordId>> combos;
      std::vector<const Record *> pick(s.pieces.size());
      std::function<void(size_t)> rec = [&](size_t j) {
        if (j == s.pieces.size()) {
          std::vector<RecordId> ids;
          for (const Record *p : pick)
            ids.push_back(p->id);
          combos.push_back(std::move(ids));
          return;
        }
        for (const Record *o : options[j]) {
          pick[j] = o;
          rec(j + 1);
        }
      };
      rec(0);
      for (const auto &ids : combos) {
        std::vector<const Record *> recs;
        for (RecordId id : ids)
          recs.push_back(&records_[id]);
        auto d = combineSplit(s, recs);
        insert(s.whole, d->first, d->second, RecordKind::splitCombine, sid,
               ids, "");
      }
    }
  }
}

long ProofGraph::propagate() {
  const size_t before = records_.size();
  // Records are processed in creation order; each derivation only ever
  // creates strictly improving records (insert()), and genus is bounded
  // below by 0 over finitely many partitions, so this terminates.
  while (!pending_.empty()) {
    std::vector<RecordId> batch;
    batch.swap(pending_);
    std::sort(batch.begin(), batch.end());
    batch.erase(std::unique(batch.begin(), batch.end()), batch.end());
    for (RecordId r : batch)
      deriveFrom(r);
  }
  return static_cast<long>(records_.size() - before);
}

// ---------------------------------------------------------------- lower bounds

int ProofGraph::lower(NodeId n, const Partition &q) const {
  const Node &node = nodes_.at(n);
  if (node.linking && !linkingAllows(q, *node.linking))
    return kNoSurface;
  int v = std::max(0, node.genusLowerBound.value_or(0));
  // Every surface refining q refines each coarser partition, so a bound
  // stored for any partition that q refines applies to q.
  for (const auto &[labels, val] : lower_.at(n))
    if (val > v && q.refines(Partition::fromLabels(labels)))
      v = val;
  return v;
}

void ProofGraph::clearLowerBounds() {
  for (Node &node : nodes_) {
    node.genusLowerBound.reset();
    node.lowerBoundSource.clear();
  }
  for (auto &m : lower_) m.clear();
  for (auto &m : lowerReason_) m.clear();
  ++lowerVersion_;
}

bool ProofGraph::raiseLower(NodeId n, const Partition &q, int value,
                            const LowerReason &why) {
  if (nodes_.at(n).components > kMaxLowerComponents)
    return false;
  if (value <= lower(n, q))
    return false;
  lower_[n][q.labels()] = std::min(value, kNoSurface);
  lowerReason_[n][q.labels()] = why;
  ++lowerVersion_;
  return true;
}

ProofGraph::LowerFact ProofGraph::lowerWhy(NodeId n, const Partition &q) const {
  // The same choice lower() makes, with its reason.
  const Node &node = nodes_.at(n);
  LowerFact f;
  f.storedFor = q;
  if (node.linking && !linkingAllows(q, *node.linking)) {
    f.value = kNoSurface;
    f.reason.kind = LowerReason::Kind::linking;
    return f;
  }
  if (node.genusLowerBound && *node.genusLowerBound > 0) {
    f.value = *node.genusLowerBound;
    f.storedFor = Partition::coarsest(node.components);
    f.reason.kind = LowerReason::Kind::literature;
  }
  for (const auto &[labels, val] : lower_.at(n))
    if (val > f.value && q.refines(Partition::fromLabels(labels))) {
      f.value = val;
      f.storedFor = Partition::fromLabels(labels);
      f.reason = lowerReason_.at(n).at(labels);
    }
  return f;
}

int ProofGraph::transportedLower(const WitnessEdge &e, bool toIsIn,
                                 const Partition &p, Transport *detail) const {
  // A surface F for the `to` end with partition p, capped onto e, gives a
  // surface for the other end whose partition we compute, of genus at most
  // genus(F) + (what glue() adds with a genus-0 cap). So genus(F) is at
  // least the other end's bound for that partition minus that addition.
  const NodeId to = toIsIn ? e.in : e.out, other = toIsIn ? e.out : e.in;
  const std::vector<int> &toMap = toIsIn ? e.inMap : e.outMap;
  const std::vector<int> &otherMap = toIsIn ? e.outMap : e.inMap;
  if (nodes_[to].linking && !linkingAllows(p, *nodes_[to].linking))
    return kNoSurface; // no surface has partition p
  std::vector<int> curveLabels(toMap.size());
  for (size_t i = 0; i < toMap.size(); ++i)
    curveLabels[i] = p.blockOf(toMap[i]);
  auto g = glue(e.shape, toIsIn ? Side::incoming : Side::outgoing,
                Partition::fromLabels(curveLabels), 0);
  if (!g)
    return 0; // the other end has no curves: nothing transported
  std::vector<int> otherLabels(otherMap.size());
  for (size_t j = 0; j < otherMap.size(); ++j)
    otherLabels[otherMap[j]] = g->partition.blockOf(j);
  const Partition op = Partition::fromLabels(otherLabels);
  const int lo = lower(other, op);
  if (detail) {
    detail->otherPartition = op;
    detail->addition = g->genus;
    detail->other = lo;
  }
  return lo >= kNoSurface ? kNoSurface : lo - g->genus;
}

std::string ProofGraph::profileFields(NodeId n) const {
  const Node &node = nodes_.at(n);
  std::ostringstream o;
  o << "\"components\":" << node.components;
  if (node.linking) o << ",\"linking\":" << json::matrix(*node.linking);
  if (node.genusLowerBound)
    o << ",\"genus_lower\":" << *node.genusLowerBound;
  std::vector<ProfileEntry> es = node.profile.entries();
  std::sort(es.begin(), es.end(), [](const ProfileEntry &a, const ProfileEntry &b) {
    return a.partition == b.partition ? a.genus < b.genus : a.partition < b.partition;
  });
  o << ",\"entries\":[";
  for (size_t i = 0; i < es.size(); ++i)
    o << (i ? "," : "") << "{\"p\":\"" << es[i].partition.str() << "\",\"g\":" << es[i].genus
      << ",\"r\":" << es[i].record << '}';
  o << ']';
  if (node.components <= kMaxLowerComponents) {
    std::vector<Partition> all = allPartitions(node.components);
    std::sort(all.begin(), all.end());
    o << ",\"lower\":[";
    bool first = true;
    for (const Partition &q : all) {
      const int lo = lower(n, q);
      if (lo <= 0)
        continue;
      o << (first ? "" : ",") << "{\"p\":\"" << q.str() << "\",";
      if (lo >= kNoSurface)
        o << "\"forbidden\":true}";
      else
        o << "\"lo\":" << lo << '}';
      first = false;
    }
    o << ']';
  }
  return o.str();
}

int ProofGraph::lowerAcross(const WitnessEdge &e, bool toIsIn,
                            const Partition &q) const {
  // Every surface refining q is bounded by the minimum of transportedLower
  // over the refinements of q, and that minimum is at q itself
  // (lem:transport-monotone): refining the cap adds vertices to the gluing
  // graph, so add = E - V + #components can only fall, the other end's
  // partition only refines, and lower() is monotone under refinement.
  return transportedLower(e, toIsIn, q);
}

std::optional<int> ProofGraph::lowerIf(const std::vector<LowerSeed> &seeds,
                                       NodeId target,
                                       const Partition &goal) const {
  ProofGraph what = *this;
  const size_t before = what.contradictions().size();
  LowerReason seed;
  seed.kind = LowerReason::Kind::seed;
  for (const LowerSeed &s : seeds)
    what.raiseLower(s.node, s.partition, s.value, seed);
  what.propagateLower();
  if (what.contradictions().size() != before)
    return std::nullopt;
  return what.lower(target, goal);
}

long ProofGraph::propagateLower() {
  long improved = 0;
  for (size_t n = 0; n < nodes_.size(); ++n)
    if (nodes_[n].genusLowerBound) {
      LowerReason why;
      why.kind = LowerReason::Kind::literature;
      improved += raiseLower(static_cast<NodeId>(n),
                             Partition::coarsest(nodes_[n].components),
                             *nodes_[n].genusLowerBound, why);
    }
  const int maxPasses = 200;
  bool changed = true;
  int pass = 0;
  for (; changed && pass < maxPasses; ++pass) {
    changed = false;
    for (const WitnessEdge &e : witnesses_)
      for (bool toIsIn : {true, false}) {
        const NodeId to = toIsIn ? e.in : e.out;
        if (nodes_[to].components > kMaxLowerComponents)
          continue;
        for (const Partition &q : allPartitions(nodes_[to].components)) {
          Transport t;
          const int v = transportedLower(e, toIsIn, q, &t);
          if (v <= lower(to, q))
            continue;
          LowerReason why;
          why.kind = LowerReason::Kind::witness;
          why.edge = e.id;
          why.toIsIn = toIsIn;
          why.fromPartition = t.otherPartition.labels();
          why.from = t.other;
          why.addition = t.addition;
          if (raiseLower(to, q, v, why)) {
            changed = true;
            ++improved;
          }
        }
      }
    for (const SplitEdge &s : splits_) {
      const Node &w = nodes_[s.whole];
      if (w.components <= kMaxLowerComponents) {
        // The whole, from its pieces, for partitions that never mix pieces.
        for (const Partition &q : allPartitions(w.components)) {
          int sum = 0;
          bool mixes = false;
          LowerReason why;
          why.kind = LowerReason::Kind::splitWhole;
          for (size_t k = 0; k < s.pieces.size() && !mixes; ++k) {
            auto r = restrictSplit(s, k, q);
            if (!r) {
              mixes = true;
              break;
            }
            const int lo = lower(s.pieces[k], *r);
            sum = lo >= kNoSurface ? kNoSurface : std::min(kNoSurface, sum + lo);
            why.pieces.push_back(r->labels());
          }
          if (!mixes && raiseLower(s.whole, q, sum, why)) {
            changed = true;
            ++improved;
          }
        }
      }
      // A piece, from the whole and PROVED surfaces for the other pieces:
      // F_k u (achieved surfaces for the rest) bounds the whole, so
      // genus(F_k) >= lower(whole, combined) - their genera.
      for (size_t k = 0; k < s.pieces.size(); ++k) {
        const Node &pk = nodes_[s.pieces[k]];
        if (pk.components > kMaxLowerComponents ||
            w.components > kMaxLowerComponents)
          continue;
        std::vector<std::vector<const ProfileEntry *>> others(s.pieces.size());
        bool anyEmpty = false;
        for (size_t j = 0; j < s.pieces.size(); ++j) {
          if (j == k) continue;
          for (const ProfileEntry &pe : nodes_[s.pieces[j]].profile.entries())
            others[j].push_back(&pe);
          if (others[j].empty()) anyEmpty = true;
        }
        if (anyEmpty) continue;
        for (const Partition &q : allPartitions(pk.components)) {
          int bestV = 0;
          LowerReason why;
          why.kind = LowerReason::Kind::splitPiece;
          std::vector<const ProfileEntry *> pick(s.pieces.size(), nullptr);
          std::function<void(size_t)> rec = [&](size_t j) {
            if (j == s.pieces.size()) {
              std::vector<int> labels(w.components, -1);
              int offset = 0, genus = 0;
              for (size_t t = 0; t < s.pieces.size(); ++t) {
                const Partition &pt = t == k ? q : pick[t]->partition;
                for (size_t c = 0; c < s.pieceMap[t].size(); ++c)
                  labels[s.pieceMap[t][c]] = offset + pt.blockOf(c);
                offset += pt.blocks();
                if (t != k) genus += pick[t]->genus;
              }
              const Partition wp = Partition::fromLabels(labels);
              const int lo = lower(s.whole, wp);
              const int v = lo >= kNoSurface ? kNoSurface : lo - genus;
              if (v > bestV) {
                bestV = v;
                why.fromPartition = wp.labels();
                why.from = lo;
                why.addition = genus;
                why.records.clear();
                for (size_t t = 0; t < s.pieces.size(); ++t)
                  if (t != k) why.records.push_back(pick[t]->record);
              }
              return;
            }
            if (j == k) {
              rec(j + 1);
              return;
            }
            for (const ProfileEntry *pe : others[j]) {
              pick[j] = pe;
              rec(j + 1);
            }
          };
          rec(0);
          if (raiseLower(s.pieces[k], q, bestV, why)) {
            changed = true;
            ++improved;
          }
        }
      }
    }
  }
  // Sums along components: the whole's CONNECTED bound from one summand's
  // (paper cor:sum-pieces): every surface for the whole gives, with the
  // other summands capped by proved connected surfaces, a surface for piece
  // i of genus at most g + sum_j (h_j + n_j - 1). Kept outside the pass
  // loop since it reads only literature-seeded and transported bounds of
  // the pieces; one extra relaxation round suffices for what it adds.
  for (int round = 0; round < 2; ++round)
    for (const SumEdge &s : sums_) {
      const Node &w = nodes_[s.whole];
      if (w.components > kMaxLowerComponents) continue;
      for (size_t i = 0; i < s.pieces.size(); ++i) {
        const int lo = lower(s.pieces[i], Partition::coarsest(nodes_[s.pieces[i]].components));
        if (lo <= 0 || lo >= kNoSurface) continue;
        int subtract = 0;
        bool known = true;
        LowerReason why;
        why.kind = LowerReason::Kind::sumPiece;
        why.edge = s.id;
        why.piece = static_cast<int>(i);
        why.from = lo;
        for (size_t j = 0; j < s.pieces.size(); ++j) {
          if (j == i) continue;
          auto b = bestConnected(s.pieces[j]);
          if (!b) { known = false; break; }
          subtract += b->genus + nodes_[s.pieces[j]].components - 1;
          why.records.push_back(b->record);
        }
        if (!known) continue;
        why.addition = subtract;
        if (raiseLower(s.whole, Partition::coarsest(w.components), lo - subtract, why))
          ++improved;
      }
    }
  if (changed)
    contradictions_.push_back(
        "lower bounds did not converge in " + std::to_string(maxPasses) +
        " passes: the facts are inconsistent");
  // Every proved surface must respect every lower bound.
  for (const Node &n : nodes_)
    for (const ProfileEntry &e : n.profile.entries()) {
      const int lo = lower(n.id, e.partition);
      if (e.genus < lo)
        contradictions_.push_back(
            n.label + ": record " + std::to_string(e.record) + " gives " +
            e.partition.str() + " genus " + std::to_string(e.genus) +
            " below the propagated lower bound " +
            (lo >= kNoSurface ? std::string("(no such surface)") : std::to_string(lo)));
    }
  return improved;
}

long ProofGraph::saturate() {
  for (const Record &r : records_)
    pending_.push_back(r.id);
  return propagate();
}

std::optional<ProfileEntry> ProofGraph::best(NodeId n,
                                             const Partition &target) const {
  return nodes_.at(n).profile.best(target);
}

std::optional<ProfileEntry> ProofGraph::bestConnected(NodeId n) const {
  return best(n, Partition::coarsest(nodes_.at(n).components));
}

std::vector<RecordId> ProofGraph::proof(RecordId r) const {
  std::set<RecordId> seen;
  std::vector<RecordId> stack{r};
  while (!stack.empty()) {
    RecordId x = stack.back();
    stack.pop_back();
    if (!seen.insert(x).second)
      continue;
    for (RecordId c : records_.at(x).children)
      stack.push_back(c);
  }
  // Children always have smaller ids, so ascending id order is topological.
  return {seen.begin(), seen.end()};
}

std::string ProofGraph::recheck(RecordId rid) const {
  const Record &r = records_.at(rid);
  std::optional<std::pair<Partition, int>> d;
  switch (r.kind) {
  case RecordKind::leaf:
    return r.children.empty() ? "" : "a leaf with children";
  case RecordKind::witnessForward:
  case RecordKind::witnessReverse: {
    if (r.children.size() != 1)
      return "a witness record needs exactly one child";
    const WitnessEdge &e = witnesses_.at(r.edge);
    const Record &c = records_.at(r.children[0]);
    const bool fwd = r.kind == RecordKind::witnessForward;
    if (c.node != (fwd ? e.out : e.in) || r.node != (fwd ? e.in : e.out))
      return "endpoints do not match the edge";
    d = throughWitness(e, fwd, c.partition, c.genus);
    break;
  }
  case RecordKind::splitCombine: {
    const SplitEdge &s = splits_.at(r.edge);
    if (r.children.size() != s.pieces.size() || r.node != s.whole)
      return "split-combine children/endpoint mismatch";
    std::vector<const Record *> recs;
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      recs.push_back(&records_.at(r.children[k]));
      if (recs.back()->node != s.pieces[k])
        return "split-combine child on the wrong piece";
    }
    d = combineSplit(s, recs);
    break;
  }
  case RecordKind::sumCombine: {
    const SumEdge &s = sums_.at(r.edge);
    if (r.children.size() != s.pieces.size() || r.node != s.whole)
      return "sum-combine children/endpoint mismatch";
    std::vector<const Record *> recs;
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      recs.push_back(&records_.at(r.children[k]));
      if (recs.back()->node != s.pieces[k])
        return "sum-combine child on the wrong piece";
    }
    d = combineSum(s, recs);
    break;
  }
  case RecordKind::splitRestrict: {
    const SplitEdge &s = splits_.at(r.edge);
    if (r.children.size() != 1)
      return "split-restrict needs one child";
    const Record &c = records_.at(r.children[0]);
    if (c.node != s.whole)
      return "split-restrict child is not the whole";
    for (size_t k = 0; k < s.pieces.size(); ++k)
      if (s.pieces[k] == r.node)
        if (auto p = restrictSplit(s, k, c.partition))
          if (*p == r.partition && c.genus == r.genus)
            return "";
    return "split-restrict does not reproduce";
  }
  }
  if (!d)
    return "the derivation gives nothing";
  if (!(d->first == r.partition) || d->second != r.genus) {
    std::ostringstream o;
    o << "re-derived " << d->first.str() << " g" << d->second << ", recorded "
      << r.partition.str() << " g" << r.genus;
    return o.str();
  }
  return "";
}

} // namespace bounds
