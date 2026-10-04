// cobordismgraph.cpp

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

const char *kindName(DerivationKind k) {
  switch (k) {
  case DerivationKind::leaf: return "leaf";
  case DerivationKind::cobordismForward: return kFrozenKindWitnessForward;
  case DerivationKind::cobordismReverse: return kFrozenKindWitnessReverse;
  case DerivationKind::splitCombine: return "split-combine";
  case DerivationKind::splitRestrict: return "split-restrict";
  case DerivationKind::sumCombine: return "sum-combine";
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

LinkId CobordismGraph::addLink(int components, std::string label,
                           std::optional<std::vector<std::vector<int>>> lk) {
  if (components < 1)
    throw std::invalid_argument("addLink: a link has at least one component");
  if (lk && static_cast<int>(lk->size()) != components)
    throw std::invalid_argument("addLink: linking matrix size");
  GraphLink n;
  n.id = static_cast<LinkId>(links_.size());
  n.components = components;
  n.label = std::move(label);
  n.linking = std::move(lk);
  n.partitionGenera = PartitionGenera(components);
  links_.push_back(std::move(n));
  lower_.emplace_back();
  lowerReason_.emplace_back();
  return links_.back().id;
}

void CobordismGraph::setGenusLowerBound(LinkId n, int lo, std::string source) {
  GraphLink &link = links_.at(n);
  link.genusLowerBound = lo;
  link.lowerBoundSource = std::move(source);
  for (const PartitionGenus &e : link.partitionGenera.entries())
    if (e.genus < lo)
      contradictions_.push_back(link.label + ": record " +
                                std::to_string(e.derivation) + " has genus " +
                                std::to_string(e.genus) +
                                " below the lower bound " + std::to_string(lo) +
                                " (" + link.lowerBoundSource + ")");
}

DerivationId CobordismGraph::addLeaf(LinkId n, const Partition &p, int genus,
                             std::string source) {
  if (p.size() != links_.at(n).components)
    throw std::invalid_argument("addLeaf: partition size");
  return insert(n, p, genus, DerivationKind::leaf, -1, {}, std::move(source));
}

RelationId CobordismGraph::addCobordism(LinkId in, LinkId out, CobordismShape shape,
                              std::vector<int> inMap, std::vector<int> outMap,
                              std::string key) {
  shape.validate();
  requireBijection(inMap, links_.at(in).components, "addCobordism inMap");
  requireBijection(outMap, links_.at(out).components, "addCobordism outMap");
  if (inMap.size() != shape.inComponent.size() ||
      outMap.size() != shape.outComponent.size())
    throw std::invalid_argument("addCobordism: map size != curve count");
  LinkCobordism e;
  e.id = static_cast<RelationId>(cobordisms_.size());
  e.in = in;
  e.out = out;
  e.shape = std::move(shape);
  e.inMap = std::move(inMap);
  e.outMap = std::move(outMap);
  e.key = std::move(key);
  cobordisms_.push_back(e);
  links_[in].cobordisms.push_back(e.id);
  if (out != in)
    links_[out].cobordisms.push_back(e.id);
  // Existing derivations at either end can now be pushed across it.
  for (LinkId x : {in, out})
    for (const PartitionGenus &pe : links_[x].partitionGenera.entries())
      pending_.push_back(pe.derivation);
  return e.id;
}

RelationId CobordismGraph::addSplit(LinkId whole, std::vector<LinkId> pieces,
                            std::vector<std::vector<int>> pieceMap) {
  if (pieces.size() != pieceMap.size() || pieces.empty())
    throw std::invalid_argument("addSplit: pieces/maps");
  std::vector<int> all;
  for (size_t k = 0; k < pieces.size(); ++k) {
    if (static_cast<int>(pieceMap[k].size()) != links_.at(pieces[k]).components)
      throw std::invalid_argument("addSplit: piece map size");
    if (pieces[k] == whole)
      throw std::invalid_argument("addSplit: a piece is the whole");
    all.insert(all.end(), pieceMap[k].begin(), pieceMap[k].end());
  }
  requireBijection(all, links_.at(whole).components, "addSplit pieceMap");
  Split s;
  s.id = static_cast<RelationId>(splits_.size());
  s.whole = whole;
  s.pieces = std::move(pieces);
  s.pieceMap = std::move(pieceMap);
  splits_.push_back(s);
  links_[whole].splits.push_back(s.id);
  std::set<LinkId> seen;
  for (LinkId p : s.pieces)
    if (seen.insert(p).second)
      links_[p].splits.push_back(s.id);
  for (LinkId x : seen)
    for (const PartitionGenus &pe : links_[x].partitionGenera.entries())
      pending_.push_back(pe.derivation);
  for (const PartitionGenus &pe : links_[whole].partitionGenera.entries())
    pending_.push_back(pe.derivation);
  return s.id;
}

RelationId CobordismGraph::addSum(LinkId whole, std::vector<LinkId> pieces,
                          std::vector<std::vector<int>> pieceMap) {
  if (pieces.size() != pieceMap.size() || pieces.size() < 2)
    throw std::invalid_argument("addSum: at least two pieces, one map each");
  const int n = links_.at(whole).components;
  std::vector<int> hits(n, 0);
  for (size_t k = 0; k < pieces.size(); ++k) {
    if (static_cast<int>(pieceMap[k].size()) != links_.at(pieces[k]).components)
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
  Sum s;
  s.id = static_cast<RelationId>(sums_.size());
  s.whole = whole;
  s.pieces = std::move(pieces);
  s.pieceMap = std::move(pieceMap);
  sums_.push_back(s);
  links_[whole].sums.push_back(s.id);
  std::set<LinkId> seen;
  for (LinkId p : s.pieces)
    if (seen.insert(p).second)
      links_[p].sums.push_back(s.id);
  for (LinkId x : seen)
    for (const PartitionGenus &pe : links_[x].partitionGenera.entries())
      pending_.push_back(pe.derivation);
  return s.id;
}

DerivationId CobordismGraph::insert(LinkId n, const Partition &p, int genus,
                            DerivationKind kind, RelationId relation,
                            std::vector<DerivationId> children,
                            std::string source) {
  GraphLink &link = links_.at(n);
  if (link.partitionGenera.implies(p, genus))
    return -1;
  const DerivationId id = static_cast<DerivationId>(derivations_.size());
  for (DerivationId c : children)
    if (c < 0 || c >= id)
      throw std::logic_error("insert: a child record is not older");
  Derivation r;
  r.id = id;
  r.link = n;
  r.partition = p;
  r.genus = genus;
  r.kind = kind;
  r.relation = relation;
  r.children = std::move(children);
  r.source = std::move(source);
  derivations_.push_back(std::move(r));
  link.partitionGenera.insert(p, genus, id);
  pending_.push_back(id);

  if (link.genusLowerBound && genus < *link.genusLowerBound)
    contradictions_.push_back(
        link.label + ": record " + std::to_string(id) + " (" +
        kindName(kind) + ") gives genus " + std::to_string(genus) +
        " below the proved lower bound " +
        std::to_string(*link.genusLowerBound) + " (" +
        link.lowerBoundSource + ")");
  if (link.linking && !linkingAllows(p, *link.linking))
    contradictions_.push_back(
        link.label + ": record " + std::to_string(id) + " (" +
        kindName(kind) + ") claims partition " + p.str() +
        ", which its linking numbers forbid: a component map is wrong");
  return id;
}

std::optional<std::pair<Partition, int>>
CobordismGraph::throughCobordism(const LinkCobordism &e, bool forward,
                           const Partition &p, int genus) const {
  // forward: p partitions link `out`; glue on the outgoing side and read the
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
  std::vector<int> linkLabels(freeMap.size());
  for (size_t i = 0; i < freeMap.size(); ++i)
    linkLabels[freeMap[i]] = g->partition.blockOf(i);
  return std::make_pair(Partition::fromLabels(linkLabels), g->genus);
}

std::optional<std::pair<Partition, int>>
CobordismGraph::combineSplit(const Split &s,
                         const std::vector<const Derivation *> &perPiece) const {
  // Surfaces for the pieces, placed in disjoint balls: blocks never merge
  // across pieces, and genera add.
  std::vector<int> labels(links_[s.whole].components, -1);
  int genus = 0;
  int offset = 0;
  for (size_t k = 0; k < s.pieces.size(); ++k) {
    const Derivation &r = *perPiece[k];
    for (size_t c = 0; c < s.pieceMap[k].size(); ++c)
      labels[s.pieceMap[k][c]] = offset + r.partition.blockOf(c);
    offset += r.partition.blocks();
    genus += r.genus;
  }
  return std::make_pair(Partition::fromLabels(labels), genus);
}

std::pair<Partition, int> CobordismGraph::combineSum(
    const Sum &s, const std::vector<const Derivation *> &perPiece) const {
  // Boundary-connected-sum the surfaces at each sum site: the two pieces
  // containing the summed components become one, their genera add and no
  // genus is created (paper lem:sum-partitions). So the whole's blocks are
  // the unions, over the summands' blocks, of their images, merged wherever
  // two piece components share an image.
  const int n = links_[s.whole].components;
  std::vector<int> parent(n);
  std::iota(parent.begin(), parent.end(), 0);
  std::function<int(int)> find = [&](int x) {
    return parent[x] == x ? x : parent[x] = find(parent[x]);
  };
  int genus = 0;
  for (size_t k = 0; k < s.pieces.size(); ++k) {
    const Derivation &r = *perPiece[k];
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

std::optional<Partition> CobordismGraph::restrictSplit(const Split &s,
                                                   size_t piece,
                                                   const Partition &whole) const {
  // Only a partition whose blocks each lie within one piece restricts: then
  // the surface's pieces bounding this piece's components form a surface for
  // it, of total genus at most the whole's.
  std::vector<int> pieceOf(links_[s.whole].components, -1);
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

void CobordismGraph::deriveFrom(DerivationId rid) {
  // Copy: insert() may reallocate derivations_ and links_ entries' vectors.
  const Derivation r = derivations_[rid];
  const GraphLink &n = links_[r.link];
  const std::vector<RelationId> cobordismIds = n.cobordisms;
  const std::vector<RelationId> splitIds = n.splits;
  const std::vector<RelationId> sumIds = n.sums;

  // As a summand: combine with every current entry of the other summands.
  for (RelationId mid : sumIds) {
    const Sum s = sums_[mid];
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      if (s.pieces[k] != r.link) continue;
      std::vector<std::vector<DerivationId>> options(s.pieces.size());
      bool anyEmpty = false;
      for (size_t j = 0; j < s.pieces.size(); ++j) {
        if (j == k) {
          options[j] = {rid};
          continue;
        }
        for (const PartitionGenus &pe : links_[s.pieces[j]].partitionGenera.entries())
          options[j].push_back(pe.derivation);
        if (options[j].empty()) anyEmpty = true;
      }
      if (anyEmpty) continue;
      std::vector<std::vector<DerivationId>> combos;
      std::vector<DerivationId> pick(s.pieces.size());
      std::function<void(size_t)> rec = [&](size_t j) {
        if (j == s.pieces.size()) {
          combos.push_back(pick);
          return;
        }
        for (DerivationId o : options[j]) {
          pick[j] = o;
          rec(j + 1);
        }
      };
      rec(0);
      for (const auto &ids : combos) {
        std::vector<const Derivation *> recs;
        for (DerivationId id : ids) recs.push_back(&derivations_[id]);
        auto d = combineSum(s, recs);
        insert(s.whole, d.first, d.second, DerivationKind::sumCombine, mid, ids, "");
      }
    }
  }

  for (RelationId eid : cobordismIds) {
    const LinkCobordism e = cobordisms_[eid];
    if (e.out == r.link)
      if (auto d = throughCobordism(e, /*forward=*/true, r.partition, r.genus))
        insert(e.in, d->first, d->second, DerivationKind::cobordismForward, eid,
               {rid}, "");
    if (e.in == r.link)
      if (auto d = throughCobordism(e, /*forward=*/false, r.partition, r.genus))
        insert(e.out, d->first, d->second, DerivationKind::cobordismReverse, eid,
               {rid}, "");
  }

  for (RelationId sid : splitIds) {
    const Split s = splits_[sid];
    if (s.whole == r.link) {
      for (size_t k = 0; k < s.pieces.size(); ++k)
        if (auto p = restrictSplit(s, k, r.partition))
          insert(s.pieces[k], *p, r.genus, DerivationKind::splitRestrict, sid,
                 {rid}, "");
    }
    // As a piece (possibly several times, if the same link appears twice):
    // combine with every current entry of the other pieces.
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      if (s.pieces[k] != r.link)
        continue;
      std::vector<std::vector<const Derivation *>> options(s.pieces.size());
      bool anyEmpty = false;
      for (size_t j = 0; j < s.pieces.size(); ++j) {
        if (j == k) {
          options[j] = {&derivations_[rid]};
          continue;
        }
        for (const PartitionGenus &pe : links_[s.pieces[j]].partitionGenera.entries())
          options[j].push_back(&derivations_[pe.derivation]);
        if (options[j].empty())
          anyEmpty = true;
      }
      if (anyEmpty)
        continue;
      // Collect every combination first: insert() can reallocate derivations_.
      std::vector<std::vector<DerivationId>> combos;
      std::vector<const Derivation *> pick(s.pieces.size());
      std::function<void(size_t)> rec = [&](size_t j) {
        if (j == s.pieces.size()) {
          std::vector<DerivationId> ids;
          for (const Derivation *p : pick)
            ids.push_back(p->id);
          combos.push_back(std::move(ids));
          return;
        }
        for (const Derivation *o : options[j]) {
          pick[j] = o;
          rec(j + 1);
        }
      };
      rec(0);
      for (const auto &ids : combos) {
        std::vector<const Derivation *> recs;
        for (DerivationId id : ids)
          recs.push_back(&derivations_[id]);
        auto d = combineSplit(s, recs);
        insert(s.whole, d->first, d->second, DerivationKind::splitCombine, sid,
               ids, "");
      }
    }
  }
}

long CobordismGraph::propagate() {
  const size_t before = derivations_.size();
  // Derivations are processed in creation order; each one only ever
  // creates strictly improving derivations (insert()), and genus is bounded
  // below by 0 over finitely many partitions, so this terminates.
  while (!pending_.empty()) {
    std::vector<DerivationId> batch;
    batch.swap(pending_);
    std::sort(batch.begin(), batch.end());
    batch.erase(std::unique(batch.begin(), batch.end()), batch.end());
    for (DerivationId r : batch)
      deriveFrom(r);
  }
  return static_cast<long>(derivations_.size() - before);
}

// ---------------------------------------------------------------- lower bounds

int CobordismGraph::lower(LinkId n, const Partition &q) const {
  const GraphLink &link = links_.at(n);
  if (link.linking && !linkingAllows(q, *link.linking))
    return kNoSurface;
  int v = std::max(0, link.genusLowerBound.value_or(0));
  // Every surface refining q refines each coarser partition, so a bound
  // stored for any partition that q refines applies to q.
  for (const auto &[labels, val] : lower_.at(n))
    if (val > v && q.refines(Partition::fromLabels(labels)))
      v = val;
  return v;
}

void CobordismGraph::clearLowerBounds() {
  for (GraphLink &link : links_) {
    link.genusLowerBound.reset();
    link.lowerBoundSource.clear();
  }
  for (auto &m : lower_) m.clear();
  for (auto &m : lowerReason_) m.clear();
  ++lowerVersion_;
}

bool CobordismGraph::raiseLower(LinkId n, const Partition &q, int value,
                            const LowerReason &why) {
  if (links_.at(n).components > kMaxLowerComponents)
    return false;
  if (value <= lower(n, q))
    return false;
  lower_[n][q.labels()] = std::min(value, kNoSurface);
  lowerReason_[n][q.labels()] = why;
  ++lowerVersion_;
  return true;
}

CobordismGraph::LowerFact CobordismGraph::lowerWhy(LinkId n, const Partition &q) const {
  // The same choice lower() makes, with its reason.
  const GraphLink &link = links_.at(n);
  LowerFact f;
  f.storedFor = q;
  if (link.linking && !linkingAllows(q, *link.linking)) {
    f.value = kNoSurface;
    f.reason.kind = LowerReason::Kind::linking;
    return f;
  }
  if (link.genusLowerBound && *link.genusLowerBound > 0) {
    f.value = *link.genusLowerBound;
    f.storedFor = Partition::coarsest(link.components);
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

int CobordismGraph::transportedLower(const LinkCobordism &e, bool toIsIn,
                                 const Partition &p, Transport *detail) const {
  // A surface F for the `to` end with partition p, capped onto e, gives a
  // surface for the other end whose partition we compute, of genus at most
  // genus(F) + (what glue() adds with a genus-0 cap). So genus(F) is at
  // least the other end's bound for that partition minus that addition.
  const LinkId to = toIsIn ? e.in : e.out, other = toIsIn ? e.out : e.in;
  const std::vector<int> &toMap = toIsIn ? e.inMap : e.outMap;
  const std::vector<int> &otherMap = toIsIn ? e.outMap : e.inMap;
  if (links_[to].linking && !linkingAllows(p, *links_[to].linking))
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

std::string CobordismGraph::partitionGeneraFields(LinkId n) const {
  const GraphLink &link = links_.at(n);
  std::ostringstream o;
  o << "\"components\":" << link.components;
  if (link.linking) o << ",\"linking\":" << json::matrix(*link.linking);
  if (link.genusLowerBound)
    o << ",\"genus_lower\":" << *link.genusLowerBound;
  std::vector<PartitionGenus> es = link.partitionGenera.entries();
  std::sort(es.begin(), es.end(), [](const PartitionGenus &a, const PartitionGenus &b) {
    return a.partition == b.partition ? a.genus < b.genus : a.partition < b.partition;
  });
  o << ",\"entries\":[";
  for (size_t i = 0; i < es.size(); ++i)
    o << (i ? "," : "") << "{\"p\":\"" << es[i].partition.str() << "\",\"g\":" << es[i].genus
      << ",\"r\":" << es[i].derivation << '}';
  o << ']';
  if (link.components <= kMaxLowerComponents) {
    std::vector<Partition> all = allPartitions(link.components);
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

int CobordismGraph::lowerAcross(const LinkCobordism &e, bool toIsIn,
                            const Partition &q) const {
  // Every surface refining q is bounded by the minimum of transportedLower
  // over the refinements of q, and that minimum is at q itself
  // (lem:transport-monotone): refining the cap adds vertices to the gluing
  // graph, so add = E - V + #components can only fall, the other end's
  // partition only refines, and lower() is monotone under refinement.
  return transportedLower(e, toIsIn, q);
}

std::optional<int> CobordismGraph::lowerIf(const std::vector<LowerSeed> &seeds,
                                       LinkId target,
                                       const Partition &goal) const {
  CobordismGraph what = *this;
  const size_t before = what.contradictions().size();
  LowerReason seed;
  seed.kind = LowerReason::Kind::seed;
  for (const LowerSeed &s : seeds)
    what.raiseLower(s.link, s.partition, s.value, seed);
  what.propagateLower();
  if (what.contradictions().size() != before)
    return std::nullopt;
  return what.lower(target, goal);
}

long CobordismGraph::propagateLower() {
  long improved = 0;
  for (size_t n = 0; n < links_.size(); ++n)
    if (links_[n].genusLowerBound) {
      LowerReason why;
      why.kind = LowerReason::Kind::literature;
      improved += raiseLower(static_cast<LinkId>(n),
                             Partition::coarsest(links_[n].components),
                             *links_[n].genusLowerBound, why);
    }
  const int maxPasses = 200;
  bool changed = true;
  int pass = 0;
  for (; changed && pass < maxPasses; ++pass) {
    changed = false;
    for (const LinkCobordism &e : cobordisms_)
      for (bool toIsIn : {true, false}) {
        const LinkId to = toIsIn ? e.in : e.out;
        if (links_[to].components > kMaxLowerComponents)
          continue;
        for (const Partition &q : allPartitions(links_[to].components)) {
          Transport t;
          const int v = transportedLower(e, toIsIn, q, &t);
          if (v <= lower(to, q))
            continue;
          LowerReason why;
          why.kind = LowerReason::Kind::cobordism;
          why.relation = e.id;
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
    for (const Split &s : splits_) {
      const GraphLink &w = links_[s.whole];
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
        const GraphLink &pk = links_[s.pieces[k]];
        if (pk.components > kMaxLowerComponents ||
            w.components > kMaxLowerComponents)
          continue;
        std::vector<std::vector<const PartitionGenus *>> others(s.pieces.size());
        bool anyEmpty = false;
        for (size_t j = 0; j < s.pieces.size(); ++j) {
          if (j == k) continue;
          for (const PartitionGenus &pe : links_[s.pieces[j]].partitionGenera.entries())
            others[j].push_back(&pe);
          if (others[j].empty()) anyEmpty = true;
        }
        if (anyEmpty) continue;
        for (const Partition &q : allPartitions(pk.components)) {
          int bestV = 0;
          LowerReason why;
          why.kind = LowerReason::Kind::splitPiece;
          std::vector<const PartitionGenus *> pick(s.pieces.size(), nullptr);
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
                why.derivations.clear();
                for (size_t t = 0; t < s.pieces.size(); ++t)
                  if (t != k) why.derivations.push_back(pick[t]->derivation);
              }
              return;
            }
            if (j == k) {
              rec(j + 1);
              return;
            }
            for (const PartitionGenus *pe : others[j]) {
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
    for (const Sum &s : sums_) {
      const GraphLink &w = links_[s.whole];
      if (w.components > kMaxLowerComponents) continue;
      for (size_t i = 0; i < s.pieces.size(); ++i) {
        const int lo = lower(s.pieces[i], Partition::coarsest(links_[s.pieces[i]].components));
        if (lo <= 0 || lo >= kNoSurface) continue;
        int subtract = 0;
        bool known = true;
        LowerReason why;
        why.kind = LowerReason::Kind::sumPiece;
        why.relation = s.id;
        why.piece = static_cast<int>(i);
        why.from = lo;
        for (size_t j = 0; j < s.pieces.size(); ++j) {
          if (j == i) continue;
          auto b = bestConnected(s.pieces[j]);
          if (!b) { known = false; break; }
          subtract += b->genus + links_[s.pieces[j]].components - 1;
          why.derivations.push_back(b->derivation);
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
  for (const GraphLink &n : links_)
    for (const PartitionGenus &e : n.partitionGenera.entries()) {
      const int lo = lower(n.id, e.partition);
      if (e.genus < lo)
        contradictions_.push_back(
            n.label + ": record " + std::to_string(e.derivation) + " gives " +
            e.partition.str() + " genus " + std::to_string(e.genus) +
            " below the propagated lower bound " +
            (lo >= kNoSurface ? std::string("(no such surface)") : std::to_string(lo)));
    }
  return improved;
}

long CobordismGraph::saturate() {
  for (const Derivation &r : derivations_)
    pending_.push_back(r.id);
  return propagate();
}

std::optional<PartitionGenus> CobordismGraph::best(LinkId n,
                                             const Partition &target) const {
  return links_.at(n).partitionGenera.best(target);
}

std::optional<PartitionGenus> CobordismGraph::bestConnected(LinkId n) const {
  return best(n, Partition::coarsest(links_.at(n).components));
}

std::vector<DerivationId> CobordismGraph::proof(DerivationId r) const {
  std::set<DerivationId> seen;
  std::vector<DerivationId> stack{r};
  while (!stack.empty()) {
    DerivationId x = stack.back();
    stack.pop_back();
    if (!seen.insert(x).second)
      continue;
    for (DerivationId c : derivations_.at(x).children)
      stack.push_back(c);
  }
  // Children always have smaller ids, so ascending id order is topological.
  return {seen.begin(), seen.end()};
}

std::string CobordismGraph::recheck(DerivationId rid) const {
  const Derivation &r = derivations_.at(rid);
  std::optional<std::pair<Partition, int>> d;
  switch (r.kind) {
  case DerivationKind::leaf:
    return r.children.empty() ? "" : "a leaf with children";
  case DerivationKind::cobordismForward:
  case DerivationKind::cobordismReverse: {
    if (r.children.size() != 1)
      return "a cobordism derivation needs exactly one child";
    const LinkCobordism &e = cobordisms_.at(r.relation);
    const Derivation &c = derivations_.at(r.children[0]);
    const bool fwd = r.kind == DerivationKind::cobordismForward;
    if (c.link != (fwd ? e.out : e.in) || r.link != (fwd ? e.in : e.out))
      return "endpoints do not match the cobordism";
    d = throughCobordism(e, fwd, c.partition, c.genus);
    break;
  }
  case DerivationKind::splitCombine: {
    const Split &s = splits_.at(r.relation);
    if (r.children.size() != s.pieces.size() || r.link != s.whole)
      return "split-combine children/endpoint mismatch";
    std::vector<const Derivation *> recs;
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      recs.push_back(&derivations_.at(r.children[k]));
      if (recs.back()->link != s.pieces[k])
        return "split-combine child on the wrong piece";
    }
    d = combineSplit(s, recs);
    break;
  }
  case DerivationKind::sumCombine: {
    const Sum &s = sums_.at(r.relation);
    if (r.children.size() != s.pieces.size() || r.link != s.whole)
      return "sum-combine children/endpoint mismatch";
    std::vector<const Derivation *> recs;
    for (size_t k = 0; k < s.pieces.size(); ++k) {
      recs.push_back(&derivations_.at(r.children[k]));
      if (recs.back()->link != s.pieces[k])
        return "sum-combine child on the wrong piece";
    }
    d = combineSum(s, recs);
    break;
  }
  case DerivationKind::splitRestrict: {
    const Split &s = splits_.at(r.relation);
    if (r.children.size() != 1)
      return "split-restrict needs one child";
    const Derivation &c = derivations_.at(r.children[0]);
    if (c.link != s.whole)
      return "split-restrict child is not the whole";
    for (size_t k = 0; k < s.pieces.size(); ++k)
      if (s.pieces[k] == r.link)
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
