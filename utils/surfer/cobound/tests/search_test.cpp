// hoprunner_test.cpp
//
// An in-process hop (hoprunner.h) against the route a child hop takes. On
// exhaustive searches at face cap 3 (the canaries' shape, so the counts are
// a function of the code alone):
//
//   1. The accounting is verifyslicegenus's: 3_1 accepts 1,752 surfaces,
//      L2a1{0} 945 with 150 turned away as the other orientation
//      (tools/orchestrate/canaries.expected).
//   2. Every kept surface gives the same edge read in process as read back
//      from its pair signature (computed later from its faces) -- the route
//      a certificate's checker replays: the same far-side curves on the same
//      surface components, the same genus, and far sides that split into the
//      same pieces.
//   3. Keys are distinct, and a stop request ends the search.

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include <link/link.h>

#include "linknaming/diagrams/diagramiso.h"
#include "cobound/bounds/searchcobordisms.h"
#include "cobound/cobordisms/pairsigner.h"
#include "cobound/search/search.h"
#include "cobound/bounds/nodes.h"
#include "linknaming/tests/check.h"

#ifndef CASCADE_TEST_DATA
#error "CASCADE_TEST_DATA must point at cascade/tests/data"
#endif

using exactnaming::GaussDiagram;
using namespace cascade;

namespace {

// A row PD as knotbuilder reads it: labels from 0 or 1.
regina::Link linkFromRowPD(const std::string &pd) {
  std::vector<std::array<int, 4>> xs;
  std::vector<int> nums;
  std::string cur;
  for (char c : pd) {
    if (std::isdigit(static_cast<unsigned char>(c))) cur += c;
    else if (!cur.empty()) { nums.push_back(std::stoi(cur)); cur.clear(); }
  }
  if (!cur.empty()) nums.push_back(std::stoi(cur));
  const int shift = *std::min_element(nums.begin(), nums.end()) == 0 ? 1 : 0;
  for (size_t i = 0; i + 3 < nums.size(); i += 4)
    xs.push_back({nums[i] + shift, nums[i + 1] + shift, nums[i + 2] + shift,
                  nums[i + 3] + shift});
  return regina::Link::fromPD(xs.begin(), xs.end());
}

HopRow makeRow(NodeRegistry &reg, const std::string &pd) {
  const regina::Link l = linkFromRowPD(pd);
  std::vector<size_t> origin(l.countComponents());
  for (size_t i = 0; i < origin.size(); ++i) origin[i] = i;
  GaussDiagram d = GaussDiagram::of(l, origin);
  NodeMatch nm = reg.intern(simplifyKeepingComponents(d), "row");
  HopRow row;
  row.node = nm.node;
  row.diagram = d;
  row.nodeMap = nm.componentMap;
  row.pd = pd;
  row.layers = 2;
  return row;
}

// What an edge says, independent of how its surface was read: per surface
// component, the row's components and the far node's components on it
// (surface components are unlabelled, so sorted), and the far side's piece
// nodes -- with the row's components relabelled by `rowPerm` and the far
// node's by `farPerm`.
//
// The pair-signature route carries the surface back by an isomorphism
// sending its incoming curve onto L x {0}, and T has automorphisms preserving
// L (3_1's has six), so its far-side curves may sit elsewhere in T; and a
// symmetric far side (the Hopf link) may be matched to its node by either of
// two component maps. Both reads are then true, and differ by a symmetry of
// the row and one of the far side: so they are compared up to those.
std::string content(const ProofGraph &g, const HopEdge &e, const std::vector<int> &rowPerm,
                    const std::vector<int> &farPerm) {
  if (!e.ok) return "not ok: " + e.why;
  if (e.direct) return "direct";
  const WitnessEdge &w = g.witness(e.edge);
  std::map<int, std::pair<std::vector<int>, std::vector<int>>> bySurface;
  for (size_t c = 0; c < w.shape.inComponent.size(); ++c)
    bySurface[w.shape.inComponent[c]].first.push_back(rowPerm[c]);
  for (size_t j = 0; j < w.outMap.size(); ++j)
    bySurface[w.shape.outComponent[j]].second.push_back(farPerm[w.outMap[j]]);
  std::vector<std::string> parts;
  for (auto &[s, io] : bySurface) {
    std::sort(io.first.begin(), io.first.end());
    std::sort(io.second.begin(), io.second.end());
    std::string p;
    for (int c : io.first) p += std::to_string(c) + '.';
    p += '/';
    for (int c : io.second) p += std::to_string(c) + '.';
    parts.push_back(p);
  }
  std::sort(parts.begin(), parts.end());
  std::vector<NodeId> pieces;
  for (const NodeMatch &m : e.pieces) pieces.push_back(m.node);
  std::sort(pieces.begin(), pieces.end());
  std::string out = "genus " + std::to_string(w.shape.genus) + ", far node " +
                    (e.pieces.size() == 1 ? std::to_string(e.farNode) : "split") +
                    ", pieces";
  for (NodeId p : pieces) out += ' ' + std::to_string(p);
  out += ", surface components";
  for (const std::string &p : parts) out += ' ' + p;
  return out;
}

// Every permutation of `d`'s components that is a symmetry of the diagram
// (up to mirror image and global reversal, the registry's equivalence).
std::vector<std::vector<int>> symmetries(const GaussDiagram &d) {
  std::vector<int> p(d.components());
  for (size_t i = 0; i < p.size(); ++i) p[i] = static_cast<int>(i);
  std::vector<std::vector<int>> out;
  do {
    if (findDiagramIsomorphism(d, d, true, true, &p)) out.push_back(p);
  } while (std::next_permutation(p.begin(), p.end()));
  return out;
}

// Whether the two reads of one surface say the same, up to symmetries.
bool sameEdge(const ProofGraph &g, const NodeRegistry &reg, NodeId rowNode, const HopEdge &a,
              const HopEdge &b) {
  const int rowN = g.node(rowNode).components;
  std::vector<int> rowId(rowN);
  for (int i = 0; i < rowN; ++i) rowId[i] = i;
  if (!a.ok || !b.ok || a.direct || b.direct)
    return content(g, a, rowId, {}) == content(g, b, rowId, {});
  const int farN = g.node(a.farNode).components;
  std::vector<int> farId(farN);
  for (int i = 0; i < farN; ++i) farId[i] = i;
  const std::string want = content(g, a, rowId, farId);
  std::vector<std::vector<int>> farSyms = {farId};
  if (a.pieces.size() == 1 && reg.known(a.farNode))
    farSyms = symmetries(reg.info(a.farNode).diagram);
  for (const auto &rp : symmetries(reg.info(rowNode).diagram))
    for (const auto &fp : farSyms)
      if (content(g, b, rp, fp) == want) return true;
  return false;
}

// The canaries' shape: exhaustive to 3 added faces, one pass, no budget.
HopShape capThree() {
  HopShape s;
  s.maxFaces = 3;
  s.iddfsIterations = 0;
  s.iddfsStart = 0;
  s.iddfsStep = 0;
  s.rootBudgetStart = 0;
  s.resolveUnlinked = false;
  return s;
}

// Every kept surface's request and its signature, signed on its own in the
// searched thickening; testBatchSigning() signs them again all at once.
std::vector<std::pair<SignRequest, std::string>> signed_;

void checkRow(const farside::SignatureTable &sigs, const std::string &name,
              const std::string &pd, long long accepted, long long otherOrientation) {
  // One graph and registry for both reads, so equal far sides are one node.
  ProofGraph g;
  NodeRegistry reg(g);
  const HopRow row = makeRow(reg, pd);
  HopAssembler inProcess(g, reg, row);
  HopAssembler byPairSig(g, reg, row);

  HopSearcher searcher(sigs, nullptr, capThree(), 4);
  HopRun run = searcher.run(inProcess.redrawer(), name, 1'000'000'000LL, 600);
  CHECK_EQ(run.accepted, accepted, name + ": accepted, as the canaries pin it");
  CHECK_EQ(run.outcome, std::string("exhausted"), name + ": exhausted at cap 3");
  CHECK_EQ(run.accountingFailure, std::string(), name + ": the accounting balances");
  CHECK(run.accounting.find("other-orientation " + std::to_string(otherOrientation) +
                            ",") != std::string::npos,
        name + ": other orientations turned away, as the canaries pin it (" +
            run.accounting + ")");
  CHECK(!run.kept.empty(), name + ": surfaces kept");

  std::set<std::string> keys;
  int same = 0, compared = 0, fromFaces = 0;
  const farside::WitnessRedrawer &redraw = inProcess.redrawer();
  const std::vector<int> seed = redraw.rowBuild().seedFaces;
  for (size_t i = 0; i < run.kept.size(); ++i) {
    const KeptSurface &k = run.kept[i];
    keys.insert(k.key);
    const std::string key = name + "#" + std::to_string(i);
    const HopEdge a = inProcess.addRead(k.link, k.genus, key);
    const std::string sig = pairSigOf(redraw.thickening(), k.faces);
    signed_.push_back({{pd, 2, k.faces}, sig});
    const HopEdge b = byPairSig.add({sig, k.genus, key});
    ++compared;
    if (a.ok && sameEdge(g, reg, row.node, a, b)) {
      ++same;
    } else {
      std::cout << "  " << name << " kept surface " << i << " (" << k.farName << "): in process ok="
                << a.ok << " '" << a.why << "', by pair signature ok=" << b.ok << " '" << b.why
                << "'; they differ beyond the row's and far side's symmetries\n";
    }
    // A certificate's route: the faces, rebuilt in the same thickening,
    // read exactly as the search read them -- same curves, same components.
    std::string why;
    auto link = redraw.outgoingLinkFromFaces(k.faces, why);
    if (link) {
      const HopEdge c = inProcess.addRead(*link, k.genus, key);
      if (c.ok && c.farCurveEdges == a.farCurveEdges && c.shape.outComponent == a.shape.outComponent &&
          c.shape.inComponent == a.shape.inComponent && c.farNode == a.farNode)
        ++fromFaces;
    } else {
      std::cout << "  " << name << " kept surface " << i << ": from faces: " << why << "\n";
    }
  }
  CHECK_EQ(static_cast<size_t>(keys.size()), run.kept.size(),
           name + ": one surface per key");
  CHECK_EQ(same, compared,
           name + ": every kept surface reads the same in process and from its "
                  "pair signature");
  CHECK_EQ(fromFaces, compared,
           name + ": every kept surface reads exactly the same rebuilt from its faces");

  // Faces that are not a surface of this row are refused, not misread.
  const std::vector<int> &faces = run.kept.front().faces;
  std::string why;
  std::vector<int> noSeed;
  for (int f : faces)
    if (f != seed.front()) noSeed.push_back(f);
  CHECK(!redraw.outgoingLinkFromFaces(noSeed, why),
        name + ": a seed triangle missing is refused (" + why + ")");
  // One triangle outside the seed dropped: its neighbours' edges become
  // boundary inside the thickening, so the rest is not properly embedded.
  auto outside = std::find_if(faces.rbegin(), faces.rend(), [&](int f) {
    return std::find(seed.begin(), seed.end(), f) == seed.end();
  });
  if (outside != faces.rend()) {
    std::vector<int> cut;
    for (int f : faces)
      if (f != *outside) cut.push_back(f);
    CHECK(!redraw.outgoingLinkFromFaces(cut, why),
          name + ": a surface with a triangle missing is refused (" + why + ")");
  }
}

// The build digest names the thickening: the same for two builds of one row,
// different for another row or another number of layers.
void testBuildChecksum() {
  const char *trefoil = "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]";
  const farside::WitnessRedrawer a(trefoil, 2), b(trefoil, 2);
  const farside::WitnessRedrawer hopf("PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]", 2);
  const farside::WitnessRedrawer oneLayer(trefoil, 1);
  CHECK_EQ(a.buildChecksum(), b.buildChecksum(), "digest: two builds of one row agree");
  CHECK(a.buildChecksum() != hopf.buildChecksum(), "digest: another row differs");
  CHECK(a.buildChecksum() != oneLayer.buildChecksum(), "digest: another layer count differs");
  CHECK_EQ(a.buildChecksum().size(), static_cast<size_t>(16), "digest: 16 hex digits");
}

// pairSigsOf(): the rows rebuilt from their PD codes and signed with one
// ambient context each, several rows at once, give exactly the signatures
// of the searched thickenings, in request order.
void testBatchSigning() {
  std::vector<SignRequest> requests;
  for (const auto &[r, sig] : signed_) requests.push_back(r);
  std::reverse(requests.begin(), requests.end()); // rows interleaved differently
  const std::vector<std::string> sigs = pairSigsOf(requests, 3);
  CHECK_EQ(sigs.size(), requests.size(), "batch: one signature per request");
  int same = 0;
  for (size_t i = 0; i < requests.size(); ++i)
    if (sigs[i] == signed_[requests.size() - 1 - i].second) ++same;
  CHECK_EQ(same, static_cast<int>(requests.size()),
           "batch: every signature equals the one taken in the searched thickening");
  CHECK(requests.size() >= 4, "batch: surfaces from both rows");
}

void testStop(const farside::SignatureTable &sigs) {
  ProofGraph g;
  NodeRegistry reg(g);
  HopAssembler hop(g, reg, makeRow(reg, "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"));
  HopSearcher searcher(sigs, nullptr, capThree(), 4);
  int asked = 0;
  HopRun run = searcher.run(hop.redrawer(), "3_1", 1'000'000'000LL, 600,
                            [&](const KeptSurface &) { return ++asked == 1; });
  CHECK_EQ(run.outcome, std::string("stopped"), "stop: the outcome says so");
  CHECK(run.accepted < 1752, "stop: the search ended early");
  CHECK_EQ(run.accountingFailure, std::string(),
           "stop: what was described is still accounted for");
}

} // namespace

int main() {
  const std::string data = CASCADE_TEST_DATA;
  const farside::SignatureTable sigs = farside::SignatureTable::fromTables(
      data + "/knots_to_6.csv", data + "/links_to_6.csv");
  checkRow(sigs, "3_1", "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", 1752, 0);
  checkRow(sigs, "L2a1{0}", "PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]", 945, 150);
  testBuildChecksum();
  testBatchSigning();
  testStop(sigs);
  return cascadetest::finish("hoprunner_test");
}
