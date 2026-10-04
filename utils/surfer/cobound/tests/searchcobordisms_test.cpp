// hopedges_test.cpp
//
// End to end on real witnesses (tests/data, from the 2026-09-28 Phase 0
// hops): a hop row is certified, every witness becomes edges, far sides
// become nodes with the right component maps, and the proof graph derives
// what the witnesses prove. README.md, "Composing hops".

#include <fstream>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

#include <link/link.h>

#include "cobound/bounds/searchcobordisms.h"
#include "cobound/bounds/links.h"
#include "linknaming/tables.h"
#include "linknaming/tests/check.h"

#ifndef COBOUND_TEST_DATA
#error "CASCADE_TEST_DATA must point at cascade/tests/data"
#endif

using linknaming::GaussDiagram;
using namespace bounds;

namespace {

struct StoredLine {
  std::string subject, other, genus, otherComponents, pairsig;
};

// Minimal CSV reader for cobordisms.csv (quoted fields may hold commas).
std::vector<StoredLine> readCobordisms(const std::string &path) {
  std::ifstream in(path);
  std::string line;
  std::getline(in, line); // header
  std::vector<StoredLine> lines;
  while (std::getline(in, line)) {
    std::vector<std::string> f;
    std::string cur;
    bool q = false;
    for (char c : line) {
      if (c == '"') q = !q;
      else if (c == ',' && !q) { f.push_back(cur); cur.clear(); }
      else cur += c;
    }
    f.push_back(cur);
    // kind,subject,subject_components,other,other_candidates,other_components,
    // genus,tubed,pairsig,...
    lines.push_back({f[1], f[3], f[6], f[5], f[8]});
  }
  return lines;
}

GaussDiagram of(const regina::Link &l) {
  std::vector<size_t> origin(l.countComponents());
  std::iota(origin.begin(), origin.end(), 0);
  return GaussDiagram::of(l, origin);
}

// The row as searched: its own (unsimplified) diagram, interned via its
// simplification (simplify keeps component indices, so the map carries).
SearchedLink makeSearchedLink(LinkRegistry &reg, const std::string &pd) {
  GaussDiagram d = of(linknaming::linkFromTablePD(pd));
  LinkMatch nm = reg.intern(simplifyKeepingComponents(d), "row");
  SearchedLink searched;
  searched.link = nm.link;
  searched.diagram = d;
  searched.linkMap = nm.componentMap;
  searched.pd = pd;
  searched.layers = 2;
  return searched;
}

const char *PD_10_3 =
    "[[2;10;3;9];[4;17;5;18];[6;15;7;16];[8;4;9;3];[10;2;11;1];[12;20;13;19];"
    "[14;7;15;8];[16;5;17;6];[18;14;19;13];[20;12;1;11]]";
const char *PD_L11n33 =
    "[[11;16;12;17];[17;10;0;11];[9;15;10;14];[8;18;9;21];[6;19;7;20];"
    "[5;12;6;13];[13;4;14;5];[20;3;21;4];[2;8;3;7];[18;2;19;1];[15;1;16;0]]";

void test10_3() {
  ProofGraph g;
  LinkRegistry reg(g);
  SearchedLink searched = makeSearchedLink(reg, PD_10_3);
  CobordismAssembler assembler(g, reg, searched); // certifies the row (throws if not)
  CHECK(true, "10_3: the row is certified");
  auto ws = readCobordisms(std::string(COBOUND_TEST_DATA) + "/hop_10_3_witnesses.csv");
  CHECK_EQ(static_cast<int>(ws.size()), 9, "10_3: nine witnesses");
  int ok = 0, selfLoops = 0, splits = 0;
  for (const StoredLine &r : ws) {
    AddedCobordism e = assembler.add({r.pairsig, std::stoi(r.genus), r.other});
    CHECK(e.ok, "10_3 witness assembles: " + r.other + " " + e.why);
    if (!e.ok) continue;
    ++ok;
    CHECK_EQ(static_cast<int>(e.shape.outComponent.size()), std::stoi(r.otherComponents),
             "far-side curve count matches the witness file: " + r.other);
    if (r.other == "10_3") {
      CHECK_EQ(e.outgoing, searched.link, "the far side 10_3 is the row's own node");
      ++selfLoops;
    }
    if (r.other.find(" u ") != std::string::npos || r.other == "2-component unlink") {
      CHECK(e.pieces.size() >= 2, "a split far side has several pieces: " + r.other);
      ++splits;
    }
  }
  CHECK_EQ(ok, 9, "every witness assembled");
  CHECK(selfLoops >= 1, "10_3's identity witness is a self-loop");
  CHECK(splits >= 3, "split far sides recognised");
  g.propagate();
  auto best = g.bestConnected(searched.link);
  CHECK(best.has_value() && best->genus == 0, "10_3 is proved slice");
  if (best) {
    for (RecordId r : g.proof(best->record))
      CHECK(g.recheck(r).empty(), "every record of the proof rechecks");
    // The proof rests only on the unknot's disc: constructive.
    bool onlyUnknot = true;
    for (RecordId r : g.proof(best->record))
      if (g.record(r).kind == RecordKind::leaf &&
          g.record(r).source.find("unknot") == std::string::npos)
        onlyUnknot = false;
    CHECK(onlyUnknot, "the slice proof uses no literature");
  }
  CHECK(g.contradictions().empty(), "no contradictions");
}

void testL11n33() {
  ProofGraph g;
  LinkRegistry reg(g);
  SearchedLink searched = makeSearchedLink(reg, PD_L11n33);
  CobordismAssembler assembler(g, reg, searched);
  // The certificate's map: knotbuilder's component order onto the node's.
  CHECK_EQ(static_cast<int>(assembler.incomingToLink().size()), 2, "two row components");
  auto ws = readCobordisms(std::string(COBOUND_TEST_DATA) + "/hop_L11n33_witnesses.csv");
  int ok = 0;
  for (const StoredLine &r : ws) {
    AddedCobordism e = assembler.add({r.pairsig, std::stoi(r.genus), r.other});
    CHECK(e.ok, "L11n33 witness assembles: " + r.other + " " + e.why);
    if (!e.ok) continue;
    ++ok;
    if (r.other == "L11n33{1}" && r.genus == "0") {
      CHECK_EQ(e.outgoing, searched.link, "the identity far side is the row's node");
      // A self-loop of two annuli maps component i to component i.
      const LinkCobordism &we = g.cobordism(e.edge);
      if (e.shape.components == 2) {
        bool product = true;
        for (int i = 0; i < 2; ++i) {
          int ci = we.shape.inComponent[i];
          for (int j = 0; j < 2; ++j)
            if (we.outMap[j] == i && we.shape.outComponent[j] != ci) product = false;
        }
        CHECK(product, "the identity witness joins each component to itself");
      }
    }
  }
  CHECK(ok >= 1, "L11n33 witnesses assembled");
  g.propagate();
  CHECK(g.contradictions().empty(), "no contradictions");
}

void testFastMatchesReference() {
  // G1: the fast read (outgoingLinkFast, no KnottedSurface rebuild) and the
  // reference read give the same cobordism shape and the same far-side
  // pieces for every real witness.
  for (auto [pd, file] : {std::pair{PD_10_3, "hop_10_3_witnesses.csv"},
                          std::pair{PD_L11n33, "hop_L11n33_witnesses.csv"}}) {
    ProofGraph gf, gr;
    LinkRegistry rf(gf), rr(gr);
    CobordismAssembler fast(gf, rf, makeSearchedLink(rf, pd), CobordismAssembler::Read::fast);
    CobordismAssembler ref(gr, rr, makeSearchedLink(rr, pd), CobordismAssembler::Read::reference);
    auto ws = readCobordisms(std::string(COBOUND_TEST_DATA) + "/" + file);
    int compared = 0;
    for (const StoredLine &r : ws) {
      AddedCobordism a = fast.add({r.pairsig, std::stoi(r.genus), r.other});
      AddedCobordism b = ref.add({r.pairsig, std::stoi(r.genus), r.other});
      CHECK(a.ok && b.ok, std::string("both reads assemble: ") + r.other);
      if (!a.ok || !b.ok) continue;
      // The two reads may list the far-side curves in different orders.
      // Match curves by their edge sets (a bijection, or the reads
      // disagree), then compare everything under that relabelling.
      const size_t m = a.outgoingCurveEdges.size();
      CHECK_EQ(m, b.outgoingCurveEdges.size(), "same number of far curves");
      std::vector<int> perm(m, -1); // a's curve j is b's curve perm[j]
      bool bijection = m == b.outgoingCurveEdges.size();
      std::vector<char> used(b.outgoingCurveEdges.size(), 0);
      for (size_t j = 0; j < m && bijection; ++j) {
        for (size_t k = 0; k < b.outgoingCurveEdges.size(); ++k)
          if (!used[k] && a.outgoingCurveEdges[j] == b.outgoingCurveEdges[k]) {
            perm[j] = static_cast<int>(k);
            used[k] = 1;
            break;
          }
        if (perm[j] < 0) bijection = false;
      }
      CHECK(bijection, std::string("far curves match by edge set: ") + r.other);
      if (!bijection) continue;
      // Shapes: incoming curves are in node order in both; outgoing curves
      // relabelled by perm; surface components compared by first appearance.
      auto canon = [](std::vector<int> v) {
        std::map<int, int> seen;
        for (int &x : v) x = seen.try_emplace(x, static_cast<int>(seen.size())).first->second;
        return v;
      };
      std::vector<int> la = a.shape.inComponent, lb = b.shape.inComponent;
      std::vector<int> bOut(m);
      for (size_t j = 0; j < m; ++j) bOut[j] = b.shape.outComponent[perm[j]];
      la.insert(la.end(), a.shape.outComponent.begin(), a.shape.outComponent.end());
      lb.insert(lb.end(), bOut.begin(), bOut.end());
      CHECK(canon(la) == canon(lb) && a.shape.components == b.shape.components &&
                a.shape.genus == b.shape.genus,
            std::string("fast and reference shapes agree up to curve order: ") + r.other);
      auto relabelled = [&](std::vector<std::vector<size_t>> v, bool apply) {
        for (auto &x : v) {
          if (apply)
            for (size_t &o : x) o = static_cast<size_t>(perm[o]);
          std::sort(x.begin(), x.end());
        }
        std::sort(v.begin(), v.end());
        return v;
      };
      CHECK(relabelled(a.pieceOrigins, true) == relabelled(b.pieceOrigins, false),
            std::string("fast and reference pieces agree up to curve order: ") + r.other);
      CHECK_EQ(a.direct, b.direct, "direct agrees");
      ++compared;
    }
    CHECK(compared == static_cast<int>(ws.size()), "every witness compared");
  }
}

} // namespace

int main() {
  testFastMatchesReference();
  test10_3();
  testL11n33();
  return checks::finish("hopedges_test");
}
