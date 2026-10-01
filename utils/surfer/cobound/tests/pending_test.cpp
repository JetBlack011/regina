// keptstore_test.cpp
//
// A hop's kept surfaces into the atlas's witness store (keptstore.h), on a
// real exhaustive row (3_1 at face cap 3, the canaries' shape):
//
//   1. kept.csv round-trips: every witness column and every face comes back.
//   2. The store gets exactly one line per witness identity, the sweep's
//      rule, each with the pair signature taken in the searched thickening,
//      other_candidates from the name table, and the hop's provenance.
//   3. Storing the same surfaces again, or against a store that already
//      holds them (--dedupe-against), appends nothing.
//   4. A torn last line of kept.csv (a run killed mid-append) is skipped.

#include <array>
#include <filesystem>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <unistd.h>

#include <link/link.h>

#include "cobound/bounds/searchcobordisms.h"
#include "cobound/cobordisms/pairsigner.h"
#include "cobound/search/search.h"
#include "cobound/cobordisms/pending.h"
#include "cobound/bounds/links.h"
#include "linknaming/tests/check.h"
#include "surfer/report/csvwriter.h"
#include "cobound/cobordisms/database.h"
#include "cobound/solver/literature.h"

#ifndef CASCADE_TEST_DATA
#error "CASCADE_TEST_DATA must point at cascade/tests/data"
#endif

using exactnaming::GaussDiagram;
using namespace cascade;
namespace fs = std::filesystem;

namespace {

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

bool sameWitness(const cobordismgraph::Witness &a, const cobordismgraph::Witness &b) {
  return witnessstore::formatWitness(a) == witnessstore::formatWitness(b);
}

} // namespace

int main() {
  const std::string data = CASCADE_TEST_DATA;
  const std::string knots = data + "/knots_to_6.csv", links = data + "/links_to_6.csv";
  const farside::SignatureTable sigs = farside::SignatureTable::fromTables(knots, links);
  const fs::path dir = fs::temp_directory_path() / ("keptstore_test." + std::to_string(::getpid()));
  const fs::path hopDir = dir / "hop_0_n0";
  fs::create_directories(hopDir);

  const std::string pd = "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]";
  ProofGraph g;
  NodeRegistry reg(g);
  const HopRow row = makeRow(reg, pd);
  HopAssembler hop(g, reg, row);
  HopSearcher searcher(sigs, nullptr, capThree(), 4);
  HopRun run = searcher.run(hop.redrawer(), "3_1", 1'000'000'000LL, 600);
  CHECK_EQ(run.accepted, 1752LL, "the canaries' count for 3_1 at cap 3");

  // What cascadesearch does with them (expand()).
  std::vector<PendingWitness> pending;
  std::map<std::string, std::string> sigOfIdentity; // first surface of each identity
  for (const KeptSurface &ks : run.kept) {
    PendingWitness p{ks.witness, row.pd, row.layers, ks.faces};
    p.witness.sourceRow = "3_1";
    p.witness.thickenLayers = row.layers;
    p.witness.maxFaces = 3;
    sigOfIdentity.emplace(cobordismgraph::witnessIdentity(p.witness),
                          pairSigOf(hop.redrawer().thickening(), ks.faces));
    pending.push_back(std::move(p));
  }
  CHECK(pending.size() >= sigOfIdentity.size(), "at least one kept surface per identity");

  // 1. kept.csv round trip (in two appends, as two hops would write it).
  const size_t half = pending.size() / 2;
  appendKept(hopDir.string(), {pending.begin(), pending.begin() + half});
  appendKept(hopDir.string(), {pending.begin() + half, pending.end()});
  std::vector<PendingWitness> back = readKept(dir.string());
  CHECK_EQ(back.size(), pending.size(), "kept.csv: every kept surface comes back");
  int same = 0;
  for (size_t i = 0; i < back.size() && i < pending.size(); ++i)
    if (sameWitness(back[i].witness, pending[i].witness) && back[i].faces == pending[i].faces &&
        back[i].rowPD == pending[i].rowPD && back[i].layers == pending[i].layers)
      ++same;
  CHECK_EQ(same, static_cast<int>(pending.size()), "kept.csv: every field and face round-trips");

  // 2. Into a store.
  cobordismgraph::NameTable names;
  witnessstore::loadNameTable(knots, names);
  witnessstore::loadNameTable(links, names);
  const std::string store = (dir / "cobordisms.csv").string();
  // Every surface offered twice: within one batch, an identity is recorded once.
  std::vector<PendingWitness> twice = back;
  twice.insert(twice.end(), back.begin(), back.end());
  StoreResult s = storeKept(twice, store, {}, names, 3);
  CHECK_EQ(s.kept, 2 * pending.size(), "store: every surface offered");
  CHECK_EQ(s.fresh, sigOfIdentity.size(), "store: one fresh surface per witness identity");
  CHECK_EQ(s.appended, sigOfIdentity.size(), "store: all of them appended");
  std::vector<cobordismgraph::Witness> stored = witnessstore::loadWitnesses(store, false);
  CHECK_EQ(stored.size(), sigOfIdentity.size(), "store: one line per identity");
  int sigOk = 0, provenanceOk = 0, candidatesOk = 0;
  {
    std::ifstream in(store);
    std::string line;
    std::getline(in, line);
    CHECK_EQ(line, std::string(witnessstore::COBORDISMS_HEADER), "store: the header");
    while (std::getline(in, line)) {
      cobordismgraph::Witness w;
      witnessstore::witnessFromFields(parseCsvLine(line), w, true, false, store);
      auto it = sigOfIdentity.find(cobordismgraph::witnessIdentity(w));
      if (it != sigOfIdentity.end() && it->second == w.pairSig) ++sigOk;
      if (w.subject == "3_1" && w.sourceRow == "3_1" && w.thickenLayers == 2 && w.maxFaces == 3)
        ++provenanceOk;
      if (w.kind == cobordismgraph::WitnessKind::direct ||
          w.otherCandidates == names.candidates(w.other, w.otherComponents))
        ++candidatesOk;
    }
  }
  CHECK_EQ(sigOk, static_cast<int>(stored.size()),
           "store: each pair signature is the one taken in the searched thickening");
  CHECK_EQ(provenanceOk, static_cast<int>(stored.size()), "store: subject and provenance");
  CHECK_EQ(candidatesOk, static_cast<int>(stored.size()),
           "store: other_candidates as the sweep fills them");

  // The row sidecar: one line per stored witness, its key a stored pair
  // signature's and its row the hop's own PD.
  {
    std::ifstream in(store + ".rows.csv");
    std::string line;
    std::getline(in, line);
    CHECK_EQ(line, std::string("witness,layers,row_pd"), "rows sidecar: header");
    std::map<std::string, std::string> rows;
    while (std::getline(in, line)) {
      std::vector<std::string> f = parseCsvLine(line);
      if (f.size() == 3) rows[f[0]] = f[2];
    }
    CHECK_EQ(rows.size(), stored.size(), "rows sidecar: one line per stored witness");
    int keyed = 0;
    for (const auto &w : witnessstore::loadWitnesses(store, true))
      if (rows.count(w.pairSigKey) && rows[w.pairSigKey] == pd) ++keyed;
    CHECK_EQ(keyed, static_cast<int>(stored.size()),
             "rows sidecar: every stored witness keyed to its hop row's PD");
  }

  // 3. Deduplication.
  s = storeKept(back, store, {}, names, 3);
  CHECK(s.fresh == 0 && s.appended == 0, "again: nothing new, nothing appended");
  const std::string other = (dir / "other.csv").string();
  s = storeKept(back, other, {store}, names, 3);
  CHECK(s.fresh == 0 && s.appended == 0 && !fs::exists(other),
        "against a store holding them: nothing appended");

  // 4. A torn last line is skipped.
  std::ofstream(hopDir / "kept.csv", std::ios::app) << "cobordism,3_1,1,Unkn";
  CHECK_EQ(readKept(dir.string()).size(), pending.size(), "a torn kept line is skipped");

  fs::remove_all(dir);
  return cascadetest::finish("keptstore_test");
}
