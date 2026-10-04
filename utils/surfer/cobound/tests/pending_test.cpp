// pending_test.cpp
//
// A search's kept surfaces into the atlas's database (cobordisms/pending.h),
// on a real exhaustive search (3_1 at face cap 3, the canaries' shape):
//
//   1. kept.csv round-trips: every cobordism column and every face comes back.
//   2. The database gets exactly one line per cobordism identity, the sweep's
//      rule, each with the pair signature taken in the searched thickening,
//      other_candidates from the name table, and the search's provenance.
//   3. Signing the same surfaces again, or against a database that already
//      holds them (dedupe_against), appends nothing.
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
#include "linknaming/tables.h"
#include "linknaming/tests/check.h"
#include "surfer/report/csvwriter.h"
#include "cobound/cobordisms/database.h"
#include "cobound/solver/literature.h"

#ifndef COBOUND_TEST_DATA
#error "CASCADE_TEST_DATA must point at cascade/tests/data"
#endif

using linknaming::GaussDiagram;
using namespace bounds;
using namespace cobordisms;
using namespace search;
namespace fs = std::filesystem;

namespace {

SearchedLink makeSearchedLink(LinkRegistry &reg, const std::string &pd) {
  const regina::Link l = linknaming::linkFromTablePD(pd);
  std::vector<size_t> origin(l.countComponents());
  for (size_t i = 0; i < origin.size(); ++i) origin[i] = i;
  GaussDiagram d = GaussDiagram::of(l, origin);
  LinkMatch nm = reg.intern(simplifyKeepingComponents(d), "row");
  SearchedLink searched;
  searched.link = nm.link;
  searched.diagram = d;
  searched.linkMap = nm.componentMap;
  searched.pd = pd;
  searched.layers = 2;
  return searched;
}

RunShape capThree() {
  RunShape s;
  s.layers = 2;
  s.maxFaces = 3;
  s.iddfsIterations = 0;
  s.iddfsStart = 0;
  s.iddfsStep = 0;
  s.rootBudgetStart = 0;
  s.rootBudgetGrowth = 2;
  s.resolveUnlinked = false;
  return s;
}

bool sameCobordism(const cobordisms::Cobordism &a, const cobordisms::Cobordism &b) {
  return cobordisms::formatCobordism(a) == cobordisms::formatCobordism(b);
}

} // namespace

int main() {
  const std::string data = COBOUND_TEST_DATA;
  const std::string knots = data + "/knots_to_6.csv", links = data + "/links_to_6.csv";
  const linknaming::SignatureTable sigs = linknaming::SignatureTable::fromTables(knots, links);
  const fs::path dir = fs::temp_directory_path() / ("keptstore_test." + std::to_string(::getpid()));
  const fs::path searchDir = dir / "hop_0_n0";
  fs::create_directories(searchDir);

  const std::string pd = "[[1;5;2;4];[3;1;4;6];[5;3;6;2]]";
  CobordismGraph g;
  LinkRegistry reg(g);
  const SearchedLink searched = makeSearchedLink(reg, pd);
  CobordismAssembler assembler(g, reg, searched);
  Searcher searcher(sigs, nullptr, capThree(), 4);
  SearchResult run = searcher.run(assembler.redrawer(), "3_1", 1'000'000'000LL, 600);
  CHECK_EQ(run.accepted, 1752LL, "the canaries' count for 3_1 at cap 3");

  // What a goal run does with them (Scheduler::expand()).
  std::vector<PendingCobordism> pending;
  std::map<std::string, std::string> sigOfIdentity; // first surface of each identity
  for (const KeptSurface &ks : run.kept) {
    PendingCobordism p{ks.cobordism, searched.pd, searched.layers, ks.faces};
    p.cobordism.sourceSearch = "3_1";
    p.cobordism.thickenLayers = searched.layers;
    p.cobordism.maxFaces = 3;
    sigOfIdentity.emplace(cobordisms::cobordismIdentity(p.cobordism),
                          pairSigOf(assembler.redrawer().thickening(), ks.faces));
    pending.push_back(std::move(p));
  }
  CHECK(pending.size() >= sigOfIdentity.size(), "at least one kept surface per identity");

  // 1. kept.csv round trip (in two appends, as two searches would write it).
  const size_t half = pending.size() / 2;
  appendKept(searchDir.string(), {pending.begin(), pending.begin() + half});
  appendKept(searchDir.string(), {pending.begin() + half, pending.end()});
  std::vector<PendingCobordism> back = readKept(dir.string());
  CHECK_EQ(back.size(), pending.size(), "kept.csv: every kept surface comes back");
  int same = 0;
  for (size_t i = 0; i < back.size() && i < pending.size(); ++i)
    if (sameCobordism(back[i].cobordism, pending[i].cobordism) && back[i].faces == pending[i].faces &&
        back[i].incomingPD == pending[i].incomingPD && back[i].layers == pending[i].layers)
      ++same;
  CHECK_EQ(same, static_cast<int>(pending.size()), "kept.csv: every field and face round-trips");

  // 2. Into a database.
  solver::NameTable names;
  solver::loadNameTable(knots, names);
  solver::loadNameTable(links, names);
  const std::string database = (dir / "cobordisms.csv").string();
  // Every surface offered twice: within one batch, an identity is recorded once.
  std::vector<PendingCobordism> twice = back;
  twice.insert(twice.end(), back.begin(), back.end());
  SignResult s = signKept(twice, database, {}, names, 3);
  CHECK_EQ(s.kept, 2 * pending.size(), "store: every surface offered");
  CHECK_EQ(s.fresh, sigOfIdentity.size(), "store: one fresh surface per witness identity");
  CHECK_EQ(s.appended, sigOfIdentity.size(), "store: all of them appended");
  std::vector<cobordisms::Cobordism> inDatabase = cobordisms::loadCobordisms(database, false);
  CHECK_EQ(inDatabase.size(), sigOfIdentity.size(), "store: one line per identity");
  int sigOk = 0, provenanceOk = 0, candidatesOk = 0;
  {
    std::ifstream in(database);
    std::string line;
    std::getline(in, line);
    CHECK_EQ(line, std::string(cobordisms::COBORDISMS_HEADER), "store: the header");
    while (std::getline(in, line)) {
      cobordisms::Cobordism w;
      cobordisms::cobordismFromFields(parseCsvLine(line), w, true, false, database);
      auto it = sigOfIdentity.find(cobordisms::cobordismIdentity(w));
      if (it != sigOfIdentity.end() && it->second == w.pairSig) ++sigOk;
      if (w.subject == "3_1" && w.sourceSearch == "3_1" && w.thickenLayers == 2 && w.maxFaces == 3)
        ++provenanceOk;
      if (w.kind == cobordisms::CobordismKind::direct ||
          w.otherCandidates == names.candidates(w.other, w.otherComponents))
        ++candidatesOk;
    }
  }
  CHECK_EQ(sigOk, static_cast<int>(inDatabase.size()),
           "store: each pair signature is the one taken in the searched thickening");
  CHECK_EQ(provenanceOk, static_cast<int>(inDatabase.size()), "store: subject and provenance");
  CHECK_EQ(candidatesOk, static_cast<int>(inDatabase.size()),
           "store: other_candidates as the sweep fills them");

  // The `.rows.csv` sidecar: one line per signed cobordism, its key a
  // signed pair signature's and its row_pd the search's own PD.
  {
    std::ifstream in(database + ".rows.csv");
    std::string line;
    std::getline(in, line);
    CHECK_EQ(line, std::string("witness,layers,row_pd"), "rows sidecar: header");
    std::map<std::string, std::string> sidecar;
    while (std::getline(in, line)) {
      std::vector<std::string> f = parseCsvLine(line);
      if (f.size() == 3) sidecar[f[0]] = f[2];
    }
    CHECK_EQ(sidecar.size(), inDatabase.size(), "rows sidecar: one line per stored witness");
    int keyed = 0;
    for (const auto &w : cobordisms::loadCobordisms(database, true))
      if (sidecar.count(w.pairSigKey) && sidecar[w.pairSigKey] == pd) ++keyed;
    CHECK_EQ(keyed, static_cast<int>(inDatabase.size()),
             "rows sidecar: every stored witness keyed to its hop row's PD");
  }

  // 3. Deduplication.
  s = signKept(back, database, {}, names, 3);
  CHECK(s.fresh == 0 && s.appended == 0, "again: nothing new, nothing appended");
  const std::string other = (dir / "other.csv").string();
  s = signKept(back, other, {database}, names, 3);
  CHECK(s.fresh == 0 && s.appended == 0 && !fs::exists(other),
        "against a store holding them: nothing appended");

  // 4. A torn last line is skipped.
  std::ofstream(searchDir / "kept.csv", std::ios::app) << "cobordism,3_1,1,Unkn";
  CHECK_EQ(readKept(dir.string()).size(), pending.size(), "a torn kept line is skipped");

  fs::remove_all(dir);
  return checks::finish("keptstore_test");
}
