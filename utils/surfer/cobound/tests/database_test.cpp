//
//  witnessstore_test.cpp
//
//  The witness file (witnessstore.h), which verifyslicegenus and
//  cascadesearch both append to and every merge and solve reads:
//
//   1. A witness written and read back is the same witness, a name with a
//      comma included (csvField() quotes it), resolved_vertices written
//      empty when 0.
//   2. An append creates the file with its header, keeps every earlier byte,
//      gives each witness its line's offset and drops its pair signature.
//   3. A torn last line (no newline) is ignored on load and truncated before
//      the next append.
//   4. A 12-column file is refused for appending, and --rewrite-witnesses
//      migrates it, each line gaining the empty 13th field.
//

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "surfer/report/csvwriter.h"
#include "cobound/cobordisms/database.h"

namespace fs = std::filesystem;
using cobordismgraph::Witness;

namespace {

int failures = 0;

void check(bool ok, const std::string &what) {
  std::cout << (ok ? "  ok: " : "  FAIL: ") << what << "\n";
  if (!ok) ++failures;
}

std::string slurp(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  std::ostringstream s;
  s << in.rdbuf();
  return s.str();
}

Witness sample(const std::string &subject, int genus, int resolved) {
  Witness w;
  w.kind = cobordismgraph::WitnessKind::cobordism;
  w.subject = subject;
  w.subjectComponents = 2;
  w.other = "#{L2a1{0},L2a1{1}}"; // a name with a comma
  w.otherCandidates = {"#{L2a1{0},L2a1{1}}"};
  w.otherComponents = 3;
  w.genus = genus;
  w.tubed = true;
  w.pairSig = "-cabcdef" + std::to_string(genus);
  w.sourceRow = subject;
  w.thickenLayers = 2;
  w.maxFaces = 5;
  w.resolvedVertices = resolved;
  return w;
}

bool same(const Witness &a, const Witness &b) {
  return a.kind == b.kind && a.subject == b.subject &&
         a.subjectComponents == b.subjectComponents && a.other == b.other &&
         a.otherCandidates == b.otherCandidates && a.otherComponents == b.otherComponents &&
         a.genus == b.genus && a.tubed == b.tubed && a.sourceRow == b.sourceRow &&
         a.thickenLayers == b.thickenLayers && a.maxFaces == b.maxFaces &&
         a.resolvedVertices == b.resolvedVertices;
}

} // namespace

int main() {
  std::cout << "witnessstore\n";
  const fs::path dir = fs::temp_directory_path() /
                       ("witnessstore_test." + std::to_string(::getpid()));
  fs::create_directories(dir);
  const fs::path store = dir / "cobordisms.csv";

  // 1. Round trip of a line.
  {
    const Witness w = sample("cascade:run/L6a5{0;1}/n3", 1, 0);
    const std::string line = witnessstore::formatWitness(w);
    Witness back;
    check(witnessstore::witnessFromFields(parseCsvLine(line), back, true, false, store) &&
              same(w, back) && back.pairSig == w.pairSig,
          "a witness round-trips through its line");
    check(line.back() == ',', "resolved_vertices 0 is written empty");
    const Witness r = sample("L6a5{0;1}", 0, 2);
    Witness backR;
    witnessstore::witnessFromFields(parseCsvLine(witnessstore::formatWitness(r)), backR, false,
                                    false, store);
    check(backR.resolvedVertices == 2 && backR.pairSig.empty(),
          "resolved_vertices 2 round-trips; the pair signature is dropped when not kept");
  }

  // 2. Appends.
  {
    std::vector<Witness> ws = {sample("a", 0, 0), sample("b", 1, 0)};
    witnessstore::appendWitnesses(store, ws, 0);
    const std::string first = slurp(store);
    check(first.rfind(std::string(witnessstore::COBORDISMS_HEADER) + "\n", 0) == 0,
          "a new store starts with the header");
    check(ws[0].pairSig.empty() && ws[0].fileOffset > 0 && ws[1].fileOffset > ws[0].fileOffset,
          "appended witnesses get their offsets and drop their pair signatures");
    std::vector<Witness> more = {sample("c", 2, 1)};
    witnessstore::appendWitnesses(store, more, 0);
    const std::string second = slurp(store);
    check(second.compare(0, first.size(), first) == 0 && second.size() > first.size(),
          "a second append keeps every earlier byte");
    const auto loaded = witnessstore::loadWitnesses(store, true);
    check(loaded.size() == 3 && loaded[2].subject == "c" && loaded[2].resolvedVertices == 1 &&
              !loaded[0].pairSigKey.empty(),
          "the store loads back, keys computed");
  }

  // 3. A torn last line.
  {
    std::ofstream(store, std::ios::app | std::ios::binary) << "cobordism,torn,1,Unknot";
    check(witnessstore::loadWitnesses(store, false).size() == 3,
          "a torn last line is ignored on load");
    std::vector<Witness> d = {sample("d", 0, 0)};
    witnessstore::appendWitnesses(store, d, 0);
    const auto loaded = witnessstore::loadWitnesses(store, false);
    check(loaded.size() == 4 && loaded.back().subject == "d" &&
              slurp(store).find("torn") == std::string::npos,
          "the next append truncates it first");
  }

  // 4. The 12-column file.
  {
    const fs::path old = dir / "old.csv";
    {
      std::ofstream out(old, std::ios::binary);
      out << "kind,subject,subject_components,other,other_candidates,other_components,"
             "genus,tubed,pairsig,source_row,thicken_layers,max_faces\n";
      out << "cobordism,3_1,1,Unknot,Unknot,1,1,false,-cafoo,3_1,1,4\n";
    }
    std::vector<Witness> e = {sample("e", 0, 0)};
    bool refused = false;
    try {
      witnessstore::appendWitnesses(old, e, 0);
    } catch (const std::exception &) {
      refused = true;
    }
    check(refused, "a 12-column file is refused for appending");
    const size_t n = witnessstore::rewriteWitnessFile(old);
    const std::string after = slurp(old);
    check(n == 1 && after.find("cobordism,3_1,1,Unknot,Unknot,1,1,false,-cafoo,3_1,1,4,\n") !=
                        std::string::npos,
          "--rewrite-witnesses migrates a 12-column line by one empty field");
  }

  // 5. witnessIdentities(): loadWitnesses()'s identities, in byte ranges.
  {
    const fs::path big = dir / "identities.csv";
    std::vector<Witness> ws;
    for (int i = 0; i < 400; ++i)
      ws.push_back(sample("row" + std::to_string(i % 97), i % 5, i % 3));
    witnessstore::appendWitnesses(big, ws, 0);
    {
      std::ofstream out(big, std::ios::app | std::ios::binary);
      out << "\n"                            // an empty line
          << "cobordism,not,a,witness\n"     // a malformed line
          << witnessstore::formatWitness(sample("late", 4, 0)) << "\n"
          << "cobordism,torn,1,Unknot";      // a torn last line
    }
    std::unordered_set<std::string> serial;
    for (const Witness &w : witnessstore::loadWitnesses(big, false))
      serial.insert(cobordismgraph::witnessIdentity(w));
    bool allSame = true;
    for (unsigned threads : {1u, 2u, 3u, 7u, 16u})
      for (std::streamoff range : {1, 100, 4096, 64 << 20})
        allSame = allSame && witnessstore::witnessIdentities(big, threads, range) == serial;
    check(allSame && serial.size() > 100 && serial.count(cobordismgraph::witnessIdentity(
                                               sample("late", 4, 0))),
          "witnessIdentities() equals loadWitnesses()'s identities at every range split");
    check(witnessstore::witnessIdentities(dir / "absent.csv", 4).empty(),
          "a missing file has no identities");
  }

  fs::remove_all(dir);
  std::cout << (failures ? "witnessstore_test: FAILED\n" : "witnessstore_test: all passed\n");
  return failures ? 1 : 0;
}
