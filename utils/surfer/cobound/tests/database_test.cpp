//
//  witnessstore_test.cpp
//
//  The database (cobordisms/database.h), which every run and `sign`
//  append to and every merge and solve reads:
//
//   1. A cobordism written and read back is the same cobordism, a name with a
//      comma included (csvField() quotes it), resolved_vertices written
//      empty when 0.
//   2. An append creates the file with its header, keeps every earlier byte,
//      gives each cobordism its line's offset and drops its pair signature.
//   3. A torn last line (no newline) is ignored on load and truncated before
//      the next append.
//   4. A 12-column file is refused for appending.
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
using cobordisms::Cobordism;

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

Cobordism sample(const std::string &subject, int genus, int resolved) {
  Cobordism w;
  w.kind = cobordisms::CobordismKind::cobordism;
  w.subject = subject;
  w.subjectComponents = 2;
  w.other = "#{L2a1{0},L2a1{1}}"; // a name with a comma
  w.otherCandidates = {"#{L2a1{0},L2a1{1}}"};
  w.otherComponents = 3;
  w.genus = genus;
  w.tubed = true;
  w.pairSig = "-cabcdef" + std::to_string(genus);
  w.sourceSearch = subject;
  w.thickenLayers = 2;
  w.maxFaces = 5;
  w.resolvedVertices = resolved;
  return w;
}

bool same(const Cobordism &a, const Cobordism &b) {
  return a.kind == b.kind && a.subject == b.subject &&
         a.subjectComponents == b.subjectComponents && a.other == b.other &&
         a.otherCandidates == b.otherCandidates && a.otherComponents == b.otherComponents &&
         a.genus == b.genus && a.tubed == b.tubed && a.sourceSearch == b.sourceSearch &&
         a.thickenLayers == b.thickenLayers && a.maxFaces == b.maxFaces &&
         a.resolvedVertices == b.resolvedVertices;
}

} // namespace

int main() {
  std::cout << "witnessstore\n";
  const fs::path dir = fs::temp_directory_path() /
                       ("witnessstore_test." + std::to_string(::getpid()));
  fs::create_directories(dir);
  const fs::path database = dir / "cobordisms.csv";

  // 1. Round trip of a line.
  {
    const Cobordism w = sample("cascade:run/L6a5{0;1}/n3", 1, 0);
    const std::string line = cobordisms::formatCobordism(w);
    Cobordism back;
    check(cobordisms::cobordismFromFields(parseCsvLine(line), back, true, false, database) &&
              same(w, back) && back.pairSig == w.pairSig,
          "a witness round-trips through its line");
    check(line.back() == ',', "resolved_vertices 0 is written empty");
    const Cobordism r = sample("L6a5{0;1}", 0, 2);
    Cobordism backR;
    cobordisms::cobordismFromFields(parseCsvLine(cobordisms::formatCobordism(r)), backR, false,
                                    false, database);
    check(backR.resolvedVertices == 2 && backR.pairSig.empty(),
          "resolved_vertices 2 round-trips; the pair signature is dropped when not kept");
  }

  // 2. Appends.
  {
    std::vector<Cobordism> ws = {sample("a", 0, 0), sample("b", 1, 0)};
    cobordisms::appendCobordisms(database, ws, 0);
    const std::string first = slurp(database);
    check(first.rfind(std::string(cobordisms::COBORDISMS_HEADER) + "\n", 0) == 0,
          "a new store starts with the header");
    check(ws[0].pairSig.empty() && ws[0].fileOffset > 0 && ws[1].fileOffset > ws[0].fileOffset,
          "appended witnesses get their offsets and drop their pair signatures");
    std::vector<Cobordism> more = {sample("c", 2, 1)};
    cobordisms::appendCobordisms(database, more, 0);
    const std::string second = slurp(database);
    check(second.compare(0, first.size(), first) == 0 && second.size() > first.size(),
          "a second append keeps every earlier byte");
    const auto loaded = cobordisms::loadCobordisms(database, true);
    check(loaded.size() == 3 && loaded[2].subject == "c" && loaded[2].resolvedVertices == 1 &&
              !loaded[0].pairSigKey.empty(),
          "the store loads back, keys computed");
  }

  // 3. A torn last line.
  {
    std::ofstream(database, std::ios::app | std::ios::binary) << "cobordism,torn,1,Unknot";
    check(cobordisms::loadCobordisms(database, false).size() == 3,
          "a torn last line is ignored on load");
    std::vector<Cobordism> d = {sample("d", 0, 0)};
    cobordisms::appendCobordisms(database, d, 0);
    const auto loaded = cobordisms::loadCobordisms(database, false);
    check(loaded.size() == 4 && loaded.back().subject == "d" &&
              slurp(database).find("torn") == std::string::npos,
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
    std::vector<Cobordism> e = {sample("e", 0, 0)};
    bool refused = false;
    try {
      cobordisms::appendCobordisms(old, e, 0);
    } catch (const std::exception &) {
      refused = true;
    }
    check(refused, "a 12-column file is refused for appending");
  }

  // 5. cobordismIdentities(): loadCobordisms()'s identities, in byte ranges.
  {
    const fs::path big = dir / "identities.csv";
    std::vector<Cobordism> ws;
    for (int i = 0; i < 400; ++i)
      ws.push_back(sample("row" + std::to_string(i % 97), i % 5, i % 3));
    cobordisms::appendCobordisms(big, ws, 0);
    {
      std::ofstream out(big, std::ios::app | std::ios::binary);
      out << "\n"                            // an empty line
          << "cobordism,not,a,witness\n"     // a malformed line
          << cobordisms::formatCobordism(sample("late", 4, 0)) << "\n"
          << "cobordism,torn,1,Unknot";      // a torn last line
    }
    std::unordered_set<std::string> serial;
    for (const Cobordism &w : cobordisms::loadCobordisms(big, false))
      serial.insert(cobordisms::cobordismIdentity(w));
    bool allSame = true;
    for (unsigned threads : {1u, 2u, 3u, 7u, 16u})
      for (std::streamoff range : {1, 100, 4096, 64 << 20})
        allSame = allSame && cobordisms::cobordismIdentities(big, threads, range) == serial;
    check(allSame && serial.size() > 100 && serial.count(cobordisms::cobordismIdentity(
                                               sample("late", 4, 0))),
          "witnessIdentities() equals loadWitnesses()'s identities at every range split");
    check(cobordisms::cobordismIdentities(dir / "absent.csv", 4).empty(),
          "a missing file has no identities");
  }

  // 6. other_candidates splits at ';' outside a tag only (phase 3's fix).
  {
    using cobordisms::splitCandidates;
    check(splitCandidates("L8a4{0;1};L8a4{1;0}") ==
              std::vector<std::string>{"L8a4{0;1}", "L8a4{1;0}"},
          "a tag's ';' does not split a candidate");
    check(splitCandidates("3_1;m3_1") == std::vector<std::string>{"3_1", "m3_1"},
          "untagged names split at ';'");
    check(splitCandidates("") .empty(), "an empty field has no candidates");
    check(splitCandidates(";L2a1{0};;") == std::vector<std::string>{"L2a1{0}"},
          "empty pieces are dropped, as before");
    Cobordism w = sample("tagged", 1, 0);
    w.otherCandidates = {"L8a4{0;1}", "L8a4{1;0}"};
    Cobordism back;
    check(cobordisms::parseCobordismLine(cobordisms::formatCobordism(w), back, false, false,
                                         "x") &&
              back.otherCandidates == w.otherCandidates,
          "a three-component candidate list round-trips");
  }

  // 7. The other readers: readCobordisms() keeps pair signatures, the index
  //    finds lines by subject and outgoing base, PairSigReader reads one back;
  //    all leave a torn last line out.
  {
    const fs::path db = dir / "indexed.csv";
    std::vector<Cobordism> ws = {sample("K", 1, 0), sample("L", 2, 0), sample("K", 3, 0)};
    ws[1].other = "m3_1";
    ws[1].otherCandidates = {"m3_1"};
    cobordisms::appendCobordisms(db, ws, 0);
    std::ofstream(db, std::ios::app | std::ios::binary) << "cobordism,K,2,torn";
    const std::vector<Cobordism> all = cobordisms::readCobordisms(db);
    check(all.size() == 3 && all[0].pairSig == "-cabcdef1" && all[2].pairSig == "-cabcdef3",
          "readWitnesses(): every complete line, pair signatures kept, the torn one left out");
    const cobordisms::DatabaseIndex index(db.string());
    check(index.subjects() == 2 && index.has("K") && !index.has("torn"),
          "the index: two subjects, the torn line left out");
    const auto k = index.ofSubject("K");
    check(k.size() == 2 && k[0].cobordism.genus == 1 && k[1].cobordism.genus == 3 &&
              k[1].cobordism.pairSig == "-cabcdef3" && k[0].incomingPD.empty(),
          "rows(): a subject's lines, in file order, pair signatures kept");
    const auto byBase = index.byOutgoing(cobordisms::DatabaseIndex::base("3_1"), 10);
    check(byBase.size() == 1 && byBase[0].cobordism.subject == "L",
          "byOutgoing(): the outgoing base, mirror mark dropped");
    cobordisms::PairSigReader reader;
    reader.setPath(db);
    check(reader.at(k[1].cobordism.fileOffset) == "-cabcdef3" && reader.at(-1).empty(),
          "PairSigReader reads a line's pair signature back by offset");
  }

  fs::remove_all(dir);
  std::cout << (failures ? "witnessstore_test: FAILED\n" : "witnessstore_test: all passed\n");
  return failures ? 1 : 0;
}
