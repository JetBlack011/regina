//
//  databasecobordisms.cpp
//

#include "cobound/bounds/databasecobordisms.h"

#include <algorithm>
#include <condition_variable>
#include <iomanip>
#include <iostream>
#include <memory>
#include <mutex>
#include <sstream>
#include <thread>

#include "cobound/bounds/searchcobordisms.h"
#include "cobound/cobordisms/cobordismkey.h"
#include "cobound/driver/timers.h"
#include "cobound/json.h"
#include "cobound/outgoing/fromdatabase.h"
#include "cobound/frozen.h"
#include "linknaming/diagrams/simplification.h"

namespace bounds {

using linknaming::GaussDiagram;
using timers::Clock;
using timers::secondsSince;

namespace {

// A database line as a goal run carries it (cobordisms/database's reader
// parses it): the outgoing name, the pair signature, genus, outgoing
// components, layers and the PD it was searched on.
struct Cobordism {
  std::string other, pairsig;
  int genus = 0, otherComponents = 0, layers = 2;
  /// The PD it was searched on, when the file's `.rows.csv` sidecar records
  /// it (a goal run's search is on its link's simplified diagram, not the
  /// table's); empty for the table's PD.
  std::string incomingPD;
};

Cobordism carried(const cobordisms::StoredCobordism &s) {
  return {s.cobordism.other, s.cobordism.pairSig, s.cobordism.genus, s.cobordism.otherComponents,
          s.cobordism.thickenLayers, s.incomingPD};
}

} // namespace

DatabaseCobordisms::DatabaseCobordisms(const std::string &database,
                                       const std::vector<std::string> &tables)
    : index_(database) {
  // name -> the table's PD string, as the atlas's searches were run on it.
  for (const std::string &file : tables)
    for (const linknaming::TableRow &row : linknaming::readTableRows(file))
      tablePD_[row.name] = row.pd;
}

bool DatabaseCobordisms::subjectsFor(LinkId n, const LinkAxioms &axioms,
                                 const linknaming::Tables &tables,
                                 std::vector<std::string> *variants) const {
  auto it = axioms.tableName.find(n);
  if (it == axioms.tableName.end()) return false;
  // Every table entry of this link's class (one oriented link up to mirror
  // and global reversal) is the same link; each is a subject of its own.
  const std::string canon = axioms.classOf(it->second);
  const linknaming::TableEntry *e = tables.entry(it->second);
  if (!e) return false;
  bool any = !index_.byOutgoing(cobordisms::DatabaseIndex::base(it->second), 1).empty();
  for (const linknaming::TableEntry *v : tables.variants(e->base))
    if (axioms.classOf(v->name) == canon && index_.has(v->name)) {
      any = true;
      if (variants) variants->push_back(v->name);
    }
  return any;
}

double DatabaseCobordisms::load(LinkId n, DatabaseLoad &ld) {
  const auto tLoad = Clock::now();
  done_.insert(n);
  std::vector<std::string> own;
  if (!subjectsFor(n, ld.axioms, ld.tables, &own)) return secondsSince(tLoad);
  CobordismGraph &g = ld.g;
  LinkRegistry &reg = ld.reg;
  const linknaming::Tables &tables = ld.tables;
  std::map<LinkId, std::string> &tableName = ld.axioms.tableName;
  double readSeconds = 0; // phase A: the assembly's waits for the readers
  const size_t linksBefore = g.linkCount();
  int assembled = 0, failed = 0, refusedSubjects = 0;
  // A cobordism can be on an upper proof only if its genus is at most the
  // goal (glue() never lowers a genus), and on a lower proof only if the
  // charge it costs, at least its genus, is affordable: at most the largest
  // special source's bound minus the lower goal. Anything above both is
  // skipped before it is read back (the genus is in the CSV line).
  const int maxGenus = ld.maxGenus;
  size_t skippedGenus = 0;
  auto keep = [&](const Cobordism &w) {
    if (w.genus <= maxGenus) return true;
    ++skippedGenus;
    return false;
  };
  // Cobordisms by subject: this link's own subjects (every cobordism), and
  // other subjects whose recorded outgoing link names this link's base (a hint: the
  // outgoing link is redrawn and named exactly like any other).
  std::map<std::string, std::vector<Cobordism>> bySubject;
  std::set<std::string> ownSubjects(own.begin(), own.end());
  for (const std::string &name : own)
    for (const cobordisms::StoredCobordism &s : index_.ofSubject(name))
      if (const Cobordism w = carried(s); keep(w)) bySubject[s.cobordism.subject].push_back(w);
  size_t reverse = 0;
  for (const cobordisms::StoredCobordism &s :
       index_.byOutgoing(cobordisms::DatabaseIndex::base(tableName[n]), 300))
    if (const Cobordism w = carried(s); !ownSubjects.count(s.cobordism.subject) && keep(w)) {
      bySubject[s.cobordism.subject].push_back(w);
      ++reverse;
    }
  // Two phases. First every subject's cobordisms are read back, subjects in
  // parallel: a read-back needs only the subject's redrawer (no graph, no
  // registry), it was most of a run's single-threaded time, and reading one
  // enumerates isomorphisms onto the subject's thickening -- kept across runs
  // in read_back_cache. Then, serially and in the same order as before, each
  // subject is interned, certified and its read-backs added: the graph is
  // exactly what the one-phase loop built.
  struct SubjectRead {
    std::string name, pd;
    int layers = 2;
    std::vector<const Cobordism *> ws;
    std::unique_ptr<outgoing::OutgoingReader> redraw;
    std::string buildError; // the redrawer could not be built
    std::vector<std::optional<outgoing::OutgoingLink>> links;
    std::vector<std::string> why;
    std::vector<char> invariant; // the read-back broke an invariant (logic_error)
    size_t cacheHits = 0;
  };
  std::vector<SubjectRead> reads;
  for (const auto &[name, ws] : bySubject) {
    const linknaming::TableEntry *e = tables.entry(name);
    auto pdIt = tablePD_.find(name);
    if (!e || pdIt == tablePD_.end()) continue;
    // One read per (searched PD, layers): the table's PD unless the sidecar
    // recorded the diagram the cobordism was really searched on.
    std::map<std::pair<std::string, int>, std::vector<const Cobordism *>> byDiagram;
    for (const Cobordism &w : ws)
      byDiagram[{w.incomingPD.empty() ? pdIt->second : w.incomingPD, w.layers}].push_back(&w);
    for (const auto &[key, group] : byDiagram) {
      SubjectRead r;
      r.name = name;
      r.pd = key.first;
      r.layers = key.second;
      r.ws = group;
      reads.push_back(std::move(r));
    }
  }
  // Read by a pool of `threads` readers that lives for the whole load, in
  // subject order, at most `window` subjects ahead of the assembly below, which
  // takes each subject as soon as it is read: so at most a window of subjects'
  // thickenings is held at once (reading every subject first held them all:
  // 5.5 GB on L10a174), and no reader waits at a batch's end for the batch's
  // slowest subject (2026-10-02: a fresh pool per batch of 2 x threads
  // subjects, joined before
  // the batch was assembled, cost L10n112's cold load 6% of its wall).
  const std::string &readBackCache = ld.readBackCache;
  auto readSubject = [&](size_t i) {
    SubjectRead &r = reads[i];
    try {
      r.redraw = std::make_unique<outgoing::OutgoingReader>(r.pd, r.layers);
    } catch (const std::exception &ex) {
      r.buildError = ex.what();
      return;
    }
    outgoing::ReadBacks cache(readBackCache, r.pd, r.layers,
                       readBackCache.empty() ? std::string() : r.redraw->buildChecksum());
    r.links.resize(r.ws.size());
    r.why.resize(r.ws.size());
    r.invariant.assign(r.ws.size(), 0);
    for (size_t k = 0; k < r.ws.size(); ++k) {
      const std::string key = cobordisms::cobordismKey(r.ws[k]->pairsig);
      if (const outgoing::CachedReadBack *c = cache.get(key)) {
        r.links[k] = c->link;
        r.why[k] = c->why;
        continue;
      }
      std::string why;
      try {
        r.links[k] = r.redraw->outgoingLinkFast(r.ws[k]->pairsig, why);
        r.why[k] = why;
        cache.put(key, {r.links[k], why});
      } catch (const std::logic_error &ex) {
        r.invariant[k] = 1; // reported, never cached
        r.why[k] = ex.what();
      } catch (const std::exception &ex) {
        r.why[k] = ex.what(); // not cached either: it may be transient
      }
    }
    r.cacheHits = cache.hits();
    try {
      cache.flush();
    } catch (const std::exception &ex) {
      std::cerr << "[!] read-back cache not written for " << r.name << ": " << ex.what()
                << "\n";
    }
  };
  const size_t window = 2 * static_cast<size_t>(std::max(ld.threads, 1u));
  std::mutex readMutex;
  std::condition_variable readCv;
  std::vector<char> readDone(reads.size(), 0);
  size_t assembling = 0;      // subjects before this one are assembled (under readMutex)
  bool stopReading = false;   // the assembly left early (under readMutex)
  size_t nextRead = 0;        // the next subject to read (under readMutex)
  auto reader = [&] {
    for (;;) {
      // A reader takes a subject only once it is inside the window, so no
      // reader sits on a subject it may not read yet while the subjects before it
      // wait for a thread.
      size_t i;
      {
        std::unique_lock<std::mutex> lock(readMutex);
        readCv.wait(lock, [&] {
          return stopReading || nextRead >= reads.size() || nextRead < assembling + window;
        });
        if (stopReading || nextRead >= reads.size()) return;
        i = nextRead++;
      }
      readSubject(i);
      {
        std::lock_guard<std::mutex> lock(readMutex);
        readDone[i] = 1;
      }
      readCv.notify_all();
    }
  };
  // Joined on every way out of the assembly, an exception included: the
  // readers stop at once, finishing only the subjects they hold.
  struct ReaderPool {
    std::vector<std::thread> threads;
    std::function<void()> release;
    ~ReaderPool() {
      release();
      for (std::thread &t : threads) t.join();
    }
  } pool;
  pool.release = [&] {
    {
      std::lock_guard<std::mutex> lock(readMutex);
      stopReading = true;
    }
    readCv.notify_all();
  };
  for (size_t t = 0, k = std::min<size_t>(std::max(ld.threads, 1u), reads.size()); t < k; ++t)
    pool.threads.emplace_back(reader);
  // Subject i, read: what the assembly waited for it is read-back time.
  auto awaitSubject = [&](size_t i) {
    const auto tWait = Clock::now();
    {
      std::unique_lock<std::mutex> lock(readMutex);
      assembling = i;
      readCv.notify_all();
      readCv.wait(lock, [&] { return readDone[i] != 0; });
    }
    readSeconds += secondsSince(tWait);
  };
  const auto tSubjects = Clock::now();
  size_t cacheHits = 0;
  // A subject is interned once per PD it was searched on: the table's own
  // diagram (simplified, then the registry's exact tests), or the recorded
  // diagram of a goal run's search, which is a link's diagram of another run:
  // already reduced, so interned as it is. An own subject must be THIS link.
  struct SubjectLink {
    LinkMatch match;
    GaussDiagram diagram; ///< the searched diagram, in its PD's component order
  };
  std::map<std::string, SubjectLink> interned; // by searched PD
  for (size_t i = 0; i < reads.size(); ++i) {
    awaitSubject(i);
    SubjectRead &r = reads[i];
    cacheHits += r.cacheHits;
    auto seen = interned.find(r.pd);
    if (seen == interned.end()) {
      SubjectLink rn;
      if (r.pd == tablePD_.at(r.name)) {
        rn.diagram = GaussDiagram::of(tables.entry(r.name)->diagram);
        rn.match = reg.intern(linknaming::simplifyKeepingComponents(rn.diagram),
                              kFrozenRowLabel + r.name);
      } else {
        rn.diagram = GaussDiagram::of(linknaming::linkFromTablePD(r.pd));
        rn.match = reg.intern(rn.diagram, kFrozenRowLabel + r.name + " (recorded diagram)");
      }
      seen = interned.emplace(r.pd, std::move(rn)).first;
      if (ownSubjects.count(r.name) && seen->second.match.link != n) {
        ++refusedSubjects;
        std::cout << "[!] master row " << r.name << " did not intern as node " << n
                  << " (got " << seen->second.match.link << "); not used\n";
      }
    }
    const LinkMatch &m = seen->second.match;
    if (ownSubjects.count(r.name) && m.link != n) continue;
    SearchedLink searched;
    searched.link = m.link;
    searched.diagram = seen->second.diagram;
    searched.linkMap = m.componentMap;
    searched.pd = r.pd;
    searched.layers = r.layers;
    std::unique_ptr<CobordismAssembler> assembler;
    try {
      if (!r.buildError.empty()) throw std::runtime_error(r.buildError);
      assembler = std::make_unique<CobordismAssembler>(g, reg, searched, std::move(r.redraw));
    } catch (const std::exception &ex) {
      ++refusedSubjects;
      std::cout << "[!] master row " << r.name << " refused: " << ex.what() << "\n";
      continue;
    }
    for (size_t k = 0; k < r.ws.size(); ++k) {
      const Cobordism *w = r.ws[k];
      const std::string key = cobordisms::cobordismKey(w->pairsig);
      AddedCobordism added;
      if (r.invariant[k]) {
        ++ld.invariantFailures;
        std::cout << "[!!] master witness " << key << ": INVARIANT: " << r.why[k] << "\n";
      } else if (!r.links[k]) {
        added.why = r.why[k];
      } else {
        try {
          added = assembler->addRead(*r.links[k], w->genus, "master:" + key);
        } catch (const std::logic_error &ex) {
          ++ld.invariantFailures;
          std::cout << "[!!] master witness " << key << ": INVARIANT: " << ex.what() << "\n";
        } catch (const std::exception &ex) {
          added.why = ex.what();
        }
      }
      if (!added.ok) { ++failed; continue; }
      ++assembled;
      CobordismSource info{"master", searched.pd, key, added, r.layers, w->pairsig, searched.linkMap};
      if (added.direct) ld.sources.direct["master:" + key] = info;
      else ld.sources.byCobordism[added.cobordism] = info;
    }
  }
  // Phase B, the serial assembly: the subject loop less its read-backs.
  const double assembleSeconds = secondsSince(tSubjects) - readSeconds;
  const auto tName = Clock::now();
  std::vector<LinkId> fresh;
  for (size_t m = linksBefore; m < g.linkCount(); ++m) fresh.push_back(static_cast<LinkId>(m));
  ld.axioms.name(fresh, ld.axioms.depth.at(n) + 1);
  const double nameSeconds = secondsSince(tName);
  const auto tProp = Clock::now();
  g.propagate();
  g.propagateLower();
  const double propagateSeconds = secondsSince(tProp);
  const double loadSeconds = secondsSince(tLoad);
  auto best = g.best(ld.target, ld.goal);
  std::ostringstream o;
  o << std::fixed << std::setprecision(1);
  o << "{\"master\":\"" << json::escape(tableName[n]) << "\",\"node\":" << n
    << ",\"wall\":" << loadSeconds << ",\"readback_s\":" << readSeconds
    << ",\"assemble_s\":" << assembleSeconds << ",\"name_s\":" << nameSeconds
    << ",\"propagate_s\":" << propagateSeconds
    << ",\"own_rows\":" << own.size() << ",\"reverse_witnesses\":" << reverse
    << ",\"subject_rows\":" << bySubject.size() << ",\"refused_rows\":" << refusedSubjects
    << ",\"skipped_genus\":" << skippedGenus
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"read_back_cache_hits\":" << cacheHits
    << ",\"nodes\":" << g.linkCount() << ",\"target_best\":"
    << (best ? std::to_string(best->genus) : "null") << "}";
  ld.log(o.str());
  std::cout << kFrozenMasterRowsLine << n << " (" << tableName[n] << "): "
            << assembled << kFrozenWitnessesAssembled << failed << " failed, " << skippedGenus
            << " skipped by genus; target best " << (best ? std::to_string(best->genus) : "none")
            << "\n";
  return loadSeconds;
}

} // namespace bounds
