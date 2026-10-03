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
#include "linknaming/diagrams/simplification.h"

namespace cascade {

using exactnaming::GaussDiagram;
using timers::Clock;
using timers::secondsSince;

namespace {

// A database line as a goal run carries it (cobordisms/database's reader
// parses it): the outgoing name, the pair signature, genus, outgoing
// components, layers and the searched row's PD.
struct Witness {
  std::string other, pairsig;
  int genus = 0, otherComponents = 0, layers = 2;
  /// The PD its row was searched on, when the file's `.rows.csv` sidecar
  /// records it (a cascade hop's row is its node's simplified diagram, not
  /// the table's); empty for the table's PD.
  std::string rowPD;
};

Witness carried(const witnessstore::StoredCobordism &s) {
  return {s.witness.other, s.witness.pairSig, s.witness.genus, s.witness.otherComponents,
          s.witness.thickenLayers, s.rowPD};
}

} // namespace

DatabaseCobordisms::DatabaseCobordisms(const std::string &database,
                                       const std::vector<std::string> &tables)
    : index_(database) {
  // name -> the table's PD string, as the atlas's rows were searched from it.
  for (const std::string &file : tables)
    for (const exactnaming::TableRow &row : exactnaming::readTableRows(file))
      tablePD_[row.name] = row.pd;
}

bool DatabaseCobordisms::rowsFor(NodeId n, const NodeAxioms &axioms,
                                 const exactnaming::ExactTables &tables,
                                 std::vector<std::string> *rows) const {
  auto it = axioms.tableName.find(n);
  if (it == axioms.tableName.end()) return false;
  // Every table entry of this link's class (one oriented link up to mirror
  // and global reversal) is the same node; each is a row of its own.
  const std::string canon = axioms.classOf(it->second);
  const exactnaming::TableEntry *e = tables.entry(it->second);
  if (!e) return false;
  bool any = !index_.byOutgoing(witnessstore::DatabaseIndex::base(it->second), 1).empty();
  for (const exactnaming::TableEntry *v : tables.variants(e->base))
    if (axioms.classOf(v->name) == canon && index_.has(v->name)) {
      any = true;
      if (rows) rows->push_back(v->name);
    }
  return any;
}

double DatabaseCobordisms::load(NodeId n, DatabaseLoad &ld) {
  const auto tLoad = Clock::now();
  done_.insert(n);
  std::vector<std::string> own;
  if (!rowsFor(n, ld.axioms, ld.tables, &own)) return secondsSince(tLoad);
  ProofGraph &g = ld.g;
  NodeRegistry &reg = ld.reg;
  const exactnaming::ExactTables &tables = ld.tables;
  std::map<NodeId, std::string> &tableName = ld.axioms.tableName;
  double readSeconds = 0; // phase A: the assembly's waits for the readers
  const size_t nodesBefore = g.nodeCount();
  int assembled = 0, failed = 0, refusedRows = 0;
  // A witness can be on an upper proof only if its genus is at most the
  // goal (glue() never lowers a genus), and on a lower proof only if the
  // charge it costs, at least its genus, is affordable: at most the largest
  // special source's bound minus the lower goal. Anything above both is
  // skipped before it is read back (the genus is in the CSV line).
  const int maxGenus = ld.maxGenus;
  size_t skippedGenus = 0;
  auto keep = [&](const Witness &w) {
    if (w.genus <= maxGenus) return true;
    ++skippedGenus;
    return false;
  };
  // Witnesses by subject row: this node's own rows (every witness), and
  // other rows whose recorded far side names this node's base (a hint: the
  // far side is redrawn and identified exactly like any other).
  std::map<std::string, std::vector<Witness>> bySubject;
  std::set<std::string> ownRows(own.begin(), own.end());
  for (const std::string &name : own)
    for (const witnessstore::StoredCobordism &s : index_.rows(name))
      if (const Witness w = carried(s); keep(w)) bySubject[s.witness.subject].push_back(w);
  size_t reverse = 0;
  for (const witnessstore::StoredCobordism &s :
       index_.byOutgoing(witnessstore::DatabaseIndex::base(tableName[n]), 300))
    if (const Witness w = carried(s); !ownRows.count(s.witness.subject) && keep(w)) {
      bySubject[s.witness.subject].push_back(w);
      ++reverse;
    }
  // Two phases. First every subject row's witnesses are read back, rows in
  // parallel: a read-back needs only the row's redrawer (no graph, no
  // registry), it was most of a run's single-threaded time, and reading one
  // enumerates isomorphisms onto the row's thickening -- kept across runs in
  // read_back_cache. Then, serially and in the same order as before, each
  // row is interned, certified and its read-backs added: the graph is
  // exactly what the one-phase loop built.
  struct RowRead {
    std::string name, pd;
    int layers = 2;
    std::vector<const Witness *> ws;
    std::unique_ptr<farside::WitnessRedrawer> redraw;
    std::string buildError; // the redrawer could not be built
    std::vector<std::optional<farside::OutgoingLink>> links;
    std::vector<std::string> why;
    std::vector<char> invariant; // the read-back broke an invariant (logic_error)
    size_t cacheHits = 0;
  };
  std::vector<RowRead> reads;
  for (const auto &[name, ws] : bySubject) {
    const exactnaming::TableEntry *e = tables.entry(name);
    auto pdIt = tablePD_.find(name);
    if (!e || pdIt == tablePD_.end()) continue;
    // One read per (row PD, layers): the table's PD unless the sidecar
    // recorded the row the witness was really searched on.
    std::map<std::pair<std::string, int>, std::vector<const Witness *>> byRow;
    for (const Witness &w : ws)
      byRow[{w.rowPD.empty() ? pdIt->second : w.rowPD, w.layers}].push_back(&w);
    for (const auto &[key, group] : byRow) {
      RowRead r;
      r.name = name;
      r.pd = key.first;
      r.layers = key.second;
      r.ws = group;
      reads.push_back(std::move(r));
    }
  }
  // Read by a pool of `threads` readers that lives for the whole load, in
  // row order, at most `window` rows ahead of the assembly below, which takes
  // each row as soon as it is read: so at most a window of rows' thickenings
  // is held at once (reading every row first held them all: 5.5 GB on
  // L10a174), and no reader waits at a batch's end for the batch's slowest
  // row (2026-10-02: a fresh pool per batch of 2 x threads rows, joined before
  // the batch was assembled, cost L10n112's cold load 6% of its wall).
  const std::string &readBackCache = ld.readBackCache;
  auto readRow = [&](size_t i) {
    RowRead &r = reads[i];
    try {
      r.redraw = std::make_unique<farside::WitnessRedrawer>(r.pd, r.layers);
    } catch (const std::exception &ex) {
      r.buildError = ex.what();
      return;
    }
    RowReadBacks cache(readBackCache, r.pd, r.layers,
                       readBackCache.empty() ? std::string() : r.redraw->buildChecksum());
    r.links.resize(r.ws.size());
    r.why.resize(r.ws.size());
    r.invariant.assign(r.ws.size(), 0);
    for (size_t k = 0; k < r.ws.size(); ++k) {
      const std::string key = witnesskey::witnessKey(r.ws[k]->pairsig);
      if (const CachedReadBack *c = cache.get(key)) {
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
  size_t assembling = 0;      // rows before this one are assembled (under readMutex)
  bool stopReading = false;   // the assembly left early (under readMutex)
  size_t nextRead = 0;        // the next row to read (under readMutex)
  auto reader = [&] {
    for (;;) {
      // A reader takes a row only once it is inside the window, so no
      // reader sits on a row it may not read yet while the rows before it
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
      readRow(i);
      {
        std::lock_guard<std::mutex> lock(readMutex);
        readDone[i] = 1;
      }
      readCv.notify_all();
    }
  };
  // Joined on every way out of the assembly, an exception included: the
  // readers stop at once, finishing only the rows they hold.
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
  // Row i, read: what the assembly waited for it is read-back time.
  auto awaitRow = [&](size_t i) {
    const auto tWait = Clock::now();
    {
      std::unique_lock<std::mutex> lock(readMutex);
      assembling = i;
      readCv.notify_all();
      readCv.wait(lock, [&] { return readDone[i] != 0; });
    }
    readSeconds += secondsSince(tWait);
  };
  const auto tRows = Clock::now();
  size_t cacheHits = 0;
  // A row is interned once per PD it was searched on: the table's own
  // diagram (simplified, then the registry's exact tests), or the recorded
  // diagram of a cascade hop, which is a node's diagram of another run:
  // already reduced, so interned as it is. An own row must be THIS node.
  struct RowNode {
    NodeMatch match;
    GaussDiagram diagram; ///< the row's diagram, in the row PD's component order
  };
  std::map<std::string, RowNode> interned; // by row PD
  for (size_t i = 0; i < reads.size(); ++i) {
    awaitRow(i);
    RowRead &r = reads[i];
    cacheHits += r.cacheHits;
    auto seen = interned.find(r.pd);
    if (seen == interned.end()) {
      RowNode rn;
      if (r.pd == tablePD_.at(r.name)) {
        rn.diagram = GaussDiagram::of(tables.entry(r.name)->diagram);
        rn.match = reg.intern(simplifyKeepingComponents(rn.diagram), "row " + r.name);
      } else {
        rn.diagram = GaussDiagram::of(exactnaming::linkFromTablePD(r.pd));
        rn.match = reg.intern(rn.diagram, "row " + r.name + " (recorded diagram)");
      }
      seen = interned.emplace(r.pd, std::move(rn)).first;
      if (ownRows.count(r.name) && seen->second.match.node != n) {
        ++refusedRows;
        std::cout << "[!] master row " << r.name << " did not intern as node " << n
                  << " (got " << seen->second.match.node << "); not used\n";
      }
    }
    const NodeMatch &m = seen->second.match;
    if (ownRows.count(r.name) && m.node != n) continue;
    HopRow row;
    row.node = m.node;
    row.diagram = seen->second.diagram;
    row.nodeMap = m.componentMap;
    row.pd = r.pd;
    row.layers = r.layers;
    std::unique_ptr<HopAssembler> hop;
    try {
      if (!r.buildError.empty()) throw std::runtime_error(r.buildError);
      hop = std::make_unique<HopAssembler>(g, reg, row, std::move(r.redraw));
    } catch (const std::exception &ex) {
      ++refusedRows;
      std::cout << "[!] master row " << r.name << " refused: " << ex.what() << "\n";
      continue;
    }
    for (size_t k = 0; k < r.ws.size(); ++k) {
      const Witness *w = r.ws[k];
      const std::string key = witnesskey::witnessKey(w->pairsig);
      HopEdge he;
      if (r.invariant[k]) {
        ++ld.invariantFailures;
        std::cout << "[!!] master witness " << key << ": INVARIANT: " << r.why[k] << "\n";
      } else if (!r.links[k]) {
        he.why = r.why[k];
      } else {
        try {
          he = hop->addRead(*r.links[k], w->genus, "master:" + key);
        } catch (const std::logic_error &ex) {
          ++ld.invariantFailures;
          std::cout << "[!!] master witness " << key << ": INVARIANT: " << ex.what() << "\n";
        } catch (const std::exception &ex) {
          he.why = ex.what();
        }
      }
      if (!he.ok) { ++failed; continue; }
      ++assembled;
      EdgeInfo info{"master", row.pd, key, he, r.layers, w->pairsig, row.nodeMap};
      if (he.direct) ld.edges.direct["master:" + key] = info;
      else ld.edges.byEdge[he.edge] = info;
    }
  }
  // Phase B, the serial assembly: the row loop less its read-backs.
  const double assembleSeconds = secondsSince(tRows) - readSeconds;
  const auto tName = Clock::now();
  std::vector<NodeId> fresh;
  for (size_t m = nodesBefore; m < g.nodeCount(); ++m) fresh.push_back(static_cast<NodeId>(m));
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
    << ",\"subject_rows\":" << bySubject.size() << ",\"refused_rows\":" << refusedRows
    << ",\"skipped_genus\":" << skippedGenus
    << ",\"assembled\":" << assembled << ",\"failed\":" << failed
    << ",\"read_back_cache_hits\":" << cacheHits
    << ",\"nodes\":" << g.nodeCount() << ",\"target_best\":"
    << (best ? std::to_string(best->genus) : "null") << "}";
  ld.log(o.str());
  std::cout << "[+] master rows of node " << n << " (" << tableName[n] << "): "
            << assembled << " witnesses assembled, " << failed << " failed, " << skippedGenus
            << " skipped by genus; target best " << (best ? std::to_string(best->genus) : "none")
            << "\n";
  return loadSeconds;
}

} // namespace cascade
