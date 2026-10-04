//
//  databasecobordisms.h
//
//  A database's cobordisms of a link about to be searched, as graph cobordisms.
//

#ifndef SURFER_COBOUND_DATABASECOBORDISMS_H
#define SURFER_COBOUND_DATABASECOBORDISMS_H

#include <functional>
#include <map>
#include <set>
#include <string>
#include <vector>

#include "cobound/bounds/axioms.h"
#include "cobound/bounds/certificate.h"
#include "cobound/bounds/cobordismgraph.h"
#include "cobound/bounds/links.h"
#include "cobound/bounds/partition.h"
#include "cobound/cobordisms/database.h"
#include "linknaming/tables.h"

/*! \file utils/surfer/cobound/bounds/databasecobordisms.h
 *  \brief A goal run's free cobordisms: those a read-only database
 *  (master_cobordisms) holds for a table link, turned into the graph's
 *  cobordisms before that link is searched (they are what its search would
 *  find again). It holds no reading code: lines are selected and
 *  offset-indexed by cobordisms/database (DatabaseIndex), and each
 *  cobordism's outgoing link is read back by outgoing/fromdatabase
 *  (OutgoingReader, ReadBacks). Without a goal nothing is loaded at all
 *  (plan, "Startup per process").
 */

namespace bounds {

/// Where loading puts what it reads: the run's graph, its links and their
/// names and outside facts, and how each new cobordism's surface is found.
struct DatabaseLoad {
  CobordismGraph &g;
  LinkRegistry &reg;
  LinkAxioms &axioms;
  const linknaming::Tables &tables;
  CobordismSources &sources;
  /// Read-backs and assemblies that broke an invariant (reported, dropped).
  int &invariantFailures;
  /// A cobordism of greater genus is on no proof the run can use, so it is
  /// skipped before it is read back (the genus is in its line).
  int maxGenus = 0;
  unsigned threads = 1;
  /// read_back_cache: read-backs kept across runs; "" for none.
  std::string readBackCache;
  /// The target and its goal partition, for the record lines' target best.
  LinkId target = -1;
  Partition goal;
  /// Appends one line to cascade.jsonl.
  std::function<void(const std::string &)> log;
};

class DatabaseCobordisms {
public:
  /// Indexes `database` (read only) and the tables' PDs.
  /// \exception std::runtime_error the database cannot be read.
  DatabaseCobordisms(const std::string &database, const std::vector<std::string> &tables);

  /// Subjects the database holds cobordisms of.
  size_t subjects() const { return index_.subjects(); }

  /// Whether the database holds cobordisms for link `n` (a table link):
  /// cobordisms of its class's table entries (named into `variants`), or of
  /// other subjects whose outgoing link has its base name.
  bool subjectsFor(LinkId n, const LinkAxioms &axioms, const linknaming::Tables &tables,
               std::vector<std::string> *variants = nullptr) const;
  /// Whether `n`'s cobordisms were loaded already.
  bool loaded(LinkId n) const { return done_.count(n) > 0; }

  /**
   * Loads link `n`'s cobordisms into the graph: its own subjects' and the
   * hinted other subjects', read back (on `load.threads` readers, at most a
   * window of subjects ahead of the serial assembly), each subject interned and
   * certified, its cobordisms added; then names the new links (at n's
   * depth + 1) and relaxes the graph. Writes the cascade.jsonl record and
   * the `[+] master rows of node` line. Returns its wall seconds.
   */
  double load(LinkId n, DatabaseLoad &load);

private:
  cobordisms::DatabaseIndex index_;
  std::map<std::string, std::string> tablePD_;
  std::set<LinkId> done_;
};

} // namespace bounds

#endif // SURFER_COBOUND_DATABASECOBORDISMS_H
