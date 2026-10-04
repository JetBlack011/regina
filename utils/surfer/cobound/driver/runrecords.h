//
//  runrecords.h
//
//  A goal run's directory: the records it writes (names and formats frozen).
//

#ifndef SURFER_COBOUND_RUNRECORDS_H
#define SURFER_COBOUND_RUNRECORDS_H

#include <functional>
#include <map>
#include <string>

#include "cobound/bounds/cobordismgraph.h"
#include "cobound/bounds/links.h"
#include "cobound/bounds/partition.h"
#include "linknaming/tables.h"

/*! \file utils/surfer/cobound/driver/runrecords.h
 *  \brief The files a goal run writes into its directory besides each
 *  search's hop_<k>_n<node>/ and the certificates (bounds/certificate.h):
 *  cascade.jsonl (one record per search, per database load, and the run's
 *  own), profiles.jsonl, node_bounds.jsonl, lower_report.jsonl and
 *  nodes.csv. Every name and format is frozen (plan, Hard constraint 2):
 *  the atlas's cascade_layer.py, cascade_record.py and search_breadth.py
 *  read them.
 *
 *  Every writer throws std::runtime_error when its file cannot be opened or
 *  a write to it fails (the stream is flushed and checked): a run whose
 *  record is incomplete halts and claims no goal (scheduler.cpp,
 *  Scheduler::recordWriteFailure()).
 */

namespace runrecords {

using bounds::LinkId;

/// What the records read of a goal run's graph.
struct GraphView {
  const bounds::CobordismGraph &g;
  const bounds::LinkRegistry &reg;
  const std::map<LinkId, std::string> &tableName;
  const std::map<LinkId, int> &depth;
  LinkId target = -1;
  bounds::Partition goal; ///< the target's goal partition
};

/// Appends one line to <work>/cascade.jsonl.
void append(const std::string &work, const std::string &line);

/**
 * <work>/profiles.jsonl: every link's Pareto (partition, genus) entries and
 * per-partition lower bounds (CobordismGraph::partitionGeneraFields()), with its name
 * (`subjectName`, or the unlink's spelling for a crossingless link other
 * than the target), table name, label, depth, crossings and whether it was
 * searched (`searched`).
 */
void writePartitionGenera(const std::string &work, const GraphView &v,
                   const std::function<std::string(LinkId)> &subjectName,
                   const std::function<bool(LinkId)> &searched);

/**
 * <work>/node_bounds.jsonl at a run's end: every link's identity, diagram,
 * proved partition genera and lower bounds (README.md, "Lower goals").
 */
void writeLinkBounds(const std::string &work, const GraphView &v);

/**
 * <work>/lower_report.jsonl (lower_report): for every tabulated link, what
 * its literature lower bound carries to the target, as what-ifs on the
 * run's `threads`; prints the `[+] lower report:` line.
 */
void writeLowerReport(const std::string &work, const GraphView &v,
                      const linknaming::Tables &tables, const std::string &targetName,
                      const std::map<std::string, bool> &special, unsigned threads);

/// <work>/nodes.csv: the `cascade:` subjects searched (`subjects`, by link),
/// with their diagrams, for the atlas's results/cascade/nodes.csv.
void writeLinksCsv(const std::string &work, const std::map<LinkId, std::string> &subjects,
                   const bounds::LinkRegistry &reg);

} // namespace runrecords

#endif // SURFER_COBOUND_RUNRECORDS_H
