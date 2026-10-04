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
 */

namespace runrecords {

using bounds::LinkId;

/// What the records read of a goal run's graph.
struct GraphView {
  const bounds::ProofGraph &g;
  const bounds::LinkRegistry &reg;
  const std::map<LinkId, std::string> &tableName;
  const std::map<LinkId, int> &depth;
  LinkId target = -1;
  bounds::Partition goal; ///< the target's goal partition
};

/// Appends one line to <work>/cascade.jsonl.
void append(const std::string &work, const std::string &line);

/**
 * <work>/profiles.jsonl: every node's Pareto (partition, genus) entries and
 * per-partition lower bounds (ProofGraph::profileFields()), with its name
 * (`subjectName`, or the unlink's spelling for a crossingless node other
 * than the target), table name, label, depth, crossings and whether it was
 * searched (`searched`).
 */
void writeProfiles(const std::string &work, const GraphView &v,
                   const std::function<std::string(LinkId)> &subjectName,
                   const std::function<bool(LinkId)> &searched);

/**
 * <work>/node_bounds.jsonl at a run's end: every node's identity, diagram,
 * proved profile entries and lower bounds (README.md, "Lower-bound mode").
 */
void writeLinkBounds(const std::string &work, const GraphView &v);

/**
 * <work>/lower_report.jsonl (lower_report): for every tabulated node, what
 * its literature lower bound carries to the target, as what-ifs on the
 * run's `threads`; prints the `[+] lower report:` line.
 */
void writeLowerReport(const std::string &work, const GraphView &v,
                      const linknaming::ExactTables &tables, const std::string &targetName,
                      const std::map<std::string, bool> &special, unsigned threads);

/// <work>/nodes.csv: the cascade: subjects searched (`subjects`, by node),
/// with their diagrams, for the atlas's results/cascade/nodes.csv.
void writeLinksCsv(const std::string &work, const std::map<LinkId, std::string> &subjects,
                   const bounds::LinkRegistry &reg);

} // namespace runrecords

#endif // SURFER_COBOUND_RUNRECORDS_H
