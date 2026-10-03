//
//  targets.h
//
//  What a run searches from: table rows, or one diagram.
//

#ifndef SURFER_COBOUND_TARGETS_H
#define SURFER_COBOUND_TARGETS_H

#include <filesystem>
#include <map>
#include <string>
#include <unordered_map>
#include <vector>

#include "cobound/solver/literature.h"
#include "cobound/solver/verdicts.h"

/*! \file utils/surfer/cobound/driver/targets.h
 *  \brief A run's targets: the rows of a table (`targets`), each with its
 *  literature interval, or one diagram (`target_pd`). Without a goal they
 *  are searched once each, in crossing order, those within max_crossings;
 *  a goal run has one.
 */

namespace targets {

/**
 * The rows of a table to search (Name,PD Notation,Genus-4D), in file order.
 * A row to search must state its literature bounds: a malformed field stops
 * the run, never becomes a bound.
 * \throws std::runtime_error a malformed 4-genus field.
 */
std::vector<cobordismgraph::InputRow> loadInputCsv(const std::filesystem::path &path);

/**
 * One diagram as a row: its name, PD and crossings, and its literature
 * interval from `tables` when the name is a table row there, else [0, 99]
 * (as the retired child hop wrote an untabulated row).
 */
cobordismgraph::InputRow oneDiagram(const std::string &name, const std::string &pd,
                                    const std::vector<std::filesystem::path> &tables);

/**
 * The rows a run without a goal searches, in order: those of at most
 * `maxCrossings` crossings, stable-sorted by crossings. Each row above it
 * that the verdicts do not hold yet is recorded there as `skipped`.
 */
std::vector<cobordismgraph::InputRow>
searchOrder(const std::vector<cobordismgraph::InputRow> &rows, int maxCrossings,
            std::unordered_map<std::string, verdicts::OutputRow> &outputRows);

/** name -> the table's PD string, as the atlas's rows were searched from it. */
std::map<std::string, std::string> tablePDs(const std::vector<std::string> &files);

} // namespace targets

#endif // SURFER_COBOUND_TARGETS_H
