//
//  commands.h
//
//  cobound's commands: run | solve | sign | draw | name | meridians.
//

#ifndef SURFER_COBOUND_COMMANDS_H
#define SURFER_COBOUND_COMMANDS_H

#include <string>
#include <vector>

/*! \file utils/surfer/cobound/driver/commands.h
 *  \brief The commands of the one program, each given its arguments after
 *  the command word and returning the process's exit code. run, solve and
 *  sign take their configuration from `--config FILE` and `--set
 *  key=value` (driver/config.h); draw and name also take their reference
 *  tools' flags (farsidediagram's, farsidename's), which the atlas and the
 *  cascade's checker pass; meridians takes peripheral_slopes's one
 *  subcommand word.
 *
 *  | command   | was                               | does |
 *  |-----------|-----------------------------------|------|
 *  | run       | verifyslicegenus (a search run), cascadesearch | searches from each target: without a goal each once (sweep.h), with one by the scheduler (scheduler.h) |
 *  | solve     | verifyslicegenus --solve-only     | re-derives every verdict from the database |
 *  | sign      | cascadesearch --sign-only         | signs a work directory's pending cobordisms into the database |
 *  | draw      | farsidediagram                    | stored cobordisms' outgoing links, drawn |
 *  | name      | farsidename                       | stored cobordisms' outgoing links, named |
 *  | meridians | peripheral_slopes                 | sig, dump, dump-link, dump-subset, slope |
 */

namespace commands {

int run(const std::vector<std::string> &args);
int solve(const std::vector<std::string> &args);
int sign(const std::vector<std::string> &args);
int draw(const std::vector<std::string> &args);
int name(const std::vector<std::string> &args);
int meridians(const std::vector<std::string> &args);

} // namespace commands

#endif // SURFER_COBOUND_COMMANDS_H
