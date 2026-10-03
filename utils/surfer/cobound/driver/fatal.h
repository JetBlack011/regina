//
//  fatal.h
//
//  A halt on a state that cannot occur (a run without a goal, and solve).
//

#ifndef SURFER_COBOUND_FATAL_H
#define SURFER_COBOUND_FATAL_H

#include <functional>
#include <string>

/*! \file utils/surfer/cobound/driver/fatal.h
 *  \brief Fatal-bug detection, as verifyslicegenus had it: a found surface
 *  implying a genus BELOW an established lower bound (a contradiction in
 *  the cobordism graph), an impossible accounting state, a broken seed
 *  invariant or a petal linking number that disagrees under the audit is a
 *  mathematical impossibility, not a data problem: it means the search
 *  itself computed something wrong, and nothing else the run would report
 *  can be trusted. So the run halts, exit code 2, once what it found is
 *  written (beforeHalt()). Goal runs halt their own way (scheduler.h).
 */

namespace fatal {

/// Records `message` as the reason for a halt, if none is recorded yet
/// (the first detector wins). Safe from any thread.
void flag(std::string message);

/// Whether flag() was called.
bool flagged();

/// What a halt writes first (a run's pending cobordisms, signed).
void beforeHalt(std::function<void()> write);

/// If flag() was called: runs beforeHalt()'s writer, prints the banner and
/// the message on stderr, and exits with code 2. Call only from the main
/// thread with no search running.
void haltIfFlagged();

} // namespace fatal

#endif // SURFER_COBOUND_FATAL_H
