//
//  progress.h
//
//  A search's live progress block on stderr.
//

/*! \file utils/surfer/surfer/report/progress.h
 *  \brief A block of status lines redrawn in place (the format of each
 *  driver's block is its own; this only redraws it).
 */

#ifndef SURFER_REPORT_PROGRESS_H
#define SURFER_REPORT_PROGRESS_H

#include <cstddef>
#include <mutex>
#include <string>

namespace report {

/**
 * A block of status lines redrawn in place on stderr: each draw() moves the
 * cursor back to the previous block's first line, erases to the end of the
 * screen, and prints the new block (ANSI `ESC[<n>F ESC[0J`).
 *
 * Thread-safe: draw() is normally called from one reporter thread, but
 * commitLine() can fire from any worker thread at any time (surfer's
 * queue-drain-pause notices); both take one mutex, so two threads never
 * interleave writes to std::cerr or race on the block's height.
 */
class RollingReport {
  public:
    /** Replaces the block on screen with `text` (its lines end in '\n'). */
    void draw(const std::string &text);

    /** Prints `text` as permanent output, then forget()s the block, so the
     *  next draw() starts underneath `text` instead of erasing it. */
    void commitLine(const std::string &text);

    /** Forgets the block on screen: the next draw() erases nothing. For a
     *  caller that has printed permanent lines of its own since. */
    void forget();

  private:
    std::mutex mutex_;
    std::size_t prevLines_ = 0;
};

} // namespace report

#endif
