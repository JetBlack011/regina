//
//  appendonly.h
//
//  Append-only files: locks, complete writes, fsyncs, torn last lines.
//

#ifndef SURFER_COBOUND_APPENDONLY_H
#define SURFER_COBOUND_APPENDONLY_H

#include <string>

/*! \file utils/surfer/cobound/cobordisms/appendonly.h
 *  \brief How cobound appends to a file that is never rewritten: the
 *  cobordism database and its `.rows.csv` sidecar, a search's pending
 *  cobordisms (kept.csv), the read-back cache. One lock, one complete write,
 *  one fsync, one rule for a torn last line. (A file replaced whole goes
 *  through surfer/report/atomicwrite.h instead.)
 */

namespace appendonly {

/** An exclusive flock(2) on `path` (created if absent), held for the
 *  object's life. */
class FileLock {
  public:
    explicit FileLock(const std::string &path);
    ~FileLock();
    FileLock(const FileLock &) = delete;
    FileLock &operator=(const FileLock &) = delete;

  private:
    int fd_ = -1;
};

/** Writes all of `bytes` to `fd` (retrying short writes and EINTR).
 *  `what` names the file in the exception. */
void writeAll(int fd, const std::string &bytes, const std::string &what);

/** fsync(2)s `fd`; `what` names the file in the exception. */
void sync(int fd, const std::string &what);

/**
 * A torn append (a crash mid-write) leaves a last line with no newline. If
 * `fd`'s file does not end in '\n', truncates it back to just after its last
 * '\n' (to empty if it has none) and says so on stderr. Returns the number of
 * bytes cut.
 */
long cutTornLine(int fd, const std::string &what);

/** Whether `append()` fsyncs. */
enum class Sync { no, yes };

/**
 * Appends `bytes` to `path` (created if absent): a torn last line is cut
 * first, then all of `bytes` written, then, with Sync::yes, the file fsynced.
 * \exception std::runtime_error any step fails.
 */
void append(const std::string &path, const std::string &bytes, Sync sync);

} // namespace appendonly

#endif // SURFER_COBOUND_APPENDONLY_H
