//
//  atomicwrite.h
//
//  Writing a whole file so that no reader, and no crash, ever sees half of it.
//

#ifndef SURFER_REPORT_ATOMICWRITE_H
#define SURFER_REPORT_ATOMICWRITE_H

#include <filesystem>
#include <functional>
#include <ostream>

/*! \file utils/surfer/surfer/report/atomicwrite.h
 *  \brief atomicWrite(): the one way a file is replaced whole.
 *
 *  The library's own files (a search's frontier, a pair-signature context)
 *  and cobound's (the verdicts table, the read-back cache, the database's one
 *  full rewrite) are all replaced through it. It lives in surfer, the lowest
 *  part that writes such a file; cobound's appends are
 *  cobound/cobordisms/appendonly.h's.
 */

namespace report {

/**
 * Replaces `path` with what `write` writes, atomically and durably: `write`
 * fills a temporary file beside `path` (named uniquely, by process and call,
 * so two writers never share one), which is flushed and fsynced, renamed
 * over `path`, and the directory fsynced. A reader sees the old file or the
 * new one, never part of either, and once this returns the new one survives
 * a crash.
 *
 * \exception std::runtime_error the temporary file cannot be written or
 * synced, or the rename fails; `path` is then untouched and the temporary
 * file removed. Whatever `write` throws propagates the same way.
 */
void atomicWrite(const std::filesystem::path &path,
                 const std::function<void(std::ostream &)> &write);

} // namespace report

#endif // SURFER_REPORT_ATOMICWRITE_H
