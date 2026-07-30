//
//  csvwriter.h
//
//  Created by John Teague on 06/19/2024.
//

#ifndef CSVWRITER_H

#define CSVWRITER_H

#include <filesystem>
#include <fstream>
#include <memory>
#include <mutex>
#include <string>
#include <string_view>
#include <vector>

/*! \file utils/surfer/csvwriter.h
 *  \brief A sharded, multi-threaded CSV writer.
 */

/**
 * Escapes a single CSV field per RFC 4180: wraps in double quotes (doubling
 * any embedded quotes) only when the field contains a comma, quote, or
 * newline -- otherwise returned unquoted.
 */
std::string csvField(std::string_view s);

/**
 * Buffers CSV rows, distributing writers across a small, bounded pool of
 * shard files (at most MAX_SHARDS, and never more than numThreads) rather
 * than one shard per thread -- a caller's own thread count can be set far
 * higher than the system's open-file-descriptor limit (ulimit -n), and one
 * persistently-open file handle per thread risks exhausting it (silently,
 * before this fix: failed opens/writes went unchecked, so hitting the
 * limit meant rows -- or the entire output -- vanished with no error, just
 * a misleading "Wrote CSV output" message). Threads sharing a shard are
 * serialized by that shard's own mutex, which only contends among threads
 * mapped to the same shard, not globally.
 *
 * finalize() concatenates every shard into the requested output path
 * (header first) and removes the shard files -- call it once,
 * single-threaded, only after every thread that might call writeRow() has
 * already joined.
 *
 * \warning No usable output exists until finalize() has been called: a
 * crash or kill mid-run leaves only orphaned, unconsolidated shard files.
 * Fine for a one-shot dump of everything a single search finds (surfer.cpp's
 * -o, or another program's own analogous debug/log output), but not
 * appropriate for output that needs to survive an interruption -- a caller
 * needing that should use a different, simpler durability pattern instead
 * (e.g. rewrite-and-atomically-rename after each unit of progress).
 *
 * shardForThisThread() caches its result in a thread_local, function-static
 * pointer -- correct as long as at most one CsvWriter of a given instance
 * is ever written to concurrently by a given thread (i.e. a thread reusing
 * itself across two different CsvWriter instances gets a fresh cache miss
 * against the second one, which is fine -- the cache is per-instance state,
 * not per-thread-global).
 */
class CsvWriter {
public:
  CsvWriter(std::filesystem::path outputPath, std::string headerLine,
           unsigned numThreads);

  void writeRow(const std::string &row);

  // Any failure here (opening the final output, reading a shard back,
  // writing to it) throws rather than limping on -- a partially-written
  // CSV that looks complete but silently isn't is worse than a loud
  // failure that makes clear the run needs to be retried (e.g. with a
  // lower thread count, or a raised ulimit -n).
  void finalize();

private:
  struct Shard {
    std::filesystem::path path;
    std::ofstream file;
    std::string buffer;
    std::mutex mutex;
  };

  static constexpr size_t FLUSH_THRESHOLD = 1 << 16; // ~64KB
  static constexpr unsigned MAX_SHARDS = 64;

  // Caller must hold shard.mutex.
  static void flush(Shard &shard);

  Shard &shardForThisThread();

  std::filesystem::path outputPath_;
  std::string headerLine_;
  unsigned maxShards_;
  std::mutex registryMutex_;
  size_t nextShard_ = 0; // guarded by registryMutex_
  std::vector<std::unique_ptr<Shard>> shards_;
};

#endif // CSVWRITER_H
