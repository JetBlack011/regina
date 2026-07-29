//
//  surfer.cpp
//
//  Created by John Teague on 06/19/2024.
//

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <memory>
#include <mutex>
#include <optional>
#include <sstream>
#include <stdexcept>
#include <string_view>
#include <thread>

#include <unistd.h>

#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "collar.h"
#include "embeddingsearch.h"
#include "surfacesearch.h"
#include "knotbuilder.h"
#include "linkcomplement.h"

namespace {

// Redraws a block of text in place, using the same ANSI cursor-rewind
// sequence the old inline reporters used, so each call's output overwrites
// the previous one instead of scrolling the terminal.
//
// Thread-safe: draw() is normally called once a second from a single
// dedicated reporter thread, but commitLine() (see below) can also fire
// from an arbitrary worker thread at an arbitrary time relative to that
// (e.g. SurfaceSearch's queue-drain-pause callbacks) -- both take the same
// mutex, so two threads never interleave writes to std::cerr or race on
// prevLines_.
class RollingReport {
  std::mutex mutex_;
  size_t prevLines_ = 0;

public:
  void draw(const std::string &text) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (prevLines_ > 0)
      std::cerr << "\x1b[" << prevLines_ << "F\x1b[0J";
    std::cerr << text;
    prevLines_ =
        static_cast<size_t>(std::count(text.begin(), text.end(), '\n'));
  }

  // Prints `text` as permanent output, then resets this report's own
  // erase-tracking so its *next* draw() call redraws fresh underneath
  // `text` instead of erasing it -- for messages that need to survive
  // independently of this report's own periodic redraw cycle (e.g. a
  // queue-drain-pause/resume notice fired from a worker thread), rather
  // than being silently overwritten by the next draw().
  void commitLine(const std::string &text) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::cerr << text;
    prevLines_ = 0;
  }
};

// Escapes a single CSV field per RFC 4180: wraps in double quotes (doubling
// any embedded quotes) only when the field contains a comma, quote, or
// newline -- otherwise returned unquoted.
std::string csvField(std::string_view s) {
  if (s.find_first_of(",\"\n") == std::string_view::npos)
    return std::string(s);
  std::string out = "\"";
  for (char c : s) {
    if (c == '"')
      out += '"';
    out += c;
  }
  out += '"';
  return out;
}

// The columns every CSV row shares regardless of BoundaryCondition.
std::string csvRow(bool orientable, int genus, int punctures,
                   long long triangles, BoundaryCondition mostRestrictive) {
  std::ostringstream row;
  row << (orientable ? "true" : "false") << ',' << genus << ',' << punctures
      << ',' << triangles << ',' << boundaryConditionName(mostRestrictive);
  return row.str();
}

// Buffers CSV rows, distributing writers across a small, bounded pool of
// shard files (at most MAX_SHARDS, and never more than numThreads) rather
// than one shard per thread -- --threads can be set far higher than the
// system's open-file-descriptor limit (ulimit -n), and one persistently-
// open file handle per thread risks exhausting it (silently, before this
// fix: failed opens/writes went unchecked, so hitting the limit meant rows
// -- or the entire output -- vanished with no error, just a misleading
// "Wrote CSV output" message). Threads sharing a shard are serialized by
// that shard's own mutex, which only contends among threads mapped to the
// same shard, not globally.
//
// finalize() concatenates every shard into the requested output path
// (header first) and removes the shard files -- call it once,
// single-threaded, only after every thread that might call writeRow() has
// already joined (i.e. after search() itself has returned).
//
// shardForThisThread() caches its result in a thread_local, function-static
// pointer -- correct as long as at most one CsvWriter is ever alive per
// process, which is exactly how surfer.cpp uses it (constructed at most
// once per run, in runSearch()).
class CsvWriter {
public:
  CsvWriter(std::filesystem::path outputPath, std::string headerLine,
           unsigned numThreads)
      : outputPath_(std::move(outputPath)),
        headerLine_(std::move(headerLine)),
        maxShards_(std::max<unsigned>(1, std::min(numThreads, MAX_SHARDS))) {}

  void writeRow(const std::string &row) {
    Shard &shard = shardForThisThread();
    std::lock_guard<std::mutex> lock(shard.mutex);
    shard.buffer += row;
    shard.buffer += '\n';
    if (shard.buffer.size() >= FLUSH_THRESHOLD)
      flush(shard);
  }

  // Any failure here (opening the final output, reading a shard back,
  // writing to it) throws rather than limping on -- a partially-written
  // CSV that looks complete but silently isn't is worse than a loud
  // failure that makes clear the run needs to be retried (e.g. with a
  // lower --threads, or a raised ulimit -n).
  void finalize() {
    std::ofstream out(outputPath_, std::ios::trunc);
    if (!out)
      throw std::runtime_error("CsvWriter: failed to open " +
                               outputPath_.string() + " for writing");
    out << headerLine_ << "\n";
    for (auto &shard : shards_) {
      std::lock_guard<std::mutex> lock(shard->mutex);
      flush(*shard);
      shard->file.close();
      std::ifstream in(shard->path);
      if (!in)
        throw std::runtime_error("CsvWriter: failed to reopen " +
                                 shard->path.string() + " for consolidation");
      out << in.rdbuf();
      in.close();
      std::filesystem::remove(shard->path);
    }
    if (!out)
      throw std::runtime_error("CsvWriter: failed while writing " +
                               outputPath_.string());
    std::cerr << "[+] Wrote CSV output to " << outputPath_.string() << "\n";
  }

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
  static void flush(Shard &shard) {
    if (shard.buffer.empty())
      return;
    shard.file << shard.buffer;
    if (!shard.file)
      throw std::runtime_error("CsvWriter: failed writing to " +
                               shard.path.string());
    shard.buffer.clear();
  }

  Shard &shardForThisThread() {
    static thread_local Shard *cached = nullptr;
    if (cached)
      return *cached;
    std::lock_guard<std::mutex> lock(registryMutex_);
    if (shards_.size() < maxShards_) {
      auto shard = std::make_unique<Shard>();
      shard->path = outputPath_.string() + ".tmp." +
                   std::to_string(getpid()) + ".shard" +
                   std::to_string(shards_.size());
      shard->file.open(shard->path, std::ios::trunc);
      if (!shard->file)
        throw std::runtime_error("CsvWriter: failed to open " +
                                 shard->path.string() + " for writing");
      shards_.push_back(std::move(shard));
    }
    cached = shards_[nextShard_++ % shards_.size()].get();
    return *cached;
  }

  std::filesystem::path outputPath_;
  std::string headerLine_;
  unsigned maxShards_;
  std::mutex registryMutex_;
  size_t nextShard_ = 0; // guarded by registryMutex_
  std::vector<std::unique_ptr<Shard>> shards_;
};

void usage(const char *progName, const std::string &error = std::string()) {
  if (!error.empty())
    std::cerr << error << "\n\n";

  std::cerr << "Usage:\n";
  std::cerr
      << "    " << progName
      << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
         "    [ --threads N ] [ --iddfs-iterations N --iddfs-step D ]\n"
         "    [ --iddfs-start N ] [ --iddfs-final-threads N ] <isosig>\n\n"
      << "    " << progName
      << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
         "    [ --threads N ] [ --thicken-layers N ] [ --cone | --no-cone ]\n"
         "    [ --collar-layers N ]\n"
         "    [ --iddfs-iterations N --iddfs-step D ] [ --iddfs-start N ]\n"
         "    [ --iddfs-final-threads N ] --pd <pdcode>\n\n"
      << "    " << progName << " [ -v, --version | -h, --help ]\n\n";
  std::cerr
      << "    -a, --all      : Find all embedded submanifolds, regardless of "
         "boundary\n"
         "                     conditions (default)\n";
  std::cerr << "    -c, --closed   : Find only closed embedded submanifolds\n";
  std::cerr << "    -p, --proper   : Find embedded submanifolds whose "
               "boundary is contained\n"
               "                     entirely in the boundary of the ambient "
               "triangulation\n";
  std::cerr << "    --connected    : Same as --proper, but additionally "
               "require at most one\n"
               "                     boundary component of the embedded "
               "submanifold per\n"
               "                     boundary component of the ambient "
               "triangulation\n\n";
  std::cerr << "    --threads N    : Number of worker threads to search "
               "with (default: the\n"
               "                     number of hardware threads available)\n"
               "\n";
  std::cerr
      << "    --iddfs-iterations N : Run N capped \"fast sweep\" passes "
         "over every\n"
         "                     search root before the final, unbounded "
         "pass, so\n"
         "                     shallow (few-face) results are found "
         "promptly even\n"
         "                     when some roots' full search trees are very "
         "deep or\n"
         "                     effectively unbounded -- pass i is capped "
         "at\n"
         "                     --iddfs-start + (i - 1) * --iddfs-step "
         "faces. Requires\n"
         "                     --iddfs-step. Default: 0 (a single "
         "unbounded pass, as\n"
         "                     before this option existed). Too small a "
         "step wastes\n"
         "                     the capped passes (they find almost "
         "nothing new each\n"
         "                     time); too large a step makes each capped "
         "pass itself\n"
         "                     slow -- a good starting point is found "
         "empirically by\n"
         "                     watching how fast \"roots completed\" "
         "climbs through\n"
         "                     the first pass at a trial value.\n";
  std::cerr << "    --iddfs-step D : Face-count increment per "
               "--iddfs-iterations pass\n"
               "                     after the first. Only used when "
               "--iddfs-iterations >\n"
               "                     0.\n";
  std::cerr
      << "    --iddfs-start N : Face-count cap for the first "
         "--iddfs-iterations pass\n"
         "                     (default: --iddfs-step, i.e. pass i is "
         "capped at i *\n"
         "                     --iddfs-step, as before this option "
         "existed). Useful\n"
         "                     to start the very first pass smaller (or "
         "larger) than\n"
         "                     the step between later passes. Only used "
         "when\n"
         "                     --iddfs-iterations > 0.\n";
  std::cerr
      << "    --iddfs-final-threads N : Number of threads to use for the "
         "final,\n"
         "                     unbounded pass after --iddfs-iterations "
         "capped\n"
         "                     passes (default: --threads). Since every "
         "capped\n"
         "                     pass is bounded and fast regardless of "
         "thread count,\n"
         "                     this only needs to be lower than --threads "
         "if the\n"
         "                     final unbounded pass is the one you want to "
         "throttle,\n"
         "                     e.g. to run longer with less memory "
         "pressure.\n\n";
  std::cerr
      << "    -o, --output <path> : Also write one CSV row per found "
         "surface to <path>,\n"
         "                     alongside the usual progress/summary "
         "output. Columns are\n"
         "                     orientable,genus,punctures,triangles,"
         "condition, plus a\n"
         "                     boundary column when condition is --proper "
         "or --connected.\n\n";
  std::cerr
      << "    --no-simplify  : Skip Triangulation<3>::simplify() when "
         "building a\n"
         "                     boundary complement (see linkcomplement.h). "
         "Yields\n"
         "                     much larger triangulations whose isomorphism "
         "signatures\n"
         "                     rarely recur, so the recognition cache and "
         "census\n"
         "                     lookups become far less effective -- a "
         "profiling/\n"
         "                     diagnostic option, trading recognition "
         "quality for a\n"
         "                     possibly much cheaper boundary-processing "
         "phase.\n\n";
  std::cerr
      << "    --pending-surface-cap N : Cap on how many found surfaces may "
         "wait\n"
         "                     for boundary-link processing at once "
         "(default:\n"
         "                     500000). Past this, DFS worker threads fully "
         "drain\n"
         "                     the backlog themselves (down to empty, not "
         "just\n"
         "                     back under cap) before continuing the "
         "search --\n"
         "                     bounds memory during a long --proper/"
         "--connected\n"
         "                     run at some throughput cost.\n";
  std::cerr
      << "    --petal-cache-limit N : Distinct-petal entry-count threshold "
         "past\n"
         "                     which the shared local-flatness/self-"
         "intersection\n"
         "                     memoization cache clears itself entirely "
         "(default:\n"
         "                     2000000). This cache grows with total DFS "
         "candidates\n"
         "                     examined, not just found surfaces, so it's "
         "usually\n"
         "                     the largest single contributor to unbounded "
         "memory\n"
         "                     growth on a long search.\n";
  std::cerr
      << "    --recognition-cache-limit N : Distinct-isoSig entry-count "
         "threshold\n"
         "                     past which the boundary-complement "
         "recognition cache\n"
         "                     (see linkcomplement.h) clears itself "
         "entirely\n"
         "                     (default: 200000).\n";
  std::cerr
      << "    --boundary-signature-cache-limit N : As above, for each "
         "ambient\n"
         "                     boundary component's pre-triangulation "
         "signature\n"
         "                     cache (default: 200000, per component).\n";
  std::cerr
      << "    --boundary-tally-cap N : Distinct-boundary-descriptor cap "
         "past which\n"
         "                     the live boundary tally stops recording "
         "brand-new\n"
         "                     descriptors (default: 1000000; surfaces with "
         "an\n"
         "                     already-seen boundary are still tallied). "
         "This is a\n"
         "                     reporting structure only, so nothing else "
         "depends on\n"
         "                     it.\n\n";
  std::cerr
      << "    <isosig>       : Isomorphism signature of a 4-manifold\n"
         "                     triangulation to search directly (default "
         "input mode)\n\n";
  std::cerr
      << "    --pd <pdcode>  : Planar diagram code of a knot/link. Builds a "
         "3-sphere\n"
         "                     triangulation containing it, thickens it "
         "into a\n"
         "                     4-dimensional cobordism, and searches that "
         "instead of\n"
         "                     a bare isosig. Mutually exclusive with "
         "<isosig>.\n";
  std::cerr << "    --thicken-layers N : Number of times to thicken() the "
               "cobordism\n"
               "                     (default: 1; only valid with --pd)\n";
  std::cerr
      << "    --cone, --no-cone      : Whether to cap the cobordism with "
         "cone()\n"
         "                     (default: --cone; only valid with --pd)\n";
  std::cerr
      << "    --collar-layers N : Number of thickening layers the search is "
         "seeded\n"
         "                     with a collar traced through the link's "
         "edges,\n"
         "                     starting from the base link (default: 1; 0 "
         "disables\n"
         "                     the collar seed entirely; only valid with "
         "--pd, and\n"
         "                     cannot exceed --thicken-layers). Layers "
         "beyond N are\n"
         "                     still thickened, just not covered by the "
         "collar, which\n"
         "                     is useful when it's not best for the collar "
         "to be\n"
         "                     contained in every layer.\n\n";
  std::cerr
      << "    -v, --version  : Show which version of Regina is being used\n";
  std::cerr << "    -h, --help     : Display this help\n";
  exit(1);
}

// Runs a SurfaceSearch over `tri`, seeded with `seedFaces` if non-empty
// (an empty seed is not a valid seeded-search argument, so this falls back
// to an unseeded search instead), and reports the same way regardless of
// which input mode built `tri`. All output happens here, via the callbacks
// passed to search() -- SurfaceSearch itself has no opinion on how (or
// whether) progress/results get displayed.
void runSearch(const regina::Triangulation<4> &tri,
              const std::vector<int> &seedFaces, BoundaryCondition cond,
              unsigned numThreads,
              const std::optional<std::string> &outputPath,
              unsigned iddfsIterations, long long iddfsStep,
              std::optional<long long> iddfsStart,
              std::optional<unsigned> iddfsFinalThreads,
              const SurfaceSearchLimits &limits) {
  std::cerr << "[+] Running with " << numThreads
            << " threads, condition = " << boundaryConditionName(cond)
            << "\n\n";

  // SurfaceSearch has no default constructor, and the seeded/unseeded
  // constructors take different arguments -- an optional, emplaced once
  // below, lets both branches share the callback wiring that follows
  // instead of duplicating it per branch.
  std::optional<SurfaceSearch> eOpt;
  if (seedFaces.empty())
    eOpt.emplace(tri);
  else
    eOpt.emplace(tri, seedFaces);
  SurfaceSearch &e = *eOpt;
  e.configureLimits(limits);

  const bool wantLinks = cond == BoundaryCondition::proper ||
                         cond == BoundaryCondition::connected;

  // Only bounded for the once-a-second progress report -- keeps its line
  // count from growing without limit over a long search that turns up many
  // distinct boundary types. The final summary (onSearchComplete) always
  // passes std::nullopt, listing every boundary found.
  static constexpr size_t ROLLING_BOUNDARY_LIMIT = 10;

  // Shared between the once-a-second progress report and the final summary:
  // the boundary-link tally when it's being tracked at all, plus the
  // surface-type tally.
  auto tallyText = [&](std::optional<size_t> maxRecentBoundaries =
                          std::nullopt) {
    std::string out;
    if (wantLinks)
      out += e.linkTally().summary(maxRecentBoundaries);
    out += e.surfaceTypeTally().summary();
    return out;
  };

  // How much recomputation the shared PetalCache (see vertexlinks.h) is
  // actually avoiding -- one cache shared across every search thread, so
  // this is already the full aggregate, not a per-thread slice.
  auto petalCacheText = [&] {
    PetalCache::Stats s = e.petalCacheStats();
    auto hitRate = [](long long hits, long long checks) {
      return checks > 0 ? 100.0 * static_cast<double>(hits) /
                               static_cast<double>(checks)
                        : 0.0;
    };
    std::ostringstream out;
    out << "[+] petal cache: isUnknot checks=" << s.unknotChecks
        << " hits=" << s.unknotCacheHits << " (" << std::fixed
        << std::setprecision(1) << hitRate(s.unknotCacheHits, s.unknotChecks)
        << "% hit rate), linking checks=" << s.linkingChecks
        << " hits=" << s.linkingCacheHits << " (" << std::fixed
        << std::setprecision(1)
        << hitRate(s.linkingCacheHits, s.linkingChecks) << "% hit rate)\n";
    out << "[+] rejected for non-local-flatness: " << s.localFlatnessRejections
        << " | rejected for a transverse self-intersection: "
        << s.transverseRejections << "\n";
    out << "[+] petal cache entries (distinct petals seen): "
        << e.petalCacheSize() << " | full resets: " << s.cacheResets << "\n";
    return out.str();
  };

  // How much recomputation the isoSig-keyed recognition cache (see
  // linkcomplement.h/.cpp) is actually avoiding -- shared process-wide, so
  // this is already the full aggregate across every search thread.
  auto recognitionCacheText = [] {
    RecognitionCacheStats s = recognitionCacheStats();
    auto hitRate = [](long long hits, long long checks) {
      return checks > 0 ? 100.0 * static_cast<double>(hits) /
                               static_cast<double>(checks)
                        : 0.0;
    };
    std::ostringstream out;
    out << "[+] recognition cache: genus checks=" << s.genusChecks
        << " hits=" << s.genusCacheHits << " (" << std::fixed
        << std::setprecision(1) << hitRate(s.genusCacheHits, s.genusChecks)
        << "% hit rate), census checks=" << s.censusChecks
        << " hits=" << s.censusCacheHits << " (" << std::fixed
        << std::setprecision(1)
        << hitRate(s.censusCacheHits, s.censusChecks) << "% hit rate)\n";
    out << "[+] recognition cache entries (distinct isoSigs seen): "
        << recognitionCacheSize() << " | full resets: " << s.cacheResets
        << "\n";
    long long misses = s.genusChecks - s.genusCacheHits;
    out << "[+] genus cache misses resolved via: group-is-Z check="
        << s.groupFastPathHits << ", SnapPea hyperbolicity check="
        << s.snapPeaFastPathHits << ", recogniseHandlebody() fallback="
        << s.recogniseHandlebodyFallbacks << " (of " << misses
        << " misses)\n";
    return out.str();
  };

  // How much recomputation the pre-triangulation boundary-signature cache
  // (see BoundarySignatureCache in linkcomplement.h/.cpp) is actually
  // avoiding -- one cache per ambient boundary component, shared across
  // every search thread, so this is already the full aggregate. A high hit
  // rate here means most boundary curves never reach
  // EdgeComplement::identify()'s buildComplement()/simplify()/isoSig() at
  // all, let alone the recognition cache or Census::lookup() above.
  auto boundarySignatureCacheText = [&] {
    BoundarySignatureCacheStats s = e.boundarySignatureCacheStats();
    double hitRate = s.checks > 0 ? 100.0 * static_cast<double>(s.hits) /
                                        static_cast<double>(s.checks)
                                  : 0.0;
    std::ostringstream out;
    out << "[+] boundary signature cache: checks=" << s.checks
        << " hits=" << s.hits << " (" << std::fixed << std::setprecision(1)
        << hitRate << "% hit rate)\n";
    out << "[+] boundary signature cache entries (distinct canonical "
           "boundary signatures seen): "
        << e.boundarySignatureCacheSize() << " | full resets: "
        << s.cacheResets << "\n";
    return out.str();
  };

  // How saturated the surface-boundary-processing queue (pendingSurfaces_)
  // is/has been -- directly shows whether the queue's cap and the "DFS
  // worker threads help drain under pressure" backpressure mechanism (see
  // SurfaceSearch::ThreadHook::onFlush()) are actually engaging on this run.
  auto pendingSurfaceQueueText = [&] {
    PendingSurfaceQueueStats s = e.pendingSurfaceQueueStats();
    std::ostringstream out;
    out << "[+] pending surface queue: size=" << s.currentSize << "/"
        << s.cap << " (peak=" << s.peakSize << "), producer-helped drains="
        << s.producerDrainEvents << " (" << s.producerDrainedEntries
        << " entries)\n";
    return out.str();
  };

  RollingReport searchReport;
  RollingReport boundaryReport;
  SurfaceSearchCallbacks callbacks;

  // Rate-tracking for the isEmbedded()-true counter: the embeddedCount and
  // elapsed-duration seen on the previous onProgress tick, diffed against
  // the current tick to compute a live "X/sec" rate. Both start at zero,
  // which correctly reports rate=0 on the very first tick.
  long long lastEmbeddedCount = 0;
  std::chrono::steady_clock::duration lastElapsed{};

  callbacks.onProgress = [&](const SearchStats &stats) {
    std::ostringstream report;
    report << tallyText(ROLLING_BOUNDARY_LIMIT);
    report << "[+] elapsed: " << formatElapsed(stats.elapsed) << "\n";
    report << "[+] roots completed: " << stats.rootsCompleted << "/"
           << stats.totalRoots << "  | candidates examined so far: "
           << stats.foundCount << "\n";
    if (stats.iddfsTotalRounds > 1) {
      report << "[+] IDDFS pass " << stats.iddfsRound << "/"
             << stats.iddfsTotalRounds;
      if (stats.iddfsCapped)
        report << " (capped at " << stats.iddfsCap << " faces)\n";
      else
        report << " (final, unbounded)\n";
    }

    std::chrono::duration<double> tickDelta = stats.elapsed - lastElapsed;
    double embeddedRate =
        tickDelta.count() > 0.0
            ? static_cast<double>(stats.embeddedCount - lastEmbeddedCount) /
                  tickDelta.count()
            : 0.0;
    report << "[+] embedded submanifolds found so far: " << stats.embeddedCount
           << " (" << std::fixed << std::setprecision(1) << embeddedRate
           << "/sec)\n";
    lastEmbeddedCount = stats.embeddedCount;
    lastElapsed = stats.elapsed;

    report << "[+] embedded submanifolds satisfying boundary condition ("
           << boundaryConditionName(cond)
           << ") found so far: " << stats.satisfyingCount << "\n";
    report << "[+] largest embedded submanifold satisfying boundary "
              "condition ("
           << boundaryConditionName(cond)
           << ") found so far: " << stats.largestSatisfying << " faces\n";
    report << "[+] average faces per embedded submanifold satisfying "
              "boundary condition ("
           << boundaryConditionName(cond) << ") found so far: " << std::fixed
           << std::setprecision(2) << stats.averageSatisfyingFaces() << "\n";
    report << petalCacheText();
    if (wantLinks)
      report << pendingSurfaceQueueText();
    searchReport.draw(report.str());
  };

  callbacks.onInterrupted = [] {
    std::cerr << "\n[!] Interrupted -- moving on to whatever comes next "
                 "with what was found so far (Ctrl+C again to quit "
                 "immediately)\n";
  };

  callbacks.onQueueDrainPause = [&](size_t queueSize, size_t cap) {
    // commitLine(), not a bare std::cerr print: this fires from a worker
    // thread at an arbitrary time relative to the reporter thread's own
    // once-a-second searchReport.draw() calls -- a plain print here would
    // (a) risk interleaving with a concurrent draw() and (b) get silently
    // erased by the very next draw(), since that call has no way to know
    // an unrelated line was printed in between its own redraws.
    std::ostringstream msg;
    msg << "\n[*] Pending surface queue over cap (" << queueSize << "/" << cap
        << "); pausing the search to fully drain it before continuing "
           "(Ctrl+C to interrupt)...\n";
    searchReport.commitLine(msg.str());
  };

  callbacks.onQueueDrainResume = [&](bool interrupted, size_t remaining) {
    std::ostringstream msg;
    if (interrupted)
      msg << "[*] Interrupted while draining -- handing off " << remaining
          << " remaining surface" << (remaining == 1 ? "" : "s")
          << " to post-search processing instead\n";
    else
      msg << "[*] Queue fully drained; resuming search.\n";
    searchReport.commitLine(msg.str());
  };

  callbacks.onSearchComplete = [&](const SearchStats &stats) {
    std::string extra = tallyText();
    if (!extra.empty())
      std::cerr << "\n" << extra;
    std::cerr << "\nTotal elapsed: " << formatElapsed(stats.elapsed) << "\n";
    std::cerr << "Total candidates examined: " << stats.foundCount
              << " (across " << numThreads << " threads)\n";

    std::chrono::duration<double> totalElapsedSecs = stats.elapsed;
    double avgEmbeddedRate = totalElapsedSecs.count() > 0.0
                                 ? static_cast<double>(stats.embeddedCount) /
                                       totalElapsedSecs.count()
                                 : 0.0;
    std::cerr << "Total embedded submanifolds found: " << stats.embeddedCount
              << " (across " << numThreads << " threads, " << std::fixed
              << std::setprecision(1) << avgEmbeddedRate << "/sec average)\n";

    std::cerr << "Total embedded submanifolds (" << boundaryConditionName(cond)
              << "): " << stats.satisfyingCount << " (across " << numThreads
              << " threads)\n";
    std::cerr << "Largest embedded submanifold (" << boundaryConditionName(cond)
              << "): " << stats.largestSatisfying << " faces\n";
    std::cerr << "Average faces per embedded submanifold ("
              << boundaryConditionName(cond) << "): " << std::fixed
              << std::setprecision(2) << stats.averageSatisfyingFaces()
              << "\n";
    std::cerr << petalCacheText();
    if (wantLinks)
      std::cerr << pendingSurfaceQueueText();
    std::cerr << recognitionCacheText();
    std::cerr << boundarySignatureCacheText();
  };

  callbacks.onBoundaryProcessingStarted = [](size_t total,
                                             unsigned workerCount) {
    std::cerr << "\n[*] Search complete; processing " << total
              << " remaining surface boundaries with " << workerCount
              << " threads...\n";
  };

  callbacks.onBoundaryProcessingProgress =
      [&](size_t processed, size_t total,
          std::chrono::steady_clock::duration elapsed) {
        long long elapsedSeconds =
            std::chrono::duration_cast<std::chrono::seconds>(elapsed)
                .count();

        std::ostringstream report;
        report << e.linkTally().summary(ROLLING_BOUNDARY_LIMIT);
        report << "[+] elapsed: " << formatElapsed(elapsed) << "\n";
        report << "[+] boundaries processed: " << processed << "/" << total
               << " (" << std::fixed << std::setprecision(2)
               << (total > 0 ? 100.0 * static_cast<double>(processed) /
                                   static_cast<double>(total)
                             : 0.0)
               << "%)\n";
        if (processed > 0 && processed < total && elapsedSeconds > 0) {
          long long etaSeconds = elapsedSeconds *
                                 static_cast<long long>(total - processed) /
                                 static_cast<long long>(processed);
          report << "[+] ETA: "
                 << formatElapsed(std::chrono::seconds(etaSeconds)) << "\n";
        }
        report << boundarySignatureCacheText();
        report << recognitionCacheText();
        boundaryReport.draw(report.str());
      };

  callbacks.onBoundaryProcessingComplete =
      [&](size_t total, std::chrono::steady_clock::duration elapsed) {
        std::chrono::duration<double> elapsedSecs = elapsed;
        double avgRate = elapsedSecs.count() > 0.0
                             ? static_cast<double>(total) / elapsedSecs.count()
                             : 0.0;
        std::cerr << "\n[+] Finished processing " << total
                  << " remaining surface boundaries in "
                  << formatElapsed(elapsed) << " (" << std::fixed
                  << std::setprecision(2) << avgRate
                  << " boundaries/sec average)\n";
        std::cerr << recognitionCacheText();
        std::cerr << boundarySignatureCacheText();
      };

  std::optional<CsvWriter> writer;
  if (outputPath) {
    writer.emplace(
        *outputPath,
        wantLinks ? "orientable,genus,punctures,triangles,condition,boundary"
                 : "orientable,genus,punctures,triangles,condition",
        numThreads);

    callbacks.onSurfaceFound = [&](const SurfaceFoundInfo &info) {
      writer->writeRow(csvRow(info.orientable, info.genus, info.punctures,
                              info.triangleCount, info.mostRestrictive));
    };
    callbacks.onSurfaceBoundaryProcessed =
        [&](const SurfaceBoundaryInfo &info) {
          writer->writeRow(csvRow(info.orientable, info.genus, info.punctures,
                                  info.triangleCount, info.mostRestrictive) +
                           "," + csvField(info.boundaryDescription));
        };
  }

  e.search(numThreads, cond, callbacks, iddfsIterations, iddfsStep,
           iddfsStart, iddfsFinalThreads);

  if (writer)
    writer->finalize();
}
} // namespace

int main(int argc, char *argv[]) {
  BoundaryCondition cond = BoundaryCondition::all;
  unsigned numThreads = std::thread::hardware_concurrency();
  if (numThreads == 0)
    numThreads = 1;

  std::string isoSig;
  bool haveIsoSig = false;

  std::string pdCode;
  bool havePD = false;

  std::optional<std::string> outputPath;

  int thickenLayers = 1;
  bool useCone = true;
  int collarLayers = 1;
  bool sawThickenLayers = false, sawCone = false, sawCollarLayers = false;

  unsigned iddfsIterations = 0;
  long long iddfsStep = 0;
  std::optional<long long> iddfsStart;
  std::optional<unsigned> iddfsFinalThreads;

  SurfaceSearchLimits limits;
  size_t recognitionCacheLimitArg = recognitionCacheLimit.load();

  for (int i = 1; i < argc; ++i) {
    std::string arg = argv[i];
    if (arg == "-h" || arg == "--help") {
      usage(argv[0]);
    } else if (arg == "-v" || arg == "--version") {
      if (argc != 2)
        usage(argv[0], "Option --version cannot be used with "
                       "any other arguments.");
      std::cout << PACKAGE_BUILD_STRING << "\n";
      return 0;
    } else if (arg == "-a" || arg == "--all") {
      cond = BoundaryCondition::all;
    } else if (arg == "-c" || arg == "--closed") {
      cond = BoundaryCondition::closed;
    } else if (arg == "-p" || arg == "--proper") {
      cond = BoundaryCondition::proper;
    } else if (arg == "--connected") {
      cond = BoundaryCondition::connected;
    } else if (arg == "--threads") {
      if (i + 1 >= argc)
        usage(argv[0], "--threads requires a value.");
      try {
        numThreads = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--threads requires an integer value.");
      }
    } else if (arg == "-o" || arg == "--output") {
      if (i + 1 >= argc)
        usage(argv[0], "-o/--output requires a value.");
      outputPath = argv[++i];
    } else if (arg == "--pd") {
      if (i + 1 >= argc)
        usage(argv[0], "--pd requires a value.");
      pdCode = argv[++i];
      havePD = true;
    } else if (arg == "--thicken-layers") {
      if (i + 1 >= argc)
        usage(argv[0], "--thicken-layers requires a value.");
      try {
        thickenLayers = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--thicken-layers requires an integer value.");
      }
      sawThickenLayers = true;
    } else if (arg == "--cone") {
      useCone = true;
      sawCone = true;
    } else if (arg == "--no-cone") {
      useCone = false;
      sawCone = true;
    } else if (arg == "--collar-layers") {
      if (i + 1 >= argc)
        usage(argv[0], "--collar-layers requires a value.");
      try {
        collarLayers = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--collar-layers requires an integer value.");
      }
      sawCollarLayers = true;
    } else if (arg == "--iddfs-iterations") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-iterations requires a value.");
      try {
        iddfsIterations = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-iterations requires an integer value.");
      }
    } else if (arg == "--iddfs-step") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-step requires a value.");
      try {
        iddfsStep = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-step requires an integer value.");
      }
    } else if (arg == "--iddfs-start") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-start requires a value.");
      try {
        iddfsStart = std::stoll(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-start requires an integer value.");
      }
    } else if (arg == "--iddfs-final-threads") {
      if (i + 1 >= argc)
        usage(argv[0], "--iddfs-final-threads requires a value.");
      try {
        iddfsFinalThreads = static_cast<unsigned>(std::stoul(argv[++i]));
      } catch (const std::exception &) {
        usage(argv[0], "--iddfs-final-threads requires an integer value.");
      }
    } else if (arg == "--no-simplify") {
      simplifyComplements = false;
    } else if (arg == "--pending-surface-cap") {
      if (i + 1 >= argc)
        usage(argv[0], "--pending-surface-cap requires a value.");
      try {
        limits.pendingSurfaceCap = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--pending-surface-cap requires an integer value.");
      }
    } else if (arg == "--petal-cache-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--petal-cache-limit requires a value.");
      try {
        limits.petalCacheLimit = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--petal-cache-limit requires an integer value.");
      }
    } else if (arg == "--recognition-cache-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--recognition-cache-limit requires a value.");
      try {
        recognitionCacheLimitArg = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--recognition-cache-limit requires an integer value.");
      }
    } else if (arg == "--boundary-signature-cache-limit") {
      if (i + 1 >= argc)
        usage(argv[0], "--boundary-signature-cache-limit requires a value.");
      try {
        limits.boundarySignatureCacheLimit = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0],
              "--boundary-signature-cache-limit requires an integer value.");
      }
    } else if (arg == "--boundary-tally-cap") {
      if (i + 1 >= argc)
        usage(argv[0], "--boundary-tally-cap requires a value.");
      try {
        limits.boundaryTallyCap = std::stoull(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--boundary-tally-cap requires an integer value.");
      }
    } else if (!arg.empty() && arg[0] == '-') {
      usage(argv[0], "Unknown option: " + arg);
    } else if (haveIsoSig) {
      usage(argv[0], "Only one isosig may be specified.");
    } else {
      isoSig = arg;
      haveIsoSig = true;
    }
  }

  if (havePD && haveIsoSig)
    usage(argv[0], "--pd and an isosig are mutually exclusive.");
  if (!havePD && !haveIsoSig)
    usage(argv[0], "Please specify an isosig or --pd <pdcode>.");
  if (!havePD && (sawThickenLayers || sawCone || sawCollarLayers))
    usage(argv[0], "--thicken-layers, --cone/--no-cone, and --collar-layers "
                   "are only valid with --pd.");
  if (havePD && collarLayers < 0)
    usage(argv[0], "--collar-layers requires a value >= 0.");
  if (havePD && collarLayers > thickenLayers)
    usage(argv[0], "--collar-layers cannot exceed --thicken-layers.");
  if (iddfsIterations > 0 && iddfsStep <= 0)
    usage(argv[0], "--iddfs-iterations > 0 requires --iddfs-step > 0.");
  if (iddfsStart && *iddfsStart <= 0)
    usage(argv[0], "--iddfs-start requires a value > 0.");
  if (limits.pendingSurfaceCap == 0)
    usage(argv[0], "--pending-surface-cap requires a value > 0.");
  if (limits.petalCacheLimit == 0)
    usage(argv[0], "--petal-cache-limit requires a value > 0.");
  if (recognitionCacheLimitArg == 0)
    usage(argv[0], "--recognition-cache-limit requires a value > 0.");
  if (limits.boundarySignatureCacheLimit == 0)
    usage(argv[0], "--boundary-signature-cache-limit requires a value > 0.");
  if (limits.boundaryTallyCap == 0)
    usage(argv[0], "--boundary-tally-cap requires a value > 0.");
  recognitionCacheLimit.store(recognitionCacheLimitArg,
                             std::memory_order_relaxed);
  if (outputPath) {
    // Fail fast, before running a potentially long search, rather than
    // discovering an unwritable path only once results are ready to flush.
    std::ofstream test(*outputPath, std::ios::app);
    if (!test)
      usage(argv[0], "Cannot open " + *outputPath + " for writing.");
  }

  std::cout << "------ SurFer (Surface Finder) \U0001F30A ------\n\n";

  if (havePD) {
    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdCode);

    knotbuilder::TriangulationWithLink link;
    try {
      link = knotbuilder::buildLink(pdcode);
    } catch (const regina::InvalidArgument &e) {
      usage(argv[0], std::string("Invalid PD code: ") + e.what());
    }
    auto &[t2, edges2] = link;

    // Indices are preserved across CobordismBuilder's internal copy/reorder
    // of t2 (see CobordismBuilder::baseTriangulation()'s doc comment), so
    // edgeIndices (taken from t2) still identify the same edges in
    // cob.baseTriangulation().
    std::vector<int> edgeIndices;
    edgeIndices.reserve(edges2.size());
    for (const regina::Edge<3> *e : edges2)
      edgeIndices.push_back(static_cast<int>(e->index()));

    CobordismBuilder<3> cob(t2);
    CollarBuilder collarBuilder(edgeIndices);
    for (int i = 0; i < thickenLayers; ++i) {
      cob.thicken();
      // Must run every layer the collar is meant to cover, not just the
      // first: CollarBuilder::addLayer only captures the most-recently-built
      // thickening layer's prisms, so tracing the collar all the way from
      // the link's original position up through --collar-layers layers
      // requires calling it once per thicken() up to that point.
      if (i < collarLayers)
        collarBuilder.addLayer(cob);
    }
    if (useCone)
      cob.cone();
    regina::Triangulation<4> tri = cob.getCobordism();

    std::vector<int> seedFaces;
    if (collarLayers > 0) {
      for (regina::Triangle<4> *t : collarBuilder.resolve())
        seedFaces.push_back(static_cast<int>(t->index()));
      std::cerr << "[+] Collar seed faces = " << seedFaces.size() << "\n";
    }

    runSearch(tri, seedFaces, cond, numThreads, outputPath, iddfsIterations,
             iddfsStep, iddfsStart, iddfsFinalThreads, limits);
  } else {
    regina::Triangulation<4> tri;
    try {
      tri = regina::Triangulation<4>(isoSig);
    } catch (const regina::InvalidArgument &e) {
      usage(argv[0], std::string("Invalid isosig: ") + e.what());
    }

    runSearch(tri, {}, cond, numThreads, outputPath, iddfsIterations,
             iddfsStep, iddfsStart, iddfsFinalThreads, limits);
  }

  return 0;
}
