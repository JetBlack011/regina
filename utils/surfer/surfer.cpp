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
#include <string_view>
#include <thread>

#include <unistd.h>

#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "collar.h"
#include "embeddingsearch.h"
#include "knotbuilder.h"

namespace {

// Redraws a block of text in place, using the same ANSI cursor-rewind
// sequence the old inline reporters used, so each call's output overwrites
// the previous one instead of scrolling the terminal.
class RollingReport {
  size_t prevLines_ = 0;

public:
  void draw(const std::string &text) {
    if (prevLines_ > 0)
      std::cerr << "\x1b[" << prevLines_ << "F\x1b[0J";
    std::cerr << text;
    prevLines_ =
        static_cast<size_t>(std::count(text.begin(), text.end(), '\n'));
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

// Buffers CSV rows per calling thread (one shard file per distinct thread
// that ever calls writeRow(), created lazily on first use) so concurrent
// writers -- the DFS workers, the background boundary-link drain thread,
// and the final boundary-batch-processing workers -- never contend on a
// lock or file handle while writing rows. finalize() concatenates every
// shard into the requested output path (header first) and removes the
// shard files -- call it once, single-threaded, only after every thread
// that might call writeRow() has already joined (i.e. after search()
// itself has returned).
//
// shardForThisThread() caches its result in a thread_local, function-static
// pointer -- correct as long as at most one CsvWriter is ever alive per
// process, which is exactly how surfer.cpp uses it (constructed at most
// once per run, in runSearch()).
class CsvWriter {
public:
  CsvWriter(std::filesystem::path outputPath, std::string headerLine)
      : outputPath_(std::move(outputPath)),
        headerLine_(std::move(headerLine)) {}

  void writeRow(const std::string &row) {
    Shard &shard = shardForThisThread();
    shard.buffer += row;
    shard.buffer += '\n';
    if (shard.buffer.size() >= FLUSH_THRESHOLD)
      flush(shard);
  }

  void finalize() {
    std::ofstream out(outputPath_, std::ios::trunc);
    out << headerLine_ << "\n";
    for (auto &shard : shards_) {
      flush(*shard);
      shard->file.close();
      std::ifstream in(shard->path);
      out << in.rdbuf();
      in.close();
      std::filesystem::remove(shard->path);
    }
    std::cerr << "[+] Wrote CSV output to " << outputPath_.string() << "\n";
  }

private:
  struct Shard {
    std::filesystem::path path;
    std::ofstream file;
    std::string buffer;
  };

  static constexpr size_t FLUSH_THRESHOLD = 1 << 16; // ~64KB

  static void flush(Shard &shard) {
    if (shard.buffer.empty())
      return;
    shard.file << shard.buffer;
    shard.buffer.clear();
  }

  Shard &shardForThisThread() {
    static thread_local Shard *cached = nullptr;
    if (cached)
      return *cached;
    std::lock_guard<std::mutex> lock(registryMutex_);
    auto shard = std::make_unique<Shard>();
    shard->path = outputPath_.string() + ".tmp." +
                 std::to_string(getpid()) + ".shard" +
                 std::to_string(shards_.size());
    shard->file.open(shard->path, std::ios::trunc);
    shards_.push_back(std::move(shard));
    cached = shards_.back().get();
    return *cached;
  }

  std::filesystem::path outputPath_;
  std::string headerLine_;
  std::mutex registryMutex_;
  std::vector<std::unique_ptr<Shard>> shards_;
};

void usage(const char *progName, const std::string &error = std::string()) {
  if (!error.empty())
    std::cerr << error << "\n\n";

  std::cerr << "Usage:\n";
  std::cerr
      << "    " << progName
      << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
         "    [ --threads N ] <isosig>\n\n"
      << "    " << progName
      << " [ -a, --all | -c, --closed | -p, --proper | --connected ]\n"
         "    [ --threads N ] [ --layers N ] [ --cone | --no-cone ]\n"
         "    [ --collar | --no-collar ] --pd <pdcode>\n\n"
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
      << "    -o, --output <path> : Also write one CSV row per found "
         "surface to <path>,\n"
         "                     alongside the usual progress/summary "
         "output. Columns are\n"
         "                     orientable,genus,punctures,triangles,"
         "condition, plus a\n"
         "                     boundary column when condition is --proper "
         "or --connected.\n\n";
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
  std::cerr << "    --layers N     : Number of times to thicken() the "
               "cobordism (default: 1;\n"
               "                     only valid with --pd)\n";
  std::cerr
      << "    --cone, --no-cone      : Whether to cap the cobordism with "
         "cone()\n"
         "                     (default: --cone; only valid with --pd)\n";
  std::cerr
      << "    --collar, --no-collar  : Whether to seed the search with a "
         "collar\n"
         "                     traced through the link's edges (default: "
         "--collar;\n"
         "                     only valid with --pd)\n\n";
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
              const std::optional<std::string> &outputPath) {
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

  const bool wantLinks = cond == BoundaryCondition::proper ||
                         cond == BoundaryCondition::connected;

  // Shared between the once-a-second progress report and the final summary:
  // the surface-type tally, plus the boundary-link tally when it's being
  // tracked at all.
  auto tallyText = [&] {
    std::string out = e.surfaceTypeTally().summary();
    if (wantLinks)
      out += e.linkTally().summary();
    return out;
  };

  RollingReport searchReport;
  RollingReport boundaryReport;
  SurfaceSearchCallbacks callbacks;

  callbacks.onProgress = [&](const SearchStats &stats) {
    std::ostringstream report;
    report << tallyText();
    report << "[+] elapsed: " << formatElapsed(stats.elapsed) << "\n";
    report << "[+] roots completed: " << stats.rootsCompleted << "/"
           << stats.totalRoots << "  | embedded submanifolds found so far: "
           << stats.foundCount << "\n";
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
    searchReport.draw(report.str());
  };

  callbacks.onInterrupted = [] {
    std::cerr << "\n[!] Interrupted -- moving on to whatever comes next "
                 "with what was found so far (Ctrl+C again to quit "
                 "immediately)\n";
  };

  callbacks.onSearchComplete = [&](const SearchStats &stats) {
    std::string extra = tallyText();
    if (!extra.empty())
      std::cerr << "\n" << extra;
    std::cerr << "\nTotal elapsed: " << formatElapsed(stats.elapsed) << "\n";
    std::cerr << "Total embedded submanifolds found: " << stats.foundCount
              << " (across " << numThreads << " threads)\n";
    std::cerr << "Total embedded submanifolds (" << boundaryConditionName(cond)
              << "): " << stats.satisfyingCount << " (across " << numThreads
              << " threads)\n";
    std::cerr << "Largest embedded submanifold (" << boundaryConditionName(cond)
              << "): " << stats.largestSatisfying << " faces\n";
    std::cerr << "Average faces per embedded submanifold ("
              << boundaryConditionName(cond) << "): " << std::fixed
              << std::setprecision(2) << stats.averageSatisfyingFaces()
              << "\n";
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
        report << e.linkTally().summary();
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
        boundaryReport.draw(report.str());
      };

  callbacks.onBoundaryProcessingComplete =
      [](std::chrono::steady_clock::duration elapsed) {
        std::cerr << "\n[+] Finished processing remaining surface "
                     "boundaries in "
                  << formatElapsed(elapsed) << "\n";
      };

  std::optional<CsvWriter> writer;
  if (outputPath) {
    writer.emplace(
        *outputPath,
        wantLinks ? "orientable,genus,punctures,triangles,condition,boundary"
                 : "orientable,genus,punctures,triangles,condition");

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

  e.search(numThreads, cond, callbacks);

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

  int layers = 1;
  bool useCone = true;
  bool useCollar = true;
  bool sawLayers = false, sawCone = false, sawCollar = false;

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
    } else if (arg == "--layers") {
      if (i + 1 >= argc)
        usage(argv[0], "--layers requires a value.");
      try {
        layers = std::stoi(argv[++i]);
      } catch (const std::exception &) {
        usage(argv[0], "--layers requires an integer value.");
      }
      sawLayers = true;
    } else if (arg == "--cone") {
      useCone = true;
      sawCone = true;
    } else if (arg == "--no-cone") {
      useCone = false;
      sawCone = true;
    } else if (arg == "--collar") {
      useCollar = true;
      sawCollar = true;
    } else if (arg == "--no-collar") {
      useCollar = false;
      sawCollar = true;
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
  if (!havePD && (sawLayers || sawCone || sawCollar))
    usage(argv[0], "--layers, --cone/--no-cone, and --collar/--no-collar "
                   "are only valid with --pd.");
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
    for (int i = 0; i < layers; ++i) {
      cob.thicken();
      // Must run every layer, not just the first: CollarBuilder::addLayer
      // only captures the most-recently-built thickening layer's prisms, so
      // tracing the collar all the way from the link's original position up
      // through every layer requires calling it once per thicken().
      if (useCollar)
        collarBuilder.addLayer(cob);
    }
    if (useCone)
      cob.cone();
    regina::Triangulation<4> tri = cob.getCobordism();

    std::vector<int> seedFaces;
    if (useCollar) {
      for (regina::Triangle<4> *t : collarBuilder.resolve())
        seedFaces.push_back(static_cast<int>(t->index()));
      std::cerr << "[+] Collar seed faces = " << seedFaces.size() << "\n";
    }

    runSearch(tri, seedFaces, cond, numThreads, outputPath);
  } else {
    regina::Triangulation<4> tri;
    try {
      tri = regina::Triangulation<4>(isoSig);
    } catch (const regina::InvalidArgument &e) {
      usage(argv[0], std::string("Invalid isosig: ") + e.what());
    }

    runSearch(tri, {}, cond, numThreads, outputPath);
  }

  return 0;
}
