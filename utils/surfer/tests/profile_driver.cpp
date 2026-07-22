//
//  profile_driver.cpp
//
//  Standalone driver for profiling EmbeddingSearch<dim,2>/SurfaceSearch
//  (embeddingsearch.h) against a range of target triangulations. Not part
//  of CTest -- build and run manually, typically under `perf record`:
//
//    ./profile_driver bary3 <subdivisions> <cond> <threads>
//    ./profile_driver bary4 <subdivisions> <cond> <threads>
//    ./profile_driver isosig4 <isosig> <cond> <threads>
//    ./profile_driver pipeline <pdcode> <thickenLayers> <cond> <threads>
//
//  cond in {all, closed, proper, connected}. The reporter thread's periodic
//  stderr output is suppressed so it doesn't pollute perf's view of steady
//  -state work.

#include <chrono>
#include <iostream>
#include <sstream>
#include <string>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>

#include "cobordismbuilder.h"
#include "embeddingsearch.h"
#include "knotbuilder.h"

namespace {

BoundaryCondition parseCond(const std::string &s) {
  if (s == "all")
    return BoundaryCondition::all;
  if (s == "closed")
    return BoundaryCondition::closed;
  if (s == "proper")
    return BoundaryCondition::proper;
  if (s == "connected")
    return BoundaryCondition::connected;
  throw std::runtime_error("Unknown condition: " + s);
}

// Silences search()'s once-a-second progress reports (written to std::cerr)
// for the duration of the call, so perf's samples aren't spent on
// terminal-escape formatting rather than the actual search.
template <typename Search>
void runSilently(Search &search, unsigned threads, BoundaryCondition cond) {
  std::ostringstream sink;
  std::streambuf *old = std::cerr.rdbuf(sink.rdbuf());
  const auto start = std::chrono::steady_clock::now();
  search.search(threads, cond);
  const auto elapsed = std::chrono::steady_clock::now() - start;
  std::cerr.rdbuf(old);
  std::cerr << "[driver] wall time: "
            << std::chrono::duration_cast<std::chrono::milliseconds>(elapsed)
                   .count()
            << " ms\n";
}

// DIAGNOSTIC: exposes EmbeddingSearch<4,2>'s protected skeleton_/graph_ and
// EmbeddednessPredicate so we can drive the DFS single-threaded ourselves and
// count, at fine grain, how many visited nodes actually satisfy a given
// BoundaryCondition (i.e. how often surfaceTypeKey()/boundaryComponentsMap-
// Injectively() actually fire) -- independent of embeddingsearch.cpp's own
// FLUSH_EVERY-batched reporting, which is too coarse-grained when the hit
// rate is low.
class DiagSearch : public EmbeddingSearch<4, 2> {
public:
  using EmbeddingSearch<4, 2>::EmbeddingSearch;

  void diagRun(BoundaryCondition cond, long long reportEvery) {
    ConnectedInducedSubgraphEnumerator enumerator(graph_.adjList.first,
                                                   graph_.adjList.second);
    EmbeddedSubmanifold<4, 2> embedding(skeleton_);
    EmbeddednessPredicate predicate(embedding, graph_.graphToSkel);

    long long visited = 0, satisfying = 0;
    const auto start = std::chrono::steady_clock::now();

    for (int s = 1; s <= graph_.adjList.first; ++s) {
      enumerator.enumerateFromRootFiltered(
          s,
          [&](const std::vector<int> &) {
            ++visited;
            if (embedding.satisfies(cond))
              ++satisfying;
            if (visited % reportEvery == 0) {
              double pct = 100.0 * static_cast<double>(satisfying) /
                           static_cast<double>(visited);
              std::cerr << "[diag] elapsed="
                        << formatElapsed(std::chrono::steady_clock::now() -
                                         start)
                        << " visited=" << visited
                        << " satisfying=" << satisfying << " (" << pct
                        << "%)\n";
            }
          },
          predicate);
    }
  }
};

} // namespace

int main(int argc, char *argv[]) {
  if (argc < 2) {
    std::cerr << "Usage:\n"
                 "  profile_driver bary3 <subdivisions> <cond> <threads>\n"
                 "  profile_driver bary4 <subdivisions> <cond> <threads>\n"
                 "  profile_driver isosig4 <isosig> <cond> <threads>\n"
                 "  profile_driver pipeline <pdcode> <thickenLayers> <cond> "
                 "<threads>\n"
                 "  profile_driver diag <pdcode> <thickenLayers> <cond> "
                 "<reportEvery>\n";
    return 1;
  }

  std::string mode = argv[1];

  if (mode == "bary3") {
    int subdivisions = std::stoi(argv[2]);
    BoundaryCondition cond = parseCond(argv[3]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[4]));

    regina::Triangulation<3> tri;
    tri.newSimplex();
    for (int i = 0; i < subdivisions; ++i)
      tri.subdivide();

    std::cerr << "[driver] dim3 tetrahedra = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    EmbeddingSearch<3, 2> search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    runSilently(search, threads, cond);
  } else if (mode == "bary4") {
    int subdivisions = std::stoi(argv[2]);
    BoundaryCondition cond = parseCond(argv[3]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[4]));

    regina::Triangulation<4> tri;
    tri.newSimplex();
    for (int i = 0; i < subdivisions; ++i)
      tri.subdivide();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    runSilently(search, threads, cond);
  } else if (mode == "isosig4") {
    std::string isosig = argv[2];
    BoundaryCondition cond = parseCond(argv[3]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[4]));

    regina::Triangulation<4> tri(isosig);
    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    runSilently(search, threads, cond);
  } else if (mode == "pipeline" || mode == "pipeline-loud") {
    std::string pdcodeStr = argv[2];
    int layers = std::stoi(argv[3]);
    BoundaryCondition cond = parseCond(argv[4]);
    unsigned threads = static_cast<unsigned>(std::stoi(argv[5]));

    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcodeStr);
    auto [t2, edges0] = knotbuilder::buildLink(pdcode);

    CobordismBuilder<3> cob(t2);
    if (layers > 0)
      cob.thicken(layers);
    regina::Triangulation<4> tri = cob.cone();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    SurfaceSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    if (mode == "pipeline-loud")
      search.search(threads, cond); // don't suppress stderr progress reports
    else
      runSilently(search, threads, cond);
  } else if (mode == "diag") {
    std::string pdcodeStr = argv[2];
    int layers = std::stoi(argv[3]);
    BoundaryCondition cond = parseCond(argv[4]);
    long long reportEvery = std::stoll(argv[5]);

    knotbuilder::PDCode pdcode = knotbuilder::parsePDCode(pdcodeStr);
    auto [t2, edges0] = knotbuilder::buildLink(pdcode);

    CobordismBuilder<3> cob(t2);
    if (layers > 0)
      cob.thicken(layers);
    regina::Triangulation<4> tri = cob.cone();

    std::cerr << "[driver] dim4 pentachora = " << tri.size()
              << ", triangles = " << tri.countTriangles() << "\n";

    DiagSearch search(tri);
    std::cerr << "[driver] embeddable faces = " << search.numEmbeddableFaces()
              << "\n";
    search.diagRun(cond, reportEvery);
  } else {
    std::cerr << "Unknown mode: " << mode << "\n";
    return 1;
  }

  return 0;
}
