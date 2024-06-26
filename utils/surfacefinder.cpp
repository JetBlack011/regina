//
//  main.cpp
//
//  Created by John Michael Teague on 06/19/24.
//

#include <gmpxx.h>
#include <link/link.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <cstring>
#include <iostream>
#include <ostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

#include "maths/perm.h"
#include "triangulation/example3.h"
#include "triangulation/example4.h"
#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"
#include "utilities/exception.h"

static const regina::Triangulation<2> GENUS_2_SURFACE =
    regina::Triangulation<2>::fromGluings(6, {{0, 2, 5, {2, 1, 0}},
                                              {0, 1, 1, {0, 2, 1}},
                                              {0, 0, 5, {1, 0, 2}},
                                              {1, 1, 2, {0, 2, 1}},
                                              {1, 0, 3, {0, 2, 1}},
                                              {2, 1, 3, {0, 2, 1}},
                                              {2, 0, 4, {0, 2, 1}},
                                              {3, 1, 4, {0, 2, 1}},
                                              {4, 1, 5, {0, 2, 1}}});

static const std::vector<int> NUM_TRIANGLES = {0, 0, 1, 4, 10, 20, 35};

class ClosedConnectedSurface {
private:
  const regina::Triangulation<2> surface_;
  bool isOrientable_;
  int genus_;

public:
  ClosedConnectedSurface(regina::Triangulation<2> &surface)
      : surface_(surface), isOrientable_(surface.isOrientable()) {
    if (!surface.isClosed() || !surface.isConnected()) {
      throw regina::InvalidArgument("Surface must be closed and connected!");
    }

    if (isOrientable_) {
      genus_ = 1 - surface.eulerChar() / 2;
    } else {
      genus_ = 2 - surface.eulerChar();
    }
  }

  friend std::ostream &operator<<(std::ostream &os,
                                  const ClosedConnectedSurface &surface) {
    if (surface.isOrientable_ && surface.genus_ <= 1) {
      if (surface.genus_ == 0) {
        os << "S^2";
      } else {
        os << "Torus";
      }
    } else if (!surface.isOrientable_ && surface.genus_ <= 2) {
      if (surface.genus_ == 1) {
        os << "RP^2";
      } else {
        os << "Klein Bottle";
      }
    } else {
      if (surface.isOrientable_) {
        os << "Orientable";
      } else {
        os << "Non-orientable";
      }
      os << " surface of genus " << surface.genus_;
    }

    return os;
  }

  friend bool operator==(const ClosedConnectedSurface &lhs,
                         const ClosedConnectedSurface &rhs) {
    return lhs.isOrientable_ == rhs.isOrientable_ && lhs.genus_ == rhs.genus_;
  }

  friend bool operator<(const ClosedConnectedSurface &lhs,
                        const ClosedConnectedSurface &rhs) {
    std::vector<int> lhsLex = {lhs.isOrientable_, lhs.genus_};
    std::vector<int> rhsLex = {rhs.isOrientable_, rhs.genus_};
    return lhsLex < rhsLex;
  }
};

template <int n> class GluingGraph {
private:
  using TriangleMap =
      std::unordered_map<const regina::Triangle<n> *, regina::Triangle<2> *>;

  struct Gluing {
    // Only included to make edge building easier
    const regina::Triangle<n> *src;
    const int srcFacet;
    const regina::Triangle<n> *dst;
    const regina::Perm<3> gluing;

    friend std::ostream &operator<<(std::ostream &os, const Gluing &gluing) {
      os << "(Triangle ";
      if (gluing.src == nullptr) {
        os << -1;
      } else {
        os << gluing.src->index();
      }
      os << ", " << gluing.srcFacet << ", Triangle " << gluing.dst->index()
         << ", " << gluing.gluing << ")";

      return os;
    }
  };

  class GluingNode {
  public:
    const Gluing gluing_;
    std::unordered_set<GluingNode *> adjList_;
    std::unordered_set<GluingNode *> invalids_; // Make sure to precompute
    bool valid_ = true;
    bool visited_ = false;

    // TODO: figure out C++
    GluingNode(const regina::Triangle<n> *src, int srcFacet,
               const regina::Triangle<n> *dst, regina::Perm<3> gluing)
        : gluing_{src, srcFacet, dst, gluing} {}

    /**
     * Comparison operators let us guarantee gluing uniqueness using
     * std::set, since I really don't want to come up with a hash function.
     */
    friend bool operator==(const GluingNode &lhs, const GluingNode &rhs) {
      // We never want two gluing nodes with the same gluing, use this to
      // distinguish the nodes
      return lhs.gluing_.src == rhs.gluing_.src &&
             lhs.gluing_.srcFacet == rhs.gluing_.srcFacet &&
             lhs.gluing_.dst == rhs.gluing_.dst &&
             lhs.gluing_.gluing == rhs.gluing_.gluing;
    }

    friend bool operator<(const GluingNode &lhs, const GluingNode &rhs) {
      // Arbitrary tie breaker. If by some miracle you manage to
      // stick in a triangulation with over LONG_MAX triangles, God help
      // you.
      std::vector<long> lhsLex = {
          lhs.gluing_.src == nullptr
              ? -1
              : static_cast<long>(lhs.gluing_.src->index()),
          static_cast<long>(lhs.gluing_.dst->index()),
          lhs.gluing_.srcFacet,
          lhs.gluing_.gluing[0],
          lhs.gluing_.gluing[1],
          lhs.gluing_.gluing[2]};
      std::vector<long> rhsLex = {
          rhs.gluing_.src == nullptr
              ? -1
              : static_cast<long>(rhs.gluing_.src->index()),
          static_cast<long>(rhs.gluing_.dst->index()),
          rhs.gluing_.srcFacet,
          rhs.gluing_.gluing[0],
          rhs.gluing_.gluing[1],
          rhs.gluing_.gluing[2]};

      return lhsLex < rhsLex;
    }

    friend std::ostream &operator<<(std::ostream &os, const GluingNode &node) {
      os << "(" << node.gluing_ << " { ";
      for (const GluingNode *adj : node.adjList_) {
        os << adj->gluing_ << " ";
      }

      os << "})";

      return os;
    }
  };

  size_t calls_ = 0;

  std::vector<GluingNode> nodes_;
  /**< Guaranteed not to be nodes with the same gluing data */

  regina::Triangulation<2> surface_;
  /**< To keep track of surfaces we built during DFS */

  std::vector<regina::Triangulation<2>> embeddedSurfaces_;
  /**< The surfaces we find embedded in the given triangulation */

  TriangleMap triangleMap_;
  /**< A mapping between the triangles in the given triangulation and those in
   * the surface we build during DFS */

  std::set<GluingNode> addTriangle_(const regina::Triangle<n> *triangle) {
    std::set<GluingNode> nodes;
    const regina::Triangulation<n> &tri = triangle->triangulation();

    // std::cout << "Adding " << *triangle << "\n";

    // For each edge of our triangle, find all of the ways adjacent
    // triangles are glued to that edge in the given triangulation
    for (int facet = 0; facet < 3; ++facet) {
      const regina::Edge<n> *edge = triangle->edge(facet);
      regina::Perm<n + 1> edgeToTriangle = triangle->edgeMapping(facet);

      for (const regina::Triangle<n> *other : tri.triangles()) {
        for (int i = 0; i < 3; ++i) {
          // Only glue identified edges and ignore identity mappings
          if (other->edge(i) != edge ||
              (other->edge(i) == edge && triangle == other))
            continue;

          regina::Perm<n + 1> edgeToOther = other->edgeMapping(i);
          // regina::Perm is immutable, use std::array instead
          std::array<int, 3> p;

          for (int i = 0; i < 3; ++i) {
            p[edgeToTriangle[i]] = edgeToOther[i];
          }

          GluingNode node = {triangle, facet, other, {p[0], p[1], p[2]}};
          // std::cout << "NEW NODE " << node;
          nodes.insert(node);
        }
      }
    }

    // std::cout << "\n";
    return nodes;
  }

  void buildGluingNodes_(const regina::Triangulation<n> &tri) {
    std::set<GluingNode> nodes;

    // TODO: Can we do better than O(|tri.triangles()|^2)?
    for (const regina::Triangle<n> *triangle : tri.triangles()) {
      std::set<GluingNode> newNodes = addTriangle_(triangle);

      for (const GluingNode &node : newNodes) {
        nodes.insert(node);
      }
    }

    for (const GluingNode &node : nodes) {
      nodes_.push_back(node);
    }
  }

  void buildGluingEdges_() {
    for (GluingNode &node : nodes_) {
      for (GluingNode &other : nodes_) {
        // std::cout << node << ", " << other;
        if (node.gluing_.dst == other.gluing_.src) {
          node.adjList_.insert(&other);
        }
      }
    }
  }

  void nspaces_(int m) {
    for (int i = 0; i < m; ++i) {
      std::cout << "  ";
    }
  }

  void dfs_(regina::Triangle<2> *triangle, GluingNode *node, int layer) {
    if (node->visited_ || !node->valid_) {
      // Don't perform the same gluing twice, and stop if adding this
      // gluing would result in an invalid triangulation
      return;
    }

    ++calls_;

    node->visited_ = true;

    regina::Triangle<2> *nextTriangle;

    auto search = triangleMap_.find(node->gluing_.dst);
    bool isNewTriangle = search == triangleMap_.end();

    if (isNewTriangle) {
      nextTriangle = surface_.newTriangle();
      triangleMap_.insert({node->gluing_.dst, nextTriangle});
    } else {
      nextTriangle = search->second;
    }

    // nspaces_(layer);
    // std::cout << "Layer " << layer << ": " << *node
    //           << ", isNewTriangle = " << isNewTriangle << "\n";

    // TODO: Is there some way to move this root check outside the
    // recursion? Probably doesn't matter
    if (triangle != nullptr) {
      try {
        triangle->join(node->gluing_.srcFacet, nextTriangle,
                       node->gluing_.gluing);
      } catch (regina::InvalidArgument &e) {
        if (isNewTriangle) {
          // If this is the first time we've seen this triangle and
          // the gluing is invalid, we're safe to remove it from
          // surface_
          surface_.removeTriangle(nextTriangle);
          // Note that find is guaranteed to return a valid iterator
          triangleMap_.erase(triangleMap_.find(node->gluing_.dst));
        }

        return;
      }
    }

    // nspaces_(layer);
    // std::cout << layer << ": Joined triangle ";
    // if (node->gluing_.src == nullptr) {
    //     std::cout << -1;
    // } else {
    //     std::cout << node->gluing_.src->index();
    // }
    // std::cout << " and " << node->gluing_.dst->index() << " along "
    //           << node->gluing_.srcFacet << " by " << node->gluing_.gluing
    //           << "\n";

    // TODO: Mark invalid gluings as invalid, and don't forget to set them
    // back at the end!
    // for (GluingNode* invalid : node->invalids_) {
    //     invalid->valid_ = false;
    // }

    // TODO: Add condition for a finished surface
    // if (surface_.isClosed()) {
    // std::cout << "New surface: " << embeddedSurfaces_.size() << "\n";
    embeddedSurfaces_.push_back(surface_);
    //}

    // Recursive step: walk along each possible next gluing
    for (GluingNode *nextNode : node->adjList_) {
      dfs_(nextTriangle, nextNode, layer + 1);
    }

    // Remember to remove triangles as we move back up the recursion and
    // mark this node as not visited
    if (isNewTriangle) {
      surface_.removeTriangle(nextTriangle);
      triangleMap_.erase(triangleMap_.find(node->gluing_.dst));
    }

    // nspaces(layer);
    // std::cout << layer << ": DONE\n";
    node->visited_ = false;
  }

public:
  GluingGraph(const regina::Triangulation<n> &tri) {
    buildGluingNodes_(tri);
    buildGluingEdges_();
    // buildGluingInvalids_();
  }

  void reset() {
    surface_ = {};
    triangleMap_.clear();
    embeddedSurfaces_.clear();

    for (GluingNode &node : nodes_) {
      node.valid_ = true;
      node.visited_ = false;
    }
  }

  std::vector<regina::Triangulation<2>> &
  findEmbeddedSurfaces(regina::Triangle<n> *triangle) {
    if (!embeddedSurfaces_.empty()) {
      return embeddedSurfaces_;
    }

    reset(); // No real need to call this here, but why not
    GluingNode root = {nullptr, -1, triangle, {}};
    for (auto &node : nodes_) {
      if (node.gluing_.src == triangle) {
        root.adjList_.insert(&node);
      }
    }
    // Call DFS with dummy root node
    dfs_(nullptr, &root, 0);

    // std::cout << "TOTAL DFS CALLS = " << calls_ << "\n";

    // std::cout << "\n";

    return embeddedSurfaces_;
  }
};

template <typename T, typename D>
std::ostream &operator<<(std::ostream &os, const std::pair<T, D> &p) {
  os << "(" << p.first << ", " << p.second << ")";
  return os;
}

template <typename T>
std::ostream &operator<<(std::ostream &os, const std::tuple<T, T, T> &p) {
  os << "(" << std::get<0>(p) << ", " << std::get<1>(p) << ", "
     << std::get<2>(p) << ")";
  return os;
}

void usage(const char *progName, const std::string &error = std::string()) {
  if (!error.empty())
    std::cerr << error << "\n\n";

  std::cerr << "Usage:\n";
  std::cerr << "    " << progName
            << " <isosig>\n"
               "    "
            << progName << " [ -v, --version | -?, --help ]\n\n";
  std::cerr << "    -v, --version : Show which version of Regina "
               "is being used\n";
  std::cerr << "    -?, --help    : Display this help\n";
  exit(1);
}

int main(int argc, char *argv[]) {
  // Check for standard arguments:
  for (int i = 1; i < argc; ++i) {
    if (strcmp(argv[i], "-?") == 0 || strcmp(argv[i], "--help") == 0)
      usage(argv[0]);
    if (strcmp(argv[i], "-v") == 0 || strcmp(argv[i], "--version") == 0) {
      if (argc != 2)
        usage(argv[0], "Option --version cannot be used with "
                       "any other arguments.");
      std::cout << PACKAGE_BUILD_STRING << "\n";
      exit(0);
    }
  }

  if (argc < 2) {
    usage(argv[0], "Please provide an isomorphism signature.");
    exit(0);
  }

  std::string isoSig = argv[1];

  // regina::Triangulation<4> tri(isoSig);

  // regina::Triangulation<3> bdry = tri.boundaryComponent(0)->build();

  // auto tri = regina::Example<4>::cp2();
  regina::Triangulation<4> tri(isoSig);
  // tri.subdivide();
  std::cout << tri.detail() << "\n\n";
  // tri.subdivide();
  //  tri.newSimplex();
  //  tri.newSimplex();
  //  tri.tetrahedron(0)->join(0, tri.tetrahedron(1), {});

  std::cout << "--- Closed Surfaces ---\n";
  std::map<ClosedConnectedSurface, int> closedSurfaces;
  int surfaceCount = 0;
  int closedCount = 0;
  GluingGraph graph(tri);
  for (int i = 0; i < tri.countTriangles(); ++i) {
    // std::cout << "Starting at Triangle " << i;
    std::vector<regina::Triangulation<2>> &ans =
        graph.findEmbeddedSurfaces(tri.triangle(i));

    // for (int i = 0; i < tri.countTriangles(); ++i) {
    //     for (auto &surface : graph.findEmbeddedSurfaces(tri.triangle(3)))
    //     {
    //         ans.push_back(surface);
    //     }
    // }

    // std::cout << "--- Surfaces ---\n";
    // for (auto &surface : ans) {
    //     if (surface.isClosed()) {
    //         ++closedCount;
    //     }
    //     std::cout << ++i << ": Euler characteristic = " <<
    //     surface.eulerChar()
    //               << ", orientable = " << surface.isOrientable()
    //               << ", and boundary components = "
    //               << surface.countBoundaryComponents() << "\n";
    // }

    // std::cout << "\n";

    for (auto &surface : ans) {
      ++surfaceCount;
      if (surface.isClosed()) {
        ++closedCount;
        ClosedConnectedSurface ccs(surface);
        auto search = closedSurfaces.find(ccs);
        if (search == closedSurfaces.end()) {
          closedSurfaces.insert({ccs, 1});
        } else {
          search->second++;
        }
        // std::cout << ++closedCount
        //           << ": Euler characteristic = " <<
        //           surface.eulerChar()
        //           << ", orientable = " << surface.isOrientable()
        //           << ", and boundary components = "
        //           << surface.countBoundaryComponents() << "\n";
      }
    }

    // std::cout << "Total surfaces = " << ans.size()
    //           << "\nTotal closed surfaces = " << closedCount << "\n\n";

    graph.reset();
  }
  std::cout << "Total surfaces = " << surfaceCount
            << "\nTotal closed surfaces = " << closedCount << "\n\n";

  if (closedSurfaces.size() == 0) {
    std::cout << "Didn't find any closed surfaces!";
  } else {
    std::cout << "Found the following surfaces:\n";
    for (auto &[ccs, count] : closedSurfaces) {
      std::cout << "- " << count << " " << (count == 1 ? "copy" : "copies")
                << " of " << ccs << "\n";
    }
  }

  return 0;
}