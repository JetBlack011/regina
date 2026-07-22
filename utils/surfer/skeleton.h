#ifndef SKELETON_H

#define SKELETON_H

#include <maths/perm.h>
#include <triangulation/forward.h>

template <int dim, int subdim>
class Skeleton {
  public:
    using Face = regina::SafeFace<dim, subdim>;
    using Facet = regina::Face<dim, subdim - 1>;

    struct Gluing {
        // Face *src;
        const size_t srcIndex;
        int srcFacet;
        // Face *dst;
        const size_t dstIndex;
        regina::Perm<subdim + 1> gluing;

        // Kept defined in-class (rather than declared here, defined in
        // skeleton.cpp) because it's a friend of a class template nested
        // inside another class template: a separately-defined function
        // template doesn't reliably bind back to the per-instantiation
        // friend declaration this would need. It's a 1-line format string,
        // so the cost of leaving it inline is negligible.
        friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
            return os << "(" << g.srcIndex << ", " << g.srcFacet << ", "
                      << g.dstIndex << ", " << g.gluing << ")";
        }
    };

  private:

    struct Node {
        Face *face;
        std::vector<Gluing> gluings;

        Node(Face *f) : face(f) {}
    };

    const regina::Triangulation<dim> *tri_;
    const size_t numVertices_;
    std::vector<Node> nodes_;

  public:
    Skeleton(const regina::Triangulation<dim> &tri);

    const regina::Triangulation<dim> &triangulation() const { return *tri_; }

    const size_t numVertices() const { return numVertices_; }

    const size_t numFaces() const { return nodes_.size(); }

    const std::vector<Node> &getNodes() const { return nodes_; }

    // Writes this graph in Graphviz DOT format: one node per subdim-face
    // and one undirected edge per shared facet. Render with e.g. `dot
    // -Tsvg`.
    //
    // Kept defined in-class for the same reason as Gluing::operator<<
    // above: it's a friend of a class template, and a separately-defined
    // function template doesn't reliably bind back to the per-instantiation
    // friend declaration this would need.
    friend std::ostream &operator<<(std::ostream &os,
                                    const Skeleton &skeleton) {
        os << "graph Skeleton {\n";

        for (const auto &node : skeleton.nodes_) {
            os << "  " << node.face->index() << ";\n";
        }

        for (const auto &node : skeleton.nodes_) {
            for (const auto &gluing : node.gluings) {
                size_t srcIndex = gluing.srcIndex;
                size_t dstIndex = gluing.dstIndex;
                int dstFacet = gluing.gluing[gluing.srcFacet];

                // Each shared facet is stored as a reciprocal pair of
                // Gluings (src/dst swapped); emit only one undirected edge
                // per pair, picking a canonical direction so self-folds
                // (srcIndex == dstIndex) are deduplicated too
                if (srcIndex > dstIndex ||
                    (srcIndex == dstIndex && gluing.srcFacet > dstFacet)) {
                    continue;
                }

                os << "  " << srcIndex << " -- " << dstIndex << " [label=\""
                   << gluing.srcFacet << ": " << gluing.gluing << "\"];\n";
            }
        }

        os << "}\n";
        return os;
    }

  private:
    void buildSkeleton_(const regina::Triangulation<dim> &tri);
};

extern template class Skeleton<3, 2>;
extern template class Skeleton<4, 2>;

#endif // SKELETON_H
