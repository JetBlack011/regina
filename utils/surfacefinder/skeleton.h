#ifndef SKELETON_H

#define SKELETON_H

#include <array>
#include <utility>

#include <maths/perm.h>
#include <triangulation/forward.h>

enum class BoundaryCondition : uint8_t { all, closed, proper, connected };

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
    Skeleton(const regina::Triangulation<dim> &tri)
        : tri_(&tri), numVertices_(tri.countVertices()) {
        static_assert(subdim >= 1 && subdim <= dim,
                      "Skeleton requires 1 <= subdim <= dim");
        buildSkeleton_(tri);
    }

    const regina::Triangulation<dim> &triangulation() const { return *tri_; }

    const size_t numVertices() const { return numVertices_; }

    const size_t numFaces() const { return nodes_.size(); }

    const std::vector<Node> &getNodes() const { return nodes_; }

    /**
     * Writes this graph in Graphviz DOT format: one node per subdim-face
     * and one undirected edge per shared facet. Render with e.g. `dot -Tsvg`
     */
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
    void buildSkeleton_(const regina::Triangulation<dim> &tri) {
        if constexpr (subdim == dim) {
            nodes_.reserve(tri.size());
            for (const auto &simplex : tri.simplices()) {
                nodes_.emplace_back(simplex);
            }
        } else {
            nodes_.reserve(tri.template countFaces<subdim>());
            for (const auto &face : tri.template faces<subdim>()) {
                nodes_.emplace_back(face);
            }
        }

        // Make an edge for each pair of subdim-faces that share a facet
        // (including self-folds)
        std::vector<std::vector<std::pair<size_t, int>>> buckets(
            tri.template countFaces<subdim - 1>());
        for (size_t fi = 0; fi < nodes_.size(); ++fi) {
            Face *face = nodes_[fi].face;
            for (int i = 0; i <= subdim; ++i) {
                Facet *facet = face->template face<subdim - 1>(i);
                buckets[facet->index()].emplace_back(fi, i);
            }
        }

        for (const auto &bucket : buckets) {
            for (const auto &[fi, i] : bucket) {
                for (const auto &[fj, j] : bucket) {
                    // Skip only comparing an occurrence to itself; a
                    // genuine self-fold (fi == fj, i != j) is kept.
                    if (fi == fj && i == j) {
                        continue;
                    }

                    Face *src = nodes_[fi].face;
                    Face *dst = nodes_[fj].face;

                    regina::Perm<dim + 1> m1 =
                        src->template faceMapping<subdim - 1>(i);
                    regina::Perm<dim + 1> m2 =
                        dst->template faceMapping<subdim - 1>(j);

                    std::array<int, subdim + 1> gluing;
                    for (int k = 0; k <= subdim; ++k) {
                        gluing[m1[k]] = m2[k];
                    }

                    nodes_[fi].gluings.emplace_back(
                        fi, i, fj, regina::Perm<subdim + 1>(gluing));
                }
            }
        }
    }
};

#endif // SKELETON_H
