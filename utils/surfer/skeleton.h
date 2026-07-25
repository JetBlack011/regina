//
//  skeleton.h
//
//  Created by John Teague on 07/15/2026.
//

#ifndef SKELETON_H

#define SKELETON_H

#include <maths/perm.h>
#include <triangulation/forward.h>

/*! \file utils/surfer/skeleton.h
 *  \brief Builds an adjacency graph over the subdim-faces of a triangulation.
 */

/**
 * The adjacency graph of a triangulation's subdim-faces: one node per
 * subdim-face, with an edge (a Gluing) for every pair of subdim-faces
 * sharing a (subdim-1)-facet.
 *
 * This is the graph that ConnectedInducedSubgraphEnumerator/EmbeddingSearch
 * search over: a connected induced subgraph corresponds to a connected set
 * of faces that may glue together into a subcomplex.
 *
 * \tparam dim the dimension of the ambient triangulation.
 * \tparam subdim the dimension of the faces used as graph nodes
 * (1 <= subdim <= dim).
 */
template <int dim, int subdim>
class Skeleton {
  public:
    using Face = regina::SafeFace<dim, subdim>; /**< A graph node: a subdim-face of the ambient triangulation. */
    using Facet = regina::Face<dim, subdim - 1>; /**< A shared facet between two Faces. */

    /**
     * One shared facet between two subdim-faces, recorded as a directed
     * (src, dst) pair together with the vertex permutation identifying
     * src's facet with dst's.
     *
     * Each shared facet produces a reciprocal pair of Gluings (src/dst
     * swapped), including self-folds where src and dst are the same face.
     */
    struct Gluing {
        const size_t srcIndex; /**< Index into Skeleton::getNodes() of the source face. */
        int srcFacet; /**< Which local facet of the source face this gluing identifies. */
        const size_t dstIndex; /**< Index into Skeleton::getNodes() of the destination face. */
        regina::Perm<subdim + 1> gluing;
            /**< The vertex permutation carrying the source face's vertex
                 labelling to the destination face's. */

        /**
         * Writes this gluing as `(srcIndex, srcFacet, dstIndex, gluing)`.
         */
        friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
            return os << "(" << g.srcIndex << ", " << g.srcFacet << ", "
                      << g.dstIndex << ", " << g.gluing << ")";
        }
    };

  private:

    /** A single subdim-face together with its Gluings to adjacent faces. */
    struct Node {
        Face *face; /**< The underlying subdim-face. */
        std::vector<Gluing> gluings; /**< This face's Gluings to every adjacent face. */

        Node(Face *f) : face(f) {}
    };

    const regina::Triangulation<dim> *tri_; /**< The ambient triangulation. */
    const size_t numVertices_; /**< tri_'s vertex count (not this graph's node count). */
    std::vector<Node> nodes_; /**< One entry per subdim-face of tri_. */

  public:
    /**
     * Builds the adjacency graph of `tri`'s subdim-faces.
     */
    Skeleton(const regina::Triangulation<dim> &tri);

    /**
     * Returns the ambient triangulation this graph was built from.
     */
    const regina::Triangulation<dim> &triangulation() const { return *tri_; }

    /**
     * Returns the ambient triangulation's vertex count.
     *
     * \note This is _not_ the number of nodes in this graph (that is
     * numFaces()) -- it is the vertex count of the underlying
     * Triangulation<dim>.
     */
    const size_t numVertices() const { return numVertices_; }

    /**
     * Returns the number of nodes (subdim-faces) in this graph.
     */
    const size_t numFaces() const { return nodes_.size(); }

    /**
     * Returns every node in this graph, indexed as in the underlying
     * subdim-face numbering.
     */
    const std::vector<Node> &getNodes() const { return nodes_; }

    /**
     * Writes this graph in Graphviz DOT format: one node per subdim-face
     * and one undirected edge per shared facet. Render with e.g.
     * `dot -Tsvg`.
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
    /** Populates nodes_ and every node's Gluings from `tri`'s skeleton. */
    void buildSkeleton_(const regina::Triangulation<dim> &tri);
};

extern template class Skeleton<3, 2>;
extern template class Skeleton<4, 2>;

#endif // SKELETON_H
