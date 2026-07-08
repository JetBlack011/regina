#ifndef GLUING_H

#define GLUING_H

#include <ostream>
#include <unordered_map>
#include <vector>

#include <maths/perm.h>
#include <triangulation/forward.h>

/** Gluing Implementation */

template <int dim, int subdim>
struct Gluing {
    regina::SafeFace<dim, subdim> *src;
    int srcFacet;
    regina::SafeFace<dim, subdim> *dst;
    regina::Perm<subdim + 1> gluing;

    Gluing() = default;

    Gluing(regina::SafeFace<dim, subdim> *src, int srcFacet,
           regina::SafeFace<dim, subdim> *dst, regina::Perm<subdim + 1> gluing)
        : src(src), srcFacet(srcFacet), dst(dst), gluing(gluing) {}

    friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
        return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
                  << g.dst->index() << ", " << g.gluing << ")";
    }
};

template <int dim>
class GluingNode {
  public:
    // Two triangles can share more than one edge with each other (common in
    // minimal/coarse triangulations, e.g. the standard one-vertex 2-triangle
    // torus, where every edge of one triangle is glued to the other) -- so
    // this must hold every gluing to a given neighbour, not just one. A map
    // keyed only by neighbour previously overwrote all but the last such
    // gluing, silently dropping the vertex identifications the others would
    // have produced and causing KnottedSurface's self-intersection tracking
    // to see incomplete adjacency (false positives on otherwise-valid,
    // embeddable surfaces).
    using AdjList =
        std::unordered_map<GluingNode *, std::vector<Gluing<dim, 2>>>;

    regina::Triangle<dim> *f;
    AdjList adjList;
    bool visited = false;
    // True while this node currently has an (undecided-or-decided) slot
    // somewhere in extend_'s frontier vector, from the moment it's pushed
    // until the inclusion decision that pushed it is backtracked -- lets
    // extend_ avoid pushing the same candidate onto the frontier more than
    // once. Safe now that self-intersection rejection is decided by
    // KnottedSurface::hasUnresolvableConflict() (order-independent) rather
    // than the old immediate hasSelfIntersection() check (order-dependent);
    // previously a candidate rejected once due to insertion order relied on
    // being re-pushed later, under different context, for a second chance.
    bool queued = false;

    GluingNode(regina::Face<dim, 2> *f) : f(f) {}

    friend std::ostream &operator<<(std::ostream &os, const GluingNode &node) {
        os << "(" << node.f->index() << " { ";
        for (const auto &[_, gluings] : node.adjList) {
            for (const auto &gluing : gluings)
                os << gluing << " ";
        }

        os << "})";

        return os;
    }
};

#endif
