#ifndef GLUING_H 

#define GLUING_H

#include "triangulation/generic/face.h"

/** Gluing Implementation */

template <int dim, int subdim>
struct Gluing {
    regina::Face<dim, subdim> *src;
    int srcFacet;
    regina::Face<dim, subdim> *dst;
    regina::Perm<subdim + 1> gluing;

    Gluing() = default;

    Gluing(regina::Face<dim, subdim> *src, int srcFacet,
           regina::Face<dim, subdim> *dst, regina::Perm<subdim + 1> gluing)
        : src(src), srcFacet(srcFacet), dst(dst), gluing(gluing) {}

    friend std::ostream &operator<<(std::ostream &os, const Gluing &g) {
        return os << "(" << g.src->index() << ", " << g.srcFacet << ", "
                  << g.dst->index() << ", " << g.gluing << ")";
    }
};

template <int dim>
class GluingNode {
   public:
    using AdjList = std::unordered_map<GluingNode *, Gluing<dim, 2>>;

    regina::Triangle<dim> *f;
    AdjList adjList;
    bool visited = false;

    GluingNode(regina::Face<dim, 2> *f) : f(f) {}

    friend std::ostream &operator<<(std::ostream &os, const GluingNode &node) {
        os << "(" << node.f->index() << " { ";
        for (const auto &[_, gluing] : node.adjList) {
            os << gluing << " ";
        }

        os << "})";

        return os;
    }
};

#endif
