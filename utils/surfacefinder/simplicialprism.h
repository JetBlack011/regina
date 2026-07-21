#ifndef SIMPLICIALPRISM_H

#define SIMPLICIALPRISM_H

#include <vector>

#include <maths/perm.h>
#include <triangulation/forward.h>

template <int dim>
class SimplicialPrism {
  private:
    std::vector<regina::Simplex<dim> *> simplices_;

  public:
    SimplicialPrism(regina::Triangulation<dim> &tri);

    // The k-th sub-simplex of this prism (see decode_/encode_ for how its
    // local vertices relate to the base simplex being thickened). Exposed
    // so callers can locate specific faces of the thickening directly --
    // e.g. the triangles swept out by a single base edge -- without needing
    // a further layer of prism-specific API.
    regina::Simplex<dim> *simplex(int k) const { return simplices_[k]; }

    // Local vertex index (within simplex(k)) holding (v, isTop), if that
    // vertex is actually present there. Public wrapper around encode_ for
    // callers that need to locate a specific base vertex's image within a
    // specific sub-simplex (see decode_'s class-level documentation).
    static int localVertex(int v, bool isTop) { return encode_(v, isTop); }

    // Glues this prism to `other` along the wall lying over the facet of
    // the "top" simplex obtained by omitting vertex `facet` (resp.
    // `otherFacet` for `other`). Because the underlying triangulation is
    // assumed to be ordered, the correspondence between the walls' two
    // sub-triangulations (and the permutation identifying them) is
    // completely determined by `facet` and `otherFacet` alone.
    void glue(int facet, SimplicialPrism<dim> &other, int otherFacet);

    // Glues this prism's top facet (the base simplex x {1}) to `next`'s
    // bottom facet (the base simplex x {0}), stacking `next` on top of this
    // prism to extend the interval factor. Since both facets triangulate
    // the same base simplex with no relabelling, the only two sub-simplices
    // involved are simplices_[dim-1] (whose sole non-top vertex is the
    // bottom copy of the top facet) and next.simplices_[0] (whose sole
    // non-bottom vertex is the top copy of the bottom facet); matching up
    // decode_/encode_ across the seam reduces to the cyclic shift
    // i -> (i + 1) mod (dim + 1).
    void stitchTop(SimplicialPrism<dim> &next);

    // Glues this prism's top facet (the base simplex x {1}) directly onto
    // `coneSimplex`, a single dim-simplex representing the base simplex
    // coned to a new apex point: `coneSimplex`'s local vertex i (for
    // i = 0,...,dim-1) is assumed to equal the base simplex's own vertex i,
    // with local vertex dim the apex (the convention used by
    // CobordismBuilder::cone()). This caps the top off in one step, using
    // the same decode_ used by glue()/stitchTop() rather than a further
    // layer of prism structure.
    void capTop(regina::Simplex<dim> *coneSimplex);

  private:
    // The unique order-preserving bijection from {0,...}\{a} to {0,...}\{b},
    // evaluated at x (where x != a).
    static int mapExcluding_(int x, int a, int b);

    // Decodes local vertex m of simplices_[k] as (base vertex, is-top-copy).
    static std::pair<int, bool> decode_(int k, int m);

    // Inverse of decode_: the local vertex index (within some simplices_[k])
    // holding (v, isTop), given that this vertex is actually present there.
    static int encode_(int v, bool isTop);

    // The local facet index (within simplices_[k]) of the facet lying in
    // wall(v), for k != v.
    static int wallFacet_(int k, int v);
};

extern template class SimplicialPrism<3>;
extern template class SimplicialPrism<4>;

#endif // SIMPLICIALPRISM_H
