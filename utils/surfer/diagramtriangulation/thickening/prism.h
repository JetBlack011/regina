//
//  simplicialprism.h
//
//  Created by John Teague on 06/10/2025.
//

#ifndef SIMPLICIALPRISM_H

#define SIMPLICIALPRISM_H

#include <vector>

#include <maths/perm.h>
#include <triangulation/forward.h>

/*! \file utils/surfer/simplicialprism.h
 *  \brief Triangulates the product of a simplex with an interval.
 */

/**
 * A triangulation of (a dim-simplex) x [0,1], built from dim glued-together
 * (dim+1)-simplices.
 *
 * This is the basic building block used to thicken a triangulation one
 * dimension up (see CobordismBuilder), one base simplex at a time.
 *
 * \tparam dim the dimension of the base simplex being thickened; the prism
 * itself lives in dimension dim+1.
 */
template <int dim>
class SimplicialPrism {
  private:
    std::vector<regina::Simplex<dim> *> simplices_;
        /**< The dim sub-simplices making up this prism, in construction
             order. */

  public:
    /**
     * Builds a new prism in `tri`, over a fresh base simplex.
     */
    SimplicialPrism(regina::Triangulation<dim> &tri);

    /**
     * Returns the k-th sub-simplex of this prism.
     *
     * See decode_()/encode_() for how a sub-simplex's local vertices relate
     * to the base simplex being thickened. Exposed so callers can locate
     * specific faces of the thickening directly (e.g. the triangles swept
     * out by a single base edge) without a further layer of prism-specific
     * API.
     */
    regina::Simplex<dim> *simplex(int k) const { return simplices_[k]; }

    /**
     * Returns the local vertex index, within whichever sub-simplex holds
     * it, of the top (\a isTop) or bottom copy of base vertex \a v.
     *
     * Public wrapper around encode_(), for callers that need to locate a
     * specific base vertex's image within a specific sub-simplex.
     */
    static int localVertex(int v, bool isTop) { return encode_(v, isTop); }

    /**
     * Glues this prism to `other` along the wall lying over the facet of
     * the "top" simplex obtained by omitting vertex `facet` (respectively
     * `otherFacet` for `other`).
     *
     * The underlying triangulation is assumed to be ordered, so the
     * correspondence between the two walls' sub-triangulations (and the
     * permutation identifying them) is completely determined by `facet`
     * and `otherFacet` alone.
     */
    void glue(int facet, SimplicialPrism<dim> &other, int otherFacet);

    /**
     * Glues this prism's top facet (the base simplex x {1}) to `next`'s
     * bottom facet (the base simplex x {0}), stacking `next` on top of this
     * prism to extend the interval factor.
     *
     * Both facets triangulate the same base simplex with no relabelling, so
     * only simplices_[dim-1] and next.simplices_[0] are involved, and
     * matching up decode_()/encode_() across the seam reduces to the cyclic
     * shift i -> (i + 1) mod (dim + 1).
     */
    void stitchTop(SimplicialPrism<dim> &next);

    /**
     * Glues this prism's top facet (the base simplex x {1}) directly onto
     * `coneSimplex`, a single dim-simplex representing the base simplex
     * coned to a new apex point.
     *
     * `coneSimplex`'s local vertex i (for i = 0,...,dim-1) is assumed to
     * equal the base simplex's own vertex i, with local vertex dim the
     * apex -- the convention used by CobordismBuilder::cone(). This caps
     * the top off in one step, using the same decode_() as glue()/
     * stitchTop() rather than a further layer of prism structure.
     */
    void capTop(regina::Simplex<dim> *coneSimplex);

  private:
    /**
     * The unique order-preserving bijection from {0,...,dim} with `a`
     * removed to {0,...,dim} with `b` removed, evaluated at `x` (where
     * `x != a`).
     */
    static int mapExcluding_(int x, int a, int b);

    /**
     * Decodes local vertex `m` of simplices_[k] as (base vertex,
     * is-top-copy).
     */
    static std::pair<int, bool> decode_(int k, int m);

    /**
     * The inverse of decode_(): the local vertex index (within some
     * simplices_[k]) holding (v, isTop), given that this vertex is
     * actually present there.
     */
    static int encode_(int v, bool isTop);

    /**
     * The local facet index (within simplices_[k]) of the facet lying in
     * wall(v), for `k != v`.
     */
    static int wallFacet_(int k, int v);
};

extern template class SimplicialPrism<3>;
extern template class SimplicialPrism<4>;

#endif // SIMPLICIALPRISM_H
