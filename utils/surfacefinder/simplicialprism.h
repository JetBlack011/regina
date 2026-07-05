#ifndef SIMPLICIALPRISM_H

#define SIMPLICIALPRISM_H

#include <unordered_set>
#include <vector>

#include <triangulation/forward.h>
#include <maths/perm.h>

template <int dim>
class SimplicialPrism {
  private:
    std::vector<regina::Simplex<dim> *> simplices_;

  public:
    SimplicialPrism(regina::Triangulation<dim> &tri) {
        for (int i = 0; i < dim; ++i) {
            simplices_.push_back(tri.newSimplex());
        }

        regina::Perm<dim + 1> id;

        for (int i = 1; i < dim; ++i) {
            simplices_[i - 1]->join(i - 1, simplices_[i], id);
        }
    }

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

    //template <int facedim>
    //std::array<const regina::Face<dim, facedim + 1> *, facedim + 1>
    //subprism(const regina::Perm<dim> &vertices) const {
    //    static_assert(
    //        facedim >= 0 && facedim < dim - 1,
    //        "SimplicialPrism::faceTimesI: facedim must be in [0, dim-2]");

    //    std::vector<int> faceVerts;
    //    for (int i = 0; i < dim; ++i) {
    //        faceVerts.push_back(vertices[i]);
    //    }
    //    std::sort(faceVerts.begin(), faceVerts.begin() + facedim + 1);

    //    std::array<const regina::Face<dim, facedim + 1> *, facedim + 1> faces;

    //    for (int i = 0; i <= facedim; ++i) {
    //        std::array<int, dim + 1> facetVerts;
    //        std::unordered_set<int> unusedVerts;
    //        for (int j = 0; j <= dim; ++j) {
    //            unusedVerts.insert(j);
    //        }

    //        int j = 0;
    //        for (; j <= i; ++j) {
    //            facetVerts[j] = faceVerts[j] == 0 ? dim : faceVerts[j] - 1;
    //            unusedVerts.erase(facetVerts[j]);
    //        }
    //        for (; j <= facedim + 1; ++j) {
    //            facetVerts[j] = faceVerts[j - 1];
    //            unusedVerts.erase(facetVerts[j]);
    //        }
    //        for (; j <= dim; ++j) {
    //            facetVerts[j] = *unusedVerts.begin();
    //            unusedVerts.erase(unusedVerts.begin());
    //        }

    //        int facetIndex =
    //            regina::Face<dim, facedim + 1>::faceNumber(facetVerts);
    //        faces[i] = simplices_[faceVerts[i]]->template face<facedim + 1>(
    //            facetIndex);
    //    }

    //    return faces;
    //}

    //inline std::array<const regina::Triangle<dim> *, 2>
    //edgeTimesI(const regina::Perm<dim> &vertices) const {
    //    return subprism<1>(vertices);
    //}

    // Glues this prism to `other` along the wall lying over the facet of
    // the "top" simplex obtained by omitting vertex `facet` (resp.
    // `otherFacet` for `other`). Because the underlying triangulation is
    // assumed to be ordered, the correspondence between the walls' two
    // sub-triangulations (and the permutation identifying them) is
    // completely determined by `facet` and `otherFacet` alone.
    void glue(int facet, SimplicialPrism<dim> &other, int otherFacet) {
        if (!(0 <= facet && facet < dim && 0 <= otherFacet &&
              otherFacet < dim))
            throw regina::InvalidArgument(
                "SimplicialPrism::glue(): Invalid faces");

        for (int k = 0; k < dim; ++k) {
            if (k == facet)
                continue;

            int kOther = mapExcluding_(k, facet, otherFacet);

            int myFacet = wallFacet_(k, facet);
            int otherLocalFacet = wallFacet_(kOther, otherFacet);

            std::array<int, dim + 1> image;
            for (int m = 0; m <= dim; ++m) {
                if (m == myFacet) {
                    image[m] = otherLocalFacet;
                    continue;
                }

                // Decode local vertex m of simplices_[k] as (base vertex,
                // is-top-copy), translate the base vertex across the
                // (ordered) base gluing, then re-encode within
                // other.simplices_[kOther]. The top/bottom level is
                // unaffected by the gluing, since it only identifies the
                // base simplices and acts as the identity on the I factor.
                auto [baseVertex, isTop] = decode_(k, m);
                int otherBaseVertex = mapExcluding_(baseVertex, facet, otherFacet);
                image[m] = encode_(otherBaseVertex, isTop);
            }

            simplices_[k]->join(myFacet, other.simplices_[kOther],
                                 regina::Perm<dim + 1>(image));
        }
    }

    // Glues this prism's top facet (the base simplex x {1}) to `next`'s
    // bottom facet (the base simplex x {0}), stacking `next` on top of this
    // prism to extend the interval factor. Since both facets triangulate
    // the same base simplex with no relabelling, the only two sub-simplices
    // involved are simplices_[dim-1] (whose sole non-top vertex is the
    // bottom copy of the top facet) and next.simplices_[0] (whose sole
    // non-bottom vertex is the top copy of the bottom facet); matching up
    // decode_/encode_ across the seam reduces to the cyclic shift
    // i -> (i + 1) mod (dim + 1).
    void stitchTop(SimplicialPrism<dim> &next) {
        std::array<int, dim + 1> image;
        for (int i = 0; i <= dim; ++i)
            image[i] = (i + 1) % (dim + 1);

        simplices_[dim - 1]->join(dim - 1, next.simplices_[0],
                                   regina::Perm<dim + 1>(image));
    }

    // Glues this prism's top facet (the base simplex x {1}) directly onto
    // `coneSimplex`, a single dim-simplex representing the base simplex
    // coned to a new apex point: `coneSimplex`'s local vertex i (for
    // i = 0,...,dim-1) is assumed to equal the base simplex's own vertex i,
    // with local vertex dim the apex (the convention used by
    // CobordismBuilder::cone()). This caps the top off in one step, using
    // the same decode_ used by glue()/stitchTop() rather than a further
    // layer of prism structure.
    void capTop(regina::Simplex<dim> *coneSimplex) {
        std::array<int, dim + 1> image;
        for (int m = 0; m <= dim; ++m) {
            if (m == dim - 1) {
                image[m] = dim;
                continue;
            }
            image[m] = decode_(dim - 1, m).first;
        }

        simplices_[dim - 1]->join(dim - 1, coneSimplex,
                                   regina::Perm<dim + 1>(image));
    }

  private:
    // The unique order-preserving bijection from {0,...}\{a} to {0,...}\{b},
    // evaluated at x (where x != a).
    static int mapExcluding_(int x, int a, int b) {
        int rank = (x < a) ? x : x - 1;
        return (rank < b) ? rank : rank + 1;
    }

    // Decodes local vertex m of simplices_[k] as (base vertex, is-top-copy).
    static std::pair<int, bool> decode_(int k, int m) {
        if (m == dim)
            return {0, true};
        if (k <= m)
            return {m, false};
        return {m + 1, true};
    }

    // Inverse of decode_: the local vertex index (within some simplices_[k])
    // holding (v, isTop), given that this vertex is actually present there.
    static int encode_(int v, bool isTop) {
        if (!isTop)
            return v;
        return v == 0 ? dim : v - 1;
    }

    // The local facet index (within simplices_[k]) of the facet lying in
    // wall(v), for k != v.
    static int wallFacet_(int k, int v) {
        if (k < v)
            return v;
        return v == 0 ? dim : v - 1;
    }
};

#endif // SIMPLICIALPRISM_H
