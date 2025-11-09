//
//  knotbuilder.h
//
//  Created by John Teague on 05/10/2024.
//
//  This is adapted from work of Srinivas Vadhiraj, Samantha Ward, Angela
//  Yuan, and Jingyuan Zhang performed for the Texas Experimental Geometry Lab
//  at UT Austin.

#ifndef KNOTBUILDER_H

#define KNOTBUILDER_H

#include <gmpxx.h>
#include <triangulation/dim2.h>
#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <unistd.h>

#include <string>

#include "triangulation/forward.h"
#include "triangulation/generic/triangulation.h"

#include "gluing.h"
#include "simplicialprism.h"

namespace knotbuilder {
using PDCode = std::vector<std::array<int, 4>>;

PDCode parsePDCode(std::string pdcode_str);

class Block {
  private:
    std::vector<regina::Tetrahedron<3> *> core_;
    std::array<std::array<regina::Tetrahedron<3> *, 3>, 4> walls_;

  public:
    Block(regina::Triangulation<3> &tri);

    void glue(size_t myWall, Block &other, size_t otherWall);

    const std::vector<regina::Edge<3> *> getEdges() const;

    friend bool operator==(const Block &lhs, const Block &rhs) {
        return lhs.core_ == rhs.core_ && lhs.walls_ == rhs.walls_;
    }
};

bool buildLink(regina::Triangulation<3> &tri, PDCode pdcode,
               std::vector<const regina::Edge<3> *> &edges);

template <int dim>
regina::Triangulation<dim> &
glueBoundaries(regina::Triangulation<dim> &tri, int bdryIndex1, int bdryIndex2,
               const regina::Isomorphism<dim - 1> &iso) {
    std::vector<Gluing<dim, dim>> gluings;
    const regina::BoundaryComponent<dim> *bdry1 =
        tri.boundaryComponent(bdryIndex1);
    const regina::BoundaryComponent<dim> *bdry2 =
        tri.boundaryComponent(bdryIndex2);

    for (int i = 0; i < bdry1->size(); ++i) {
        // Since these are boundary faces, they belong to exactly 1 simplex
        const auto &emb1 = bdry1->facet(i)->front();
        const auto &emb2 = bdry2->facet(iso.simpImage(i))->front();
        regina::Simplex<dim> *s1 = emb1.simplex();
        regina::Simplex<dim> *s2 = emb2.simplex();

        regina::Perm<dim + 1> i1 = emb1.vertices();
        regina::Perm<dim + 1> i2 = emb2.vertices();
        regina::Perm<dim> p = iso.facetPerm(i);
        std::array<int, dim + 1> isoPerm;

        for (int j = 0; j < dim; ++j) {
            isoPerm[j] = p[j];
        }
        isoPerm[dim] = dim;

        gluings.push_back({s1, emb1.face(), s2, i2 * isoPerm * i1.inverse()});
    }

    for (const auto &gluing : gluings) {
        gluing.src->join(gluing.srcFacet, gluing.dst, gluing.gluing);
    }

    return tri;
}

template <int dim>
regina::Triangulation<dim>
glueTriangulations(const regina::Triangulation<dim> &tri1, int bdryIndex1,
                   const regina::Triangulation<dim> &tri2, int bdryIndex2,
                   const regina::Isomorphism<dim - 1> &iso) {
    regina::Triangulation<dim> tri;
    tri.insertTriangulation(tri1);
    tri.insertTriangulation(tri2);

    return glueBoundaries(tri, bdryIndex1,
                          tri1.countBoundaryComponents() + bdryIndex2, iso);
}

template <int dim>
class CobordismBuilder {
  private:
    using PrismMap = std::unordered_map<const regina::Simplex<dim> *,
                                        SimplicialPrism<dim + 1>>;

    PrismMap topPrisms_;

    const regina::Triangulation<dim> &tri_;
    regina::Triangulation<dim + 1> cob_;

  public:
    CobordismBuilder(const regina::Triangulation<dim> &tri) : tri_(tri) {}

    regina::Triangulation<dim + 1> &cone() {
        regina::Triangulation<dim + 1> coneTri;

        for (auto s : tri_.simplices()) {
            coneTri.newSimplex();
        }

        for (int i = 0; i < tri_.simplices().size(); ++i) {
            auto bdrySimplex = tri_.simplex(i);
            auto coneSimplex = coneTri.simplex(i);

            for (int f = 0; f <= dim; ++f) {
                if (bdrySimplex->adjacentSimplex(f) == nullptr ||
                    coneSimplex->adjacentSimplex(f) != nullptr)
                    continue;

                regina::Perm<dim + 1> bdryGluing =
                    bdrySimplex->adjacentGluing(f);
                regina::Simplex<dim + 1> *adjConeSimplex =
                    coneTri.simplex(bdrySimplex->adjacentSimplex(f)->index());
                std::array<int, dim + 2> coneGluing;
                for (int i = 0; i < dim + 1; ++i) {
                    coneGluing[i] = bdryGluing[i];
                }
                coneGluing[dim + 1] = dim + 1;

                coneSimplex->join(f, adjConeSimplex, coneGluing);
            }
        }

        if (cob_.isConnected()) {
            return cob_ = coneTri;
        }

        cob_.insertTriangulation(coneTri);
        // Glue the cone triangulation to the cobordism triangulation
        auto iso = cob_.boundaryComponent(0)->build().isIsomorphicTo(
            cob_.boundaryComponent(2)->build());
        if (!iso.has_value()) {
            throw regina::InvalidArgument(
                "CobordismBuilder::cone(): Cannot glue cone triangulation to "
                "cobordism triangulation; boundary components are not "
                "isomorphic.");
        }
        glueBoundaries(cob_, 0, 2, iso.value());

        if (!cob_.isValid()) {
            throw regina::InvalidArgument(
                "CobordismBuilder::cone(): Resulting triangulation is not "
                "valid after gluing cone triangulation to cobordism "
                "triangulation.");
        }

        return cob_;
    }

    inline regina::Triangulation<dim + 1> &thicken() { return thicken_(); }

    inline regina::Triangulation<dim + 1> &thicken(int layers) {
        for (int i = 0; i < layers; ++i) {
            thicken_();
        }

        return cob_;
    }

    template <int facedim>
    std::vector<const regina::Face<dim + 1, facedim + 1> *> facesTimesI(
        const std::vector<const regina::Face<dim, facedim> *> &faces) const {
        std::vector<const regina::Face<dim + 1, facedim + 1> *> thickenedFaces;
        thickenedFaces.reserve(faces.size() * (facedim + 1));

        for (const auto &face : faces) {
            const regina::FaceEmbedding<dim, facedim> &emb = face->front();
            const SimplicialPrism<dim + 1> &prism =
                topPrisms_.at(emb.simplex());
            const auto facesTimesI =
                prism.template subprism<facedim>(emb.vertices());
            thickenedFaces.insert(thickenedFaces.end(), facesTimesI.begin(),
                                  facesTimesI.end());
        }

        return thickenedFaces;
    }

    inline std::vector<const regina::Triangle<dim + 1> *>
    edgesTimesI(const std::vector<const regina::Edge<dim> *> &edges) const {
        return facesTimesI(edges);
    }

  private:
    regina::Triangulation<dim + 1> &thicken_() {
        topPrisms_.clear();
        topPrisms_.reserve(tri_.size());
        for (const auto *s : tri_.simplices()) {
            topPrisms_.emplace(s, cob_);
        }
        std::set<std::pair<const regina::Simplex<dim> *, int>> visited;

        for (int i = 0; i < topPrisms_.size(); ++i) {
            for (int facet = 0; facet < dim + 1; ++facet) {
                const regina::Simplex<dim> *s = tri_.simplex(i);
                const regina::Simplex<dim> *adj = s->adjacentSimplex(facet);

                if (adj == nullptr || visited.contains({s, facet}))
                    continue;

                int adjFacet = s->adjacentFacet(facet);

                std::cout << "[*] Gluing thickened triangulation: "
                          << s->index() << " facet " << facet << " to "
                          << adj->index() << " facet " << adjFacet << "\n";
                topPrisms_.at(s).glue(facet, topPrisms_.at(adj), adjFacet);

                visited.insert({s, facet});
                visited.insert({adj, adjFacet});
            }
        }

        if (cob_.isConnected()) {
            return cob_;
        }

        // Glue the thickened triangulation to the cobordism triangulation
        auto iso = cob_.boundaryComponent(1)->build().isIsomorphicTo(
            cob_.boundaryComponent(2)->build());
        if (!iso.has_value()) {
            throw regina::InvalidArgument(
                "CobordismBuilder::thicken(): Cannot glue thickened "
                "triangulation to cobordism triangulation; boundary components "
                "are not isomorphic.");
        }

        std::cout << "[*] Gluing thickened triangulation\n";

        glueBoundaries(cob_, 1, 2, iso.value());
        std::cout << "[*] Done gluing thickened triangulation\n";

        if (!cob_.isValid()) {
            throw regina::InvalidArgument(
                "CobordismBuilder::thicken(): Resulting triangulation is not "
                "valid after gluing thickened triangulation to cobordism "
                "triangulation.");
        }

        return cob_;
    }
};

template <>
inline regina::Triangulation<4> &CobordismBuilder<3>::thicken() {
    // if (!tri_.isOrdered() && !tri_.order())
    //     throw regina::InvalidArgument(
    //         "CobordismBuilder<3>::thicken(): Could not thicken! Triangulation
    //         " "is unorderable.");
    return thicken_();
}

template <>
inline regina::Triangulation<4> &CobordismBuilder<3>::thicken(int layers) {
    // if (!tri_.isOrdered() && !tri_.order())
    //     throw regina::InvalidArgument("CobordismBuilder<3>::thicken(int): "
    //                                   "Could not thicken! Triangulation "
    //                                   "is unorderable.");
    for (int i = 0; i < layers; ++i) {
        thicken_();
    }

    return cob_;
}

} // namespace knotbuilder

#endif
