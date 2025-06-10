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

bool buildLink(regina::Triangulation<3> &tri, std::string &pdcode_str,
               std::vector<const regina::Edge<3> *> &edges);

template <int dim>
regina::Triangulation<dim>
glue(const regina::Triangulation<dim> &tri1, int bdryIndex1,
     const regina::Triangulation<dim> &tri2, int bdryIndex2,
     const regina::Isomorphism<dim - 1> &iso) {
    regina::Triangulation<dim> tri;
    tri.insertTriangulation(tri1);
    tri.insertTriangulation(tri2);
    std::vector<Gluing<dim, dim>> gluings;
    const regina::BoundaryComponent<dim> *bdry1 =
        tri.boundaryComponent(bdryIndex1);
    const regina::BoundaryComponent<dim> *bdry2 =
        tri.boundaryComponent(tri1.countBoundaryComponents() + bdryIndex2);

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

        // std::cout << "bdryIndex = " << i << " to " << iso.simpImage(i) << ",
        // IsoPerm = " << p
        //           << ", i1 = " << i1 << ", i1^-1 = " << i1.inverse()
        //           << ", i2 = " << i2 << "\n";
    }

    for (const auto &gluing : gluings) {
        // std::cout << "Gluing " << gluing.src->index() << " and "
        //           << gluing.dst->index() << " along " << gluing.srcFacet <<
        //           ", "
        //           << gluing.gluing << "\n";
        gluing.src->join(gluing.srcFacet, gluing.dst, gluing.gluing);
    }

    return tri;
}

template <int dim>
regina::Triangulation<dim + 1> cone(regina::Triangulation<dim> &bdryTri) {
    regina::Triangulation<dim + 1> coneTri;

    for (auto s : bdryTri.simplices()) {
        coneTri.newSimplex();
    }

    for (int i = 0; i < bdryTri.simplices().size(); ++i) {
        auto bdrySimplex = bdryTri.simplex(i);
        auto coneSimplex = coneTri.simplex(i);

        for (int f = 0; f <= dim; ++f) {
            if (bdrySimplex->adjacentSimplex(f) == nullptr ||
                coneSimplex->adjacentSimplex(f) != nullptr)
                continue;

            regina::Perm<dim + 1> bdryGluing = bdrySimplex->adjacentGluing(f);
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

    return coneTri;
}

template <int dim>
std::optional<regina::Triangulation<dim + 1>>
thickenHelper(regina::Triangulation<dim> &tri) {
    regina::Triangulation<dim + 1> thickenedTri;
    std::vector<SimplicialPrism<dim + 1>> prisms;
    prisms.reserve(tri.size());
    for (int i = 0; i < tri.size(); ++i) {
        prisms.emplace_back(thickenedTri);
    }
    std::set<std::pair<const regina::Simplex<dim> *, int>> visited;

    int counter = 0;
    for (int i = 0; i < prisms.size(); ++i) {
        for (int f = 0; f < dim + 1; ++f) {
            const regina::Simplex<dim> *s = tri.simplex(i);
            const regina::Simplex<dim> *adj = s->adjacentSimplex(f);

            if (adj == nullptr || visited.contains({s, f}))
                continue;

            int adjFacet = s->adjacentFacet(f);

            prisms[i].glue(f, prisms[adj->index()], adjFacet);
            ++counter;

            visited.insert({s, f});
            visited.insert({adj, adjFacet});
        }
    }

    return thickenedTri;
}

template <int dim>
std::optional<regina::Triangulation<dim + 1>> inline thicken(
    regina::Triangulation<dim> &tri) {
    thickenHelper<dim>(tri);
}

// template <>
// inline std::optional<regina::Triangulation<3>>
// thicken(regina::Triangulation<2> &tri) {
//     // std::cout << "Orienting...\n";
//     // tri.orient();
//     for (int i = 0; i < tri.simplices().size(); ++i) {
//         for (int f = 0; f < 3; ++f) {
//             if (tri.simplex(i)->adjacentSimplex(f) == nullptr)
//                 continue;
//             std::cout << "Simplex " << i << " glued along face " << f << " by
//             "
//                       << tri.simplex(i)->adjacentGluing(f) << " to "
//                       << tri.simplex(i)->adjacentSimplex(f)->index() << "\n";
//         }
//     }
//     return thickenHelper(tri);
// }

template <>
inline std::optional<regina::Triangulation<4>>
thicken(regina::Triangulation<3> &tri) {
    if (!tri.order())
        return std::nullopt;
    return thickenHelper(tri);
}

template <int dim>
std::optional<regina::Triangulation<dim + 1>>
thicken(regina::Triangulation<dim> &tri, int layers) {
    auto layeredTriOpt = thicken(tri);
    if (!layeredTriOpt.has_value())
        return std::nullopt;
    regina::Triangulation<dim + 1> layeredTri = layeredTriOpt.value();

    for (int i = 0; i < layers - 1; ++i) {
        auto thickenedTri = thicken(tri);
        if (!thickenedTri.has_value())
            return std::nullopt;

        auto iso = layeredTri.boundaryComponent(0)->build().isIsomorphicTo(
            thickenedTri.value().boundaryComponent(0)->build());

        if (iso.has_value()) {
            layeredTri = glue(layeredTri, 0, thickenedTri.value(), 0, iso.value());
        } else {
            throw regina::InvalidArgument(
                "Cannot thicken triangulation: boundary components are not "
                "isomorphic; something went wrong while layering.");
        }
    }

    return layeredTri;
}
} // namespace knotbuilder

#endif
