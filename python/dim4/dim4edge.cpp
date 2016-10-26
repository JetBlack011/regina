
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  As an exception, when this program is distributed through (i) the     *
 *  App Store by Apple Inc.; (ii) the Mac App Store by Apple Inc.; or     *
 *  (iii) Google Play by Google Inc., then that store may impose any      *
 *  digital rights management, device limits and/or redistribution        *
 *  restrictions that are required by its terms of service.               *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

#include <boost/python.hpp>
#include "triangulation/dim2.h"
#include "triangulation/dim4.h"
#include "../globalarray.h"
#include "../helpers.h"
#include "../safeheldtype.h"
#include "../generic/facehelper.h"

using namespace boost::python;
using namespace regina::python;
using regina::Dim4Edge;
using regina::EdgeEmbedding;
using regina::Face;
using regina::FaceEmbedding;
using regina::python::GlobalArray;
using regina::python::GlobalArray2D;

namespace {
    GlobalArray2D<int> Dim4Edge_edgeNumber(Dim4Edge::edgeNumber, 5);
    GlobalArray2D<int> Dim4Edge_edgeVertex(Dim4Edge::edgeVertex, 10);

    boost::python::list Dim4Edge_embeddings_list(const Dim4Edge* e) {
        boost::python::list ans;
        for (auto& emb : *e)
            ans.append(emb);
        return ans;
    }

    regina::Triangulation<2>* edge_buildLink(const Dim4Edge* e) {
        return new regina::Triangulation<2>(*(e->buildLink()));
    }

    boost::python::tuple edge_buildLinkDetail_bool(const Dim4Edge* e,
            bool labels = true) {
        regina::Dim4Isomorphism* iso;
        regina::Triangulation<2>* link = new regina::Triangulation<2>(
            *(e->buildLinkDetail(labels, &iso)));
        return boost::python::make_tuple(
            boost::python::object(boost::python::handle<>(
                regina::python::to_held_type<>::
                apply<regina::Triangulation<2>*>::type()(link))),
            boost::python::object(boost::python::handle<>(
                boost::python::manage_new_object::
                apply<regina::Dim4Isomorphism*>::type()(iso))));
    }

    boost::python::tuple edge_buildLinkDetail_void(const Dim4Edge* e) {
        return edge_buildLinkDetail_bool(e);
    }
}

void addDim4Edge() {
    class_<FaceEmbedding<4, 1>>("FaceEmbedding4_1",
            init<regina::Pentachoron<4>*, int>())
        .def(init<const EdgeEmbedding<4>&>())
        .def("simplex", &EdgeEmbedding<4>::simplex,
            return_value_policy<reference_existing_object>())
        .def("pentachoron", &EdgeEmbedding<4>::pentachoron,
            return_value_policy<reference_existing_object>())
        .def("face", &EdgeEmbedding<4>::face)
        .def("edge", &EdgeEmbedding<4>::edge)
        .def("vertices", &EdgeEmbedding<4>::vertices)
        .def(regina::python::add_output())
        .def(regina::python::add_eq_operators())
    ;

    {
        scope s = class_<Face<4, 1>, std::auto_ptr<Face<4, 1>>,
                boost::noncopyable>("Face4_1", no_init)
            .def("index", &Dim4Edge::index)
            .def("embeddings", Dim4Edge_embeddings_list)
            .def("embedding", &Dim4Edge::embedding,
                return_internal_reference<>())
            .def("front", &Dim4Edge::front,
                return_internal_reference<>())
            .def("back", &Dim4Edge::back,
                return_internal_reference<>())
            .def("triangulation", &Dim4Edge::triangulation,
                return_value_policy<to_held_type<>>())
            .def("component", &Dim4Edge::component,
                return_value_policy<reference_existing_object>())
            .def("boundaryComponent", &Dim4Edge::boundaryComponent,
                return_value_policy<reference_existing_object>())
            .def("face", &regina::python::face<Dim4Edge, 1, int>)
            .def("vertex", &Dim4Edge::vertex,
                return_value_policy<reference_existing_object>())
            .def("faceMapping", &regina::python::faceMapping<Dim4Edge, 1, 5>)
            .def("vertexMapping", &Dim4Edge::vertexMapping)
            .def("degree", &Dim4Edge::degree)
            .def("isBoundary", &Dim4Edge::isBoundary)
            .def("isLinkOrientable", &Dim4Edge::isLinkOrientable)
            .def("isValid", &Dim4Edge::isValid)
            .def("hasBadIdentification", &Dim4Edge::hasBadIdentification)
            .def("hasBadLink", &Dim4Edge::hasBadLink)
            .def("buildLink", &edge_buildLink,
                return_value_policy<to_held_type<>>())
            .def("buildLinkDetail", edge_buildLinkDetail_void)
            .def("buildLinkDetail", edge_buildLinkDetail_bool)
            .def("ordering", &Dim4Edge::ordering)
            .def("faceNumber", &Dim4Edge::faceNumber)
            .def("containsVertex", &Dim4Edge::containsVertex)
            .def(regina::python::add_output())
            .def(regina::python::add_eq_operators())
            .staticmethod("ordering")
            .staticmethod("faceNumber")
            .staticmethod("containsVertex")
        ;

        s.attr("edgeNumber") = &Dim4Edge_edgeNumber;
        s.attr("edgeVertex") = &Dim4Edge_edgeVertex;
    }

    scope().attr("Dim4EdgeEmbedding") = scope().attr("FaceEmbedding4_1");
    scope().attr("EdgeEmbedding4") = scope().attr("FaceEmbedding4_1");
    scope().attr("Dim4Edge") = scope().attr("Face4_1");
    scope().attr("Edge4") = scope().attr("Face4_1");
}

