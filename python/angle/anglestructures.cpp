
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2021, Ben Burton                                   *
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

#include "../pybind11/pybind11.h"
#include "angle/anglestructures.h"
#include "progress/progresstracker.h"
#include "triangulation/dim3.h"
#include "utilities/safeptr.h"
#include "../helpers.h"

using namespace regina::python;
using regina::AngleStructures;
using regina::ProgressTracker;
using regina::Triangulation;

void addAngleStructures(pybind11::module_& m) {
    m.def("makeAngleEquations", regina::makeAngleEquations);

    SafeIterator<AngleStructures>::addBindings(m, "AngleStructureIterator");

    pybind11::class_<AngleStructures, regina::Packet,
            regina::SafePtr<AngleStructures>>(m, "AngleStructures")
        .def(pybind11::init<const Triangulation<3>&, bool, regina::AngleAlg,
                ProgressTracker*>(),
            pybind11::arg(),
            pybind11::arg("tautOnly") = false,
            pybind11::arg("algHints") = regina::AS_ALG_DEFAULT,
            pybind11::arg("tracker") = nullptr)
        .def(pybind11::init<const AngleStructures&>())
        .def("triangulation", &AngleStructures::triangulation,
            pybind11::return_value_policy::reference_internal)
        .def("isTautOnly", &AngleStructures::isTautOnly)
        .def("algorithm", &AngleStructures::algorithm)
        .def("size", &AngleStructures::size)
        .def("structure", &AngleStructures::structure,
            pybind11::return_value_policy::reference_internal)
        .def("__iter__", [](const AngleStructures& list) {
            return SafeIterator(list);
        })
        .def("spansStrict", &AngleStructures::spansStrict)
        .def("spansTaut", &AngleStructures::spansTaut)
        .def_static("enumerate", [](Triangulation<3>& owner, bool tautOnly) {
            // This is deprecated, so we reimplement it here ourselves.
            // This means we can't use the progress tracker variant, which
            // requires threading code internal to the AngleStructures class.
            AngleStructures* ans = new AngleStructures(owner, tautOnly);
            if (regina::Packet* p = owner.packet())
                p->insertChildLast(ans);
            return ans;
        }, pybind11::arg(),
            pybind11::arg("tautOnly") = false)
        .def_static("enumerateTautDD", [](Triangulation<3>& owner) {
            // This is deprecated, so we reimplement it here ourselves.
            AngleStructures* ans = new AngleStructures(owner, true,
                regina::AS_ALG_DD);
            if (regina::Packet* p = owner.packet())
                p->insertChildLast(ans);
            return ans;
        })
        .def_property_readonly_static("typeID", [](pybind11::object) {
            return AngleStructures::typeID;
        })
    ;
}

