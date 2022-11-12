
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2022, Ben Burton                                   *
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
#include "maths/permgroup.h"
#include "../helpers.h"

using regina::Perm;
using regina::PermGroup;

template <int n>
void addPermGroup(pybind11::module_& m, const char* name) {
    auto c = pybind11::class_<PermGroup<n>>(m, name, RDOC_TODO)
        .def(pybind11::init<>(), RDOC_TODO)
        .def(pybind11::init<const PermGroup<n>&>(), RDOC_TODO)
        .def(pybind11::init<regina::NamedPermGroup>(), RDOC_TODO)
        .def(pybind11::init<int>(), RDOC_TODO)
        .def("size", &PermGroup<n>::size, RDOC_TODO)
        .def("contains", &PermGroup<n>::contains, RDOC_TODO)
        .def("__iter__", [](const PermGroup<n>& g) {
            return g.begin();
        }, pybind11::keep_alive<0, 1>(), // iterator keeps group alive
            RDOC_TODO)
    ;
    // TODO regina::python::add_output_basic(c, RDOC_TODO);
    // TODO regina::python::add_eq_operators(c, RDOC_TODO, RDOC_TODO);

    using iterator = typename PermGroup<n>::iterator;
    auto it = pybind11::class_<iterator>(c, "iterator", RDOC_TODO)
        .def("__next__", [](iterator& it) {
            if (it)
                return *it++;
            else
                throw pybind11::stop_iteration();
        }, RDOC_TODO);
    regina::python::add_eq_operators(it, RDOC_TODO, RDOC_TODO);

    c.attr("const_iterator") = c.attr("iterator");
}

void addPermGroup(pybind11::module_& m) {
    pybind11::enum_<regina::NamedPermGroup>(m, "NamedPermGroup", RDOC_TODO)
        .value("PERM_GROUP_TRIVIAL", regina::PERM_GROUP_TRIVIAL, RDOC_TODO)
        .value("PERM_GROUP_SYMMETRIC", regina::PERM_GROUP_SYMMETRIC, RDOC_TODO)
        .value("PERM_GROUP_ALTERNATING", regina::PERM_GROUP_ALTERNATING,
            RDOC_TODO)
        .export_values()
        ;

    addPermGroup<2>(m, "PermGroup2");
    addPermGroup<3>(m, "PermGroup3");
    addPermGroup<4>(m, "PermGroup4");
    addPermGroup<5>(m, "PermGroup5");
    addPermGroup<6>(m, "PermGroup6");
    addPermGroup<7>(m, "PermGroup7");
    addPermGroup<8>(m, "PermGroup8");
    addPermGroup<9>(m, "PermGroup9");
    addPermGroup<10>(m, "PermGroup10");
    addPermGroup<11>(m, "PermGroup11");
    addPermGroup<12>(m, "PermGroup12");
    addPermGroup<13>(m, "PermGroup13");
    addPermGroup<14>(m, "PermGroup14");
    addPermGroup<15>(m, "PermGroup15");
    addPermGroup<16>(m, "PermGroup16");
}

