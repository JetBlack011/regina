
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2018, Ben Burton                                   *
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
#include "surfaces/prism.h"
#include "surfaces/normalsurface.h"
#include "../helpers.h"

using regina::PrismSpec;

void addPrism(pybind11::module& m) {
    // TODO: Should we use a by-value holder type?
    auto c = pybind11::class_<PrismSpec>(m, "PrismSpec")
        .def(pybind11::init<>())
        .def(pybind11::init<unsigned long, int>())
        .def(pybind11::init<const PrismSpec&>())
        .def_readwrite("tetIndex", &PrismSpec::tetIndex)
        .def_readwrite("edge", &PrismSpec::edge)
    ;
    regina::python::add_output_ostream(c);
    regina::python::add_eq_operators(c);

    m.attr("NPrismSpec") = m.attr("PrismSpec");
}

