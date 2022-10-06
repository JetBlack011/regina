
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
#include "maths/matrix.h"
#include "maths/vector.h"
#include "../helpers.h"

using pybind11::overload_cast;
using regina::Matrix;

void addMatrixBool(pybind11::module_& m) {
    auto c = pybind11::class_<Matrix<bool>>(m, "MatrixBool")
        .def(pybind11::init<size_t>())
        .def(pybind11::init<size_t, size_t>())
        .def(pybind11::init<const Matrix<bool>&>())
        .def(pybind11::init([](pybind11::list l) {
            size_t rows = l.size();
            if (rows == 0)
                throw regina::InvalidArgument(
                    "The number of rows must be strictly positive");

            Matrix<bool>* m = nullptr;
            size_t cols = 0; // zero is unnecessary but silences warnings

            pybind11::list row;
            for (size_t i = 0; i < rows; ++i) {
                try {
                    row = l[i].cast<pybind11::list>();
                } catch (const pybind11::cast_error&) {
                    delete m;
                    throw regina::InvalidArgument(
                        "Each row must be given as a separate Python list");
                }
                if (i == 0) {
                    cols = row.size();
                    if (cols == 0)
                        throw regina::InvalidArgument(
                            "The number of columns must be strictly positive");
                    m = new Matrix<bool>(rows, cols);
                } else if (row.size() != cols) {
                    delete m;
                    throw regina::InvalidArgument(
                        "All rows must be given as lists of the same size");
                }
                for (size_t j = 0; j < cols; ++j) {
                    try {
                        m->entry(i, j) = row[j].cast<bool>();
                    } catch (const pybind11::cast_error&) {
                        delete m;
                        throw regina::InvalidArgument(
                            "Matrix element not convertible to a boolean");
                    }
                }
            }

            return m;
        }))
        .def("initialise", &Matrix<bool>::initialise)
        .def("swap", &Matrix<bool>::swap)
        .def("rows", &Matrix<bool>::rows)
        .def("columns", &Matrix<bool>::columns)
        .def("entry", overload_cast<size_t, size_t>(&Matrix<bool>::entry),
            pybind11::return_value_policy::reference_internal)
        .def("set", [](Matrix<bool>& m, size_t row, size_t col, bool value){
            m.entry(row, col) = value;
        })
        .def("transpose", &Matrix<bool>::transpose)
        .def("swapRows", &Matrix<bool>::swapRows)
        .def("swapCols", &Matrix<bool>::swapCols)
    ;
    regina::python::add_output(c);
    regina::python::add_eq_operators(c);

    regina::python::add_global_swap<Matrix<bool>>(m);
}

