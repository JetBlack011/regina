
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
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

#include "triangulation/facenumbering.h"

namespace regina::detail {

    const int FaceNumberingImpl<3, 1>::edgeNumber[4][4] = {
        { -1, 0, 1, 2 },
        {  0,-1, 3, 4 },
        {  1, 3,-1, 5 },
        {  2, 4, 5,-1 }
    };

    const int FaceNumberingImpl<3, 1>::edgeVertex[6][2] = {
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 1, 2 },
        { 1, 3 },
        { 2, 3 }
    };

    const int FaceNumberingImpl<4, 1>::edgeNumber[5][5] = {
        { -1, 0, 1, 2, 3 },
        { 0, -1, 4, 5, 6 },
        { 1, 4, -1, 7, 8 },
        { 2, 5, 7, -1, 9 },
        { 3, 6, 8, 9, -1 }};

    const int FaceNumberingImpl<4, 1>::edgeVertex[10][2] = {
        { 0, 1 },
        { 0, 2 },
        { 0, 3 },
        { 0, 4 },
        { 1, 2 },
        { 1, 3 },
        { 1, 4 },
        { 2, 3 },
        { 2, 4 },
        { 3, 4 }};

    const int FaceNumberingImpl<4, 2>::triangleNumber[5][5][5] = {
        { { -1,-1,-1,-1,-1 },
          { -1,-1, 9, 8, 7 },
          { -1, 9,-1, 6, 5 },
          { -1, 8, 6,-1, 4 },
          { -1, 7, 5, 4,-1 } },

        { { -1,-1, 9, 8, 7 },
          { -1,-1,-1,-1,-1 },
          {  9,-1,-1, 3, 2 },
          {  8,-1, 3,-1, 1 },
          {  7,-1, 2, 1,-1 } },

        { { -1, 9,-1, 6, 5 },
          {  9,-1,-1, 3, 2 },
          { -1,-1,-1,-1,-1 },
          {  6, 3,-1,-1, 0 },
          {  5, 2,-1, 0,-1 } },

        { { -1, 8, 6,-1, 4 },
          {  8,-1, 3,-1, 1 },
          {  6, 3,-1,-1, 0 },
          { -1,-1,-1,-1,-1 },
          {  4, 1, 0,-1,-1 } },

        { { -1, 7, 5, 4,-1 },
          {  7,-1, 2, 1,-1 },
          {  5, 2,-1, 0,-1 },
          {  4, 1, 0,-1,-1 },
          { -1,-1,-1,-1,-1 } } };

    const int FaceNumberingImpl<4, 2>::triangleVertex[10][3] = {
        { 2, 3, 4 },
        { 1, 3, 4 },
        { 1, 2, 4 },
        { 1, 2, 3 },
        { 0, 3, 4 },
        { 0, 2, 4 },
        { 0, 2, 3 },
        { 0, 1, 4 },
        { 0, 1, 3 },
        { 0, 1, 2 } };

} // namespace regina::detail

