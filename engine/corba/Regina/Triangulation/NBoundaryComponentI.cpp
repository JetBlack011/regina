
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2002, Ben Burton                                   *
 *  For further details contact Ben Burton (benb@acm.org).                *
 *                                                                        *
 *  This program is free software; you can redistribute it and/or         *
 *  modify it under the terms of the GNU General Public License as        *
 *  published by the Free Software Foundation; either version 2 of the    *
 *  License, or (at your option) any later version.                       *
 *                                                                        *
 *  This program is distributed in the hope that it will be useful, but   *
 *  WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU     *
 *  General Public License for more details.                              *
 *                                                                        *
 *  You should have received a copy of the GNU General Public             *
 *  License along with this program; if not, write to the Free            *
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,        *
 *  MA 02111-1307, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

#include "NBoundaryComponentI.h"
#include "NVertexI.h"
#include "NEdgeI.h"
#include "NFaceI.h"

CORBA::Long NBoundaryComponent_i::getEulerCharacteristic() {
    return MY_ENGINE_OBJECT->getEulerCharacteristic();
}
CORBA::Boolean NBoundaryComponent_i::isIdeal() {
    return MY_ENGINE_OBJECT->isIdeal();
}
CORBA::Boolean NBoundaryComponent_i::isOrientable() {
    return MY_ENGINE_OBJECT->isOrientable();
}
CORBA::Long NBoundaryComponent_i::getNumberOfFaces() {
    return MY_ENGINE_OBJECT->getNumberOfFaces();
}
CORBA::Long NBoundaryComponent_i::getNumberOfEdges() {
    return MY_ENGINE_OBJECT->getNumberOfEdges();
}
CORBA::Long NBoundaryComponent_i::getNumberOfVertices() {
    return MY_ENGINE_OBJECT->getNumberOfVertices();
}
Regina::Triangulation::NFace_ptr NBoundaryComponent_i::getFace(
        CORBA::Long index) {
    return NFace_i::newWrapper(MY_ENGINE_OBJECT->getFace(index));
}
Regina::Triangulation::NEdge_ptr NBoundaryComponent_i::getEdge(
        CORBA::Long index) {
    return NEdge_i::newWrapper(MY_ENGINE_OBJECT->getEdge(index));
}
Regina::Triangulation::NVertex_ptr NBoundaryComponent_i::getVertex(
        CORBA::Long index) {
    return NVertex_i::newWrapper(MY_ENGINE_OBJECT->getVertex(index));
}

