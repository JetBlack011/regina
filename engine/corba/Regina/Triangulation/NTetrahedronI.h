
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

#ifndef __NTETRAHEDRONI_H
#define __NTETRAHEDRONI_H

#include "triangulation/ntetrahedron.h"

#include "NTetrahedronIDL.h"
#include "ShareableObjectI.h"

class NTetrahedron_i : public virtual POA_Regina::Triangulation::NTetrahedron,
        public ShareableObject_i {
    STANDARD_ENGINE_TYPEDEFS(NTetrahedron_i, NTetrahedron,
            Regina::Triangulation::NTetrahedron)

    protected:
        NTetrahedron_i(::NTetrahedron* newCppPtr) :
                ShareableObject_i(newCppPtr) {
        }
    public:
        STANDARD_NEW_WRAPPER

        virtual Regina::Triangulation::NTetrahedron_ptr
            getAdjacentTetrahedron(CORBA::Long face);
        virtual CORBA::Char getAdjacentTetrahedronGluing(CORBA::Long face);
        virtual CORBA::Long getAdjacentFace(CORBA::Long face);
        virtual CORBA::Boolean hasBoundary();
        virtual Regina::Triangulation::NComponent_ptr getComponent();
        virtual Regina::Triangulation::NVertex_ptr getVertex(
            CORBA::Long vertex);
        virtual Regina::Triangulation::NEdge_ptr getEdge(CORBA::Long edge);
        virtual Regina::Triangulation::NFace_ptr getFace(CORBA::Long face);
        virtual CORBA::Char getEdgeMapping(CORBA::Long edge);
        virtual CORBA::Char getFaceMapping(CORBA::Long face);
        virtual void joinTo(CORBA::Long myFace,
            Regina::Triangulation::NTetrahedron_ptr you,
            CORBA::Char gluing);
        virtual Regina::Triangulation::NTetrahedron_ptr unjoin(
            CORBA::Long myFace);
        virtual void isolate();
        virtual char* getDescription();
        virtual void setDescription(const char* description);
};

#endif

