
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

#ifndef __NGROUPPRESENTATIONI_H
#define __NGROUPPRESENTATIONI_H

#include "algebra/ngrouppresentation.h"

#include "NGroupPresentationIDL.h"
#include "ShareableObjectI.h"

class ::NGroupExpression;
class ::NGroupPresentation;

class NGroupPresentation_i :
        public virtual POA_Regina::Algebra::NGroupPresentation,
        public ShareableObject_i {
    STANDARD_ENGINE_TYPEDEFS(NGroupPresentation_i, NGroupPresentation,
            Regina::Algebra::NGroupPresentation)

    protected:
        NGroupPresentation_i(::NGroupPresentation* newCppPtr) :
                ShareableObject_i(newCppPtr) {
        }
    public:
        STANDARD_NEW_WRAPPER

        virtual CORBA::Long addGenerator(CORBA::Long numToAdd);
        virtual void addRelation(Regina::Algebra::NGroupExpression_ptr rel);
        virtual CORBA::Long getNumberOfGenerators();
        virtual CORBA::Long getNumberOfRelations();
        virtual Regina::Algebra::NGroupExpression_ptr getRelation(
            CORBA::Long index);
        virtual CORBA::Boolean intelligentSimplify();
        virtual char* recogniseGroup();
};

#endif

