
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

#include "NGroupExpressionI.h"

CORBA::Long NGroupExpression_i::getNumberOfTerms() {
    return MY_ENGINE_OBJECT->getNumberOfTerms();
}
CORBA::Long NGroupExpression_i::getGenerator(CORBA::Long index) {
    return MY_ENGINE_OBJECT->getGenerator(index);
}
CORBA::Long NGroupExpression_i::getExponent(CORBA::Long index) {
    return MY_ENGINE_OBJECT->getExponent(index);
}
void NGroupExpression_i::addTermFirst(CORBA::Long gen, CORBA::Long exp) {
    MY_ENGINE_OBJECT->addTermFirst(gen, exp);
}
void NGroupExpression_i::addTermLast(CORBA::Long gen, CORBA::Long exp) {
    MY_ENGINE_OBJECT->addTermLast(gen, exp);
}
Regina::Algebra::NGroupExpression_ptr NGroupExpression_i::inverse() {
    return NGroupExpression_i::newWrapper(MY_ENGINE_OBJECT->inverse());
}
Regina::Algebra::NGroupExpression_ptr NGroupExpression_i::power(
        CORBA::Long exp) {
    return NGroupExpression_i::newWrapper(MY_ENGINE_OBJECT->power(exp));
}
CORBA::Boolean NGroupExpression_i::simplify(CORBA::Boolean cyclic) {
    return MY_ENGINE_OBJECT->simplify(cyclic);
}
CORBA::Boolean NGroupExpression_i::substitute(CORBA::Long gen,
        Regina::Algebra::NGroupExpression_ptr exp, CORBA::Boolean cyclic) {
    return MY_ENGINE_OBJECT->substitute(gen,
        *GET_ENGINE_OBJECT(NGroupExpression, exp), cyclic);
}

