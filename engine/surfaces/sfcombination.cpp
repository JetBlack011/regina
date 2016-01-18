
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2014, Ben Burton                                   *
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

/* end stub */

#include "surfaces/sfcombination.h"

#define TYPE_AND 1
#define TYPE_OR 2

namespace regina {

bool NSurfaceFilterCombination::accept(const NNormalSurface& surface) const {
    if (usesAnd_) {
        // Combine all child filters using AND.
        for (NPacket* child = firstChild(); child; child = child->nextSibling())
            if (child->type() == NSurfaceFilter::packetType)
                if (! (dynamic_cast<NSurfaceFilter*>(child)->accept(surface)))
                    return false;
        return true;
    } else {
        // Combine all child filters using OR.
        for (NPacket* child = firstChild(); child; child = child->nextSibling())
            if (child->type() == NSurfaceFilter::packetType)
                if (dynamic_cast<NSurfaceFilter*>(child)->accept(surface))
                    return true;
        return false;
    }
}

void NSurfaceFilterCombination::writeXMLFilterData(std::ostream& out) const {
    out << "    <op type=\"" << (usesAnd_ ? "and" : "or") << "\"/>\n";
}

} // namespace regina

