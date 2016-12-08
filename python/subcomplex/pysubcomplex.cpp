
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Python Interface                                                      *
 *                                                                        *
 *  Copyright (c) 1999-2016, Ben Burton                                   *
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

void addAugTriSolidTorus();
void addBlockedSFS();
void addBlockedSFSLoop();
void addBlockedSFSPair();
void addBlockedSFSTriple();
void addL31Pillow();
void addLayeredChain();
void addLayeredChainPair();
void addLayeredLensSpace();
void addLayeredLoop();
void addLayeredSolidTorus();
void addLayeredSurfaceBundle();
void addLayering();
void addPillowTwoSphere();
void addPluggedTorusBundle();
void addPlugTriSolidTorus();
void addSatAnnulus();
void addSatBlock();
void addSatBlockTypes();
void addSatRegion();
void addSnapPeaCensusTri();
void addSnappedBall();
void addSnappedTwoSphere();
void addNSpiralSolidTorus();
void addNStandardTriangulation();
void addNTriSolidTorus();
void addNTrivialTri();
void addNTxICore();

void addSubcomplexClasses() {
    addNStandardTriangulation();
    addAugTriSolidTorus();
    addL31Pillow();
    addLayeredChain();
    addLayeredChainPair();
    addLayeredLensSpace();
    addLayeredLoop();
    addLayeredSolidTorus();
    addLayeredSurfaceBundle();
    addLayering();
    addPillowTwoSphere();
    addPlugTriSolidTorus();
    addSnapPeaCensusTri();
    addSnappedBall();
    addSnappedTwoSphere();
    addNSpiralSolidTorus();
    addNTriSolidTorus();
    addNTrivialTri();
    addNTxICore();
    addSatAnnulus();
    addSatBlock();
    addSatBlockTypes();
    addSatRegion();
    addBlockedSFS();
    addBlockedSFSLoop();
    addBlockedSFSPair();
    addBlockedSFSTriple();
    addPluggedTorusBundle();
}

