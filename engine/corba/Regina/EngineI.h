
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

#ifndef __ENGINEI_H
#define __ENGINEI_H

#include "EngineIDL.h"

class Engine_i : public virtual POA_Regina::Engine,
        public PortableServer::RefCountServantBase {
    private:
        CORBA::ORB_var orb;
            /**< The ORB used by this engine. */

    public:
        /**
         * Create a new CORBA calculation engine.
         *
         * @param newORB the ORB that this engine should use.
         */
        Engine_i(CORBA::ORB_var newORB) : orb(newORB) {
        }
        virtual ~Engine_i() {
        }

        virtual Regina::Algebra::NAbelianGroup_ptr newNAbelianGroup_();
        virtual Regina::Algebra::NAbelianGroup_ptr
            newNAbelianGroup_NAbelianGroup(
            Regina::Algebra::NAbelianGroup_ptr cloneMe);
        virtual Regina::Algebra::NGroupExpression_ptr newNGroupExpression_();
        virtual Regina::Algebra::NGroupExpression_ptr
            newNGroupExpression_NGroupExpression(
            Regina::Algebra::NGroupExpression_ptr cloneMe);
        virtual Regina::Algebra::NGroupPresentation_ptr
            newNGroupPresentation_();
        virtual Regina::Algebra::NGroupPresentation_ptr
            newNGroupPresentation_NGroupPresentation(
            Regina::Algebra::NGroupPresentation_ptr cloneMe);

        virtual Regina::File::NFile_ptr newNFile();

        virtual Regina::Maths::NMatrixInt_ptr newNMatrixInt_long_long(
            CORBA::Long rows, CORBA::Long columns);
        virtual Regina::Maths::NMatrixInt_ptr newNMatrixInt_NMatrixInt(
            Regina::Maths::NMatrixInt_ptr cloneMe);

        virtual Regina::Packet::NContainer_ptr newNContainer();
        virtual Regina::Packet::NScript_ptr newNScript();
        virtual Regina::Packet::NText_ptr newNText_();
        virtual Regina::Packet::NText_ptr newNText_string(const char* text);

        virtual Regina::Progress::NProgressManager_ptr newNProgressManager();

        virtual Regina::Subcomplex::NLayeredChain_ptr
            newNLayeredChain_NTetrahedron_NPerm(
            Regina::Triangulation::NTetrahedron_ptr tet, CORBA::Char roles);
        virtual Regina::Subcomplex::NLayeredChain_ptr
            newNLayeredChain_NLayeredChain(
            Regina::Subcomplex::NLayeredChain_ptr cloneMe);
        virtual Regina::Subcomplex::NLensSpace_ptr newNLensSpace_long_long(
            CORBA::Long p, CORBA::Long q);
        virtual Regina::Subcomplex::NLensSpace_ptr newNLensSpace_NLensSpace(
            Regina::Subcomplex::NLensSpace_ptr cloneMe);
        virtual Regina::Subcomplex::NSFS_ptr newNSFS_();
        virtual Regina::Subcomplex::NSFS_ptr newNSFS_long_boolean_long(
            CORBA::Long genus, CORBA::Boolean orient, CORBA::Long punctures);
        virtual Regina::Subcomplex::NSFS_ptr newNSFS_NSFS(
            Regina::Subcomplex::NSFS_ptr cloneMe);

        virtual Regina::Surfaces::NNormalSurfaceList_ptr newNNormalSurfaceList(
            Regina::Triangulation::NTriangulation_ptr owner,
            CORBA::Long flavour, CORBA::Boolean isEmbeddedOnly);
        virtual Regina::Surfaces::NSurfaceFilter_ptr newNSurfaceFilter_();
        virtual Regina::Surfaces::NSurfaceFilter_ptr
            newNSurfaceFilter_NSurfaceFilter(
            Regina::Surfaces::NSurfaceFilter_ptr cloneMe);
        virtual Regina::Surfaces::NSurfaceFilterCombination_ptr
            newNSurfaceFilterCombination_();
        virtual Regina::Surfaces::NSurfaceFilterCombination_ptr
            newNSurfaceFilterCombination_NSurfaceFilterCombination(
            Regina::Surfaces::NSurfaceFilterCombination_ptr cloneMe);
        virtual Regina::Surfaces::NSurfaceFilterProperties_ptr
            newNSurfaceFilterProperties_();
        virtual Regina::Surfaces::NSurfaceFilterProperties_ptr
            newNSurfaceFilterProperties_NSurfaceFilterProperties(
            Regina::Surfaces::NSurfaceFilterProperties_ptr cloneMe);
        virtual Regina::Surfaces::NSurfaceSubset_ptr newNSurfaceSubset(
            Regina::Surfaces::NSurfaceSet_ptr set,
            Regina::Surfaces::NSurfaceFilter_ptr filter);

        virtual Regina::Triangulation::NTetrahedron_ptr newNTetrahedron_();
        virtual Regina::Triangulation::NTetrahedron_ptr newNTetrahedron_string(
            const char* desc);
        virtual Regina::Triangulation::NTriangulation_ptr newNTriangulation_();
        virtual Regina::Triangulation::NTriangulation_ptr
            newNTriangulation_NTriangulation(
            Regina::Triangulation::NTriangulation_ptr cloneMe);

        virtual CORBA::Long formCensus(Regina::Packet::NPacket_ptr parent,
            CORBA::Long nTetrahedra, CORBA::Char finiteness,
            CORBA::Char orientability, CORBA::Char boundary,
            CORBA::Long nBdryFaces,
            Regina::Progress::NProgressManager_ptr manager);
        virtual Regina::Subcomplex::NPillowTwoSphere_ptr
            formsPillowTwoSphere(Regina::Triangulation::NFace_ptr face1,
            Regina::Triangulation::NFace_ptr face2);
        virtual Regina::Subcomplex::NSnappedTwoSphere_ptr
            formsSnappedTwoSphere_NSnappedBall(
            Regina::Subcomplex::NSnappedBall_ptr p1,
            Regina::Subcomplex::NSnappedBall_ptr p2);
        virtual Regina::Subcomplex::NSnappedTwoSphere_ptr
            formsSnappedTwoSphere_NTetrahedron(
            Regina::Triangulation::NTetrahedron_ptr p1,
            Regina::Triangulation::NTetrahedron_ptr p2);
        virtual CORBA::Long getVersionMajor();
        virtual CORBA::Long getVersionMinor();
        virtual char* getVersionString();
        virtual Regina::Subcomplex::NAugTriSolidTorus_ptr
            isAugTriSolidTorus(Regina::Triangulation::NComponent_ptr comp);
        virtual Regina::Subcomplex::NLayeredLensSpace_ptr
            isLayeredLensSpace(Regina::Triangulation::NComponent_ptr comp);
        virtual Regina::Subcomplex::NLayeredLoop_ptr
            isLayeredLoop(Regina::Triangulation::NComponent_ptr comp);
        virtual Regina::Subcomplex::NLayeredSolidTorus_ptr
            isLayeredSolidTorusBase(Regina::Triangulation::NTetrahedron_ptr tet);
        virtual Regina::Subcomplex::NSnappedBall_ptr
            isSnappedBall(Regina::Triangulation::NTetrahedron_ptr tet);
        virtual Regina::Subcomplex::NSpiralSolidTorus_ptr
            isSpiralSolidTorus(Regina::Triangulation::NTetrahedron_ptr tet,
            CORBA::Char vertexRoles);
        virtual Regina::Subcomplex::NTriSolidTorus_ptr
            isTriSolidTorus(Regina::Triangulation::NTetrahedron_ptr tet,
            CORBA::Char vertexRoles);
        virtual Regina::Maths::NMatrixInt_ptr makeMatchingEquations(
            Regina::Triangulation::NTriangulation_ptr triangulation,
            CORBA::Long flavour);
        virtual Regina::Packet::NPacket_ptr readFromFile(
            const char* fileName);
        virtual Regina::Triangulation::NTriangulation_ptr readSnapPea(
            const char* fileName);
        virtual void smithNormalForm(Regina::Maths::NMatrixInt_ptr matrix);
        virtual CORBA::Long testEngine(CORBA::Long value);
        virtual CORBA::Boolean writeToFile(const char* fileName,
            Regina::Packet::NPacket_ptr packet);
        virtual CORBA::Boolean writeSnapPea(const char* fileName,
            Regina::Triangulation::NTriangulation_ptr tri);
};

#endif

