
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

#include "NFileI.h"
#include "NPacketI.h"

CORBA::Boolean NFile_i::open(const char* fileName, CORBA::Long openMode) {
    // Don't allow remote clients to access the server's filesystem!
    return false;
    /* return MY_ENGINE_OBJECT->open(fileName, (::NFile::mode)openMode); */
}
void NFile_i::close() {
    MY_ENGINE_OBJECT->close();
}
CORBA::Long NFile_i::getOpenMode() {
    return MY_ENGINE_OBJECT->getOpenMode();
}
CORBA::Long NFile_i::getMajorVersion() {
    return MY_ENGINE_OBJECT->getMajorVersion();
}
CORBA::Long NFile_i::getMinorVersion() {
    return MY_ENGINE_OBJECT->getMinorVersion();
}
CORBA::Boolean NFile_i::versionEarlierThan(CORBA::Long major,
        CORBA::Long minor) {
    return MY_ENGINE_OBJECT->versionEarlierThan(major, minor);
}
void NFile_i::writePacketTree(Regina::Packet::NPacket_ptr packet) {
    MY_ENGINE_OBJECT->writePacketTree(GET_ENGINE_OBJECT(NPacket, packet));
}
Regina::Packet::NPacket_ptr NFile_i::readPacketTree_() {
    return NPacket_i::newWrapper(MY_ENGINE_OBJECT->readPacketTree());
}
Regina::Packet::NPacket_ptr NFile_i::readPacketTree_NPacket(
        Regina::Packet::NPacket_ptr parent) {
    return NPacket_i::newWrapper(MY_ENGINE_OBJECT->
        readPacketTree(GET_ENGINE_OBJECT(NPacket, parent)));
}

