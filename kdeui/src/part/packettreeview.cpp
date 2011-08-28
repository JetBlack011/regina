
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2009, Ben Burton                                   *
 *  For further details contact Ben Burton (bab@debian.org).              *
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
 *  Software Foundation, Inc., 51 Franklin St, Fifth Floor, Boston,       *
 *  MA 02110-1301, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

// Regina core includes:
#include "packet/npacket.h"

// UI includes:
#include "eventids.h"
#include "packetmanager.h"
#include "packettreeview.h"
#include "regevents.h"
#include "reginapart.h"

#include <qapplication.h>
#include <qevent.h>
#include <QHeaderView>
#include <QTreeWidget>
#include <kdebug.h>
#include <klocale.h>

using regina::NPacket;

PacketTreeItem::PacketTreeItem(PacketTreeView* parent, NPacket* realPacket) :
        QTreeWidgetItem(parent), packet(realPacket), tree(parent) {
    init();
}

PacketTreeItem::PacketTreeItem(PacketTreeItem* parent,
        NPacket* realPacket) :
        QTreeWidgetItem(parent), packet(realPacket), tree(parent->tree) {
    init();
}

PacketTreeItem::PacketTreeItem(PacketTreeView* parent,
        QTreeWidgetItem* after, NPacket* realPacket) :
        QTreeWidgetItem(parent, after), packet(realPacket),
        tree(parent) {
    init();
}

PacketTreeItem::PacketTreeItem(PacketTreeItem* parent,
        QTreeWidgetItem* after, NPacket* realPacket) :
        QTreeWidgetItem(parent, after), packet(realPacket), tree(parent->tree) {
    init();
}

void PacketTreeItem::init() {
    packet->listen(this);
    refreshLabel();
    setIcon(0, PacketManager::iconSmall(packet, true));
    isEditable = packet->isPacketEditable();
}

void PacketTreeItem::fill() {
    PacketTreeItem* childTree = 0;
    for (NPacket* p = packet->getFirstTreeChild(); p;
            p = p->getNextTreeSibling()) {
        if (childTree)
            childTree = new PacketTreeItem(this, childTree, p);
        else
            childTree = new PacketTreeItem(this, p);
        childTree->fill();
    }
}

void PacketTreeItem::refreshSubtree() {
    // Is this a stale node in the tree?
    if (! packet) {
        // Yes, it's a stale node.  Delete all of its children.
        while (childCount())
            delete child(0);
        return;
    }

    // No, we're looking at a real packet.
    // Run through the child packets and child nodes and ensure they
    // match up.
    NPacket* p = packet->getFirstTreeChild();
    int itemCounter = 0;
    PacketTreeItem* item = static_cast<PacketTreeItem*>(child(itemCounter));
    PacketTreeItem* prev = 0;
    PacketTreeItem* other;
    for ( ; p; ++itemCounter, p = p->getNextTreeSibling()) {
        // INV: itemCounter is the current index of p and item.
        if (! item) {
            // We've already run out of child nodes.  Add a new one.
            if (prev)
                prev = new PacketTreeItem(this, prev, p);
            else
                prev = new PacketTreeItem(this, 0, p);
            prev->fill();

            // Variables prev and item are already correct.
        } else if (item->getPacket() == p) {
            // They match up.
            item->refreshSubtree();

            // Update our variables.
            prev = item;
            item = static_cast<PacketTreeItem*>(child(itemCounter + 1));
        } else {
            int otherCounter;
            // They both exist but they don't match up.  Hmmm.
            // Do we have a node for this packet later in the tree?
            for (otherCounter = itemCounter + 1; otherCounter < childCount();
                    ++otherCounter) {
                other = static_cast<PacketTreeItem*>(child(otherCounter));
                if (other->getPacket() == p) {
                    // We've found a node for this packet.
                    // Move it to the correct place.
                    insertChild(itemCounter, takeChild(otherCounter));
                    other->refreshSubtree();

                    // Update our variables.
                    // Note that item is already correct.
                    prev = other;
                    break;
                }
            }

            if (otherCounter == childCount() ) {
                // We couldn't find a node for this packet anywhere.
                // Insert a new one.
                if (prev)
                    prev = new PacketTreeItem(this, prev, p);
                else
                    prev = new PacketTreeItem(this, 0, p);
                prev->fill();

                // Variables prev and item are already correct.
            }
        }
    }

    // Were there any child nodes left over?
    // Note that childCount() will decrease as we delete children here.
    while (itemCounter < childCount())
        delete child(itemCounter);
}

void PacketTreeItem::refreshLabel() {
    if (packet) {
        QString newLabel = packet->getPacketLabel().c_str();
        if (packet->hasTags())
            newLabel += " (+)";
        if (text(0) != newLabel)
            setText(0, newLabel);
    } else
        setText(0, i18n("<Deleted>"));
}

void PacketTreeItem::updateEditable() {
    if (packet && packet->isPacketEditable() != isEditable) {
        // We need updating.
        isEditable = ! isEditable;
        setIcon(0, PacketManager::iconSmall(packet, true));
    }
}

void PacketTreeItem::packetWasChanged(regina::NPacket*) {
    getPart()->setModified(true);
}

void PacketTreeItem::packetWasRenamed(regina::NPacket*) {
    refreshLabel();
    getPart()->setModified(true);
}

void PacketTreeItem::packetToBeDestroyed(regina::NPacket*) {
    packet = 0;
    refreshLabel();
    getPart()->setModified(true);

    // I'm a bit worried about this line, but I understand it will
    // behave correctly. :/
    delete this;
}

void PacketTreeItem::childWasAdded(regina::NPacket*, regina::NPacket*) {
    // Be careful.  We might not be in the GUI thread.
    QApplication::postEvent(tree, new PacketTreeItemEvent(
        static_cast<QEvent::Type>(EVT_TREE_CHILD_ADDED), this));
}

void PacketTreeItem::childWasRemoved(regina::NPacket*, regina::NPacket*,
        bool inParentDestructor) {
    // If we're in the parent destructor, it's all going to be done in
    // this->packetToBeDestroyed() anyway.
    if (! inParentDestructor) {
        refreshSubtree();
        updateEditable();
        getPart()->setModified(true);
    }
}

void PacketTreeItem::childrenWereReordered(regina::NPacket*) {
    refreshSubtree();
    getPart()->setModified(true);
}

PacketTreeView::PacketTreeView(ReginaPart* newPart, QWidget* parent) 
          : QTreeWidget(parent), part(newPart) {
    setRootIsDecorated(true);
    header()->hide();
    setAlternatingRowColors(false);
    setSelectionMode(QAbstractItemView::SingleSelection);

    // Currently we use the platform default activation method (which is
    // often double-click).  To make this single-click always, change
    // itemActivated() to itemClicked().
    connect(this, SIGNAL(itemActivated(QTreeWidgetItem*, int)), this,
        SLOT(packetView(QTreeWidgetItem*)));
}

void PacketTreeView::fill(NPacket* topPacket) {
    clear();
    (new PacketTreeItem(this, topPacket))->fill();
}

PacketTreeItem* PacketTreeView::find(regina::NPacket* packet) {
    if (! packet)
        return 0;

    // Start at the root of the tree and work down.
    // Note that the invisible root item might not be a PacketTreeItem,
    // and we should not try to cast it as such.
    QTreeWidgetItem* root = invisibleRootItem();

    int itemCount = 0;
    PacketTreeItem* item;
    regina::NPacket* current;
    while (itemCount < root->childCount()) {
        item = dynamic_cast<PacketTreeItem*>(root->child(itemCount++));
        current = item->getPacket();

        if (current == packet)
            return item;
        if (current && current->isGrandparentOf(packet)) {
            root = item;
            itemCount = 0;
        }
    }

    return 0;
}

void PacketTreeView::packetView(QTreeWidgetItem* packet) {
    if (packet)
        part->packetView(dynamic_cast<PacketTreeItem*>(packet)->getPacket());
}

void PacketTreeView::refresh(NPacket* topPacket) {
    if (invisibleRootItem()->childCount() != 1)
        fill(topPacket);
    else if (((PacketTreeItem*)invisibleRootItem()->child(0))->
            getPacket() != topPacket)
        fill(topPacket);
    else
        ((PacketTreeItem*)invisibleRootItem()->child(0))->refreshSubtree();
}

void PacketTreeView::customEvent(QEvent* evt) {
    switch (evt->type()) {
        case EVT_TREE_CHILD_ADDED:
            {
                PacketTreeItem* item = static_cast<PacketTreeItemEvent*>(evt)->
                    getItem();

                item->refreshSubtree();
                item->updateEditable();
                part->setModified(true);
            }
            break;
        default:
            break;
    }
}

