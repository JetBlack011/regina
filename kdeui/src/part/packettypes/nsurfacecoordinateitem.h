
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
 *                                                                        *
 *  Copyright (c) 1999-2003, Ben Burton                                   *
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
 *  Software Foundation, Inc., 59 Temple Place, Suite 330, Boston,        *
 *  MA 02111-1307, USA.                                                   *
 *                                                                        *
 **************************************************************************/

/* end stub */

/*! \file nsurfacecoordinateitem.h
 *  \brief Provides a list view item describing a single normal surface.
 */

#ifndef __NSURFACECOORDINATEITEM_H
#define __NSURFACECOORDINATEITEM_H

#include <klistview.h>

namespace regina {
    class NNormalSurface;
    class NNormalSurfaceList;
};

/**
 * A list view item describing a single normal surface.
 */
class NSurfaceCoordinateItem : public KListViewItem {
    private:
        /**
         * The underlying normal surface.
         * Note that the surface name is a direct reference to the table
         * of local modifications in NSurfaceCoordinateUI.
         */
        const regina::NNormalSurface* surface;
        QString& name;
        unsigned long surfaceIndex;
        regina::NNormalSurfaceList* surfaces;
        int coordSystem;

    public:
        /**
         * Constructor.
         */
        NSurfaceCoordinateItem(QListView* parent,
            regina::NNormalSurfaceList* fromSurfaces,
            unsigned long newSurfaceIndex, QString& newName,
            int useCoordSystem);

        /**
         * Query this list view item.
         */
        const regina::NNormalSurface* getSurface();
        unsigned long getSurfaceIndex();

        /**
         * Query the property columns of the coordinate viewer.
         */
        static unsigned propertyColCount(bool embeddedOnly);
        static QString propertyColName(int whichCol, bool embeddedOnly);
        static QString propertyColDesc(int whichCol, bool embeddedOnly);

        /**
         * QListItem overrides.
         */
        QString text(int column) const;
        void setText(int column, const QString& str);
        int width(const QFontMetrics& fm, const QListView* lv, int c) const;
        void paintCell(QPainter* p, const QColorGroup& cg, int column,
            int width, int align);
};

inline NSurfaceCoordinateItem::NSurfaceCoordinateItem(QListView* parent,
        regina::NNormalSurfaceList* fromSurfaces,
        unsigned long newSurfaceIndex, QString& newName, int useCoordSystem) :
        KListViewItem(parent),
        surface(fromSurfaces->getSurface(newSurfaceIndex)),
        name(newName),
        surfaceIndex(newSurfaceIndex),
        surfaces(fromSurfaces),
        coordSystem(useCoordSystem) {
}

inline const regina::NNormalSurface* NSurfaceCoordinateItem::getSurface() {
    return surface;
}

inline unsigned long NSurfaceCoordinateItem::getSurfaceIndex() {
    return surfaceIndex;
}

#endif
