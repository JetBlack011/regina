
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  KDE User Interface                                                    *
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

/*! \file eltmovedialog.h
 *  \brief Provides a dialog through which the user can select and
 *  perform an elementary move.
 */

#ifndef __ELTMOVEDIALOG_H
#define __ELTMOVEDIALOG_H

#include "packet/npacketlistener.h"

#include <QDialog>
#include <vector>

class QAbstractButton;
class QButtonGroup;
class QDialogButtonBox;
class QLabel;
class QRadioButton;
class TetrahedronChooser;
class EdgeIntChooser;
template <int, int> class FaceChooser;

namespace regina {
    template <int> class Triangulation;
    typedef Triangulation<3> NTriangulation;
};

/**
 * A dialog used to select and perform an elementary move on a
 * triangulation.
 */
class EltMoveDialog : public QDialog, public regina::NPacketListener {
    Q_OBJECT

    private:
        /**
         * Internal components:
         */
        QLabel* name;
        QLabel* overview;

        FaceChooser<3, 1>* box32;
        FaceChooser<3, 2>* box23;
        TetrahedronChooser* box14;
        EdgeIntChooser* box44;
        FaceChooser<3, 1>* box20e;
        FaceChooser<3, 0>* box20v;
        EdgeIntChooser* box21;
        FaceChooser<3, 2>* boxOpenBook;
        FaceChooser<3, 1>* boxCloseBook;
        TetrahedronChooser* boxShellBdry;
        FaceChooser<3, 1>* boxCollapseEdge;
        QRadioButton* use32;
        QRadioButton* use23;
        QRadioButton* use14;
        QRadioButton* use44;
        QRadioButton* use20e;
        QRadioButton* use20v;
        QRadioButton* use21;
        QRadioButton* useOpenBook;
        QRadioButton* useCloseBook;
        QRadioButton* useShellBdry;
        QRadioButton* useCollapseEdge;
        QButtonGroup* moveTypes;
        QDialogButtonBox* buttons;

        /**
         * Packet tree structure:
         */
        regina::NTriangulation* tri;

    public:
        /**
         * Constructor and destructor.
         */
        EltMoveDialog(QWidget* parent, regina::NTriangulation* useTri);
        ~EltMoveDialog();

        /**
         * Update the overview information if the triangulation changes.
         */
        void packetWasRenamed(regina::NPacket*);
        void packetWasChanged(regina::NPacket*);
        void packetToBeDestroyed(regina::NPacket*);

    protected slots:
        /**
         * A button at the bottom of the dialog has been clicked.
         */
        virtual void clicked(QAbstractButton*);

        /**
         * One of the move types has been selected.
         */
        void updateApply();

    private:
        template <typename ChooserClass>
        void updateStates(ChooserClass* chooser, QRadioButton* button);
};

#endif
