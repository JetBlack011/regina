
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Java User Interface                                                   *
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

package normal.engine.triangulation;

import normal.engine.*;
import java.util.*;

public interface NEdge extends ShareableObject {
    public static final int edgeStart[] = { 0, 0, 0, 1, 1, 2 };
    public static final int edgeEnd[] = { 1, 2, 3, 2, 3, 3 };
    public static final int edgeNumber[][] = {
        { -1, 0, 1, 2 },
        {  0,-1, 3, 4 },
        {  1, 3,-1, 5 },
        {  2, 4, 5,-1 }};

    public NComponent getComponent();
    public NBoundaryComponent getBoundaryComponent();
    public NVertex getVertex(int vertex);

    public boolean isBoundary();
    public boolean isValid();

    public long getNumberOfEmbeddings();
    public NEdgeEmbedding getEmbedding(long index);
}

