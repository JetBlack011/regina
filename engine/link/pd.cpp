
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2021, Ben Burton                                   *
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

#include "link/link.h"
#include <cctype>
#include <climits>
#include <cstdlib>
#include <sstream>

namespace regina {

Link* Link::fromPD(const std::string& s) {
    std::vector<std::array<long, 4>> tuples;

    std::array<long, 4> tuple;
    int index = 0;

    const char* begin = s.c_str();
    char* end;
    while (true) {
        // Find the next integer.
        while (*begin && ! ::isdigit(*begin))
            ++begin;
        if (! *begin)
            break;

        long next = ::strtol(begin, &end, 10);
        if (next <= 0 || next == LONG_MAX /* overflow */)
            return nullptr;

        tuple[index++] = next;
        if (index == 4) {
            tuples.push_back(tuple);
            index = 0;
        }

        begin = end;
    }

    if (index != 0)
        return nullptr;

    return fromPD(tuples.begin(), tuples.end());
}

std::vector<std::array<int, 4>> Link::pdData() const {
    const int n = crossings_.size();

    std::vector<std::array<int, 4>> ans;
    ans.reserve(n);

    // Build a lookup table from StrandRef::id() -> PD strand number:
    int* strand = new int[2 * n];
    int next = 1;
    for (auto start : components_) {
        if (! start)
            continue;

        StrandRef s = start;
        do {
            strand[s.id()] = next++;
            ++s;
        } while (s != start);
    }

    // Now process each crossing in turn:
    for (auto start : components_) {
        if (! start)
            continue;

        StrandRef s = start;
        do {
            const StrandRef next = s.next();
            if (next.strand() == 0) {
                // Strand s enters the next crossing from beneath.
                // Identify the other two strands involved.
                const StrandRef upperOut(next.crossing(), 1);
                const StrandRef upperIn = upperOut.prev();

                if (next.crossing()->sign() > 0) {
                    // lower in, upper out, lower out, upper in
                    ans.push_back({ strand[s.id()], strand[upperOut.id()],
                        strand[next.id()], strand[upperIn.id()] });
                } else {
                    // lower in, upper in, lower out, upper out
                    ans.push_back({ strand[s.id()], strand[upperIn.id()],
                        strand[next.id()], strand[upperOut.id()] });
                }
            }

            s = next;
        } while (s != start);
    }

    delete[] strand;
    return ans;
}

std::string Link::pd() const {
    std::ostringstream out;
    pd(out);
    return out.str();
}

void Link::pd(std::ostream& out) const {
    const int n = crossings_.size();

    out << "PD[";
    bool nonEmpty = false;

    // Build a lookup table from StrandRef::id() -> PD strand number:
    int* strand = new int[2 * n];
    int next = 1;
    for (auto start : components_) {
        if (! start)
            continue;

        StrandRef s = start;
        do {
            strand[s.id()] = next++;
            ++s;
        } while (s != start);
    }

    // Now process each crossing in turn:
    for (auto start : components_) {
        if (! start)
            continue;

        StrandRef s = start;
        do {
            const StrandRef next = s.next();
            if (next.strand() == 0) {
                // Strand s enters the next crossing from beneath.
                // Identify the other two strands involved.
                const StrandRef upperOut(next.crossing(), 1);
                const StrandRef upperIn = upperOut.prev();

                if (nonEmpty)
                    out << ", ";
                else
                    nonEmpty = true;

                out << "X[" << strand[s.id()] << ", ";
                if (next.crossing()->sign() > 0) {
                    // lower in, upper out, lower out, upper in
                    out << strand[upperOut.id()] << ", "
                        << strand[next.id()] << ", "
                        << strand[upperIn.id()];
                } else {
                    // lower in, upper in, lower out, upper out
                    out << strand[upperIn.id()] << ", "
                        << strand[next.id()] << ", "
                        << strand[upperOut.id()];
                }
                out << ']';
            }

            s = next;
        } while (s != start);
    }

    delete[] strand;

    out << ']';
}

} // namespace regina

