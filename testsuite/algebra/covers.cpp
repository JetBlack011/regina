
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Test Suite                                                            *
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

#include <sstream>
#include <cppunit/extensions/HelperMacros.h>

#include "algebra/abeliangroup.h"
#include "algebra/grouppresentation.h"
#include "link/examplelink.h"
#include "link/link.h"
#include "snappea/snappeatriangulation.h"
#include "triangulation/example3.h"

#include "testsuite/algebra/testalgebra.h"

using regina::AbelianGroup;
using regina::GroupPresentation;
using regina::Link;
using regina::SnapPeaTriangulation;
using regina::Triangulation;

class CoversTest : public CppUnit::TestFixture {
    CPPUNIT_TEST_SUITE(CoversTest);

    CPPUNIT_TEST(trivial);
    CPPUNIT_TEST(manifolds);
    CPPUNIT_TEST(knots);
    CPPUNIT_TEST(freeAbelian);
    CPPUNIT_TEST(free);
    CPPUNIT_TEST(cyclic);

    CPPUNIT_TEST_SUITE_END();

    public:
        void setUp() {
        }

        void tearDown() {
        }

        template <int degree>
        std::vector<std::string> viaSnapPea(const Triangulation<3>& tri) {
            std::vector<std::string> covers;

            SnapPeaTriangulation s(tri);
            if (s.isNull()) {
                covers.push_back("Null_SnapPea");
                return covers;
            }

            s.enumerateCovers(degree, SnapPeaTriangulation::all_covers,
                    [&](const SnapPeaTriangulation& cover,
                    SnapPeaTriangulation::CoverType type) {
                const AbelianGroup* ab = cover.homologyFilled();
                if (ab)
                    covers.push_back(ab->str());
                else
                    covers.push_back("Null_AbelianGroup");
            });
            std::sort(covers.begin(), covers.end());
            return covers;
        }

        template <int degree>
        std::vector<std::string> viaRegina(const Triangulation<3>& tri) {
            std::vector<std::string> covers;

            tri.fundamentalGroup().enumerateCovers<degree>([&](
                    GroupPresentation& g) {
                covers.push_back(g.abelianisation().str());
            });
            std::sort(covers.begin(), covers.end());
            return covers;
        }

        void dumpCovers(std::ostream& out, const char* source,
                const std::vector<std::string>& covers) {
            out << source << " : ";
            bool first = true;
            for (const auto& c : covers) {
                if (first)
                    first = false;
                else
                    out << " | ";
                out << c;
            }
        }

        template <int maxDegree>
        void compareResults(Triangulation<3>* tri,
                bool deleteAfter = true, const char* name = nullptr) {
            auto invSnapPea = viaSnapPea<maxDegree>(*tri);
            auto invRegina = viaRegina<maxDegree>(*tri);
            // dumpCovers(std::cerr, tri->label().c_str(), invRegina);
            if (invSnapPea != invRegina) {
                std::ostringstream msg;
                msg << "Invariants differ for "
                    << (name ? name : tri->label()) << " : ";
                dumpCovers(msg, "SnapPea", invSnapPea);
                msg << " ; ";
                dumpCovers(msg, "Regina", invRegina);
                CPPUNIT_FAIL(msg.str());
            }

            if constexpr (maxDegree > 2) {
                compareResults<maxDegree - 1>(tri, false, name);
            }

            if (deleteAfter)
                delete tri;
        }

        template <int degree>
        void verifyResults(Triangulation<3>* tri,
                const std::vector<std::string>& expected,
                bool deleteAfter = true, const char* name = nullptr) {
            auto invRegina = viaRegina<degree>(*tri);
            if (invRegina != expected) {
                std::ostringstream msg;
                msg << "Invariants differ for "
                    << (name ? name : tri->label()) << " : ";
                dumpCovers(msg, "Expected", expected);
                msg << " ; ";
                dumpCovers(msg, "Regina", invRegina);
                CPPUNIT_FAIL(msg.str());
            }
            if (deleteAfter)
                delete tri;
        }

        template <int maxDegree>
        void compareResults(Link* link, bool deleteAfter = true,
                const char* name = nullptr) {
            compareResults<maxDegree>(link->complement(), true,
                name ? name : link->label().c_str());
            if (deleteAfter)
                delete link;
        }

        template <int maxDegree>
        void verifyResults(Link* link,
                const std::vector<std::string>& expected,
                bool deleteAfter = true, const char* name = nullptr) {
            verifyResults<maxDegree>(link->complement(), expected,
                true, name ? name : link->label().c_str());
            if (deleteAfter)
                delete link;
        }

        void trivial() {
            // No covers:
            compareResults<6>(regina::Example<3>::sphere());
        }

        void manifolds() {
            // Cover only for degree 5:
            compareResults<6>(regina::Example<3>::poincareHomologySphere());

            // Cover (which is trivial) only for degree 3:
            compareResults<6>(regina::Example<3>::lens(3, 1));

            // Many covers for degree 5:
            compareResults<6>(regina::Example<3>::weeks());

            // Many, many covers for degree 5 (and too slow to put degree 6
            // in the test suite):
            Triangulation<3>* ws = regina::Example<3>::weberSeifert();
            verifyResults<2>(ws, { }, false);
            verifyResults<3>(ws, { }, false);
            verifyResults<4>(ws, { }, false);
            verifyResults<5>(ws, {
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25", "2 Z_5 + 2 Z_25",
                "2 Z_5 + 2 Z_25",
                "4 Z + 2 Z_3",
                "6 Z_5 + Z_25", "6 Z_5 + Z_25",
                "6 Z_5 + Z_25", "6 Z_5 + Z_25",
                "6 Z_5 + Z_25", "6 Z_5 + Z_25",
                "Z_5 + 2 Z_25 + Z_75", "Z_5 + 2 Z_25 + Z_75",
                "Z_5 + 2 Z_25 + Z_75", "Z_5 + 2 Z_25 + Z_75",
                "Z_5 + 2 Z_25 + Z_75", "Z_5 + 2 Z_25 + Z_75" }, false);
            delete ws;
        }

        void knots() {
            compareResults<6>(regina::ExampleLink::trefoilRight());
            compareResults<6>(regina::ExampleLink::conway());

            // Each of the following invariants have been verified with SnapPea
            // and/or GAP.  However, SnapPea is slow to compute them for large
            // indices, so for index >= 4 we hard-code the expected results.
            //
            // Actually, we are commenting out the index 6 tests
            // entirely for link19 and link20, since Regina is also too
            // slow with these to use in the test suite.  Note that
            // almost all the running time is being spent computing
            // abelianisations (i.e., Smith normal forms of matrices).

            Link* link19 = regina::Link::fromKnotSig(
                "tabcadefghdijklmnoipkjplmefqrghbcsonqrsvvvvvvb-VzgZBa");
            link19->setLabel("19-crossing knot");
            compareResults<3>(link19, false);
            verifyResults<4>(link19,
                { "2 Z + Z_2" , "Z + Z_9 + Z_39411" }, false);
            verifyResults<5>(link19,
                { "2 Z", "3 Z", "Z + 2 Z_6691" }, false);
#if 0
            verifyResults<6>(link19,
                { "2 Z + Z_157 + Z_628", "2 Z + Z_2 + Z_12", "2 Z + Z_4379",
                  "2 Z + Z_8758", "Z + 2 Z_2 + Z_314 + Z_1375006",
                  "Z + 2 Z_4 + Z_1819388", "Z + Z_628 + Z_324048" }, false);
#endif
            delete link19;

            Link* link20 = regina::Link::fromKnotSig(
                "uabcdbefgecdhifgjklmnhijopqlkqrsaoprtnmtsRktvvvfFyWJTFl");
            link20->setLabel("20-crossing knot");
            compareResults<3>(link20, false);
            verifyResults<4>(link20,
                { "2 Z + Z_140", "Z + Z_25 + Z_91675" }, false);
            verifyResults<5>(link20,
                { "2 Z", "2 Z + Z_2", "3 Z + 2 Z_2", "3 Z + 2 Z_2",
                  "3 Z + Z_2 + Z_6", "3 Z + Z_2 + Z_6", "3 Z + Z_3",
                  "4 Z + Z_4", "Z + 2 Z_15061", "Z + Z_6 + Z_8638440" }, false);
#if 0
            verifyResults<6>(link20,
                { "2 Z + 2 Z_12", "2 Z + 2 Z_3", "2 Z + 2 Z_3", "2 Z + 2 Z_3",
                  "2 Z + 3 Z_2 + Z_36", "2 Z + 3 Z_2 + Z_36",
                  "2 Z + Z_2 + Z_114", "2 Z + Z_2 + Z_12", "2 Z + Z_2 + Z_12",
                  "2 Z + Z_2 + Z_132", "2 Z + Z_2 + Z_248 + Z_8680",
                  "2 Z + Z_2 + Z_36", "2 Z + Z_2 + Z_4 + Z_8",
                  "2 Z + Z_2 + Z_4 + Z_8", "2 Z + Z_3", "2 Z + Z_3 + Z_36",
                  "2 Z + Z_4", "2 Z + Z_456", "2 Z + Z_6 + Z_48",
                  "2 Z + Z_6 + Z_48", "3 Z", "4 Z", "4 Z + Z_2", "4 Z + Z_2",
                  "4 Z + Z_2", "4 Z + Z_2", "4 Z + Z_2", "5 Z", "5 Z",
                  "Z + 2 Z_2 + Z_248 + Z_40176",
                  "Z + 2 Z_2 + Z_496 + Z_1818832", "Z + Z_56161980" }, false);
#endif
            delete link20;
        }

        void freeAbelian() {
            // Free abelian (and indeed just free) on one generator:
            compareResults<6>(regina::Example<3>::s2xs1());
            compareResults<6>(regina::ExampleLink::unknot());

            // Free abelian on two generators:
            // Build a torus, coned in both directions, to give an ideal
            // triangulation of TxI.
            compareResults<6>(regina::Triangulation<3>::fromIsoSig(
                "eLMkbbdddpuapu"), "T x I");

            // Free abelian on three generators:
            compareResults<6>(regina::Triangulation<3>::fromIsoSig(
                "gvLQQedfedffrwawrhh"), "T x S1");
        }

        template <int index>
        void verifyFree(int rank) {
            // The Nielsen–Schreier formula says that every subgroup of g is
            // free of rank 1 + index * (rank - 1).

            // The number of conjugacy classes (at least for all our tests)
            // can be found in the OEIS (sequence A057004):
            //
            // For ranks 1, 2, 3, ...:
            //     index 2: 1, 3, 7, 15, 31, 63, 127, 255, 511
            //     index 3: 1, 7, 41, 235, 1361, 7987, 47321, 281995
            //     index 4: 1, 26, 604, 14120, 334576, 7987616, 191318464
            //     index 5: 1, 97, 13753, 1712845, 207009649
            //     index 6: 1, 624, 504243, 371515454
            //     index 7: 1, 4163, 24824785
            //     index 8: 1, 34470

            // The formulae in terms of the rank r:
            //     index 2:           2^r - 1
            //     index 3 (A057009): 6^(r-1) + 3^(r-1) - 2^(r-1)
            //     index 4 (A057010): 24^(r-1) + 8^(r-1) - 6^(r-1)
            //     index 5 (A057011): 120^(r-1) - 24^(r-1) - 12^(r-1) +
            //                        6^(r-1) + 5^(r-1) + 4^(r-1) - 2^(r-1)
            //     index 6 (A057012): OEIS has PARI code online

            // For now we just hard-code the expected results in an array,
            // which is indexed as expected[index - 2][rank - 1].
            // An array value of 0 means the result is beyond the end of
            // our hard-coded list of results for that particular index.
            static constexpr size_t maxIndex = 6;
            static constexpr size_t maxRank = 9;
            static constexpr size_t expected[maxIndex - 1][maxRank] = {
                { 1, 3, 7, 15, 31, 63, 127, 255, 511 },
                { 1, 7, 41, 235, 1361, 7987, 47321, 281995, 0 },
                { 1, 26, 604, 14120, 334576, 7987616, 191318464, 0, 0 },
                { 1, 97, 13753, 1712845, 207009649, 0, 0, 0, 0 },
                { 1, 624, 504243, 371515454, 0, 0, 0, 0, 0 }
            };

            if (index < 2 || index > maxIndex) {
                std::ostringstream msg;
                msg << "Index " << index << " is out of range for this test.";
                CPPUNIT_FAIL(msg.str());
                // We cannot continue below, but CPPUNIT_FAIL throws an
                // exception so this is fine.
            }
            if (rank > maxRank ||
                    (rank > 0 && expected[index - 2][rank - 1] == 0)) {
                std::ostringstream msg;
                msg << "Rank " << rank
                    << " is out of range for this test for index "
                    << index << ".";
                CPPUNIT_FAIL(msg.str());
                // We cannot continue below, but CPPUNIT_FAIL throws an
                // exception so this is fine.
            }

            // Now: go compute.

            size_t expectedRank = (rank > 0 ? (1 + index * (rank - 1)) : 0);
            size_t expectedCount = (rank > 0 ? expected[index - 2][rank - 1] :
                0);

            bool badCover = false;
            size_t nFound = 0;
            size_t ans = GroupPresentation(rank).enumerateCovers<index>([&](
                    GroupPresentation& g) {
                g.intelligentSimplify();

                if (g.countRelations() != 0 ||
                        g.countGenerators() != expectedRank)
                    badCover = true;

                ++nFound;
            });

            if (badCover) {
                std::ostringstream msg;
                msg << "Free(" << rank << ") yielded an index " << index
                    << " subgroup that is not free of rank "
                    << expectedRank << ".";
                CPPUNIT_FAIL(msg.str());
            }
            if (ans != nFound) {
                std::ostringstream msg;
                msg << "Free(" << rank << ") yielded inconsistent "
                    "numbers of subgroups for index " << index << ".";
                CPPUNIT_FAIL(msg.str());
            }
            if (ans != expectedCount) {
                std::ostringstream msg;
                msg << "Free(" << rank << ") yielded " << ans
                    << "subgroups for index " << index
                    << " instead of the expected "
                    << expectedCount << ".";
                CPPUNIT_FAIL(msg.str());
            }
        }

        void free() {
            // The upper bounds on the ranks below were chosen according to
            // what would finish quickly enough to be part of the test suite.
            for (int rank = 0; rank <= 9; ++rank)
                verifyFree<2>(rank);
            for (int rank = 0; rank <= 6; ++rank)
                verifyFree<3>(rank);
            for (int rank = 0; rank <= 4; ++rank)
                verifyFree<4>(rank);
            for (int rank = 0; rank <= 3; ++rank)
                verifyFree<5>(rank);
            for (int rank = 0; rank <= 2; ++rank)
                verifyFree<6>(rank);
        }

        template <int index>
        void verifyCyclic(long order) {
            // If index divides order then we should have exactly one
            // result, which must be Z_{order/index}.
            // Otherwise we should have no results at all.

            long expectedOrder = (order % index == 0 ? order / index : 0);
            bool badCover = false;
            size_t nFound = 0;

            GroupPresentation src(1);
            src.addRelation(regina::GroupExpression(0, order));
            size_t ans = src.enumerateCovers<index>([&](GroupPresentation& g) {
                g.intelligentSimplify();

                if (expectedOrder == 1) {
                    if (g.countGenerators() != 0)
                        badCover = true;
                } else {
                    if (g.countGenerators() != 1)
                        badCover = true;
                    else if (g.countRelations() != 1)
                        badCover = true;
                    else {
                        long exp = g.relation(0).terms().front().exponent;
                        if (exp != expectedOrder && exp != -expectedOrder)
                            badCover = true;
                    }
                }

                ++nFound;
            });

            if (order % index == 0) {
                if (badCover) {
                    std::ostringstream msg;
                    msg << "Z_" << order << " yielded an index " << index
                        << " subgroup that was not the expected Z_"
                        << expectedOrder << ".";
                    CPPUNIT_FAIL(msg.str());
                }
                if (ans != nFound) {
                    std::ostringstream msg;
                    msg << "Z_" << order << " yielded inconsistent "
                        "numbers of subgroups for index " << index << ".";
                    CPPUNIT_FAIL(msg.str());
                }
                if (ans != 1) {
                    std::ostringstream msg;
                    msg << "Z_" << order << " yielded " << ans
                        << "subgroups for index " << index
                        << " instead of the expected 1.";
                    CPPUNIT_FAIL(msg.str());
                }
            } else {
                if (ans != 0) {
                    std::ostringstream msg;
                    msg << "Z_" << order << " yielded an index " << index
                        << " subgroup, which should not exist.";
                    CPPUNIT_FAIL(msg.str());
                }
            }
        }

        void cyclic() {
            for (int order = 1; order <= 13; ++order)
                verifyCyclic<2>(order);
            for (int order = 1; order <= 13; ++order)
                verifyCyclic<3>(order);
            for (int order = 1; order <= 13; ++order)
                verifyCyclic<4>(order);
            for (int order = 1; order <= 13; ++order)
                verifyCyclic<5>(order);
            for (int order = 1; order <= 13; ++order)
                verifyCyclic<6>(order);
        }
};

void addCovers(CppUnit::TextUi::TestRunner& runner) {
    runner.addTest(CoversTest::suite());
}
