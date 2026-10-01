//
//  snappeaisometry_test.cpp
//
//  KernelLink (../snappeaisometry.h): far sides from the master that the
//  diagram searches could not reach, against their table entries, and two
//  kinds of pair that must NOT match -- a HOMFLY coincidence, and three far
//  sides whose complement is L11n353's while the links are not (an isometry
//  exists; none carries meridians to meridians).
//

#include <iostream>
#include <string>

#include <link/link.h>

#include "linknaming/tables.h"
#include "linknaming/isometry/isometry.h"

using exactnaming::KernelLink;
using exactnaming::linkFromTablePD;

static int passed = 0;
static int failed = 0;

static void expect(bool ok, const std::string &what) {
    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "  FAIL: " << what << "\n";
    }
}

namespace {

const char *PD_8_16 = "[[2;7;3;8];[4;10;5;9];[6;1;7;2];[8;14;9;13];[10;15;11;16];[12;6;13;5];[14;3;15;4];[16;11;1;12]]";
const char *PD_10_156 = "[[1;13;2;12];[3;8;4;9];[5;14;6;15];[7;18;8;19];[10;15;11;16];[11;1;12;20];[13;6;14;7];[16;9;17;10];[17;4;18;5];[19;3;20;2]]";
const char *PD_10_151 = "[[2;15;3;16];[4;19;5;20];[6;3;7;4];[8;14;9;13];[10;8;11;7];[11;18;12;19];[14;10;15;9];[16;1;17;2];[17;12;18;13];[20;5;1;6]]";
const char *PD_L8a1 = "PD[X[6; 1; 7; 2]; X[14; 7; 15; 8]; X[4; 15; 1; 16]; X[12; 10; 13; 9]; X[8; 4; 9; 3]; X[10; 5; 11; 6]; X[16; 11; 5; 12]; X[2; 14; 3; 13]]";
const char *PD_L11n281 = "PD[X[6; 1; 7; 2]; X[5; 12; 6; 13]; X[8; 4; 9; 3]; X[2; 14; 3; 13]; X[14; 7; 15; 8]; X[9; 18; 10; 19]; X[22; 17; 11; 18]; X[20; 11; 21; 12]; X[16; 21; 17; 22]; X[4; 15; 1; 16]; X[19; 10; 20; 5]]";
const char *PD_L7a6_0 = "PD[X[8; 1; 9; 2]; X[10; 4; 11; 3]; X[14; 10; 7; 9]; X[12; 6; 13; 5]; X[2; 7; 3; 8]; X[4; 12; 5; 11]; X[6; 14; 1; 13]]";
const char *PD_L7a6_1 = "PD[X[12; 2; 13; 1]; X[10; 3; 11; 4]; X[14; 12; 7; 11]; X[8; 5; 9; 6]; X[2; 14; 3; 13]; X[4; 9; 5; 10]; X[6; 7; 1; 8]]";
const char *PD_L11n353 ="PD[X[6; 1; 7; 2]; X[12; 7; 13; 8]; X[4; 13; 1; 14]; X[5; 16; 6; 17]; X[8; 4; 9; 3]; X[9; 21; 10; 20]; X[19; 11; 20; 10]; X[17; 22; 18; 15]; X[21; 18; 22; 19]; X[15; 14; 16; 5]; X[2; 12; 3; 11]]";

KernelLink fromSig(const std::string &sig) { return KernelLink(regina::Link::fromSig(sig)); }
KernelLink fromTable(const char *pd) { return KernelLink(linkFromTablePD(pd)); }

} // namespace

int main() {
    const KernelLink k816 = fromTable(PD_8_16), k10156 = fromTable(PD_10_156),
                     k10151 = fromTable(PD_10_151), kL8a1 = fromTable(PD_L8a1),
                     kL11n281 = fromTable(PD_L11n281), kL11n353 = fromTable(PD_L11n353);
    expect(k816.hyperbolic() && k10151.hyperbolic() && kL8a1.hyperbolic() &&
               kL11n281.hyperbolic() && kL11n353.hyperbolic(),
           "the table entries are hyperbolic");

    // Found: a 10-crossing drawing of 8_16 and of L8a1 (stuck 2 above minimal),
    // another minimal diagram of 10_151, and one of L11n281.
    expect(fromSig("k-ygSLpCidvGDyc").sameLinkAs(k816), "a 10-crossing drawing of 8_16 is 8_16");
    expect(fromSig("k-TaSUnheJKMIxb").sameLinkAs(kL8a1), "a 10-crossing drawing of L8a1 is L8a1");
    expect(fromSig("k-LbSTpoqLnsCyc").sameLinkAs(k10151), "another minimal diagram of 10_151 is 10_151");
    expect(fromSig("l-pcWQ5Gx2TknlMc").sameLinkAs(kL11n281), "a far side is L11n281");

    // Not found: 8_16 shares 10_156's HOMFLY polynomial.
    expect(!fromSig("k-ygSLpCidvGDyc").sameLinkAs(k10156), "8_16's drawing is not 10_156");
    // Not found: the complement is L11n353's, the link is not.
    for (const char *sig : {"l-ZbWs+peITLElLi", "l-JdWw6HhYSiDBQg", "l-BbWs9peKTLE7Ki"})
        expect(!fromSig(sig).sameLinkAs(kL11n353),
               std::string(sig) + ": L11n353's complement, but not L11n353");

    // A torus knot is not hyperbolic and matches nothing.
    const KernelLink trefoil(linkFromTablePD("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]"));
    expect(!trefoil.hyperbolic() && !trefoil.sameLinkAs(trefoil), "3_1 is left to diagrams");

    // Diagrams no drawing should produce, but which must not crash or match:
    // a virtual knot (the kernel refuses it), a zero-crossing unknot, a split
    // link and a composite knot (neither complement is hyperbolic).
    {
        const std::vector<int> signs{1, 1};
        const std::vector<std::vector<int>> virtualTrefoil{{1, 2, -1, -2}};
        const regina::Link v = regina::Link::fromData(signs.begin(), signs.end(),
                                                      virtualTrefoil.begin(), virtualTrefoil.end());
        expect(!v.isClassical(), "the test diagram is virtual");
        expect(!KernelLink(v).hyperbolic(), "a virtual diagram is not hyperbolic (and does not throw)");
    }
    expect(!KernelLink(regina::Link(1)).hyperbolic(), "a zero-crossing unknot is not hyperbolic");
    {
        regina::Link split = linkFromTablePD(PD_8_16);
        split.insertLink(linkFromTablePD(PD_10_151));
        expect(!KernelLink(split).hyperbolic(), "a split link (8_16 u 10_151) is not hyperbolic");
        regina::Link sum = linkFromTablePD(PD_8_16);
        sum.composeWith(linkFromTablePD(PD_10_151));
        expect(!KernelLink(sum).hyperbolic(), "a composite knot (8_16 # 10_151) is not hyperbolic");
    }
    // Orientations, read from the isometries' action on the meridians.
    {
        const regina::Link l0 = linkFromTablePD(PD_L7a6_0), l1 = linkFromTablePD(PD_L7a6_1);
        const KernelLink k0(l0), k1(l1);
        auto some = [](const std::vector<KernelLink::Meridional> &isos, bool reflects, int sign) {
            for (const auto &m : isos)
                if (m.uniform() && m.reflects == reflects && m.sign.front() == sign)
                    return true;
            return false;
        };
        expect(some(k0.meridionalIsometriesTo(k0), false, 1),
               "L7a6{0} to itself: the identity keeps S^3's and every component's orientation");
        regina::Link rev(l0);
        rev.reverse();
        expect(some(k0.meridionalIsometriesTo(KernelLink(rev)), false, -1),
               "L7a6{0} to its global reverse: every component reversed");
        regina::Link mir(l0);
        mir.changeAll();
        expect(!some(k0.meridionalIsometriesTo(KernelLink(mir)), false, 1) &&
                   !some(k0.meridionalIsometriesTo(KernelLink(mir)), false, -1),
               "L7a6{0} to its mirror: no isometry keeps S^3's orientation (L7a6 is chiral)");
        expect(k0.sameOrientedLinkAs(KernelLink(mir)), "L7a6{0} is its mirror up to mirror");
        // {0} and {1} differ by one component's orientation, and in g4 (1, 0).
        expect(k0.sameLinkAs(k1), "L7a6{0} and L7a6{1} are one unoriented link");
        expect(!k0.sameOrientedLinkAs(k1), "L7a6{0} and L7a6{1} are different oriented links");
        regina::Link flip(l0);
        flip.reverse(flip.component(1));
        expect(KernelLink(flip).sameOrientedLinkAs(k1),
               "L7a6{0} with component 1 reversed is L7a6{1}");
    }

    // A large diagram: 10_151 with forty kinks is still 10_151.
    {
        regina::Link big = linkFromTablePD(PD_10_151);
        for (int k = 0; k < 40; ++k)
            big.r1(big.crossing(0)->upper(), k % 2, k % 3 ? 1 : -1);
        expect(big.size() == 50, "the kinked diagram has 50 crossings");
        expect(KernelLink(big).sameLinkAs(k10151), "10_151 with forty kinks is 10_151");
    }

    std::cout << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
