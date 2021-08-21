// The point of this program is simply to see how large some of our data
// structures are, to help decide when to implement move semantics, and
// when it is reasonable to keep things on the stack vs allocate on the heap.

#include <iostream>

#include "manifold/graphloop.h"
#include "manifold/graphpair.h"
#include "manifold/graphtriple.h"
#include "manifold/handlebody.h"
#include "manifold/lensspace.h"
#include "manifold/sfs.h"
#include "manifold/simplesurfacebundle.h"
#include "manifold/snappeacensusmfd.h"
#include "manifold/torusbundle.h"

#include "subcomplex/augtrisolidtorus.h"
#include "subcomplex/blockedsfs.h"
#include "subcomplex/l31pillow.h"
#include "subcomplex/layeredchain.h"
#include "subcomplex/layeredchainpair.h"
#include "subcomplex/layeredlensspace.h"
#include "subcomplex/layeredloop.h"
#include "subcomplex/layeredsolidtorus.h"
#include "subcomplex/layeredtorusbundle.h"
#include "subcomplex/pillowtwosphere.h"
#include "subcomplex/satregion.h"
#include "subcomplex/snappeacensustri.h"
#include "subcomplex/snappedball.h"
#include "subcomplex/snappedtwosphere.h"
#include "subcomplex/spiralsolidtorus.h"
#include "subcomplex/trisolidtorus.h"
#include "subcomplex/trivialtri.h"

#define SHOW_SIZE(T) std::cout << #T << ": " \
    << sizeof(T) << " byte(s)" << std::endl;

using namespace regina;

int main(int argc, char* argv[]) {
    std::cout << "Manifolds:" << std::endl;
    SHOW_SIZE(GraphLoop);
    SHOW_SIZE(GraphPair);
    SHOW_SIZE(GraphTriple);
    SHOW_SIZE(Handlebody);
    SHOW_SIZE(LensSpace);
    SHOW_SIZE(SFSpace);
    SHOW_SIZE(SimpleSurfaceBundle);
    SHOW_SIZE(SnapPeaCensusManifold);
    SHOW_SIZE(TorusBundle);
    std::cout << std::endl;

    std::cout << "Standard Triangulations:" << std::endl;
    SHOW_SIZE(AugTriSolidTorus);
    SHOW_SIZE(BlockedSFS);
    SHOW_SIZE(L31Pillow);
    SHOW_SIZE(LayeredChain);
    SHOW_SIZE(LayeredChainPair);
    SHOW_SIZE(LayeredLensSpace);
    SHOW_SIZE(LayeredLoop);
    SHOW_SIZE(LayeredSolidTorus);
    SHOW_SIZE(LayeredTorusBundle);
    SHOW_SIZE(PillowTwoSphere);
    SHOW_SIZE(SatRegion);
    SHOW_SIZE(SnapPeaCensusTri);
    SHOW_SIZE(SnappedBall);
    SHOW_SIZE(SnappedTwoSphere);
    SHOW_SIZE(SpiralSolidTorus);
    SHOW_SIZE(TriSolidTorus);
    SHOW_SIZE(TrivialTri);
    std::cout << std::endl;

    return 0;
}