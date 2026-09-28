//
//  snappeakernel_isometry.cpp
//
//  The SnapPea kernel's isometry routines -- compute_isometries() and the
//  IsometryList accessors, which is what SnapPy's Manifold.is_isometric_to()
//  calls. Regina ships them in engine/snappea/kernel/unused/ but does not
//  build them; they are compiled here, against the kernel Regina does build
//  (its functions are exported from the engine library), so the engine is
//  untouched. This file is compiled with engine/snappea/kernel on its include
//  path (CMakeLists.txt), which the kernel sources expect.
//
//  isometry.c sends closed manifolds to isometry_closed.c, which needs the
//  (also unused) Dirichlet domain code. A link complement always has cusps,
//  so that branch is never taken; the stub below keeps the link complete and
//  fails safely if it ever were.
//

#include "unused/isometry.c"

#include "kernel.h"
#include "kernel_namespace.h"

FuncResult compute_closed_isometry(Triangulation * /*manifold0*/,
                                   Triangulation * /*manifold1*/,
                                   Boolean * /*are_isometric*/) {
    return func_failed;
}

#include "end_namespace.h"
