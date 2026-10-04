//
//  isometry.cpp
//

#include "linknaming/isometry/isometry.h"

#include <mutex>
#include <string>

#include <link/link.h>
#include <snappea/snappeatriangulation.h>
#include <snappea/kernel/kernel_prototypes.h>
#include <snappea/kernel/unix_file_io.h>

namespace linknaming {

namespace {

std::mutex &kernelMutex() {
    static std::mutex m;
    return m;
}

// SnapPy retries a non-geometric solution on a randomised triangulation;
// this many tries, as for its Manifold.is_isometric_to().
constexpr int RANDOMISE_TRIES = 16;

} // namespace

bool KernelLink::Meridional::uniform() const {
    for (int s : sign)
        if (s != sign.front())
            return false;
    return true;
}

KernelLink::KernelLink(const regina::Link &link, int randomisations) {
    namespace sp = regina::snappea;
    std::lock_guard<std::mutex> lock(kernelMutex());
    try {
        // Triangulated by the kernel from the diagram, so the peripheral
        // curves are the diagram's meridians and longitudes; the SnapPea
        // file text carries them into a triangulation of our own.
        const std::string file = regina::SnapPeaTriangulation(link).snapPea();
        tri_ = sp::read_triangulation_from_string(file.c_str());
        if (!tri_)
            return;
        for (int r = 0; r < randomisations; ++r)
            sp::randomize_triangulation(tri_);
        sp::SolutionType type = sp::find_complete_hyperbolic_structure(tri_);
        for (int t = 0; type != sp::geometric_solution && t < RANDOMISE_TRIES; ++t) {
            sp::randomize_triangulation(tri_);
            type = sp::find_complete_hyperbolic_structure(tri_);
        }
        hyperbolic_ = (type == sp::geometric_solution);
        if (hyperbolic_)
            volume_ = sp::volume(tri_, nullptr);
    } catch (const std::exception &) {
        hyperbolic_ = false; // regina::SnapPeaFatalError, an unusable diagram, ...
    }
}

KernelLink::~KernelLink() {
    if (tri_) {
        std::lock_guard<std::mutex> lock(kernelMutex());
        regina::snappea::free_triangulation(tri_);
    }
}

std::pair<bool, bool> KernelLink::isometries(const KernelLink &other,
                                             std::vector<Meridional> *meridional) const {
    namespace sp = regina::snappea;
    if (!hyperbolic_ || !other.hyperbolic_ || !tri_ || !other.tri_)
        return {false, false};
    std::lock_guard<std::mutex> lock(kernelMutex());
    sp::Boolean isometric = FALSE;
    sp::IsometryList *all = nullptr, *ofLinks = nullptr;
    bool any = false, extends = false;
    try {
        // Works on copies: neither triangulation is changed.
        if (sp::compute_isometries(tri_, other.tri_, &isometric, &all, &ofLinks) == sp::func_OK) {
            any = isometric;
            const int n = (isometric && ofLinks) ? sp::isometry_list_size(ofLinks) : 0;
            extends = n > 0;
            if (meridional && n > 0) {
                const int cusps = sp::isometry_list_num_cusps(ofLinks);
                for (int i = 0; i < n; ++i) {
                    Meridional m;
                    for (int c = 0; c < cusps; ++c) {
                        int image = 0;
                        int map[2][2];
                        sp::isometry_list_cusp_action(ofLinks, i, c, &image, map);
                        // map[coordinate][basis curve], M = 0, L = 1: the
                        // meridian goes to map[M][M] * meridian (map[L][M] is
                        // 0: the isometry extends to the link), and the sign
                        // of the determinant says whether the isometry keeps
                        // the orientation of the cusp torus, hence of S^3.
                        const int det = map[0][0] * map[1][1] - map[0][1] * map[1][0];
                        if (c == 0)
                            m.reflects = det < 0;
                        m.image.push_back(image);
                        m.sign.push_back((map[0][0] > 0 ? 1 : -1) * (det > 0 ? 1 : -1));
                    }
                    meridional->push_back(std::move(m));
                }
            }
        }
    } catch (const std::exception &) {
        any = extends = false;
        if (meridional)
            meridional->clear();
    }
    if (all)
        sp::free_isometry_list(all);
    if (ofLinks)
        sp::free_isometry_list(ofLinks);
    return {any, extends};
}

bool KernelLink::sameLinkAs(const KernelLink &other) const {
    return isometries(other).second;
}

bool KernelLink::sameComplementAs(const KernelLink &other) const {
    return isometries(other).first;
}

std::vector<KernelLink::Meridional> KernelLink::meridionalIsometriesTo(const KernelLink &other) const {
    std::vector<Meridional> out;
    isometries(other, &out);
    return out;
}

bool KernelLink::sameOrientedLinkAs(const KernelLink &other) const {
    for (const Meridional &m : meridionalIsometriesTo(other))
        if (m.uniform())
            return true;
    return false;
}

} // namespace linknaming
