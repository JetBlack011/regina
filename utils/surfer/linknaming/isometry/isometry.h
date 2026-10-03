//
//  snappeaisometry.h
//
//  Link complements in Regina's copy of the SnapPea kernel, compared as
//  SnapPy's Manifold.is_isometric_to() compares them.
//

/*! \file utils/surfer/linknaming/isometry/isometry.h
 *  \brief Whether two link diagrams are the same link, by an isometry of
 *  their complements that carries meridians to meridians.
 *
 *  The complement is built by the SnapPea kernel straight from the diagram
 *  (regina::SnapPeaTriangulation(const Link&)), so its peripheral curves are
 *  the diagram's own meridians and longitudes. sameLinkAs() is the kernel's
 *  compute_isometries(): it canonises both complements -- the Epstein-Penner
 *  canonical cell decomposition, found numerically -- and then tries every
 *  combinatorial isomorphism of the two canonical triangulations, computing
 *  each one's action on the peripheral curves in integers.
 *
 *  A positive answer is therefore exact: a combinatorial isomorphism between
 *  triangulations of the two complements, taking each meridian to a meridian
 *  (the kernel's isometry_extends_to_link()). Filling every cusp along its
 *  meridian extends it to a homeomorphism of S^3 carrying one link onto the
 *  other, so the links are equal up to mirror image and the orientations of
 *  their components. Floating point only steers the canonisation: if it went
 *  wrong, the two triangulations would fail to match and the answer would be
 *  "not found", never a false match.
 *
 *  Only hyperbolic complements are compared (a complete structure with every
 *  tetrahedron positively oriented, after random retriangulations if need
 *  be); the rest -- torus links, satellites -- are left to diagrams.
 *
 *  The kernel is not documented as thread-safe, so every call into it here is
 *  serialised by one mutex. Each takes milliseconds.
 */

#ifndef SURFER_EXACTNAMING_SNAPPEAISOMETRY_H
#define SURFER_EXACTNAMING_SNAPPEAISOMETRY_H

#include <utility>
#include <vector>

namespace regina {
class Link;
namespace snappea {
struct Triangulation;
} // namespace snappea
} // namespace regina

namespace linknaming {

class KernelLink {
  public:
    /** The complement of a classical diagram, with its meridians. Never
     *  throws: a diagram the kernel cannot handle gives a non-hyperbolic
     *  object, which matches nothing. `randomisations` random
     *  retriangulations are made first: a second try after a miss, since
     *  canonisation occasionally retriangulates at random and a miss is not
     *  an answer. */
    explicit KernelLink(const regina::Link &link, int randomisations = 0);
    ~KernelLink();
    KernelLink(const KernelLink &) = delete;
    KernelLink &operator=(const KernelLink &) = delete;

    /** A complete hyperbolic structure, every tetrahedron positively oriented. */
    bool hyperbolic() const { return hyperbolic_; }
    double volume() const { return volume_; }

    /** Some isometry of the complements carries every meridian to a
     *  meridian: exact when true (see the file comment); false means none
     *  was found. Both must be hyperbolic. */
    bool sameLinkAs(const KernelLink &other) const;

    /** Some isometry of the complements, whatever it does to the meridians.
     *  Never used to name anything (a complement does not determine a
     *  link); validation uses it to find the pairs sameLinkAs() must tell
     *  apart. */
    bool sameComplementAs(const KernelLink &other) const;

    /** An isometry carrying every meridian to a meridian, read with the
     *  orientations of the two diagrams. It extends to a homeomorphism of
     *  S^3 carrying this link onto the other; that homeomorphism reverses
     *  the orientation of S^3 when `reflects`, and it carries component i of
     *  this link onto the other's component `image[i]`, with its orientation
     *  kept (`sign[i]` = +1) or reversed (-1). The sign is that of the
     *  meridian's image, times -1 when the isometry reflects: the meridian
     *  is fixed by a component's orientation and that of S^3. */
    struct Meridional {
        bool reflects = false;
        std::vector<int> image;
        std::vector<int> sign;
        /** Every component kept, or every one reversed: then this diagram's
         *  oriented link is the other's, up to mirror and global reversal. */
        bool uniform() const;
    };
    /** Every isometry to `other` carrying meridians to meridians (empty when
     *  there is none, or when either is not hyperbolic). Exact, like
     *  sameLinkAs(). */
    std::vector<Meridional> meridionalIsometriesTo(const KernelLink &other) const;
    /** Some isometry carrying meridians to meridians with a uniform sign:
     *  the same ORIENTED link up to mirror and global reversal. */
    bool sameOrientedLinkAs(const KernelLink &other) const;

  private:
    /** compute_isometries(): whether any isometry exists, and every one that
     *  carries meridians to meridians (`meridional` may be null). */
    std::pair<bool, bool> isometries(const KernelLink &other,
                                     std::vector<Meridional> *meridional = nullptr) const;

    regina::snappea::Triangulation *tri_ = nullptr;
    bool hyperbolic_ = false;
    double volume_ = 0;
};

} // namespace linknaming

#endif
