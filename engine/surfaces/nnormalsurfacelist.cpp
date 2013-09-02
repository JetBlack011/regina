
/**************************************************************************
 *                                                                        *
 *  Regina - A Normal Surface Theory Calculator                           *
 *  Computational Engine                                                  *
 *                                                                        *
 *  Copyright (c) 1999-2013, Ben Burton                                   *
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

#include <iterator>
#include <list>
#include "enumerate/ndoubledescription.h"
#include "enumerate/nhilbertcd.h"
#include "enumerate/nhilbertdual.h"
#ifndef EXCLUDE_NORMALIZ
#include "enumerate/nhilbertprimal.h"
#endif
#include "maths/matrixops.h"
#include "maths/nmatrixint.h"
#include "progress/nprogressmanager.h"
#include "progress/nprogresstypes.h"
#include "surfaces/nnormalsurfacelist.h"
#include "surfaces/flavourregistry.h"
#include "triangulation/ntriangulation.h"
#include "utilities/xmlutils.h"

namespace regina {

namespace {
    // Since legacy coordinate systems don't appear in the flavour
    // registry, give them a consistent name here.
    const char* AN_LEGACY_NAME =
        "Legacy standard almost normal (pruned tri-quad-oct)";
}

const NormalCoords NNormalSurfaceList::STANDARD = NS_STANDARD;;
const NormalCoords NNormalSurfaceList::QUAD = NS_QUAD;
const NormalCoords NNormalSurfaceList::AN_LEGACY = NS_AN_LEGACY;
const NormalCoords NNormalSurfaceList::AN_QUAD_OCT = NS_AN_QUAD_OCT;
const NormalCoords NNormalSurfaceList::AN_STANDARD = NS_AN_STANDARD;
const NormalCoords NNormalSurfaceList::EDGE_WEIGHT = NS_EDGE_WEIGHT;
const NormalCoords NNormalSurfaceList::FACE_ARCS = NS_FACE_ARCS;
const NormalCoords NNormalSurfaceList::ORIENTED = NS_ORIENTED;
const NormalCoords NNormalSurfaceList::ORIENTED_QUAD = NS_ORIENTED_QUAD;

void NNormalSurfaceList::writeAllSurfaces(std::ostream& out) const {
    unsigned long n = getNumberOfSurfaces();
    out << "Number of surfaces is " << n << '\n';
    for (unsigned long i = 0; i < n; i++) {
        getSurface(i)->writeTextShort(out);
        out << '\n';
    }
}

namespace {
    struct ZeroVector : public Returns<NNormalSurfaceVector*> {
        const NTriangulation* tri_;

        ZeroVector(const NTriangulation* tri) : tri_(tri) {}

        template <typename Flavour>
        inline NNormalSurfaceVector* operator() (Flavour f) {
            return Flavour::Vector::makeZeroVector(tri_);
        }
    };
}

NNormalSurfaceVector* makeZeroVector(const NTriangulation* triangulation,
        NormalCoords flavour) {
    return forFlavour(flavour, ZeroVector(triangulation), 0);
}

namespace {
    struct MatchingEquations : public Returns<NMatrixInt*> {
        NTriangulation* tri_;

        MatchingEquations(NTriangulation* tri) : tri_(tri) {}

        template <typename Flavour>
        inline NMatrixInt* operator() (Flavour f) {
            return Flavour::Vector::makeMatchingEquations(tri_);
        }
    };
}

NMatrixInt* makeMatchingEquations(NTriangulation* triangulation,
        NormalCoords flavour) {
    return forFlavour(flavour, MatchingEquations(triangulation), 0);
}

namespace {
    struct EmbeddedConstraints : public Returns<NEnumConstraintList*> {
        NTriangulation* tri_;

        EmbeddedConstraints(NTriangulation* tri) : tri_(tri) {}

        template <typename Flavour>
        inline NEnumConstraintList* operator() (Flavour f) {
            return Flavour::Vector::makeEmbeddedConstraints(tri_);
        }
    };
}

NEnumConstraintList* makeEmbeddedConstraints(NTriangulation* triangulation,
        NormalCoords flavour) {
    return forFlavour(flavour, EmbeddedConstraints(triangulation), 0);
}

namespace {
    struct EnumeratorBase {
        NMatrixInt* eqns_;
        NEnumConstraintList* constraints_;
    };

    // Template on the surface output iterator type, since the type we
    // want to use (NNormalSurfaceList::SurfaceIterator) is protected.
    template <typename OutputIterator>
    struct DDEnumerate : public EnumeratorBase {
        OutputIterator out_;
        NProgressNumber* progress_;

        DDEnumerate(OutputIterator out) : out_(out) {}

        template <typename Flavour>
        inline void operator() (Flavour f) {
            NDoubleDescription::enumerateExtremalRays<typename Flavour::Vector>(
                out_, *eqns_, constraints_, progress_);
        }
    };
}

void* NNormalSurfaceList::VertexEnumerator::run(void*) {
    NProgressNumber* progress = 0;
    if (manager) {
        progress = new NProgressNumber(0, 1);
        manager->setProgress(progress);
    }

    // Choose the most appropriate algorithm for the job.
    if (list->flavour == NS_STANDARD && list->which_.has(NS_EMBEDDED_ONLY) &&
            triang->isValid() && ! triang->isIdeal()) {
        // Enumerate solutions in standard space by going via quad space.
        list->enumerateStandardViaReduced<NormalSpec>(triang, progress);
    } else if (list->flavour == NS_AN_STANDARD &&
            list->which_.has(NS_EMBEDDED_ONLY) && triang->isValid() &&
            ! triang->isIdeal()) {
        // Enumerate solutions in standard almost normal space by going
        // via quad-oct space.
        list->enumerateStandardViaReduced<AlmostNormalSpec>(triang, progress);
    } else {
        // The catch-all double description method.
        DDEnumerate<SurfaceInserter> dd(SurfaceInserter(*list, triang));
        dd.progress_ = progress;

        // Fetch any necessary validity constraints.
        if (list->which_.has(NS_EMBEDDED_ONLY))
            dd.constraints_ = makeEmbeddedConstraints(triang, list->flavour);
        else
            dd.constraints_ = 0;

        // Form the matching equations and starting cone.
        dd.eqns_ = makeMatchingEquations(triang, list->flavour);

        // Find the normal surfaces.
        forFlavour(list->flavour, dd);

        delete dd.eqns_;
        delete dd.constraints_;
    }

    // All done!
    triang->insertChildLast(list);

    if (progress) {
        progress->incCompleted();
        progress->setFinished();
    }

    return 0;
}

#ifndef EXCLUDE_NORMALIZ

namespace {
    // Template on the surface output iterator type, since the type we
    // want to use (NNormalSurfaceList::SurfaceIterator) is protected.
    template <typename OutputIterator>
    struct HPrimalEnumerate : public EnumeratorBase {
        OutputIterator out_;
        NNormalSurfaceList* vtx_;
        NProgressMessage* progress_;

        HPrimalEnumerate(OutputIterator out) : out_(out) {}

        template <typename Flavour>
        inline void operator() (Flavour f) {
            NHilbertPrimal::enumerateHilbertBasis<typename Flavour::Vector>(
                out_, vtx_->beginVectors(), vtx_->endVectors(),
                constraints_, progress_);
        }
    };
}

void* NNormalSurfaceList::FundPrimalEnumerator::run(void*) {
    HPrimalEnumerate<SurfaceInserter> hp(SurfaceInserter(*list, triang));

    if (manager) {
        hp.progress_ = new NProgressMessage("Initialising enumeration");
        manager->setProgress(hp.progress_);
    } else
        hp.progress_ = 0;

    // Fetch validity constraints from the registry.
    if (list->which_.has(NS_EMBEDDED_ONLY))
        hp.constraints_ = makeEmbeddedConstraints(triang, list->flavour);
    else
        hp.constraints_ = 0;

    hp.vtx_ = vtxSurfaces;
    if (! hp.vtx_) {
        // Enumerate all vertex normal surfaces using the default (and
        // hopefully best possible) algorithm.
        if (hp.progress_)
            hp.progress_->setMessage("Enumerating extremal rays");

        hp.vtx_ = new NNormalSurfaceList(list->flavour,
            NS_VERTEX | (list->which_.has(NS_EMBEDDED_ONLY) ?
                NS_EMBEDDED_ONLY : NS_IMMERSED_SINGULAR),
            NS_ALG_DEFAULT);
        VertexEnumerator e(hp.vtx_, triang, 0);
        e.run(0);
    }

    if (hp.progress_)
        hp.progress_->setMessage("Enumerating Hilbert basis");

    // Find all fundamental normal surfaces.
    forFlavour(list->flavour, hp);

    delete hp.constraints_;
    if (! vtxSurfaces)
        delete hp.vtx_;

    // All done!
    triang->insertChildLast(list);

    if (hp.progress_) {
        hp.progress_->setMessage("Finished enumeration");
        hp.progress_->setFinished();
    }

    return 0;
}

#endif // EXCLUDE_NORMALIZ

namespace {
    // Template on the surface output iterator type, since the type we
    // want to use (NNormalSurfaceList::SurfaceIterator) is protected.
    template <typename OutputIterator>
    struct HDualEnumerate : public EnumeratorBase {
        OutputIterator out_;
        NProgressNumber* progress_;

        HDualEnumerate(OutputIterator out) : out_(out) {}

        template <typename Flavour>
        inline void operator() (Flavour f) {
            NHilbertDual::enumerateHilbertBasis<typename Flavour::Vector>(
                out_, *eqns_, constraints_, progress_);
        }
    };
}

void* NNormalSurfaceList::FundDualEnumerator::run(void*) {
    HDualEnumerate<SurfaceInserter> hd(SurfaceInserter(*list, triang));

    if (manager) {
        hd.progress_ = new NProgressNumber(0, 1);
        manager->setProgress(hd.progress_);
    } else
        hd.progress_ = 0;

    // Fetch validity constraints from the registry.
    if (list->which_.has(NS_EMBEDDED_ONLY))
        hd.constraints_ = makeEmbeddedConstraints(triang, list->flavour);
    else
        hd.constraints_ = 0;

    // Form the matching equations and starting cone.
    hd.eqns_ = makeMatchingEquations(triang, list->flavour);

    // Find the normal surfaces.
    forFlavour(list->flavour, hd);

    delete hd.eqns_;
    delete hd.constraints_;

    // All done!
    triang->insertChildLast(list);

    if (hd.progress_) {
        hd.progress_->incCompleted();
        hd.progress_->setFinished();
    }

    return 0;
}

NNormalSurfaceList* NNormalSurfaceList::enumerate(NTriangulation* owner,
        NormalCoords newFlavour, bool embeddedOnly, NProgressManager* manager) {
    NNormalSurfaceList* ans = new NNormalSurfaceList(newFlavour,
        NS_VERTEX | (embeddedOnly ? NS_EMBEDDED_ONLY : NS_IMMERSED_SINGULAR),
        NS_ALG_DEFAULT);
    VertexEnumerator* e = new VertexEnumerator(ans, owner, manager);

    if (manager) {
        if (! e->start(0, true)) {
            delete ans;
            return 0;
        }
        return ans;
    } else {
        e->run(0);
        delete e;
        return ans;
    }
}

#ifndef EXCLUDE_NORMALIZ

NNormalSurfaceList* NNormalSurfaceList::enumerateFundPrimal(
        NTriangulation* owner, NormalCoords newFlavour, bool embeddedOnly,
        NNormalSurfaceList* vtxSurfaces, NProgressManager* manager) {
    NNormalSurfaceList* ans = new NNormalSurfaceList(newFlavour,
        NS_FUNDAMENTAL |
            (embeddedOnly ? NS_EMBEDDED_ONLY : NS_IMMERSED_SINGULAR),
        NS_HILBERT_PRIMAL);
    FundPrimalEnumerator* e = new FundPrimalEnumerator(ans, owner,
        vtxSurfaces, manager);

    if (manager) {
        if (! e->start(0, true)) {
            delete ans;
            return 0;
        }
        return ans;
    } else {
        e->run(0);
        delete e;
        return ans;
    }
}

#endif // EXCLUDE_NORMALIZ

NNormalSurfaceList* NNormalSurfaceList::enumerateFundDual(
        NTriangulation* owner, NormalCoords newFlavour, bool embeddedOnly,
        NProgressManager* manager) {
    NNormalSurfaceList* ans = new NNormalSurfaceList(newFlavour,
        NS_FUNDAMENTAL |
            (embeddedOnly ? NS_EMBEDDED_ONLY : NS_IMMERSED_SINGULAR),
        NS_HILBERT_DUAL);
    FundDualEnumerator* e = new FundDualEnumerator(ans, owner, manager);

    if (manager) {
        if (! e->start(0, true)) {
            delete ans;
            return 0;
        }
        return ans;
    } else {
        e->run(0);
        delete e;
        return ans;
    }
}

NNormalSurfaceList* NNormalSurfaceList::enumerateStandardDirect(
        NTriangulation* owner) {
    NNormalSurfaceList* list = new NNormalSurfaceList(NS_STANDARD,
        NS_VERTEX | NS_EMBEDDED_ONLY, NS_VERTEX_STD_DIRECT);

    // Run a vanilla enumeration in standard coordinates.
    NEnumConstraintList* constraints =
        NNormalSurfaceVectorStandard::makeEmbeddedConstraints(owner);
    NMatrixInt* eqns = makeMatchingEquations(owner, NS_STANDARD);

    NDoubleDescription::enumerateExtremalRays<NNormalSurfaceVectorStandard>(
        SurfaceInserter(*list, owner), *eqns, constraints);

    delete eqns;
    delete constraints;

    // All done!
    owner->insertChildLast(list);
    return list;
}

NNormalSurfaceList* NNormalSurfaceList::enumerateStandardANDirect(
        NTriangulation* owner) {
    NNormalSurfaceList* list = new NNormalSurfaceList(NS_AN_STANDARD,
        NS_VERTEX | NS_EMBEDDED_ONLY, NS_VERTEX_STD_DIRECT);

    // Run a vanilla enumeration in standard almost normal coordinates.
    NEnumConstraintList* constraints =
        NNormalSurfaceVectorANStandard::makeEmbeddedConstraints(owner);
    NMatrixInt* eqns = makeMatchingEquations(owner, NS_AN_STANDARD);

    NDoubleDescription::enumerateExtremalRays<NNormalSurfaceVectorANStandard>(
        SurfaceInserter(*list, owner), *eqns, constraints);

    delete eqns;
    delete constraints;

    // All done!
    owner->insertChildLast(list);
    return list;
}

#ifndef EXCLUDE_NORMALIZ

NNormalSurfaceList* NNormalSurfaceList::enumerateFundFullCone(
        NTriangulation* owner, NormalCoords newFlavour, bool embeddedOnly) {
    NNormalSurfaceList* ans = new NNormalSurfaceList(newFlavour,
        NS_FUNDAMENTAL |
            (embeddedOnly ? NS_EMBEDDED_ONLY : NS_IMMERSED_SINGULAR),
        NS_HILBERT_FULLCONE);
    NMatrixInt* eqns = makeMatchingEquations(owner, newFlavour);

    unsigned rank = rowBasis(*eqns);
    unsigned long dim = eqns->columns();

    std::list<std::vector<mpz_class> > input;
    unsigned r, c;
    for (r = 0; r < rank; ++r) {
        input.push_back(std::vector<mpz_class>());
        std::vector<mpz_class>& v(input.back());
        v.reserve(eqns->columns());
        for (c = 0; c < eqns->columns(); ++c) {
            NLargeInteger& entry(eqns->entry(r, c));
            if (entry.isNative())
                v.push_back(mpz_class(entry.longValue()));
            else
                v.push_back(mpz_class(entry.rawData()));
        }
    }
    delete eqns;

    libnormaliz::Cone<mpz_class> cone(std::list<std::vector<mpz_class> >(),
        input /* equalities */, std::list<std::vector<mpz_class> >());
    libnormaliz::ConeProperties wanted(
        libnormaliz::ConeProperty::HilbertBasis);
    cone.compute(wanted);

    if (! cone.isComputed(libnormaliz::ConeProperty::HilbertBasis)) {
        // TODO: Bail properly.
        delete ans;
        return 0;
    }

    // Fetch validity constraints from the registry.
    NEnumConstraintList* constraints = 0;
    if (embeddedOnly)
        constraints = makeEmbeddedConstraints(owner, newFlavour);

    bool broken;
    int nonZero;
    int i;
    std::list<std::vector<mpz_class> >::const_iterator hlit;
    NEnumConstraintList::const_iterator eit;
    std::set<unsigned>::const_iterator sit;
    NNormalSurfaceVector* v;
    NLargeInteger tmpInt;
    NewNormalSurfaceVector newVec(dim);
    for (hlit = cone.getHilbertBasis().begin();
            hlit != cone.getHilbertBasis().end(); ++hlit) {
        broken = false;
        if (constraints) {
            for (eit = constraints->begin(); eit != constraints->end(); ++eit) {
                nonZero = 0;
                for (sit = eit->begin(); sit != eit->end(); ++sit) {
                    if ((*hlit)[*sit] != 0) {
                        if (++nonZero > 1)
                            break;
                    }
                }
                if (nonZero > 1) {
                    broken = true;
                    break;
                }
            }
        }
        if (! broken) {
            // Insert a new surface.
            v = forFlavour(newFlavour, newVec, 0);
            if (! v) {
                // Coordinate system not recognised.
                delete ans;
                return 0;
            }
            for (i = 0; i < dim; ++i) {
                // We make two copies of the GMP integer instead of one,
                // since NVector/NRay does not give us direct non-const
                // access to its elements.
                tmpInt.setRaw((*hlit)[i].get_mpz_t());
                v->setElement(i, tmpInt);
            }
            ans->surfaces.push_back(new NNormalSurface(owner, v));
        }
    }

    // All done!
    owner->insertChildLast(ans);
    return ans;
}

#endif // EXCLUDE_NORMALIZ

namespace {
    // Template on the surface output iterator type, since the type we
    // want to use (NNormalSurfaceList::SurfaceIterator) is protected.
    template <typename OutputIterator>
    struct HCDEnumerate : public EnumeratorBase {
        OutputIterator out_;

        HCDEnumerate(OutputIterator out) : out_(out) {}

        template <typename Flavour>
        inline void operator() (Flavour f) {
            NHilbertCD::enumerateHilbertBasis<typename Flavour::Vector>(
                out_, *eqns_, constraints_);
        }
    };
}

NNormalSurfaceList* NNormalSurfaceList::enumerateFundCD(
        NTriangulation* owner, NormalCoords newFlavour, bool embeddedOnly) {
    NNormalSurfaceList* ans = new NNormalSurfaceList(newFlavour,
        NS_FUNDAMENTAL |
            (embeddedOnly ? NS_EMBEDDED_ONLY : NS_IMMERSED_SINGULAR),
        NS_HILBERT_CD);

    HCDEnumerate<SurfaceInserter> hcd(SurfaceInserter(*ans, owner));

    hcd.eqns_ = makeMatchingEquations(owner, newFlavour);
    if (embeddedOnly)
        hcd.constraints_ = makeEmbeddedConstraints(owner, newFlavour);
    else
        hcd.constraints_ = 0;

    forFlavour(newFlavour, hcd);

    delete hcd.constraints_;
    delete hcd.eqns_;

    // All done!
    owner->insertChildLast(ans);
    return ans;
}

NTriangulation* NNormalSurfaceList::getTriangulation() const {
    return dynamic_cast<NTriangulation*>(getTreeParent());
}

namespace {
    struct AlmostNormalFunction : public Returns<bool> {
        template <typename Flavour>
        inline bool operator() (Flavour f) { return f.almostNormal; }
    };
}

bool NNormalSurfaceList::allowsAlmostNormal() const {
    if (flavour == NS_AN_LEGACY)
        return true;
    else
        return forFlavour(flavour, AlmostNormalFunction(), false);
}

namespace {
    struct SpunFunction : public Returns<bool> {
        template <typename Flavour>
        inline bool operator() (Flavour f) { return f.spun; }
    };
}

bool NNormalSurfaceList::allowsSpun() const {
    // Both the default and the NS_AN_LEGACY cases should return false.
    return forFlavour(flavour, SpunFunction(), false);
}

namespace {
    struct OrientedFunction : public Returns<bool> {
        template <typename Flavour>
        inline bool operator() (Flavour f) { return f.oriented; }
    };
}

bool NNormalSurfaceList::allowsOriented() const {
    // Both the default and the NS_AN_LEGACY cases should return false.
    return forFlavour(flavour, OrientedFunction(), false);
}

namespace {
    struct NameFunction : public Returns<const char*> {
        template <typename Flavour>
        inline const char* operator() (Flavour f) { return f.name(); }
    };
}

void NNormalSurfaceList::writeTextShort(std::ostream& out) const {
    out << surfaces.size() << " vertex normal surface";
    if (surfaces.size() != 1)
        out << 's';
    out << " (";
    if (flavour == NS_AN_LEGACY)
        out << AN_LEGACY_NAME;
    else
        out << forFlavour(flavour, NameFunction(), "Unknown");
    out << ')';
}

void NNormalSurfaceList::writeTextLong(std::ostream& out) const {
    if (isEmbeddedOnly())
        out << "Embedded ";
    else
        out << "Embedded, immersed & singular ";
    out << "vertex normal surfaces\n";
    out << "Coordinates: ";
    if (flavour == NS_AN_LEGACY)
        out << AN_LEGACY_NAME << '\n';
    else
        out << forFlavour(flavour, NameFunction(), "Unknown") << '\n';
    writeAllSurfaces(out);
}

void NNormalSurfaceList::writeXMLPacketData(std::ostream& out) const {
    // Write the surface list parameters.
    out << "  <params "
        << "type=\"" << which_.intValue() << "\" "
        << "algorithm=\"" << algorithm_.intValue() << "\" "
        << "flavourid=\"" << flavour << "\"\n";
    out << "\tflavour=\"";
    if (flavour == NS_AN_LEGACY)
        out << regina::xml::xmlEncodeSpecialChars(AN_LEGACY_NAME);
    else
        out << regina::xml::xmlEncodeSpecialChars(forFlavour(
            flavour, NameFunction(), "Unknown"));
    out << "\"/>\n";

    // Write the individual surfaces.
    std::vector<NNormalSurface*>::const_iterator it;
    for (it = surfaces.begin(); it != surfaces.end(); it++)
        (*it)->writeXMLData(out);
}

NNormalSurfaceList* NNormalSurfaceList::filterForLocallyCompatiblePairs()
        const {
    // Sanity check:
    if (! isEmbeddedOnly())
        return 0;

    NNormalSurfaceList* ans = new NNormalSurfaceList(
        flavour, NS_CUSTOM | NS_EMBEDDED_ONLY, NS_ALG_CUSTOM);

    // Find all surfaces that have a compatible partner.
    std::vector<NNormalSurface*>::const_iterator first, second;
    for (first = surfaces.begin(); first != surfaces.end(); ++first) {
        for (second = surfaces.begin(); second != surfaces.end();
                ++second) {
            if (second == first)
                continue;

            if ((*first)->locallyCompatible(**second)) {
                ans->surfaces.push_back((*first)->clone());
                break;
            }
        }
    }

    getTriangulation()->insertChildLast(ans);
    return ans;
}

NNormalSurfaceList* NNormalSurfaceList::filterForDisjointPairs() const {
    // Sanity check:
    if (! isEmbeddedOnly())
        return 0;

    NNormalSurfaceList* ans = new NNormalSurfaceList(
        flavour, NS_CUSTOM | NS_EMBEDDED_ONLY, NS_ALG_CUSTOM);

    // Collect all the surfaces that we might care about.
    // This means non-empty, connected and compact.
    std::vector<NNormalSurface*> interesting;
    for (std::vector<NNormalSurface*>::const_iterator it = surfaces.begin();
            it != surfaces.end(); ++it) {
        if ((*it)->isEmpty())
            continue;
        if (! (*it)->isCompact())
            continue;
        if (! (*it)->isConnected())
            continue;

        interesting.push_back(*it);
    }

    // Find all surfaces that have a disjoint partner.
    std::vector<NNormalSurface*>::iterator first, second;
    for (first = interesting.begin(); first != interesting.end(); ++first) {
        for (second = interesting.begin(); second != interesting.end();
                ++second) {
            if (second == first)
                continue;

            if ((*first)->disjoint(**second)) {
                ans->surfaces.push_back((*first)->clone());
                break;
            }
        }
    }

    getTriangulation()->insertChildLast(ans);
    return ans;
}

NNormalSurfaceList* NNormalSurfaceList::filterForPotentiallyIncompressible()
        const {
    // Sanity check:
    if (! isEmbeddedOnly())
        return 0;

    NNormalSurfaceList* ans = new NNormalSurfaceList(
        flavour, NS_CUSTOM | NS_EMBEDDED_ONLY, NS_ALG_CUSTOM);

    NTriangulation* t;
#ifdef DEBUG
    int which = 0;
#endif
    for (std::vector<NNormalSurface*>::const_iterator it = surfaces.begin();
            it != surfaces.end(); ++it) {
#ifdef DEBUG
        std::cout << "Processing surface " << which++ << "..." << std::endl;
#endif
        if ((*it)->isVertexLinking())
            continue;
        if ((*it)->isThinEdgeLink().first)
            continue;

        // If we have a one-sided surface, don't worry about taking the
        // two-sided double cover.  If the complement of the one-sided
        // surface has a compressing disc, then the complement of the
        // double cover has the same compressing disc, and this surface
        // can happily be tossed away.
        t = (*it)->cutAlong();
        if (! t->hasSimpleCompressingDisc())
            ans->surfaces.push_back((*it)->clone());
        delete t;
    }

    getTriangulation()->insertChildLast(ans);
    return ans;
}

NPacket* NNormalSurfaceList::internalClonePacket(NPacket* /* parent */) const {
    NNormalSurfaceList* ans = new NNormalSurfaceList(
        flavour, which_, algorithm_);
    transform(surfaces.begin(), surfaces.end(), back_inserter(ans->surfaces),
        FuncNewClonePtr<NNormalSurface>());
    return ans;
}

} // namespace regina

