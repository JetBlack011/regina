//
//  pairsig.cpp
//
//  Created by John Teague on 07/26/2026.
//

#include "pairsig.h"

#include <algorithm>
#include <cassert>
#include <sstream>
#include <utility>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <utilities/exception.h>
#include <utilities/sigutils.h>

namespace {

// A delimiter separating a plain isoSig from the marked-face index list
// appended after it. Base64Encoder::spare[0] ('_') is documented to never
// occur amongst the base64 characters isoSigDetail()'s default encoding
// uses (`a..zA..Z0..9+-`, see utilities/sigutils.h), and is distinct from
// spare[1] ('.'), which Regina's own isoSig implementation already uses
// internally to guard an optional simplex/facet-lock suffix -- so this
// delimiter can never collide with anything isoSigDetail() itself produces.
constexpr char delimiter = regina::Base64Encoder::spare[0];

// The image of every marked face under psi0 alone (see pairSig()), computed
// once, independent of any automorphism of `canon`.
template <int dim, int subdim>
std::vector<FaceDescriptor<dim>> mapFacesThroughPsi0(
        const regina::Triangulation<dim> &ambient,
        const regina::Isomorphism<dim> &psi0,
        const std::vector<int> &markedFaces) {
    std::vector<FaceDescriptor<dim>> result;
    result.reserve(markedFaces.size());
    for (int f : markedFaces)
        result.push_back(
            applyIsomorphism(psi0, faceDescriptor<dim, subdim>(ambient, f)));
    return result;
}

// The sorted list of image face indices (in `canon`'s own numbering) that
// `underPsi0` maps to under the composed isomorphism alpha ∘ psi0.
template <int dim, int subdim>
std::vector<size_t> imageUnderAlpha(
        const regina::Triangulation<dim> &canon,
        const regina::Isomorphism<dim> &alpha,
        const std::vector<FaceDescriptor<dim>> &underPsi0) {
    std::vector<size_t> image;
    image.reserve(underPsi0.size());
    for (const auto &fu : underPsi0)
        image.push_back(resolveFaceIndex<dim, subdim>(
            canon, applyIsomorphism(alpha, fu)));
    std::ranges::sort(image);
    return image;
}

// Splits a pair signature into its isoSig prefix and marked-face index
// list, then reconstructs the ambient triangulation. Shared by
// fromPairSig() and fromKnottedSurfaceSig().
//
// The marked-face suffix is a sequence of fixed-width base64 fields (see
// pairSig()): each entry uses w = Base64Encoder::integerWidth(M - 1)
// characters, where M is the ambient's own subdim-face count -- the same
// scheme IsoSigPrintable uses for isoSig()'s own gluing data, and why no
// separators or explicit entry count are needed: w is recomputed here from
// the already-reconstructed `ambient`, and the entry count falls out of
// the suffix's total length, since this suffix is always the final
// component of the string.
template <int dim, int subdim>
std::pair<std::unique_ptr<regina::Triangulation<dim>>, std::vector<int>>
decodePairSigParts(const std::string &sigStr) {
    auto pos = sigStr.find(delimiter);
    if (pos == std::string::npos)
        throw regina::InvalidArgument(
            "fromPairSig(): missing delimiter in signature string");

    std::string sig = sigStr.substr(0, pos);
    std::string suffix = sigStr.substr(pos + 1);

    auto ambient = std::make_unique<regina::Triangulation<dim>>(
        regina::Triangulation<dim>::fromSig(sig));

    std::vector<int> markedFaces;
    if (!suffix.empty()) {
        size_t M = ambient->template countFaces<subdim>();
        int w = regina::Base64Encoder::integerWidth(M == 0 ? 0 : M - 1);

        if (suffix.size() % w != 0)
            throw regina::InvalidArgument(
                "fromPairSig(): malformed marked-face index list "
                "(suffix length is not a multiple of the per-index width)");

        try {
            regina::Base64Decoder decoder(suffix.begin(), suffix.end());
            size_t count = suffix.size() / w;
            markedFaces.reserve(count);
            for (size_t i = 0; i < count; ++i)
                markedFaces.push_back(decoder.template decodeInt<int>(w));
        } catch (const regina::InvalidInput &) {
            throw regina::InvalidArgument(
                "fromPairSig(): malformed marked-face index list "
                "(invalid base64 character)");
        }
    }

    return {std::move(ambient), std::move(markedFaces)};
}

} // namespace

template <int dim, int subdim>
typename PairSigContext<dim, subdim>::Detail
PairSigContext<dim, subdim>::detailFor(
        const regina::Triangulation<dim> &ambient) {
    if (ambient.isEmpty() || !ambient.isConnected())
        throw regina::FailedPrecondition(
            "pairSig(): ambient must be non-empty and connected");
    return ambient.isoSigDetail();
}

template <int dim, int subdim>
PairSigContext<dim, subdim>::PairSigContext(
        const regina::Triangulation<dim> &ambient, Detail detail)
    : ambient_(&ambient),
      sig_(std::move(detail.first)),
      psi0_(std::move(detail.second)),
      canon_(regina::Triangulation<dim>::fromSig(sig_)),
      numFaces_(canon_.template countFaces<subdim>()),
      width_(regina::Base64Encoder::integerWidth(
          numFaces_ == 0 ? 0 : numFaces_ - 1)) {
    // THREAD SAFETY, and the reason this is spelled out rather than left to
    // fall out of the initialisers above.
    //
    // sig() is called concurrently by every drain thread, and reads both
    // *ambient_ and canon_ through skeletal queries (face<subdim>(),
    // simplex()->face<subdim>()). Regina computes a triangulation's skeleton
    // LAZILY, on first skeletal query, mutating the object -- so if either
    // triangulation reached those threads with an uncomputed skeleton, the
    // first concurrent queries would race inside Regina.
    //
    // Both happen to be forced already -- detailFor() calls isConnected() on
    // the ambient, and the numFaces_ initialiser calls countFaces<subdim>()
    // on canon_ -- but only incidentally, and an assert would be no help:
    // this is built -O3 -DNDEBUG, so any assert here is compiled out of
    // exactly the binary that runs the concurrent drain.
    //
    // So ESTABLISH the invariant rather than checking it. These two calls are
    // near-free once the skeleton exists (and are what computes it if some
    // future edit removes the incidental triggers above), they run in every
    // build configuration, and they execute here while still single-threaded,
    // before call_once publishes this object to the drain threads. Regina
    // computes the skeleton as a whole, so one query per triangulation is
    // enough. ensureSkeleton() would say this more directly but is protected.
    static_cast<void>(ambient_->isConnected());
    static_cast<void>(canon_.template countFaces<subdim>());

    // Every automorphism of `canon_`, collected once. sig() minimises the
    // marked-face image over these, and which automorphism wins depends on
    // the surface -- but the SET of them does not, so enumerating them here
    // is the single largest saving after isoSigDetail() itself.
    canon_.findAllIsomorphisms(
        canon_, [this](const regina::Isomorphism<dim> &alpha) {
            autos_.push_back(alpha);
            return false; // keep enumerating every automorphism
        });
}

template <int dim, int subdim>
std::string PairSigContext<dim, subdim>::sig(
        const std::vector<int> &markedFaces) const {
    std::ostringstream out;
    out << sig_ << delimiter;

    if (markedFaces.empty())
        return out.str();

    auto underPsi0 =
        mapFacesThroughPsi0<dim, subdim>(*ambient_, psi0_, markedFaces);

    // Minimize the sorted image-index list over every automorphism of
    // `canon_`, so that two isomorphic (ambient, markedFaces) pairs always
    // settle on the same encoding, regardless of which arbitrary relabeling
    // isoSigDetail() happened to return as psi0.
    std::vector<size_t> best;
    bool haveBest = false;
    for (const auto &alpha : autos_) {
        auto candidate = imageUnderAlpha<dim, subdim>(canon_, alpha, underPsi0);
        if (!haveBest || candidate < best) {
            best = std::move(candidate);
            haveBest = true;
        }
    }

    // Encode the minimized index list the same way isoSig() itself encodes
    // its gluing data: fixed-width base64 fields (IsoSigPrintable's own
    // scheme, see utilities/sigutils.h), with no separators -- the decoder
    // recomputes the same width w from the reconstructed ambient
    // triangulation alone, and the entry count from the suffix's length.
    regina::Base64Encoder enc;
    enc.encodeInts(best, width_);
    out << enc.str();
    return out.str();
}

template <int dim, int subdim>
std::string pairSig(const regina::Triangulation<dim> &ambient,
                     const std::vector<int> &markedFaces) {
    return PairSigContext<dim, subdim>(ambient).sig(markedFaces);
}

template <int dim, int subdim>
DecodedPairSig<dim, subdim> fromPairSig(const std::string &sigStr) {
    auto [ambient, markedFaces] = decodePairSigParts<dim, subdim>(sigStr);
    auto skeleton = std::make_unique<Skeleton<dim, subdim>>(*ambient);
    auto submanifold = std::make_unique<EmbeddedSubmanifold<dim, subdim>>(
        *skeleton, markedFaces);
    return DecodedPairSig<dim, subdim>{
        .ambient = std::move(ambient), .skeleton = std::move(skeleton),
        .submanifold = std::move(submanifold)};
}

DecodedKnottedSurfaceSig fromKnottedSurfaceSig(const std::string &sigStr) {
    auto [ambient, markedFaces] = decodePairSigParts<4, 2>(sigStr);
    auto skeleton = std::make_unique<Skeleton<4, 2>>(*ambient);
    auto surface = std::make_unique<KnottedSurface>(*skeleton, markedFaces);
    return DecodedKnottedSurfaceSig{
        .ambient = std::move(ambient), .skeleton = std::move(skeleton),
        .surface = std::move(surface)};
}

template class PairSigContext<3, 2>;
template class PairSigContext<4, 2>;

template std::string pairSig<3, 2>(
    const regina::Triangulation<3> &, const std::vector<int> &);
template std::string pairSig<4, 2>(
    const regina::Triangulation<4> &, const std::vector<int> &);

template DecodedPairSig<3, 2> fromPairSig<3, 2>(const std::string &);
template DecodedPairSig<4, 2> fromPairSig<4, 2>(const std::string &);
