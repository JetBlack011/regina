//
//  pairsig.cpp
//
//  Created by John Teague on 07/26/2026.
//

#include "pairsig.h"

#include <algorithm>
#include <cctype>
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
template <int dim>
std::pair<std::unique_ptr<regina::Triangulation<dim>>, std::vector<int>>
decodePairSigParts(const std::string &sigStr) {
    auto pos = sigStr.find(delimiter);
    if (pos == std::string::npos)
        throw regina::InvalidArgument(
            "fromPairSig(): missing delimiter in signature string");

    std::string sig = sigStr.substr(0, pos);
    std::string indexList = sigStr.substr(pos + 1);

    auto ambient = std::make_unique<regina::Triangulation<dim>>(
        regina::Triangulation<dim>::fromSig(sig));

    std::vector<int> markedFaces;
    if (!indexList.empty()) {
        std::istringstream in(indexList);
        std::string token;
        while (std::getline(in, token, ',')) {
            if (token.empty() ||
                    !std::all_of(token.begin(), token.end(),
                        [](unsigned char c) { return std::isdigit(c); }))
                throw regina::InvalidArgument(
                    "fromPairSig(): malformed face index in signature "
                    "string");
            markedFaces.push_back(std::stoi(token));
        }
    }

    return {std::move(ambient), std::move(markedFaces)};
}

} // namespace

template <int dim, int subdim>
std::string pairSig(const regina::Triangulation<dim> &ambient,
                     const std::vector<int> &markedFaces) {
    if (ambient.isEmpty() || !ambient.isConnected())
        throw regina::FailedPrecondition(
            "pairSig(): ambient must be non-empty and connected");

    auto [sig, psi0] = ambient.isoSigDetail();

    std::ostringstream out;
    out << sig << delimiter;

    if (markedFaces.empty())
        return out.str();

    regina::Triangulation<dim> canon = regina::Triangulation<dim>::fromSig(sig);
    auto underPsi0 = mapFacesThroughPsi0<dim, subdim>(ambient, psi0, markedFaces);

    // Minimize the sorted image-index list over every automorphism of
    // `canon`, so that two isomorphic (ambient, markedFaces) pairs always
    // settle on the same encoding, regardless of which arbitrary relabeling
    // isoSigDetail() happened to return as psi0.
    std::vector<size_t> best;
    bool haveBest = false;
    canon.findAllIsomorphisms(canon,
        [&](const regina::Isomorphism<dim> &alpha) {
            auto candidate = imageUnderAlpha<dim, subdim>(canon, alpha, underPsi0);
            if (!haveBest || candidate < best) {
                best = std::move(candidate);
                haveBest = true;
            }
            return false; // keep enumerating every automorphism
        });

    for (size_t i = 0; i < best.size(); ++i) {
        if (i)
            out << ',';
        out << best[i];
    }
    return out.str();
}

template <int dim, int subdim>
DecodedPairSig<dim, subdim> fromPairSig(const std::string &sigStr) {
    auto [ambient, markedFaces] = decodePairSigParts<dim>(sigStr);
    auto skeleton = std::make_unique<Skeleton<dim, subdim>>(*ambient);
    auto submanifold = std::make_unique<EmbeddedSubmanifold<dim, subdim>>(
        *skeleton, markedFaces);
    return DecodedPairSig<dim, subdim>{
        .ambient = std::move(ambient), .skeleton = std::move(skeleton),
        .submanifold = std::move(submanifold)};
}

DecodedKnottedSurfaceSig fromKnottedSurfaceSig(const std::string &sigStr) {
    auto [ambient, markedFaces] = decodePairSigParts<4>(sigStr);
    auto skeleton = std::make_unique<Skeleton<4, 2>>(*ambient);
    auto surface = std::make_unique<KnottedSurface>(*skeleton, markedFaces);
    return DecodedKnottedSurfaceSig{
        .ambient = std::move(ambient), .skeleton = std::move(skeleton),
        .surface = std::move(surface)};
}

template std::string pairSig<3, 2>(
    const regina::Triangulation<3> &, const std::vector<int> &);
template std::string pairSig<4, 2>(
    const regina::Triangulation<4> &, const std::vector<int> &);

template DecodedPairSig<3, 2> fromPairSig<3, 2>(const std::string &);
template DecodedPairSig<4, 2> fromPairSig<4, 2>(const std::string &);
