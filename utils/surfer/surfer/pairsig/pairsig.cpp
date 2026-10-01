//
//  pairsig.cpp
//
//  Created by John Teague on 07/26/2026.
//

#include "surfer/pairsig/pairsig.h"

#include "surfer/report/atomicwrite.h"

#include <algorithm>
#include <cassert>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include <triangulation/dim3.h>
#include <triangulation/dim4.h>
#include <utilities/exception.h>
#include <utilities/sigutils.h>

#include "surfer/pairsig/parallelisosig.h"
#include "surfer/pairsig/sha1.h"

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
        const regina::Triangulation<dim> &ambient, unsigned threads) {
    if (ambient.isEmpty() || !ambient.isConnected())
        throw regina::FailedPrecondition(
            "pairSig(): ambient must be non-empty and connected");
    return parallelIsoSigDetail(ambient, threads);
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

template <int dim, int subdim>
PairSigContext<dim, subdim>::PairSigContext(
        const regina::Triangulation<dim> &ambient, Detail detail,
        std::vector<regina::Isomorphism<dim>> autos)
    : ambient_(&ambient),
      sig_(std::move(detail.first)),
      psi0_(std::move(detail.second)),
      canon_(regina::Triangulation<dim>::fromSig(sig_)),
      autos_(std::move(autos)),
      numFaces_(canon_.template countFaces<subdim>()),
      width_(regina::Base64Encoder::integerWidth(
          numFaces_ == 0 ? 0 : numFaces_ - 1)) {
    // The same invariant as the building constructor: both skeletons exist
    // before any concurrent sig().
    static_cast<void>(ambient_->isConnected());
    static_cast<void>(canon_.template countFaces<subdim>());
}

template <int dim, int subdim>
bool PairSigContext<dim, subdim>::verifies() const {
    if (psi0_.size() != ambient_->size() || autos_.empty())
        return false;
    if (psi0_(*ambient_) != canon_)
        return false;
    for (const auto &alpha : autos_)
        if (alpha.size() != canon_.size() || alpha(canon_) != canon_)
            return false;
    return true;
}

template <int dim, int subdim>
std::string PairSigContext<dim, subdim>::ambientKey(
        const regina::Triangulation<dim> &ambient) {
    std::ostringstream s;
    s << "dim " << dim << " size " << ambient.size();
    for (size_t i = 0; i < ambient.size(); ++i) {
        const auto *simplex = ambient.simplex(i);
        for (int f = 0; f <= dim; ++f) {
            const auto *adj = simplex->adjacentSimplex(f);
            s << ' ';
            if (adj)
                s << adj->index() << ':'
                  << simplex->adjacentGluing(f).SnIndex();
            else
                s << '-';
        }
    }
    return pairsig::sha1Hex(s.str());
}

namespace {

// One isomorphism as "image:SnIndex" per simplex.
template <int dim>
void writeIso(std::ostream &out, const regina::Isomorphism<dim> &iso) {
    out << "iso";
    for (size_t i = 0; i < iso.size(); ++i)
        out << ' ' << iso.simpImage(i) << ':' << iso.facetPerm(i).SnIndex();
    out << '\n';
}

template <int dim>
regina::Isomorphism<dim> readIso(std::istream &in, size_t n) {
    std::string line;
    if (!std::getline(in, line))
        throw std::runtime_error("truncated");
    std::istringstream fields(line);
    std::string tag;
    fields >> tag;
    if (tag != "iso")
        throw std::runtime_error("expected an isomorphism");
    regina::Isomorphism<dim> iso(n);
    for (size_t i = 0; i < n; ++i) {
        long image = 0;
        int index = 0;
        char colon = 0;
        if (!(fields >> image >> colon >> index) || colon != ':' || image < 0 ||
            static_cast<size_t>(image) >= n || index < 0 ||
            index >= static_cast<int>(regina::Perm<dim + 1>::nPerms))
            throw std::runtime_error("malformed isomorphism");
        iso.simpImage(i) = image;
        iso.facetPerm(i) = regina::Perm<dim + 1>::Sn[index];
    }
    return iso;
}

} // namespace

template <int dim, int subdim>
void PairSigContext<dim, subdim>::save_(const std::string &path) const {
    report::atomicWrite(path, [&](std::ostream &out) {
        out << "surfer-pairsig-context 1\n"
            << "dim " << dim << " subdim " << subdim << " key "
            << ambientKey(*ambient_) << " size " << ambient_->size()
            << " autos " << autos_.size() << '\n'
            << "sig " << sig_ << '\n';
        writeIso<dim>(out, psi0_);
        for (const auto &alpha : autos_)
            writeIso<dim>(out, alpha);
        out << "end\n";
    });
}

template <int dim, int subdim>
std::unique_ptr<PairSigContext<dim, subdim>>
PairSigContext<dim, subdim>::cached(const regina::Triangulation<dim> &ambient,
                                    const std::string &cacheDir,
                                    bool *loaded, unsigned threads) {
    if (loaded)
        *loaded = false;
    const std::string key = ambientKey(ambient);
    const std::string path = cacheDir + "/" + key + "." +
                             std::to_string(dim) + "-" +
                             std::to_string(subdim) + ".pairsigctx";
    if (std::ifstream in(path); in) {
        try {
            std::string line, fileKey;
            int fileDim = 0, fileSubdim = 0;
            size_t size = 0, autoCount = 0;
            std::getline(in, line);
            if (line != "surfer-pairsig-context 1")
                throw std::runtime_error("unknown format");
            std::getline(in, line);
            std::istringstream head(line);
            std::string k1, k2, k3, k4, k5;
            head >> k1 >> fileDim >> k2 >> fileSubdim >> k3 >> fileKey >> k4 >>
                size >> k5 >> autoCount;
            if (!head || fileDim != dim || fileSubdim != subdim ||
                fileKey != key || size != ambient.size())
                throw std::runtime_error("not this ambient's");
            std::getline(in, line);
            if (line.rfind("sig ", 0) != 0)
                throw std::runtime_error("expected the signature");
            std::string sig = line.substr(4);
            regina::Isomorphism<dim> psi0 = readIso<dim>(in, size);
            std::vector<regina::Isomorphism<dim>> autos;
            for (size_t a = 0; a < autoCount; ++a)
                autos.push_back(readIso<dim>(in, size));
            if (!std::getline(in, line) || line != "end")
                throw std::runtime_error("truncated");
            std::unique_ptr<PairSigContext> ctx(new PairSigContext(
                ambient, Detail(std::move(sig), std::move(psi0)),
                std::move(autos)));
            if (ctx->verifies()) {
                if (loaded)
                    *loaded = true;
                return ctx;
            }
        } catch (const std::exception &) {
            // Unreadable or not this ambient's: rebuilt and replaced below.
        }
    }
    auto ctx = std::make_unique<PairSigContext>(ambient, threads);
    try {
        std::filesystem::create_directories(cacheDir);
        ctx->save_(path);
    } catch (const std::exception &) {
        // A cache that cannot be written costs speed only.
    }
    return ctx;
}

template class PairSigContext<3, 2>;
template class PairSigContext<4, 2>;

template std::string pairSig<3, 2>(
    const regina::Triangulation<3> &, const std::vector<int> &);
template std::string pairSig<4, 2>(
    const regina::Triangulation<4> &, const std::vector<int> &);

template DecodedPairSig<3, 2> fromPairSig<3, 2>(const std::string &);
template DecodedPairSig<4, 2> fromPairSig<4, 2>(const std::string &);
