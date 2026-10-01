//
//  namecache.cpp
//

#include "surfer/enumeration/namecache.h"

#include <algorithm>
#include <sstream>

#include "surfer/pairsig/pairsig.h"

namespace identify {

BoundarySignatureCache::BoundarySignatureCache(
    const regina::Triangulation<3> &boundary, size_t clearThreshold)
    : boundary_(&boundary), clearThreshold_(clearThreshold) {}

void BoundarySignatureCache::ensureAutomorphismGroup_() {
    std::call_once(groupOnce_, [this] {
        size_t numEdges = boundary_->countEdges();
        boundary_->findAllIsomorphisms(*boundary_,
            [&](const regina::Isomorphism<3> &alpha) {
                std::vector<size_t> image(numEdges);
                for (size_t e = 0; e < numEdges; ++e)
                    image[e] = resolveFaceIndex<3, 1>(
                        *boundary_,
                        applyIsomorphism(alpha,
                            faceDescriptor<3, 1>(*boundary_,
                                static_cast<int>(e))));
                automorphismEdgeImages_.push_back(std::move(image));
                return false; // keep enumerating every automorphism
            });
    });
}

std::string BoundarySignatureCache::canonicalKey_(
        const std::vector<size_t> &edgeIndices) {
    ensureAutomorphismGroup_();

    // Starting from edgeIndices itself (already sorted, per this method's
    // precondition) means correctness never depends on the identity
    // automorphism actually being among automorphismEdgeImages_ -- only on
    // finding it, or something at least as good, being harmless.
    std::vector<size_t> best = edgeIndices;
    std::vector<size_t> candidate;
    for (const auto &perm : automorphismEdgeImages_) {
        candidate.clear();
        candidate.reserve(edgeIndices.size());
        for (size_t e : edgeIndices)
            candidate.push_back(perm[e]);
        std::ranges::sort(candidate);
        if (candidate < best)
            best = candidate;
    }

    std::ostringstream out;
    for (size_t i = 0; i < best.size(); ++i) {
        if (i)
            out << ',';
        out << best[i];
    }
    return out.str();
}

std::string BoundarySignatureCache::identifyCached(
        const std::vector<size_t> &edgeIndices,
        const std::function<std::string()> &compute) {
    std::string key = canonicalKey_(edgeIndices);

    {
        std::lock_guard<std::mutex> lock(cacheMutex_);
        ++stats_.checks;
        auto it = cache_.find(key);
        if (it != cache_.end()) {
            ++stats_.hits;
            return it->second;
        }
    }

    std::string result = compute();

    {
        std::lock_guard<std::mutex> lock(cacheMutex_);
        if (cache_.find(key) == cache_.end() &&
                cache_.size() >= clearThreshold_) {
            cache_.clear();
            ++stats_.cacheResets;
        }
        cache_.emplace(std::move(key), result);
    }
    return result;
}

void BoundarySignatureCache::prime(const std::vector<size_t> &edgeIndices,
                                   const std::string &name) {
    std::string key = canonicalKey_(edgeIndices);
    std::lock_guard<std::mutex> lock(cacheMutex_);
    cache_.insert_or_assign(std::move(key), name);
}

BoundarySignatureCacheStats BoundarySignatureCache::stats() const {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    return stats_;
}

size_t BoundarySignatureCache::size() const {
    std::lock_guard<std::mutex> lock(cacheMutex_);
    return cache_.size();
}

} // namespace identify
