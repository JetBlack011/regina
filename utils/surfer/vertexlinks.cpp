//
//  vertexlinks.cpp
//
//  Created by John Teague on 07/26/2026.
//

#include "vertexlinks.h"

#include <algorithm>

size_t
PetalCache::CornersHash::operator()(const std::vector<Corner> &corners) const {
    size_t h = corners.size();
    for (const auto &[f, local] : corners) {
        size_t seed =
            (static_cast<size_t>(f) << 2) | static_cast<size_t>(local);
        // boost::hash_combine's mixing formula -- doesn't need to be collision-
        // free (unordered_map's equality fallback handles that correctly
        // regardless), just well-distributed.
        h ^= seed + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2);
    }
    return h;
}

uint64_t PetalCache::linkKey_(int a, int b) {
    auto lo = static_cast<uint32_t>(std::min(a, b));
    auto hi = static_cast<uint32_t>(std::max(a, b));
    return (static_cast<uint64_t>(hi) << 32) | lo;
}

PetalCache::PetalId PetalCache::internPetal(std::vector<Corner> corners) {
    std::ranges::sort(corners);

    std::lock_guard<std::mutex> lock(mutex_);
    auto it = byCorners_.find(corners);
    if (it != byCorners_.end())
        return makeId_(epoch_, it->second);

    if (unknotById_.size() >= clearThreshold_) {
        byCorners_.clear();
        unknotById_.clear();
        linkingCache_.clear();
        petalSetCache_.clear();
        ++epoch_;
        ++stats_.cacheResets;
    }

    int id = static_cast<int>(unknotById_.size());
    unknotById_.emplace_back(std::nullopt);
    byCorners_.emplace(std::move(corners), id);
    return makeId_(epoch_, id);
}

std::optional<bool> PetalCache::lookupUnknot(PetalId id) const {
    std::lock_guard<std::mutex> lock(mutex_);
    ++stats_.unknotChecks;
    if (epochOf_(id) != epoch_)
        return std::nullopt; // stale: a clear happened since id was interned
    const auto &cached = unknotById_[static_cast<size_t>(localIdOf_(id))];
    if (cached)
        ++stats_.unknotCacheHits;
    return cached;
}

void PetalCache::recordUnknot(PetalId id, bool isUnknot) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (epochOf_(id) != epoch_)
        return; // stale: safe to drop -- the petal is simply re-interned fresh
                // next time
    unknotById_[static_cast<size_t>(localIdOf_(id))] = isUnknot;
}

std::optional<bool> PetalCache::lookupLinksNonzero(PetalId a, PetalId b) const {
    std::lock_guard<std::mutex> lock(mutex_);
    ++stats_.linkingChecks;
    if (epochOf_(a) != epoch_ || epochOf_(b) != epoch_)
        return std::nullopt;
    auto it = linkingCache_.find(linkKey_(localIdOf_(a), localIdOf_(b)));
    if (it == linkingCache_.end())
        return std::nullopt;
    ++stats_.linkingCacheHits;
    return it->second;
}

void PetalCache::recordLinksNonzero(PetalId a, PetalId b, bool nonzero) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (epochOf_(a) != epoch_ || epochOf_(b) != epoch_)
        return;
    linkingCache_[linkKey_(localIdOf_(a), localIdOf_(b))] = nonzero;
}

size_t
PetalCache::PetalSetHash::operator()(const std::vector<int> &key) const {
    size_t h = key.size();
    for (int x : key)
        h ^= static_cast<size_t>(x) + 0x9e3779b97f4a7c15ULL + (h << 6) +
             (h >> 2);
    return h;
}

std::optional<std::vector<int>>
PetalCache::petalSetKey_(SetQuery query,
                         const std::vector<PetalId> &ids) const {
    std::vector<int> key;
    key.reserve(ids.size() + 1);
    key.push_back(static_cast<int>(query));
    std::vector<int> local;
    local.reserve(ids.size());
    for (PetalId id : ids) {
        if (epochOf_(id) != epoch_)
            return std::nullopt;
        local.push_back(localIdOf_(id));
    }
    std::ranges::sort(local);
    key.insert(key.end(), local.begin(), local.end());
    return key;
}

std::optional<bool>
PetalCache::lookupPetalSet(SetQuery query,
                           const std::vector<PetalId> &ids) const {
    std::lock_guard<std::mutex> lock(mutex_);
    ++stats_.petalSetChecks;
    auto key = petalSetKey_(query, ids);
    if (!key)
        return std::nullopt;
    auto it = petalSetCache_.find(*key);
    if (it == petalSetCache_.end())
        return std::nullopt;
    ++stats_.petalSetCacheHits;
    return it->second;
}

void PetalCache::recordPetalSet(SetQuery query,
                                const std::vector<PetalId> &ids,
                                bool answer) {
    std::lock_guard<std::mutex> lock(mutex_);
    // Bounded like the rest of the cache: past the same threshold, clear
    // everything and start a new epoch (see internPetal()). `ids` then belong
    // to the old epoch, so petalSetKey_() declines them and nothing is
    // recorded -- a harmless miss next time.
    if (petalSetCache_.size() >= clearThreshold_) {
        byCorners_.clear();
        unknotById_.clear();
        linkingCache_.clear();
        petalSetCache_.clear();
        ++epoch_;
        ++stats_.cacheResets;
    }
    auto key = petalSetKey_(query, ids);
    if (key)
        petalSetCache_[std::move(*key)] = answer;
}

void PetalCache::setClearThreshold(size_t threshold) {
    std::lock_guard<std::mutex> lock(mutex_);
    clearThreshold_ = threshold;
}

size_t PetalCache::size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return byCorners_.size();
}

void PetalCache::recordLocalFlatnessRejection() {
    std::lock_guard<std::mutex> lock(mutex_);
    ++stats_.localFlatnessRejections;
}

void PetalCache::recordTransverseRejection() {
    std::lock_guard<std::mutex> lock(mutex_);
    ++stats_.transverseRejections;
}

PetalCache::Stats PetalCache::stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return stats_;
}
