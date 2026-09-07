#include "surfacesearch.h"

#include "identifycomplement.h"

namespace {

/** SurfaceFoundInfo's `tubedGenus`/`closedComponents` pair; see below. */
struct TubedFields {
    int genus;
    int closedComponents;
};

/**
 * Computes SurfaceFoundInfo::tubedGenus/closedComponents from a surface
 * whose (orientable, genus, punctures) classification is already known.
 */
TubedFields tubedFieldsFor(const regina::Triangulation<2> &surface, int genus,
                           int punctures) {
    if (surface.isConnected())
        return punctures > 0 ? TubedFields{genus, 0} : TubedFields{0, 1};
    auto tubed = KnottedSurface::tubedSurfaceType(surface);
    return {tubed.genus, tubed.closedComponents};
}

} // namespace

SurfaceSearch::SurfaceSearch(
    const regina::Triangulation<4> &tri, const std::vector<int> &seedFaces,
    std::optional<size_t> protectedBoundaryComponent)
    : EmbeddingSearch<4, 2>(tri, seedFaces, protectedBoundaryComponent) {
    KnottedSurface(skeleton_, petalCache_, seedFaces);
}

void SurfaceSearch::configureLimits(const SurfaceSearchLimits &limits) {
    limits_ = limits;
    pendingSurfaces_.setCap(limits.pendingSurfaceCap);
    petalCache_.setClearThreshold(limits.petalCacheLimit);
    linkTally_.setCap(limits.boundaryTallyCap);
    // boundarySigCaches_ isn't built yet. ensureBoundarySigCaches_() reads
    // limits_.boundarySignatureCacheLimit when constructing each entry.
}

void SurfaceSearch::SurfaceTypeTally::merge(
    std::map<SurfaceTypeKey, long long> &local) {
    if (local.empty())
        return;
    std::lock_guard<std::mutex> lock(mutex_);
    for (const auto &[key, n] : local)
        counts_[key] += n;
    local.clear();
}

std::string SurfaceSearch::SurfaceTypeTally::summary() const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream out;
    out << "Number of surfaces found:\n";
    if (counts_.empty()) {
        out << "  (none yet)\n";
    } else {
        for (const auto &[key, n] : counts_)
            out << "  " << KnottedSurface::formatSurfaceType(key) << " = " << n
                << "\n";
    }
    return out.str();
}

void SurfaceSearch::PendingSurfaceBatch::merge(
    std::vector<std::vector<int>> &local) {
    if (local.empty())
        return;
    std::lock_guard<std::mutex> lock(mutex_);
    for (auto &faceIndices : local)
        pending_.push_back(std::move(faceIndices));
    local.clear();
    peakSize_ = std::max(peakSize_, pending_.size());
}

std::vector<std::vector<int>> SurfaceSearch::PendingSurfaceBatch::drain() {
    std::lock_guard<std::mutex> lock(mutex_);
    std::vector<std::vector<int>> out;
    std::swap(out, pending_);
    return out;
}

std::vector<std::vector<int>>
SurfaceSearch::PendingSurfaceBatch::popSome(size_t maxCount) {
    std::lock_guard<std::mutex> lock(mutex_);
    size_t count = std::min(maxCount, pending_.size());
    std::vector<std::vector<int>> out;
    out.reserve(count);
    for (size_t i = 0; i < count; ++i)
        out.push_back(std::move(pending_[pending_.size() - count + i]));
    pending_.resize(pending_.size() - count);
    return out;
}

size_t SurfaceSearch::PendingSurfaceBatch::size() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return pending_.size();
}

size_t SurfaceSearch::PendingSurfaceBatch::cap() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return cap_;
}

void SurfaceSearch::PendingSurfaceBatch::setCap(size_t cap) {
    std::lock_guard<std::mutex> lock(mutex_);
    cap_ = cap;
}

void SurfaceSearch::PendingSurfaceBatch::recordProducerDrain(
    size_t entriesDrained) {
    std::lock_guard<std::mutex> lock(mutex_);
    ++producerDrainEvents_;
    producerDrainedEntries_ += static_cast<long long>(entriesDrained);
}

PendingSurfaceQueueStats SurfaceSearch::PendingSurfaceBatch::stats() const {
    std::lock_guard<std::mutex> lock(mutex_);
    return PendingSurfaceQueueStats{
        .currentSize = pending_.size(),
        .cap = cap_,
        .peakSize = peakSize_,
        .producerDrainEvents = producerDrainEvents_,
        .producerDrainedEntries = producerDrainedEntries_};
}

void SurfaceSearch::LinkBoundaryTally::record(const std::string &descriptor,
                                              const SurfaceTypeKey &type) {
    std::lock_guard<std::mutex> lock(mutex_);
    auto it = descriptorSurfaceTypes_.find(descriptor);
    if (it == descriptorSurfaceTypes_.end()) {
        if (descriptorSurfaceTypes_.size() >= descriptorCap_) {
            ++refusedCount_;
            return;
        }
        order_.push_back(descriptor);
    }
    ++descriptorSurfaceTypes_[descriptor][type];
}

void SurfaceSearch::LinkBoundaryTally::setCap(size_t cap) {
    std::lock_guard<std::mutex> lock(mutex_);
    descriptorCap_ = cap;
}

std::string
SurfaceSearch::LinkBoundaryTally::summary(std::optional<size_t> maxRecent) const {
    std::lock_guard<std::mutex> lock(mutex_);
    std::ostringstream out;
    out << "Surface boundaries found:\n";
    if (descriptorSurfaceTypes_.empty()) {
        out << "  (none yet)\n";
    } else {
        std::vector<std::string> descriptors;
        if (maxRecent && *maxRecent < order_.size()) {
            descriptors.assign(order_.end() - static_cast<long>(*maxRecent),
                               order_.end());
            out << "  (showing the " << *maxRecent << " most recently "
                << "encountered of " << order_.size() << " distinct "
                << "boundaries)\n";
        } else {
            descriptors = order_;
            std::sort(descriptors.begin(), descriptors.end());
        }

        for (const std::string &descriptor : descriptors) {
            out << "  " << descriptor << " - ";
            bool first = true;
            for (const auto &[type, count] :
                 descriptorSurfaceTypes_.at(descriptor)) {
                if (!first)
                    out << ", ";
                out << KnottedSurface::formatSurfaceType(type) << " (" << count
                    << ")";
                first = false;
            }
            out << "\n";
        }
    }
    if (refusedCount_ > 0)
        out << "  (" << refusedCount_ << " distinct boundaries beyond the "
            << "cap of " << descriptorCap_ << " were not recorded; "
            << "surfaces with an already-seen boundary are still "
            << "tallied)\n";
    return out.str();
}

void SurfaceSearch::ensureBoundarySigCaches_() const {
    std::call_once(boundaryCachesOnce_, [this] {
        const auto &tri = skeleton_.triangulation();
        boundaryComponentTris_.reserve(tri.countBoundaryComponents());
        for (size_t c = 0; c < tri.countBoundaryComponents(); ++c)
            boundaryComponentTris_.push_back(
                tri.boundaryComponent(c)->build());

        boundarySigCaches_.reserve(boundaryComponentTris_.size());
        for (const auto &bc : boundaryComponentTris_)
            boundarySigCaches_.push_back(
                std::make_unique<identify::BoundarySignatureCache>(
                    bc, limits_.boundarySignatureCacheLimit));
    });
}

std::pair<std::string, std::vector<BoundaryComponentNames>>
SurfaceSearch::describeBoundary_(
    const std::vector<std::pair<size_t, Link>> &links) {
    ensureBoundarySigCaches_();

    std::ostringstream out;
    std::vector<BoundaryComponentNames> structured;
    structured.reserve(links.size());
    bool firstComponent = true;
    for (const auto &[component, link] : links) {
        if (!firstComponent)
            out << ", ";
        firstComponent = false;

        identify::BoundarySignatureCache &cache = *boundarySigCaches_[component];

        out << (component + 1) << ": ";
        std::vector<std::string> curveNames;
        curveNames.reserve(link.comps_.size());
        bool firstCurve = true;
        for (const Knot &curve : link.comps_) {
            if (!firstCurve)
                out << ", ";
            firstCurve = false;
            std::string name = cache.identifyCached(
                curve.edgeIndices(),
                [&curve] { return identify::identify(curve); });
            out << name;
            curveNames.push_back(std::move(name));
        }
        std::optional<std::string> linkName;
        if (link.comps_.size() > 1) {
            linkName = cache.identifyCached(
                link.edgeIndices(),
                [&link] { return identify::identify(link); });
            out << " (" << *linkName << ")";
        }
        structured.push_back(BoundaryComponentNames{
            component, std::move(curveNames), std::move(linkName)});
    }
    return {out.str(), std::move(structured)};
}

identify::BoundarySignatureCacheStats
SurfaceSearch::boundarySignatureCacheStats() const {
    ensureBoundarySigCaches_();
    identify::BoundarySignatureCacheStats total;
    for (const auto &cache : boundarySigCaches_) {
        auto s = cache->stats();
        total.checks += s.checks;
        total.hits += s.hits;
        total.cacheResets += s.cacheResets;
    }
    return total;
}

size_t SurfaceSearch::boundarySignatureCacheSize() const {
    ensureBoundarySigCaches_();
    size_t total = 0;
    for (const auto &cache : boundarySigCaches_)
        total += cache->size();
    return total;
}

BoundaryCondition SurfaceSearch::classifyByLinks_(
    const std::vector<std::pair<size_t, Link>> &links) {
    if (links.empty())
        return BoundaryCondition::closed;
    bool connected =
        std::all_of(links.begin(), links.end(), [](const auto &entry) {
            return entry.second.comps_.size() <= 1;
        });
    return connected ? BoundaryCondition::connected : BoundaryCondition::proper;
}

BoundaryCondition
SurfaceSearch::classifyCheaply_(const EmbeddedSubmanifold<4, 2> &embedding) {
    if (embedding.isClosed())
        return BoundaryCondition::closed;
    if (embedding.isProper())
        return BoundaryCondition::proper;
    return BoundaryCondition::all;
}

void SurfaceSearch::ThreadHook::onFound(EmbeddedSubmanifold<4, 2> &embedding,
                                        const std::vector<int> &U,
                                        long long faceCount) {
    auto type = KnottedSurface::surfaceTypeKey(embedding.triangulation());
    ++localTypeCounts_[type];
    if (wantLinks_) {
        std::vector<int> faceIndices;
        for (int v : U)
            for (int f : owner_.graph_.graphToSkel[v - 1])
                faceIndices.push_back(f);
        localPending_.push_back(std::move(faceIndices));
    } else if (callbacks_.onSurfaceFound) {
        // Idempotent, so this costs nothing after the first call.
        embedding.usePairSigContext(&owner_.pairSigCtx_);
        auto [orientable, genus, punctures] = type;
        auto tubed =
            tubedFieldsFor(embedding.triangulation(), genus, punctures);
        callbacks_.onSurfaceFound(SurfaceFoundInfo{
            .orientable = orientable,
            .genus = genus,
            .tubedGenus = tubed.genus,
            .closedComponents = tubed.closedComponents,
            .punctures = punctures,
            .connected = embedding.triangulation().isConnected(),
            .triangleCount = faceCount,
            .mostRestrictive = classifyCheaply_(embedding),
            .capturePairSig =
                owner_.limits_.capturePairSig
                    ? std::function<std::string()>(
                          [&embedding] { return embedding.pairSig(); })
                    : std::function<std::string()>{}});
    }
}

void SurfaceSearch::ThreadHook::onFlush() {
    tally_.merge(localTypeCounts_);
    if (!wantLinks_)
        return;
    owner_.pendingSurfaces_.merge(localPending_);

    if (owner_.pendingSurfaces_.size() <= owner_.pendingSurfaces_.cap())
        return;

    bool expected = false;
    if (!owner_.pauseRequested_.compare_exchange_strong(
            expected, true, std::memory_order_relaxed))
        return;

    if (callbacks_.onQueueDrainPause)
        callbacks_.onQueueDrainPause(owner_.pendingSurfaces_.size(),
                                    owner_.pendingSurfaces_.cap());

    constexpr size_t HELPER_DRAIN_BATCH = 64;
    bool interrupted = false;
    while (true) {
        if (owner_.stopRequested_.load(std::memory_order_relaxed)) {
            interrupted = true;
            break;
        }
        auto batch = owner_.pendingSurfaces_.popSome(HELPER_DRAIN_BATCH);
        if (batch.empty())
            break; // fully drained
        if (!helperEmbedding_) {
            helperEmbedding_.emplace(owner_.skeleton_, owner_.petalCache_);
            helperEmbedding_->usePairSigContext(&owner_.pairSigCtx_);
        }
        for (const auto &faceIndices : batch)
            owner_.processEntry_(*helperEmbedding_, faceIndices, callbacks_);
        owner_.pendingSurfaces_.recordProducerDrain(batch.size());
    }

    owner_.pauseRequested_.store(false, std::memory_order_relaxed);

    if (callbacks_.onQueueDrainResume)
        callbacks_.onQueueDrainResume(interrupted,
                                      owner_.pendingSurfaces_.size());
}

std::thread SurfaceSearch::AuxHooks::spawn(std::atomic<bool> &workersFinished) {
    if (!wantLinks_)
        return {};
    return std::thread([this, &workersFinished]() {
        owner_.backgroundDrainLoop_(workersFinished, callbacks_);
    });
}

void SurfaceSearch::AuxHooks::afterJoin() {
    owner_.processRemainingSurfaceBoundaries(numThreads_, callbacks_);
}

void SurfaceSearch::backgroundDrainLoop_(
    const std::atomic<bool> &workersFinished,
    const SurfaceSearchCallbacks &callbacks) {
    using namespace std::chrono_literals;
    constexpr size_t POP_BATCH = 64;
    KnottedSurface embedding(skeleton_, petalCache_);
    embedding.usePairSigContext(&pairSigCtx_);
    while (!workersFinished.load(std::memory_order_relaxed)) {
        auto items = pendingSurfaces_.popSome(POP_BATCH);
        if (items.empty()) {
            std::this_thread::sleep_for(20ms);
            continue;
        }
        for (size_t i = 0; i < items.size(); ++i) {
            if (workersFinished.load(std::memory_order_relaxed)) {
                std::vector<std::vector<int>> remainder(
                    std::make_move_iterator(items.begin() +
                                            static_cast<ptrdiff_t>(i)),
                    std::make_move_iterator(items.end()));
                pendingSurfaces_.merge(remainder);
                return;
            }
            processEntry_(embedding, items[i], callbacks);
        }
    }
}

void SurfaceSearch::processRemainingSurfaceBoundaries(
    unsigned numThreads, const SurfaceSearchCallbacks &callbacks) {
    processBatchParallel_(pendingSurfaces_.drain(), numThreads, callbacks);
}

void SurfaceSearch::processBatchParallel_(
    std::vector<std::vector<int>> batch, unsigned numThreads,
    const SurfaceSearchCallbacks &callbacks) {
    if (batch.empty())
        return;

    const auto phaseStart = std::chrono::steady_clock::now();
    const size_t total = batch.size();
    const unsigned workerCount =
        static_cast<unsigned>(std::min<size_t>(numThreads, total));

    if (callbacks.onBoundaryProcessingStarted)
        callbacks.onBoundaryProcessingStarted(total, workerCount);

    constexpr size_t CHUNK = 64;
    std::atomic<size_t> nextIndex{0};
    std::atomic<size_t> processedCount{0};
    std::atomic<bool> done{false};

    std::thread reporter([&]() {
        using namespace std::chrono_literals;
        while (!done.load(std::memory_order_relaxed)) {
            std::this_thread::sleep_for(1s);
            if (!callbacks.onBoundaryProcessingProgress)
                continue;
            callbacks.onBoundaryProcessingProgress(
                processedCount.load(std::memory_order_relaxed), total,
                std::chrono::steady_clock::now() - phaseStart);
        }
    });

    auto worker = [&]() {
        KnottedSurface embedding(skeleton_, petalCache_);
        embedding.usePairSigContext(&pairSigCtx_);
        while (true) {
            if (skipRemainingDrain_.load(std::memory_order_relaxed))
                break;
            size_t begin =
                nextIndex.fetch_add(CHUNK, std::memory_order_relaxed);
            if (begin >= total)
                break;
            size_t end = std::min(begin + CHUNK, total);
            processBatchRange_(embedding, batch, begin, end, callbacks,
                               &processedCount);
        }
    };

    std::vector<std::thread> threads;
    threads.reserve(workerCount);
    for (unsigned t = 0; t < workerCount; ++t)
        threads.emplace_back(worker);
    for (auto &th : threads)
        th.join();

    done.store(true, std::memory_order_relaxed);
    reporter.join();
    if (callbacks.onBoundaryProcessingComplete)
        callbacks.onBoundaryProcessingComplete(
            total, std::chrono::steady_clock::now() - phaseStart);
}

void SurfaceSearch::processEntry_(KnottedSurface &embedding,
                                  const std::vector<int> &faceIndices,
                                  const SurfaceSearchCallbacks &callbacks) {
    for (int idx : faceIndices)
        embedding.addFace(idx);

    SurfaceTypeKey type = embedding.surfaceType();
    auto links = embedding.boundaryLinks();
    std::string descriptor;
    std::vector<BoundaryComponentNames> boundaryComponents;
    if (!links.empty()) {
        std::tie(descriptor, boundaryComponents) = describeBoundary_(links);
        linkTally_.record(descriptor, type);
    }

    if (callbacks.onSurfaceBoundaryProcessed) {
        auto [orientable, genus, punctures] = type;
        auto tubed =
            tubedFieldsFor(embedding.triangulation(), genus, punctures);
        callbacks.onSurfaceBoundaryProcessed(SurfaceBoundaryInfo{
            SurfaceFoundInfo{
                .orientable = orientable,
                .genus = genus,
                .tubedGenus = tubed.genus,
                .closedComponents = tubed.closedComponents,
                .punctures = punctures,
                .connected = embedding.triangulation().isConnected(),
                .triangleCount = static_cast<long long>(faceIndices.size()),
                .mostRestrictive = classifyByLinks_(links),
                .capturePairSig =
                    limits_.capturePairSig
                        ? std::function<std::string()>(
                              [&embedding] { return embedding.pairSig(); })
                        : std::function<std::string()>{}},
            descriptor, boundaryComponents,
            [&embedding] { return embedding.orientedBoundaryLinks(); }});
    }

    // Reverse order, mirroring how the DFS itself would back out --
    // resets embedding to empty for the next entry.
    for (auto it = faceIndices.rbegin(); it != faceIndices.rend(); ++it)
        embedding.removeFace(*it);
}

void SurfaceSearch::processBatchRange_(
    KnottedSurface &embedding, const std::vector<std::vector<int>> &batch,
    size_t begin, size_t end, const SurfaceSearchCallbacks &callbacks,
    std::atomic<size_t> *processedCounter) {
    for (size_t i = begin; i < end; ++i) {
        processEntry_(embedding, batch[i], callbacks);
        if (processedCounter)
            processedCounter->fetch_add(1, std::memory_order_relaxed);
    }
}

SearchStats SurfaceSearch::search(unsigned numThreads, BoundaryCondition cond,
                                  const SurfaceSearchCallbacks &callbacks,
                                  unsigned iddfsIterations, long long iddfsStep,
                                  std::optional<long long> iddfsStart,
                                  std::optional<unsigned> finalThreads,
                                  bool orientableOnly,
                                  std::optional<long long> hardFaceCap,
                                  long long rootBudgetStart,
                                  long long rootBudgetGrowth) {
    const bool wantLinks = cond == BoundaryCondition::proper ||
                           cond == BoundaryCondition::connected;

    // Thread safety (and performance): Vertex<4>::buildLink() caches
    // its result as a plain, unsynchronized lazily-constructed pointer
    // (see engine/triangulation/dim4/vertex4.cpp)
    for (auto v : skeleton_.triangulation().vertices())
        v->buildLink().isValid();

    // Same reasoning, for BoundaryComponent<4>
    for (size_t c = 0; c < skeleton_.triangulation().countBoundaryComponents();
         ++c)
        skeleton_.triangulation().boundaryComponent(c)->build();

    AuxHooks auxHooks(*this, numThreads, wantLinks, callbacks);
    return runSearch_<KnottedSurface>(
        numThreads, cond,
        // Must stay a prvalue: KnottedSurface owns a PetalCache holding a
        // std::mutex, so it is not movable, and only guaranteed copy elision
        // makes this compile. The context is attached in ThreadHook::onFound
        // instead, where the embedding is actually used.
        [this] { return KnottedSurface(skeleton_, petalCache_); },
        [this, wantLinks, &callbacks] {
            return std::make_unique<ThreadHook>(*this, surfaceTypeTally_,
                                                wantLinks, callbacks);
        },
        [this, wantLinks, &callbacks](const std::vector<int> &seedFaces) {
            KnottedSurface probe(skeleton_, petalCache_, seedFaces);
            probe.usePairSigContext(&pairSigCtx_);
            SurfaceTypeKey type = probe.surfaceType();
            std::map<SurfaceTypeKey, long long> seedTypeCounts{{type, 1}};
            surfaceTypeTally_.merge(seedTypeCounts);
            auto [orientable, genus, punctures] = type;
            auto triangleCount = static_cast<long long>(seedFaces.size());
            std::function<std::string()> capturePairSig =
                limits_.capturePairSig
                    ? std::function<std::string()>(
                          [&probe] { return probe.pairSig(); })
                    : std::function<std::string()>{};
            if (wantLinks) {
                auto links = probe.boundaryLinks();
                std::string descriptor;
                std::vector<BoundaryComponentNames> boundaryComponents;
                if (!links.empty()) {
                    std::tie(descriptor, boundaryComponents) =
                        describeBoundary_(links);
                    linkTally_.record(descriptor, type);
                }
                auto tubed = tubedFieldsFor(probe.triangulation(), genus,
                                            punctures);
                if (callbacks.onSurfaceBoundaryProcessed)
                    callbacks.onSurfaceBoundaryProcessed(SurfaceBoundaryInfo{
                        SurfaceFoundInfo{.orientable = orientable,
                                        .genus = genus,
                                        .tubedGenus = tubed.genus,
                                        .closedComponents =
                                            tubed.closedComponents,
                                        .punctures = punctures,
                                        .connected =
                                            probe.triangulation().isConnected(),
                                        .triangleCount = triangleCount,
                                        .mostRestrictive =
                                            classifyByLinks_(links),
                                        .capturePairSig = capturePairSig},
                        descriptor, boundaryComponents,
                        [&probe] { return probe.orientedBoundaryLinks(); }});
            } else if (callbacks.onSurfaceFound) {
                auto tubed = tubedFieldsFor(probe.triangulation(), genus,
                                            punctures);
                callbacks.onSurfaceFound(SurfaceFoundInfo{
                    .orientable = orientable,
                    .genus = genus,
                    .tubedGenus = tubed.genus,
                    .closedComponents = tubed.closedComponents,
                    .punctures = punctures,
                    .connected = probe.triangulation().isConnected(),
                    .triangleCount = triangleCount,
                    .mostRestrictive = classifyCheaply_(probe),
                    .capturePairSig = capturePairSig});
            }
        },
        callbacks, auxHooks,
        iddfsIterations, iddfsStep, iddfsStart, finalThreads, orientableOnly,
        hardFaceCap, rootBudgetStart, rootBudgetGrowth);
}

