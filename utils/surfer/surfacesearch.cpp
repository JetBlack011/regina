#include "surfacesearch.h"

SurfaceSearch::SurfaceSearch(const regina::Triangulation<4> &tri,
                             const std::vector<int> &seedFaces)
    : EmbeddingSearch<4, 2>(tri, seedFaces) {
    // Additional validation beyond the base class's own (facet-level-only)
    // seeded-constructor check above: KnottedSurface's seeded constructor
    // runs the enhanced addFace()/addFaces() (transverse-self-intersection/
    // local-flatness checks included), throwing regina::InvalidArgument if
    // any face fails them -- this temporary exists purely for that side
    // effect.
    KnottedSurface(skeleton_, petalCache_, seedFaces);
}

void SurfaceSearch::configureLimits(const SurfaceSearchLimits &limits) {
    limits_ = limits;
    pendingSurfaces_.setCap(limits.pendingSurfaceCap);
    petalCache_.setClearThreshold(limits.petalCacheLimit);
    linkTally_.setCap(limits.boundaryTallyCap);
    // boundarySigCaches_ isn't built yet -- ensureBoundarySigCaches_() reads
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
            // Most recently *first-encountered* descriptors, oldest of the
            // shown ones first -- discovery order, not alphabetical, so a
            // live progress report reads as a log rather than reshuffling
            // every time a new descriptor enters/leaves the window.
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
            boundarySigCaches_.push_back(std::make_unique<BoundarySignatureCache>(
                bc, limits_.boundarySignatureCacheLimit));
    });
}

std::string SurfaceSearch::describeBoundary_(
    const std::vector<std::pair<size_t, Link>> &links) {
    ensureBoundarySigCaches_();

    std::ostringstream out;
    bool firstComponent = true;
    for (const auto &[component, link] : links) {
        if (!firstComponent)
            out << ", ";
        firstComponent = false;

        BoundarySignatureCache &cache = *boundarySigCaches_[component];

        out << (component + 1) << ": ";
        bool firstCurve = true;
        for (const Knot &curve : link.comps_) {
            if (!firstCurve)
                out << ", ";
            firstCurve = false;
            out << cache.identifyCached(curve.edgeIndices(),
                                        [&curve] { return curve.identify(); });
        }
        if (link.comps_.size() > 1)
            out << " (" << cache.identifyCached(
                               link.edgeIndices(),
                               [&link] { return link.identify(); })
                << ")";
    }
    return out.str();
}

BoundarySignatureCacheStats SurfaceSearch::boundarySignatureCacheStats() const {
    ensureBoundarySigCaches_();
    BoundarySignatureCacheStats total;
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
        auto [orientable, genus, punctures] = type;
        callbacks_.onSurfaceFound(SurfaceFoundInfo{
            .orientable = orientable,
            .genus = genus,
            .punctures = punctures,
            .triangleCount = faceCount,
            .mostRestrictive = classifyCheaply_(embedding)});
    }
}

void SurfaceSearch::ThreadHook::onFlush() {
    tally_.merge(localTypeCounts_);
    if (!wantLinks_)
        return;
    owner_.pendingSurfaces_.merge(localPending_);

    if (owner_.pendingSurfaces_.size() <= owner_.pendingSurfaces_.cap())
        return;

    // Backpressure: once the shared queue is over its cap, one thread
    // claims exclusive responsibility for a full drain-to-empty pass and
    // PAUSES THE WHOLE SEARCH while it runs -- setting owner_.pauseRequested_
    // makes every other worker's tryAdd() block (see InterruptiblePredicate)
    // rather than reject, so no in-progress DFS branch is lost and no other
    // thread can push more work into the queue while this drain is in
    // progress. That's what makes this a genuine full drain (down to
    // empty), not just "back under cap": nothing new can arrive mid-drain.
    //
    // The compare-exchange means only the first thread to notice the queue
    // over cap actually pauses/drains; any other thread whose onFlush()
    // also notices it around the same time just returns here -- it's about
    // to block in its own tryAdd() anyway, within one DFS step, once
    // pauseRequested_ takes effect.
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
        // Checked before popping another batch (not mid-batch), so nothing
        // popped is ever abandoned unprocessed -- whatever's still merged
        // into pendingSurfaces_ when this breaks stays there for
        // processRemainingSurfaceBoundaries()'s single final pass.
        if (owner_.stopRequested_.load(std::memory_order_relaxed)) {
            interrupted = true;
            break;
        }
        auto batch = owner_.pendingSurfaces_.popSome(HELPER_DRAIN_BATCH);
        if (batch.empty())
            break; // fully drained
        if (!helperEmbedding_)
            helperEmbedding_.emplace(owner_.skeleton_, owner_.petalCache_);
        for (const auto &faceIndices : batch)
            owner_.processEntry_(*helperEmbedding_, faceIndices, callbacks_);
        owner_.pendingSurfaces_.recordProducerDrain(batch.size());
    }

    // Unpause before the resume callback, not after: a caller watching for
    // "search has resumed" shouldn't be able to observe pauseRequested_
    // still set once notified.
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
    // Matches processBatchParallel_'s own CHUNK size -- not load-bearing
    // that they're equal, just a reasonable shared "small batch" constant.
    constexpr size_t POP_BATCH = 64;
    KnottedSurface embedding(skeleton_, petalCache_);
    while (!workersFinished.load(std::memory_order_relaxed)) {
        auto items = pendingSurfaces_.popSome(POP_BATCH);
        if (items.empty()) {
            std::this_thread::sleep_for(20ms);
            continue;
        }
        for (size_t i = 0; i < items.size(); ++i) {
            if (workersFinished.load(std::memory_order_relaxed)) {
                // The search ended partway through this small batch: hand
                // the unprocessed remainder back to the shared queue --
                // popSome() already removed it from pending_, so
                // processRemainingSurfaceBoundaries()'s drain() would
                // never otherwise see it.
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
        while (true) {
            size_t begin =
                nextIndex.fetch_add(CHUNK, std::memory_order_relaxed);
            if (begin >= total)
                break;
            size_t end = std::min(begin + CHUNK, total);
            processBatchRange_(embedding, batch, begin, end, callbacks);
            processedCount.fetch_add(end - begin, std::memory_order_relaxed);
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
    if (!links.empty()) {
        descriptor = describeBoundary_(links);
        linkTally_.record(descriptor, type);
    }

    if (callbacks.onSurfaceBoundaryProcessed) {
        auto [orientable, genus, punctures] = type;
        callbacks.onSurfaceBoundaryProcessed(SurfaceBoundaryInfo{
            SurfaceFoundInfo{
                .orientable = orientable,
                .genus = genus,
                .punctures = punctures,
                .triangleCount = static_cast<long long>(faceIndices.size()),
                .mostRestrictive = classifyByLinks_(links)},
            descriptor});
    }

    // Reverse order, mirroring how the DFS itself would back out --
    // resets embedding to empty for the next entry.
    for (auto it = faceIndices.rbegin(); it != faceIndices.rend(); ++it)
        embedding.removeFace(*it);
}

void SurfaceSearch::processBatchRange_(
    KnottedSurface &embedding, const std::vector<std::vector<int>> &batch,
    size_t begin, size_t end, const SurfaceSearchCallbacks &callbacks) {
    for (size_t i = begin; i < end; ++i)
        processEntry_(embedding, batch[i], callbacks);
}

SearchStats SurfaceSearch::search(unsigned numThreads, BoundaryCondition cond,
                                  const SurfaceSearchCallbacks &callbacks,
                                  unsigned iddfsIterations, long long iddfsStep,
                                  std::optional<long long> iddfsStart,
                                  std::optional<unsigned> finalThreads) {
    const bool wantLinks = cond == BoundaryCondition::proper ||
                           cond == BoundaryCondition::connected;

    // Thread safety, not just performance: Vertex<4>::buildLink() caches
    // its result as a plain, unsynchronized lazily-constructed pointer
    // (see engine/triangulation/dim4/vertex4.cpp) -- since every worker
    // thread below shares the same ambient Triangulation<4> (and hence the
    // same Vertex<4> objects) even though each has its own KnottedSurface,
    // two workers concurrently reaching faces around the same ambient
    // vertex could otherwise race on that cache. Forcing every vertex's
    // link to build once, single-threaded, here -- before runSearch_
    // spawns any worker thread -- means every later call (even
    // concurrent) only ever takes the already-built, read-only path.
    for (auto v : skeleton_.triangulation().vertices())
        v->buildLink();

    // Same reasoning, for BoundaryComponent<4>::build(): every
    // KnottedSurface constructor independently calls
    // tri.boundaryComponent(c)->build() to populate its own bdryComponents_
    // (embeddedsubmanifold.cpp), but build() itself lazily caches its
    // result on the shared BoundaryComponent<4> object (common to every
    // KnottedSurface, since it belongs to the one ambient Triangulation<4>)
    // -- so two threads racing to build the SAME boundary component for
    // the first time (e.g. backgroundDrainLoop_'s aux thread against a DFS
    // worker's own embedding, or one worker's helperEmbedding_ under
    // backpressure against another's primary embedding) can otherwise race
    // on that cache. TSan caught exactly this (a write race inside
    // Triangulation<3>'s copy constructor, reached via two concurrent
    // build() calls) once backpressure started constructing KnottedSurface
    // objects from more threads at once. Pre-building every boundary
    // component once here, before any thread spawns, closes it the same
    // way the vertex-link loop above already does.
    for (size_t c = 0; c < skeleton_.triangulation().countBoundaryComponents();
         ++c)
        skeleton_.triangulation().boundaryComponent(c)->build();

    AuxHooks auxHooks(*this, numThreads, wantLinks, callbacks);
    return runSearch_<KnottedSurface>(
        numThreads, cond,
        [this] { return KnottedSurface(skeleton_, petalCache_); },
        [this, wantLinks, &callbacks] {
            return std::make_unique<ThreadHook>(*this, surfaceTypeTally_,
                                                wantLinks, callbacks);
        },
        [this, wantLinks, &callbacks](const std::vector<int> &seedFaces) {
            // One-off, not hot-path: rebuilds the seed as a KnottedSurface
            // (rather than reusing runSearch_'s generic proto embedding) to
            // get at surfaceType()/boundaryLinks(), which only KnottedSurface
            // exposes. seedFaces is added in the same order used to validate
            // it originally (see buildSeededGraph_), so this always embeds.
            KnottedSurface probe(skeleton_, petalCache_, seedFaces);
            SurfaceTypeKey type = probe.surfaceType();
            std::map<SurfaceTypeKey, long long> seedTypeCounts{{type, 1}};
            surfaceTypeTally_.merge(seedTypeCounts);
            auto [orientable, genus, punctures] = type;
            auto triangleCount = static_cast<long long>(seedFaces.size());
            if (wantLinks) {
                auto links = probe.boundaryLinks();
                std::string descriptor;
                if (!links.empty()) {
                    descriptor = describeBoundary_(links);
                    linkTally_.record(descriptor, type);
                }
                if (callbacks.onSurfaceBoundaryProcessed)
                    callbacks.onSurfaceBoundaryProcessed(SurfaceBoundaryInfo{
                        SurfaceFoundInfo{.orientable = orientable,
                                        .genus = genus,
                                        .punctures = punctures,
                                        .triangleCount = triangleCount,
                                        .mostRestrictive =
                                            classifyByLinks_(links)},
                        descriptor});
            } else if (callbacks.onSurfaceFound) {
                callbacks.onSurfaceFound(SurfaceFoundInfo{
                    .orientable = orientable,
                    .genus = genus,
                    .punctures = punctures,
                    .triangleCount = triangleCount,
                    .mostRestrictive = classifyCheaply_(probe)});
            }
        },
        callbacks, auxHooks,
        iddfsIterations, iddfsStep, iddfsStart, finalThreads);
}

