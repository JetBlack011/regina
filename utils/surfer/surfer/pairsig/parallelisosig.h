//
//  parallelisosig.h
//
//  Regina's isoSigDetail() on several threads, byte for byte.
//

#pragma once

#include <algorithm>
#include <cstddef>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include <triangulation/dim4.h>
#include <triangulation/isosig.h>
#include <utilities/exception.h>

/**
 * Exactly `tri.isoSigDetail()` -- the same signature and the same
 * isomorphism onto its canonical form -- computed on `threads` threads.
 *
 * isoSigDetail() (engine/triangulation/detail/isosig-impl.h) runs through
 * every start (simplex, permutation) in IsoSigClassic's order, fills a
 * canonical labelling from each, encodes it with IsoSigPrintable, and keeps
 * the first start whose encoding is least (a strict `<`). Each start is
 * independent of every other, so start i goes to thread i mod `threads`,
 * which keeps its own least (encoding, i); the answer is the least over the
 * threads by (encoding, i) -- the first least start in the serial order,
 * hence the same isomorphism as well as the same signature.
 *
 * It calls the engine's own fillFrom() and IsoSigPrintable::encode(), which
 * the engine instantiates, so nothing of the algorithm is reimplemented.
 *
 * For a 12-crossing link's thickening this is ~99.8% of a pair-signature
 * context (pairsig.h), which a campaign's searches used to build on one thread
 * (README.md, "Performance").
 *
 * \pre `tri` is non-empty and connected, with its skeleton computed (any
 * skeletal query does that; this function makes one before it starts
 * threads, since Regina computes a skeleton lazily and would race).
 */
template <int dim>
std::pair<std::string, regina::Isomorphism<dim>>
parallelIsoSigDetail(const regina::Triangulation<dim> &tri, unsigned threads) {
    if (tri.isEmpty())
        throw regina::FailedPrecondition(
            "parallelIsoSigDetail() requires a non-empty triangulation");
    if (!tri.isConnected()) // also computes the skeleton, single-threaded
        throw regina::FailedPrecondition(
            "parallelIsoSigDetail() requires a connected triangulation");
    const size_t starts = tri.size() * regina::Perm<dim + 1>::nPerms;
    threads = static_cast<unsigned>(
        std::min<size_t>(std::max(threads, 1u), starts));
    if (threads == 1) return tri.isoSigDetail();

    regina::Component<dim> *c = tri.component(0);
    struct Best {
        bool have = false;
        size_t start = 0; ///< its index in IsoSigClassic's order
        std::string sig;
        regina::Isomorphism<dim> iso;
        explicit Best(size_t n) : iso(n) {}
    };
    std::vector<Best> best;
    best.reserve(threads);
    for (unsigned t = 0; t < threads; ++t) best.emplace_back(tri.size());

    auto work = [&](unsigned t) {
        Best &b = best[t];
        regina::IsoSigClassic<dim> type(*c, false /* first-generation sigs are unoriented */);
        regina::IsoSigData<1, dim> data(c);
        regina::Isomorphism<dim> curr(tri.size());
        size_t i = 0;
        do {
            if (i % threads == t) {
                data.fillFrom(c->simplex(type.simplex()), type.perm(), &curr);
                std::string s = regina::IsoSigPrintable::encode(data);
                if (!b.have || s < b.sig) {
                    b.sig.swap(s);
                    b.iso.swap(curr);
                    b.start = i;
                    b.have = true;
                }
            }
            ++i;
        } while (type.next());
    };
    std::vector<std::thread> pool;
    for (unsigned t = 1; t < threads; ++t) pool.emplace_back(work, t);
    work(0);
    for (auto &th : pool) th.join();

    size_t win = 0;
    for (size_t t = 1; t < best.size(); ++t) {
        const Best &b = best[t], &w = best[win];
        if (!b.have) continue;
        if (!w.have || b.sig < w.sig || (b.sig == w.sig && b.start < w.start)) win = t;
    }
    return {std::move(best[win].sig), std::move(best[win].iso)};
}
