//
//  parallelfor.h
//
//  The one thread loop: every index of a range, on a few threads.
//

#ifndef SURFER_COBOUND_PARALLELFOR_H
#define SURFER_COBOUND_PARALLELFOR_H

#include <algorithm>
#include <atomic>
#include <cstddef>
#include <thread>
#include <vector>

/*! \file utils/surfer/cobound/parallelfor.h
 *  \brief parallelFor(): independent jobs over an index range, on a few
 *  threads, each thread taking the next index as it finishes one.
 */

/**
 * Calls `work(i)` exactly once for every `i` in `[0, n)`, on
 * `min(max(threads, 1), n)` threads -- the calling thread one of them --
 * each taking the next unclaimed index as it finishes the last, and returns
 * when every call has. The calls run in no particular order: `work` must
 * write only to its own index's slots (or lock), and must not throw (catch
 * inside, as every caller records its failures).
 */
template <typename Work>
void parallelFor(size_t n, unsigned threads, Work &&work) {
    if (n == 0)
        return;
    const size_t k = std::min<size_t>(std::max(threads, 1u), n);
    std::atomic<size_t> next{0};
    auto loop = [&] {
        for (size_t i; (i = next.fetch_add(1)) < n;)
            work(i);
    };
    std::vector<std::thread> pool;
    pool.reserve(k - 1);
    for (size_t t = 1; t < k; ++t)
        pool.emplace_back(loop);
    loop();
    for (std::thread &t : pool)
        t.join();
}

#endif // SURFER_COBOUND_PARALLELFOR_H
