//
//  pairsigner.h
//
//  Pair signatures for cobordisms found by a search.
//

/*! \file utils/surfer/cobound/cobordisms/pairsigner.h
 *  \brief Pair signatures for the surfaces a search kept, from their faces
 *  and the thickening they were found in: one at a time, or many at once
 *  with each thickening's ambient part (PairSigContext) computed once.
 */

#pragma once

#include <atomic>
#include <chrono>
#include <condition_variable>
#include <deque>
#include <functional>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include <triangulation/dim4.h>

#include "cobound/cobordisms/cobordism.h"
#include "surfer/pairsig/pairsig.h"

namespace cascade {

/// The pair signature of a kept surface, from its faces and the thickening
/// it was found in (rebuilt from the row's PD by rowsearch::buildRow(), or
/// the searched one itself).
std::string pairSigOf(const regina::Triangulation<4> &thickening,
                      const std::vector<int> &faces);

/// A kept surface to sign: the row it was found in, and its faces there.
struct SignRequest {
  std::string rowPD;
  int layers = 2;
  std::vector<int> faces;
};

/// pairSigOf() for many surfaces at once. Almost all of a signature's cost is
/// the ambient's own (~50 s for a 10-crossing row, 2026-09-28), so each
/// distinct row is rebuilt (rowsearch::buildRow(), deterministic, so faces
/// index it as they did the searched one) and its ambient part computed once
/// (PairSigContext), up to `threads` rows at a time. In request order.
/// With a `cacheDir`, each row's context is read from there when stored and
/// stored when built (PairSigContext::cached()), so a row signed again in a
/// later run -- a node the cascade meets often -- costs no rebuild.
std::vector<std::string> pairSigsOf(const std::vector<SignRequest> &requests,
                                    unsigned threads, const std::string &cacheDir = "");

/// A thickening's pair-signature context (its ambient part, nearly all of a
/// signature's cost) on `threads` threads: built in memory, or with a
/// `cacheDir`, through the one cache format (PairSigContext::cached(): a
/// `.pairsigctx` file per ambient, checked on every use). `loaded` (if
/// given) says whether the cache held it.
std::unique_ptr<PairSigContext<4, 2>>
pairSigContextFor(const regina::Triangulation<4> &thickening, const std::string &cacheDir,
                  unsigned threads = 1, bool *loaded = nullptr);

// Signs a row's new witnesses off the drain threads (2026-09-29).
//
// A witness's pair signature needs the row's pair-signature context (the
// whole ambient's isoSigDetail: ~27 s at 10 crossings on halcyon, 2-3 min on
// a loaded yoga), built by the first signature. Signed where it was found, the
// first new witness -- the seed, on a fresh row, taken by the aux drain thread
// -- held up the whole drain for that long, and any drain thread that met
// another new witness waited too.
//
// So the drain only queues a witness's faces (captureFaces(), exactly what
// pairSig() signs), and one signer thread signs them in order, publishing
// each as it is signed, so checkpoints still persist them mid-row. The
// signer builds the context when the first witness arrives -- at once on a
// fresh row, whose seed is new -- so it overlaps the search, and a row that
// finds nothing new never builds it (LazyPairSigContext's reason for being
// lazy). After the drain, finish() signs what is left with every thread the
// row had and waits, before the row's witnesses are written.
class WitnessSigner {
public:
  using Context = std::function<const PairSigContext<4, 2> &()>;
  using Publish = std::function<void(cobordismgraph::Witness &&)>;

  WitnessSigner(Context context, Publish publish)
      : context_(std::move(context)), publish_(std::move(publish)),
        start_(std::chrono::steady_clock::now()) {
    threads_.emplace_back([this] { loop(); });
  }

  WitnessSigner(const WitnessSigner &) = delete;
  WitnessSigner &operator=(const WitnessSigner &) = delete;

  // On any exit, the threads are joined (what is queued is still signed).
  ~WitnessSigner() {
    try {
      finish(1);
    } catch (...) {
    }
  }

  void add(cobordismgraph::Witness w, std::vector<int> faces) {
    {
      std::lock_guard<std::mutex> lock(mutex_);
      queue_.push_back({std::move(w), std::move(faces)});
    }
    cv_.notify_one();
  }

  // Signs everything still queued with `threads` threads in all, and waits.
  // Rethrows a failure to build the context or to sign.
  void finish(unsigned threads) {
    if (finished_)
      return;
    const auto t0 = std::chrono::steady_clock::now();
    {
      std::lock_guard<std::mutex> lock(mutex_);
      closing_ = true;
    }
    cv_.notify_all();
    for (unsigned t = 1; t < threads; ++t)
      threads_.emplace_back([this] { loop(); });
    for (std::thread &t : threads_)
      t.join();
    threads_.clear();
    finished_ = true;
    finishSeconds_ = std::chrono::duration<double>(
                         std::chrono::steady_clock::now() - t0)
                         .count();
    if (!error_.empty())
      throw std::runtime_error("signing a witness: " + error_);
  }

  long long signedCount() const { return signed_.load(); }
  long long signMillis() const { return signMillis_.load(); }
  // From the row's start until the context was ready; 0 if never needed.
  double contextSeconds() const { return contextSeconds_; }
  // What finish() added to the row: the wait for the last signatures.
  double finishSeconds() const { return finishSeconds_; }

private:
  struct Pending {
    cobordismgraph::Witness w;
    std::vector<int> faces;
  };

  void loop() {
    while (true) {
      Pending p;
      {
        std::unique_lock<std::mutex> lock(mutex_);
        cv_.wait(lock, [this] { return !queue_.empty() || closing_; });
        if (queue_.empty() || !error_.empty())
          return;
        p = std::move(queue_.front());
        queue_.pop_front();
      }
      try {
        // The first call builds the context; the others wait for it.
        const PairSigContext<4, 2> &context = context_();
        if (!contextReady_.exchange(true))
          contextSeconds_ = std::chrono::duration<double>(
                                std::chrono::steady_clock::now() - start_)
                                .count();
        const auto t0 = std::chrono::steady_clock::now();
        p.w.pairSig = context.sig(p.faces);
        signMillis_.fetch_add(
            std::chrono::duration_cast<std::chrono::milliseconds>(
                std::chrono::steady_clock::now() - t0)
                .count());
        signed_.fetch_add(1);
        publish_(std::move(p.w));
      } catch (const std::exception &e) {
        std::lock_guard<std::mutex> lock(mutex_);
        if (error_.empty())
          error_ = e.what();
        return;
      }
    }
  }

  Context context_;
  Publish publish_;
  std::chrono::steady_clock::time_point start_;
  std::mutex mutex_;
  std::condition_variable cv_;
  std::deque<Pending> queue_;
  bool closing_ = false;
  bool finished_ = false;
  std::string error_;
  std::vector<std::thread> threads_;
  std::atomic<long long> signed_{0};
  std::atomic<long long> signMillis_{0};
  std::atomic<bool> contextReady_{false};
  double contextSeconds_ = 0;
  double finishSeconds_ = 0;
};

} // namespace cascade
