// readbackcache.h
//
// Stored witnesses' read-backs, kept across runs.
//
// Reading a master witness back (WitnessRedrawer::outgoingLinkFast())
// enumerates the isomorphisms from its pair signature's ambient onto the
// row's thickening -- most of a run's single-threaded driver time on
// 2026-09-30, and the same witnesses of the same popular rows (L6a5{0;1},
// L7a7{0;0}, L10a147{0;0}, ...) were read again by nearly every target's run.
// The result is fixed by the witness and the row's exact thickening, so it is
// kept here: one file per row, one line per witness.
//
//   <dir>/<rowkey>.readback       rowkey = witnessKey(row PD + "|" + layers)
//     readback 1 <build digest>   the row's WitnessRedrawer::buildChecksum()
//     <witness key>\tok\t<link>   a read-back (serialiseLink())
//     <witness key>\tfail\t<why>  a read-back that fails (it always will)
//
// A file whose digest is not the row's as built now is ignored and replaced:
// knotbuilder or the thickening changed, so its edge numbers mean nothing.
// Appends hold an flock on the file; a torn last line is ignored.
#ifndef CASCADE_READBACKCACHE_H
#define CASCADE_READBACKCACHE_H

#include <optional>
#include <string>
#include <unordered_map>

#include "farsidecurves.h"

namespace cascade {

/// One cached read-back: the link, or why reading it back fails.
struct CachedReadBack {
  std::optional<farside::OutgoingLink> link;
  std::string why;
};

std::string serialiseLink(const farside::OutgoingLink &link);
/// Inverse of serialiseLink(); nullopt if the text is malformed.
std::optional<farside::OutgoingLink> parseLink(const std::string &text);

/// One row's read-backs: loaded from its file, appended as new ones are
/// computed. Not shared between threads (one per row being read).
class RowReadBacks {
public:
  /// An empty dir means no cache (get() finds nothing, put() keeps nothing).
  RowReadBacks(const std::string &dir, const std::string &rowPD, int layers,
               const std::string &buildDigest);

  const CachedReadBack *get(const std::string &witnessKey) const;
  /// Records one read-back, in memory and (buffered) for the file.
  void put(const std::string &witnessKey, const CachedReadBack &r);
  /// Appends what put() recorded since the last flush, under the file's lock.
  void flush();

  size_t hits() const { return hits_; }
  size_t loaded() const { return entries_.size(); }

private:
  std::string path_, digest_;
  bool rewrite_ = false; // the file's digest was stale (or absent): start it afresh
  std::unordered_map<std::string, CachedReadBack> entries_;
  std::string pending_;
  mutable size_t hits_ = 0;
};

} // namespace cascade

#endif
