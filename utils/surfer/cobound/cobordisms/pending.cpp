// keptstore.cpp

#include "cobound/cobordisms/pending.h"

#include <algorithm>
#include <cassert>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

#include "surfer/report/atomicwrite.h"
#include "surfer/report/csvwriter.h"
#include "cobound/cobordisms/appendonly.h"
#include "cobound/cobordisms/pairsigner.h"
#include "cobound/cobordisms/database.h"

namespace fs = std::filesystem;

namespace cobordisms {

namespace {

// An exclusive flock(2) on the store's sidecar lock file, released on
// destruction.
struct StoreLock : appendonly::FileLock {
  explicit StoreLock(const std::string &store) : appendonly::FileLock(store + ".lock") {}
};

void identitiesOf(const std::string &path, std::unordered_set<std::string> &into,
                  unsigned threads = 1) {
  if (!fs::exists(path)) return;
  std::unordered_set<std::string> ids = cobordisms::witnessIdentities(path, threads);
  if (into.empty()) into.swap(ids);
  else into.merge(ids);
}

// The identities of the witness lines appended to `path` after its first
// `from` bytes (complete lines only, as loadWitnesses() reads them).
void identitiesSince(const std::string &path, std::uintmax_t from,
                     std::unordered_set<std::string> &into) {
  if (!fs::exists(path) || fs::file_size(path) <= from) return;
  std::ifstream in(path, std::ios::binary);
  in.seekg(static_cast<std::streamoff>(from));
  std::string line;
  if (from == 0) std::getline(in, line); // the header
  while (std::getline(in, line)) {
    if (in.eof()) break; // no newline: a torn last line
    if (line.empty()) continue;
    cobordisms::Witness w;
    if (cobordisms::parseWitnessLine(line, w, false, false, path))
      into.insert(cobordisms::witnessIdentity(w));
  }
}

int hopNumber(const fs::path &dir) {
  // hop_<k>_n<node>
  const std::string name = dir.filename().string();
  try {
    return std::stoi(name.substr(4));
  } catch (const std::exception &) {
    return 1 << 30;
  }
}

} // namespace

std::string formatKept(const PendingWitness &p) {
  std::ostringstream faces;
  for (size_t i = 0; i < p.faces.size(); ++i) faces << (i ? " " : "") << p.faces[i];
  return cobordisms::formatWitness(p.witness) + ',' + csvField(faces.str()) + ',' +
         csvField(p.rowPD) + ',' + std::to_string(p.layers) + '\n';
}

void appendKept(const std::string &hopDir, const std::vector<PendingWitness> &kept) {
  if (kept.empty()) return;
  std::string buffer;
  for (const PendingWitness &p : kept) buffer += formatKept(p);
  appendonly::append(hopDir + "/kept.csv", buffer, appendonly::Sync::yes);
}

long long signedThrough(const std::string &path) {
  std::ifstream in(path + ".signed");
  long long bytes = 0;
  if (in >> bytes) return bytes;
  return 0;
}

std::vector<PendingWitness> readKept(const std::string &work,
                                     std::vector<std::pair<std::string, long long>> *readTo) {
  std::vector<fs::path> dirs;
  if (fs::exists(work))
    for (const auto &e : fs::directory_iterator(work))
      if (e.is_directory() && e.path().filename().string().rfind("hop_", 0) == 0 &&
          fs::exists(e.path() / "kept.csv"))
        dirs.push_back(e.path());
  std::sort(dirs.begin(), dirs.end(), [](const fs::path &a, const fs::path &b) {
    return hopNumber(a) < hopNumber(b);
  });
  std::vector<PendingWitness> out;
  for (const fs::path &dir : dirs) {
    const fs::path path = dir / "kept.csv";
    std::ifstream in(path, std::ios::binary);
    const long long from = signedThrough(path.string());
    long long to = from;
    in.seekg(from);
    std::string line;
    while (std::getline(in, line)) {
      if (in.eof()) break; // no newline: a torn last line
      to += static_cast<long long>(line.size()) + 1;
      if (line.empty()) continue;
      std::vector<std::string> f = parseCsvLine(line);
      if (f.size() != 16)
        throw std::runtime_error(path.string() + ": a kept line has " +
                                 std::to_string(f.size()) + " fields, not 16");
      PendingWitness p;
      if (!cobordisms::witnessFromFields(f, p.witness, false, false, path))
        throw std::runtime_error(path.string() + ": a malformed kept witness");
      std::istringstream faces(f[13]);
      for (int t; faces >> t;) p.faces.push_back(t);
      p.rowPD = f[14];
      p.layers = std::stoi(f[15]);
      out.push_back(std::move(p));
    }
    if (readTo) readTo->emplace_back(path.string(), to);
  }
  return out;
}

StoreResult signPending(const std::string &work, const std::string &store,
                        const std::vector<std::string> &dedupeAgainst,
                        const solver::NameTable &names, unsigned threads,
                        const std::string &pairSigCache, const LoadedStore &loaded,
                        const std::function<bool(const PendingWitness &)> &sidecarLine) {
  std::vector<std::pair<std::string, long long>> readTo;
  std::vector<PendingWitness> pending = readKept(work, &readTo);
  StoreResult r = storeKept(std::move(pending), store, dedupeAgainst, names, threads,
                            pairSigCache, loaded, sidecarLine);
  // Signed (or already in the store): recorded only once the store holds
  // them, so a failure leaves the file to be signed again.
  for (const auto &[path, bytes] : readTo)
    if (bytes > signedThrough(path))
      report::atomicWrite(path + ".signed",
                          [&](std::ostream &out) { out << bytes << '\n'; });
  return r;
}

StoreResult storeKept(std::vector<PendingWitness> pending, const std::string &store,
                      const std::vector<std::string> &dedupeAgainst,
                      const solver::NameTable &names, unsigned threads,
                      const std::string &pairSigCache, const LoadedStore &loaded,
                      const std::function<bool(const PendingWitness &)> &sidecarLine) {
  StoreResult r;
  r.kept = pending.size();
  if (pending.empty()) return r;

  // The store's identities as they stand: the run's own copy of what it
  // loaded and what was appended since, or the whole file.
  auto storeIdentities = [&](std::unordered_set<std::string> &into) {
    if (!loaded.identities) {
      identitiesOf(store, into);
      return;
    }
    for (const std::string &id : *loaded.identities) into.insert(id);
    identitiesSince(store, loaded.bytes, into);
  };

  // Fresh against the read-only stores (once) and the store as it stands.
  const auto tDedupe = std::chrono::steady_clock::now();
  std::unordered_set<std::string> seen;
  for (const std::string &path : dedupeAgainst) identitiesOf(path, seen, threads);
  {
    StoreLock lock(store);
    storeIdentities(seen);
  }
  r.dedupeSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - tDedupe).count();
  std::vector<PendingWitness> fresh;
  for (PendingWitness &p : pending)
    if (seen.insert(cobordisms::witnessIdentity(p.witness)).second)
      fresh.push_back(std::move(p));
  r.fresh = fresh.size();
  if (fresh.empty()) return r;

  std::vector<SignRequest> requests;
  requests.reserve(fresh.size());
  for (const PendingWitness &p : fresh) requests.push_back({p.rowPD, p.layers, p.faces});
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<std::string> sigs = pairSigsOf(requests, threads, pairSigCache);
  r.signSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  std::vector<cobordisms::Witness> out;
  // parallel to out: row PD, layers, and whether the sidecar records it
  std::vector<std::pair<std::string, int>> rowOf;
  std::vector<char> sidecar;
  out.reserve(fresh.size());
  for (size_t i = 0; i < fresh.size(); ++i) {
    rowOf.emplace_back(fresh[i].rowPD, fresh[i].layers);
    sidecar.push_back(!sidecarLine || sidecarLine(fresh[i]));
    cobordisms::Witness w = std::move(fresh[i].witness);
    w.pairSig = std::move(sigs[i]);
    if (w.pairSig.empty())
      throw std::runtime_error("storeKept: an empty pair signature for a " + w.subject +
                               " witness");
    w.otherCandidates.clear();
    if (w.kind == cobordisms::WitnessKind::cobordism)
      w.otherCandidates = names.candidates(w.other, w.otherComponents);
    out.push_back(std::move(w));
  }

  // Under the lock: re-read what the store holds now (another run may have
  // appended since), and append only what is still new.
  StoreLock lock(store);
  std::unordered_set<std::string> now;
  storeIdentities(now);
  std::vector<cobordisms::Witness> append;
  std::vector<std::pair<std::string, int>> appendRows;
  std::vector<char> appendSidecar;
  for (size_t i = 0; i < out.size(); ++i)
    if (!now.count(cobordisms::witnessIdentity(out[i]))) {
      append.push_back(std::move(out[i]));
      appendRows.push_back(rowOf[i]);
      appendSidecar.push_back(sidecar[i]);
    }
  cobordisms::appendWitnesses(store, append, 0);
  r.appended = append.size();
  // Which diagram each pair signature's ambient was built from: a hop's row
  // is a node's own diagram, not a table PD, so the atlas's farsidename
  // (which rebuilds the row to redraw a far side) needs it. witness key
  // (sha1(pairsig)[:12]), layers, row PD; appended beside the store.
  std::string lines;
  for (size_t i = 0; i < append.size(); ++i)
    if (appendSidecar[i])
      lines += append[i].pairSigKey + ',' + std::to_string(appendRows[i].second) + ',' +
               csvField(appendRows[i].first) + '\n';
  if (!lines.empty()) {
    const std::string sidecarPath = store + ".rows.csv";
    if (!fs::exists(sidecarPath)) lines = "witness,layers,row_pd\n" + lines;
    // fsynced like the store it describes (it was not, before phase 3).
    appendonly::append(sidecarPath, lines, appendonly::Sync::yes);
  }
  return r;
}

namespace {

constexpr long long CHECKPOINT_INTERVAL_MS = 60'000;

long long tickNow() {
  return std::chrono::duration_cast<std::chrono::milliseconds>(
             std::chrono::steady_clock::now().time_since_epoch())
      .count();
}

} // namespace

PendingWriter::PendingWriter(fs::path path) : path_(std::move(path)) {
  std::error_code ec;
  const auto size = fs::file_size(path_, ec);
  synced_ = ec ? 0 : static_cast<long long>(size);
  lastWriteTick_.store(tickNow(), std::memory_order_relaxed);
}

void PendingWriter::add(const PendingWitness &p) {
  std::string line = formatKept(p);
  std::lock_guard<std::mutex> lock(mutex_);
  queued_ += line;
}

void PendingWriter::write_() {
  if (queued_.empty()) return;
  appendonly::append(path_.string(), queued_, appendonly::Sync::yes);
  queued_.clear();
  synced_ = static_cast<long long>(fs::file_size(path_));
}

void PendingWriter::checkpoint(bool force) {
  const long long now = tickNow();
  if (!force &&
      now - lastWriteTick_.load(std::memory_order_relaxed) < CHECKPOINT_INTERVAL_MS)
    return;
  std::lock_guard<std::mutex> lock(mutex_);
  if (!force &&
      now - lastWriteTick_.load(std::memory_order_relaxed) < CHECKPOINT_INTERVAL_MS)
    return;
  lastWriteTick_.store(now, std::memory_order_relaxed);
  try {
    write_();
  } catch (const std::exception &e) {
    // A failed checkpoint must not kill a running search: what is queued
    // stays queued, and the next write (the search's end, at the latest)
    // tries again and may throw.
    std::cerr << "[!] pending checkpoint failed: " << e.what() << "\n";
  }
}

void PendingWriter::flush() {
  std::lock_guard<std::mutex> lock(mutex_);
  write_();
}

long long PendingWriter::syncedBytes() const {
  std::lock_guard<std::mutex> lock(mutex_);
  return synced_;
}

namespace {

// The length of `path` up to the end of its last complete line: a torn last
// line is not loaded (loadWitnesses()), and the next append cuts it.
std::uintmax_t completeBytes(const fs::path &path) {
  std::error_code ec;
  const std::uintmax_t size = fs::file_size(path, ec);
  if (ec || size == 0) return 0;
  std::ifstream in(path, std::ios::binary);
  std::uintmax_t end = size;
  char c = 0;
  while (end > 0) {
    in.seekg(static_cast<std::streamoff>(end - 1));
    if (!in.get(c)) return 0;
    if (c == '\n') return end;
    --end;
  }
  return 0;
}

} // namespace

RecordedWitnesses::RecordedWitnesses(fs::path path,
                                     std::vector<cobordisms::Witness> loaded)
    : path_(std::move(path)), witnesses_(std::move(loaded)) {
  bytes_ = completeBytes(path_);
  identities_.reserve(witnesses_.size() * 2 + 1024);
  for (const cobordisms::Witness &w : witnesses_)
    identities_.insert(cobordisms::witnessIdentity(w));
}

} // namespace cobordisms
