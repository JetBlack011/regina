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
#include "cobound/frozen.h"

namespace fs = std::filesystem;

namespace cobordisms {

namespace {

// An exclusive flock(2) on the store's sidecar lock file, released on
// destruction.
struct DatabaseLock : appendonly::FileLock {
  explicit DatabaseLock(const std::string &database) : appendonly::FileLock(database + ".lock") {}
};

void identitiesOf(const std::string &path, std::unordered_set<std::string> &into,
                  unsigned threads = 1) {
  if (!fs::exists(path)) return;
  std::unordered_set<std::string> ids = cobordisms::cobordismIdentities(path, threads);
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
    cobordisms::Cobordism w;
    if (cobordisms::parseCobordismLine(line, w, false, false, path))
      into.insert(cobordisms::cobordismIdentity(w));
  }
}

int searchNumber(const fs::path &dir) {
  // hop_<k>_n<node>
  const std::string name = dir.filename().string();
  try {
    return std::stoi(name.substr(sizeof kFrozenHopDirPrefix - 1));
  } catch (const std::exception &) {
    return 1 << 30;
  }
}

} // namespace

std::string formatKept(const PendingCobordism &p) {
  std::ostringstream faces;
  for (size_t i = 0; i < p.faces.size(); ++i) faces << (i ? " " : "") << p.faces[i];
  return cobordisms::formatCobordism(p.cobordism) + ',' + csvField(faces.str()) + ',' +
         csvField(p.incomingPD) + ',' + std::to_string(p.layers) + '\n';
}

void appendKept(const std::string &searchDir, const std::vector<PendingCobordism> &kept) {
  if (kept.empty()) return;
  std::string buffer;
  for (const PendingCobordism &p : kept) buffer += formatKept(p);
  appendonly::append(searchDir + "/kept.csv", buffer, appendonly::Sync::yes);
}

long long signedThrough(const std::string &path) {
  std::ifstream in(path + ".signed");
  long long bytes = 0;
  if (in >> bytes) return bytes;
  return 0;
}

std::vector<PendingCobordism> readKept(const std::string &work,
                                     std::vector<std::pair<std::string, long long>> *readTo) {
  std::vector<fs::path> dirs;
  if (fs::exists(work))
    for (const auto &e : fs::directory_iterator(work))
      if (e.is_directory() && e.path().filename().string().rfind(kFrozenHopDirPrefix, 0) == 0 &&
          fs::exists(e.path() / "kept.csv"))
        dirs.push_back(e.path());
  std::sort(dirs.begin(), dirs.end(), [](const fs::path &a, const fs::path &b) {
    return searchNumber(a) < searchNumber(b);
  });
  std::vector<PendingCobordism> out;
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
      PendingCobordism p;
      if (!cobordisms::cobordismFromFields(f, p.cobordism, false, false, path))
        throw std::runtime_error(path.string() + ": a malformed kept witness");
      std::istringstream faces(f[13]);
      for (int t; faces >> t;) p.faces.push_back(t);
      p.incomingPD = f[14];
      p.layers = std::stoi(f[15]);
      out.push_back(std::move(p));
    }
    if (readTo) readTo->emplace_back(path.string(), to);
  }
  return out;
}

SignResult signPending(const std::string &work, const std::string &database,
                        const std::vector<std::string> &dedupeAgainst,
                        const solver::NameTable &names, unsigned threads,
                        const std::string &pairSigCache, const LoadedPrefix &loaded,
                        const std::function<bool(const PendingCobordism &)> &sidecarLine) {
  std::vector<std::pair<std::string, long long>> readTo;
  std::vector<PendingCobordism> pending = readKept(work, &readTo);
  SignResult r = signKept(std::move(pending), database, dedupeAgainst, names, threads,
                            pairSigCache, loaded, sidecarLine);
  // Signed (or already in the store): recorded only once the store holds
  // them, so a failure leaves the file to be signed again.
  for (const auto &[path, bytes] : readTo)
    if (bytes > signedThrough(path))
      report::atomicWrite(path + ".signed",
                          [&](std::ostream &out) { out << bytes << '\n'; });
  return r;
}

SignResult signKept(std::vector<PendingCobordism> pending, const std::string &database,
                      const std::vector<std::string> &dedupeAgainst,
                      const solver::NameTable &names, unsigned threads,
                      const std::string &pairSigCache, const LoadedPrefix &loaded,
                      const std::function<bool(const PendingCobordism &)> &sidecarLine) {
  SignResult r;
  r.kept = pending.size();
  if (pending.empty()) return r;

  // The store's identities as they stand: the run's own copy of what it
  // loaded and what was appended since, or the whole file.
  auto databaseIdentities = [&](std::unordered_set<std::string> &into) {
    if (!loaded.identities) {
      identitiesOf(database, into);
      return;
    }
    for (const std::string &id : *loaded.identities) into.insert(id);
    identitiesSince(database, loaded.bytes, into);
  };

  // Fresh against the read-only stores (once) and the store as it stands.
  const auto tDedupe = std::chrono::steady_clock::now();
  std::unordered_set<std::string> seen;
  for (const std::string &path : dedupeAgainst) identitiesOf(path, seen, threads);
  {
    DatabaseLock lock(database);
    databaseIdentities(seen);
  }
  r.dedupeSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - tDedupe).count();
  std::vector<PendingCobordism> fresh;
  for (PendingCobordism &p : pending)
    if (seen.insert(cobordisms::cobordismIdentity(p.cobordism)).second)
      fresh.push_back(std::move(p));
  r.fresh = fresh.size();
  if (fresh.empty()) return r;

  std::vector<SignRequest> requests;
  requests.reserve(fresh.size());
  for (const PendingCobordism &p : fresh) requests.push_back({p.incomingPD, p.layers, p.faces});
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<std::string> sigs = pairSigsOf(requests, threads, pairSigCache);
  r.signSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  std::vector<cobordisms::Cobordism> out;
  // parallel to out: row PD, layers, and whether the sidecar records it
  std::vector<std::pair<std::string, int>> diagramOf;
  std::vector<char> sidecar;
  out.reserve(fresh.size());
  for (size_t i = 0; i < fresh.size(); ++i) {
    diagramOf.emplace_back(fresh[i].incomingPD, fresh[i].layers);
    sidecar.push_back(!sidecarLine || sidecarLine(fresh[i]));
    cobordisms::Cobordism w = std::move(fresh[i].cobordism);
    w.pairSig = std::move(sigs[i]);
    if (w.pairSig.empty())
      throw std::runtime_error("storeKept: an empty pair signature for a " + w.subject +
                               " witness");
    w.otherCandidates.clear();
    if (w.kind == cobordisms::CobordismKind::cobordism)
      w.otherCandidates = names.candidates(w.other, w.otherComponents);
    out.push_back(std::move(w));
  }

  // Under the lock: re-read what the store holds now (another run may have
  // appended since), and append only what is still new.
  DatabaseLock lock(database);
  std::unordered_set<std::string> now;
  databaseIdentities(now);
  std::vector<cobordisms::Cobordism> append;
  std::vector<std::pair<std::string, int>> appendDiagrams;
  std::vector<char> appendSidecar;
  for (size_t i = 0; i < out.size(); ++i)
    if (!now.count(cobordisms::cobordismIdentity(out[i]))) {
      append.push_back(std::move(out[i]));
      appendDiagrams.push_back(diagramOf[i]);
      appendSidecar.push_back(sidecar[i]);
    }
  cobordisms::appendCobordisms(database, append, 0);
  r.appended = append.size();
  // Which diagram each pair signature's ambient was built from: a hop's row
  // is a node's own diagram, not a table PD, so the atlas's farsidename
  // (which rebuilds the row to redraw a far side) needs it. witness key
  // (sha1(pairsig)[:12]), layers, row PD; appended beside the store.
  std::string lines;
  for (size_t i = 0; i < append.size(); ++i)
    if (appendSidecar[i])
      lines += append[i].pairSigKey + ',' + std::to_string(appendDiagrams[i].second) + ',' +
               csvField(appendDiagrams[i].first) + '\n';
  if (!lines.empty()) {
    const std::string sidecarPath = database + kFrozenRowsSidecarSuffix;
    if (!fs::exists(sidecarPath)) lines = kFrozenRowsSidecarHeader + lines;
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

void PendingWriter::add(const PendingCobordism &p) {
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

LoadedDatabase::LoadedDatabase(fs::path path,
                                     std::vector<cobordisms::Cobordism> loaded)
    : path_(std::move(path)), cobordisms_(std::move(loaded)) {
  bytes_ = completeBytes(path_);
  identities_.reserve(cobordisms_.size() * 2 + 1024);
  for (const cobordisms::Cobordism &w : cobordisms_)
    identities_.insert(cobordisms::cobordismIdentity(w));
}

} // namespace cobordisms
