// keptstore.cpp

#include "keptstore.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <unordered_set>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

#include "csvwriter.h"
#include "hoprunner.h"
#include "witnessstore.h"

namespace fs = std::filesystem;

namespace cascade {

namespace {

// An exclusive flock(2) on a sidecar lock file, released on destruction.
class StoreLock {
public:
  explicit StoreLock(const std::string &store) {
    const std::string path = store + ".lock";
    fd_ = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
    if (fd_ < 0)
      throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
    while (::flock(fd_, LOCK_EX) != 0)
      if (errno != EINTR)
        throw std::runtime_error("cannot lock " + path + ": " + std::strerror(errno));
  }
  ~StoreLock() {
    if (fd_ >= 0) {
      ::flock(fd_, LOCK_UN);
      ::close(fd_);
    }
  }
  StoreLock(const StoreLock &) = delete;
  StoreLock &operator=(const StoreLock &) = delete;

private:
  int fd_ = -1;
};

void identitiesOf(const std::string &path, std::unordered_set<std::string> &into,
                  unsigned threads = 1) {
  if (!fs::exists(path)) return;
  std::unordered_set<std::string> ids = witnessstore::witnessIdentities(path, threads);
  if (into.empty()) into.swap(ids);
  else into.merge(ids);
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

void appendKept(const std::string &hopDir, const std::vector<PendingWitness> &kept) {
  if (kept.empty()) return;
  std::string buffer;
  for (const PendingWitness &p : kept) {
    std::ostringstream faces;
    for (size_t i = 0; i < p.faces.size(); ++i) faces << (i ? " " : "") << p.faces[i];
    buffer += witnessstore::formatWitness(p.witness) + ',' + csvField(faces.str()) + ',' +
              csvField(p.rowPD) + ',' + std::to_string(p.layers) + '\n';
  }
  const std::string path = hopDir + "/kept.csv";
  const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_APPEND | O_CLOEXEC, 0644);
  if (fd < 0) throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
  const char *p = buffer.data();
  size_t left = buffer.size();
  while (left > 0) {
    const ssize_t n = ::write(fd, p, left);
    if (n < 0) {
      if (errno == EINTR) continue;
      ::close(fd);
      throw std::runtime_error("write to " + path + " failed: " + std::strerror(errno));
    }
    p += n;
    left -= static_cast<size_t>(n);
  }
  const int synced = ::fsync(fd);
  ::close(fd);
  if (synced != 0) throw std::runtime_error("fsync of " + path + " failed");
}

std::vector<PendingWitness> readKept(const std::string &work) {
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
    std::string line;
    while (std::getline(in, line)) {
      if (in.eof()) break; // no newline: a torn last line
      if (line.empty()) continue;
      std::vector<std::string> f = parseCsvLine(line);
      if (f.size() != 16)
        throw std::runtime_error(path.string() + ": a kept line has " +
                                 std::to_string(f.size()) + " fields, not 16");
      PendingWitness p;
      if (!witnessstore::witnessFromFields(f, p.witness, false, false, path))
        throw std::runtime_error(path.string() + ": a malformed kept witness");
      std::istringstream faces(f[13]);
      for (int t; faces >> t;) p.faces.push_back(t);
      p.rowPD = f[14];
      p.layers = std::stoi(f[15]);
      out.push_back(std::move(p));
    }
  }
  return out;
}

StoreResult storeKept(std::vector<PendingWitness> pending, const std::string &store,
                      const std::vector<std::string> &dedupeAgainst,
                      const cobordismgraph::NameTable &names, unsigned threads,
                      const std::string &pairSigCache) {
  StoreResult r;
  r.kept = pending.size();
  if (pending.empty()) return r;

  // Fresh against the read-only stores (once) and the store as it stands.
  const auto tDedupe = std::chrono::steady_clock::now();
  std::unordered_set<std::string> seen;
  for (const std::string &path : dedupeAgainst) identitiesOf(path, seen, threads);
  {
    StoreLock lock(store);
    identitiesOf(store, seen);
  }
  r.dedupeSeconds =
      std::chrono::duration<double>(std::chrono::steady_clock::now() - tDedupe).count();
  std::vector<PendingWitness> fresh;
  for (PendingWitness &p : pending)
    if (seen.insert(cobordismgraph::witnessIdentity(p.witness)).second)
      fresh.push_back(std::move(p));
  r.fresh = fresh.size();
  if (fresh.empty()) return r;

  std::vector<SignRequest> requests;
  requests.reserve(fresh.size());
  for (const PendingWitness &p : fresh) requests.push_back({p.rowPD, p.layers, p.faces});
  const auto t0 = std::chrono::steady_clock::now();
  std::vector<std::string> sigs = pairSigsOf(requests, threads, pairSigCache);
  r.signSeconds = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();

  std::vector<cobordismgraph::Witness> out;
  std::vector<std::pair<std::string, int>> rowOf; // parallel to out: row PD, layers
  out.reserve(fresh.size());
  for (size_t i = 0; i < fresh.size(); ++i) {
    rowOf.emplace_back(fresh[i].rowPD, fresh[i].layers);
    cobordismgraph::Witness w = std::move(fresh[i].witness);
    w.pairSig = std::move(sigs[i]);
    if (w.pairSig.empty())
      throw std::runtime_error("storeKept: an empty pair signature for a " + w.subject +
                               " witness");
    w.otherCandidates.clear();
    if (w.kind == cobordismgraph::WitnessKind::cobordism)
      w.otherCandidates = names.candidates(w.other, w.otherComponents);
    out.push_back(std::move(w));
  }

  // Under the lock: re-read what the store holds now (another run may have
  // appended since), and append only what is still new.
  StoreLock lock(store);
  std::unordered_set<std::string> now;
  identitiesOf(store, now);
  std::vector<cobordismgraph::Witness> append;
  std::vector<std::pair<std::string, int>> appendRows;
  for (size_t i = 0; i < out.size(); ++i)
    if (!now.count(cobordismgraph::witnessIdentity(out[i]))) {
      append.push_back(std::move(out[i]));
      appendRows.push_back(rowOf[i]);
    }
  witnessstore::appendWitnesses(store, append, 0);
  r.appended = append.size();
  // Which diagram each pair signature's ambient was built from: a hop's row
  // is a node's own diagram, not a table PD, so the atlas's farsidename
  // (which rebuilds the row to redraw a far side) needs it. witness key
  // (sha1(pairsig)[:12]), layers, row PD; appended beside the store.
  if (!append.empty()) {
    std::string buffer;
    const std::string sidecar = store + ".rows.csv";
    if (!fs::exists(sidecar)) buffer += "witness,layers,row_pd\n";
    for (size_t i = 0; i < append.size(); ++i)
      buffer += append[i].pairSigKey + ',' + std::to_string(appendRows[i].second) + ',' +
                csvField(appendRows[i].first) + '\n';
    std::ofstream(sidecar, std::ios::app) << buffer;
  }
  return r;
}

} // namespace cascade
