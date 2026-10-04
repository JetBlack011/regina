//
//  database.cpp
//
//  The append-only cobordism database (cobordisms.csv), shared by
//  verifyslicegenus and cascadesearch. Moved verbatim from
//  verifyslicegenus.cpp (2026-09-28), as witnessstore.cpp.
//

#include "cobound/cobordisms/database.h"

#include "cobound/cobordisms/appendonly.h"
#include "cobound/parallelfor.h"


#include <cctype>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <algorithm>
#include <thread>
#include <iostream>
#include <sstream>
#include <stdexcept>

#include <fcntl.h>
#include <unistd.h>

#include "surfer/report/csvwriter.h"
#include "cobound/cobordisms/cobordismkey.h"
#include "cobound/frozen.h"
#include "linknaming/names.h"

namespace cobordisms {
// ─────────────────────────────────────────────────────────────────────────
// Witness file (--cobordisms) I/O
// ─────────────────────────────────────────────────────────────────────────
//
// The point of persisting these separately from --output is that a witness
// is a fact ("a surface with this boundary and this genus exists") while an
// --output row is a conclusion. Conclusions get better whenever the solver,
// the literature tables, or the identification improves; facts do not. So
// every fact a search paid for is written here once and never re-searched,
// and `--solve-only` re-derives all the conclusions from them in seconds.

// resolved_vertices (Witness::resolvedVertices) is written EMPTY when 0, which
// is every witness but those found under --resolve-unlinked. That keeps a row
// written before the column existed, or merged in by a Python DictWriter
// (which fills a missing field with ""), byte-identical after a --solve-only
// round trip -- the invariant merge_cobordisms.py relies on.
std::string formatCobordism(const cobordisms::Cobordism &w) {
  std::ostringstream candidates;
  for (size_t i = 0; i < w.otherCandidates.size(); ++i) {
    if (i)
      candidates << ';';
    candidates << w.otherCandidates[i];
  }
  std::ostringstream out;
  out << (w.kind == cobordisms::CobordismKind::direct ? "direct"
                                                        : "cobordism")
      << ',' << csvField(w.subject) << ',' << w.subjectComponents << ','
      << csvField(w.other) << ',' << csvField(candidates.str()) << ','
      << w.otherComponents << ',' << w.genus << ','
      << (w.tubed ? "true" : "false") << ',' << csvField(w.pairSig) << ','
      << csvField(w.sourceSearch) << ',' << w.thickenLayers << ',' << w.maxFaces
      << ',';
  if (w.resolvedVertices != 0)
    out << w.resolvedVertices;
  return out.str();
}

// The header of a witness file written before resolved_vertices existed.
constexpr const char *COBORDISMS_HEADER_12 =
    "kind,subject,subject_components,other,other_candidates,other_components,"
    "genus,tubed,pairsig,source_row,thicken_layers,max_faces";

// Reads the first 12 or 13 fields of a witness line into `w`, keeping the
// pair signature only if `keepPairSig`. Returns false for a malformed line.
bool cobordismFromFields(std::vector<std::string> f, cobordisms::Cobordism &w,
                       bool keepPairSig, bool wantPairSigKey,
                       const std::filesystem::path &path) {
  if (f.size() < 12)
    return false;
  w.kind = f[0] == "direct" ? cobordisms::CobordismKind::direct
                            : cobordisms::CobordismKind::cobordism;
  w.subject = f[1];
  try {
    w.subjectComponents = std::stoi(f[2]);
    w.otherComponents = std::stoi(f[5]);
    w.genus = std::stoi(f[6]);
    w.thickenLayers = std::stoi(f[10]);
    w.maxFaces = std::stoll(f[11]);
  } catch (const std::exception &) {
    return false;
  }
  w.other = f[3];
  w.otherCandidates = splitCandidates(f[4]);
  w.tubed = f[7] == "true";
  if (wantPairSigKey && !f[8].empty())
    w.pairSigKey = cobordisms::cobordismKey(f[8]);
  if (keepPairSig)
    w.pairSig = std::move(f[8]);
  w.sourceSearch = f[9];
  // Optional 13th field (absent in files written before it existed).
  // Informational only, so an unreadable value never costs the witness
  // itself.
  if (f.size() > 12 && !f[12].empty()) {
    try {
      w.resolvedVertices = std::stoi(f[12]);
    } catch (const std::exception &) {
      std::cerr << "[!] " << path.string() << ": unreadable resolved_vertices '"
                << f[12] << "' for a " << w.subject
                << " witness; keeping the witness, recording 0\n";
    }
  }
  return true;
}

std::vector<std::string> splitCandidates(const std::string &field) {
  // Only at depth 0: "L8a4{0;1};L8a4{1;0}" is two names, not four (before
  // phase 3 every ';' split, which was inert: no recorded candidate list
  // held a tagged link of three or more components).
  std::vector<std::string> out;
  std::string cur;
  int depth = 0;
  for (char c : field) {
    if (c == '{') ++depth;
    else if (c == '}' && depth > 0) --depth;
    if (c == ';' && depth == 0) {
      if (!cur.empty()) out.push_back(std::move(cur));
      cur.clear();
    } else {
      cur += c;
    }
  }
  if (!cur.empty()) out.push_back(std::move(cur));
  return out;
}

// Parses one witness line (12 or 13 fields) into `w`, keeping the pair
// signature only if `keepPairSig`. Returns false for a malformed line.
bool parseCobordismLine(const std::string &line, cobordisms::Cobordism &w,
                      bool keepPairSig, bool wantPairSigKey,
                      const std::filesystem::path &path) {
  return cobordismFromFields(parseCsvLine(line), w, keepPairSig, wantPairSigKey,
                           path);
}

namespace {

// Every complete witness line of `path` (a torn last line ignored, malformed
// lines counted and skipped), each with its line's byte offset.
std::vector<cobordisms::Cobordism>
readLines(const std::filesystem::path &path, bool keepPairSigs, bool wantPairSigKeys) {
  std::vector<cobordisms::Cobordism> result;
  std::ifstream in(path, std::ios::binary);
  if (!in)
    return result;

  std::string line;
  std::getline(in, line); // header
  if (line != COBORDISMS_HEADER && line != COBORDISMS_HEADER_12)
    std::cerr << "[!] " << path.string()
              << ": unexpected witness-file header; reading it anyway\n";
  size_t malformed = 0;
  while (true) {
    const std::streamoff offset = in.tellg();
    if (!std::getline(in, line))
      break;
    if (in.eof()) {
      if (!line.empty())
        std::cerr << "[!] " << path.string() << ": ignoring a torn last line ("
                  << line.size() << " bytes, no newline)\n";
      break;
    }
    if (line.empty())
      continue;
    cobordisms::Cobordism w;
    if (!parseCobordismLine(line, w, keepPairSigs, wantPairSigKeys, path)) {
      ++malformed;
      continue;
    }
    w.fileOffset = static_cast<long long>(offset);
    result.push_back(std::move(w));
  }
  if (malformed > 0)
    std::cerr << "[!] " << path.string() << ": skipped " << malformed
              << " malformed witness lines\n";
  return result;
}

} // namespace

// Loads every complete witness line, WITHOUT its pair signature: each
// witness keeps only the byte offset of its line (Witness::fileOffset), and
// anything that needs the signature reads it back from there. That is what
// keeps a solve's memory proportional to the number of witnesses rather
// than to the ~11 KB signature each one carries.
//
// A final line with no terminating newline is a torn append (the process
// died mid-write) and is ignored here; appendWitnesses() truncates it away
// before it next appends.
std::vector<cobordisms::Cobordism>
loadCobordisms(const std::filesystem::path &path, bool wantPairSigKeys) {
  return readLines(path, /*keepPairSigs=*/false, wantPairSigKeys);
}

std::vector<cobordisms::Cobordism> readCobordisms(const std::filesystem::path &path) {
  return readLines(path, /*keepPairSigs=*/true, /*wantPairSigKeys=*/false);
}

void PairSigReader::setPath(std::filesystem::path path) {
  std::lock_guard<std::mutex> lock(mutex_);
  path_ = std::move(path);
  cache_.clear();
}

std::string PairSigReader::at(long long offset) {
  if (offset < 0)
    return {};
  std::lock_guard<std::mutex> lock(mutex_);
  if (auto it = cache_.find(offset); it != cache_.end())
    return it->second;
  std::ifstream in(path_, std::ios::binary);
  std::string line;
  if (!in || !in.seekg(offset) || !std::getline(in, line))
    throw std::runtime_error("cannot read the witness at byte " +
                             std::to_string(offset) + " of " + path_.string());
  auto f = parseCsvLine(line);
  if (f.size() < 12)
    throw std::runtime_error("no witness line at byte " + std::to_string(offset) +
                             " of " + path_.string());
  return cache_.emplace(offset, std::move(f[8])).first->second;
}

std::string DatabaseIndex::base(std::string name) {
  name = linknaming::stripOrientationTag(name);
  if (name.size() > 1 && name[0] == 'm' && std::isdigit(static_cast<unsigned char>(name[1])))
    name.erase(0, 1);
  return name;
}

DatabaseIndex::DatabaseIndex(const std::string &path) : path_(path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot read " + path);
  std::string line;
  // The row-PD sidecar (pending.cpp: witness,layers,row_pd), if any: a
  // cobordism recorded by a cascade hop was searched on its node's
  // simplified diagram, and can only be read back on that row.
  if (std::ifstream side(path + kFrozenRowsSidecarSuffix); side) {
    std::getline(side, line);
    while (std::getline(side, line)) {
      auto f = parseCsvLine(line);
      if (f.size() >= 3 && !f[2].empty())
        incomingPD_[f[0]] = f[2];
    }
  }
  std::getline(in, line);
  std::streamoff at = in.tellg();
  while (std::getline(in, line)) {
    if (in.eof())
      break; // a torn last line, as loadWitnesses() leaves it out
    const auto a = line.find(','), b = line.find(',', a + 1);
    if (a != std::string::npos && b != std::string::npos) {
      offsets_[line.substr(a + 1, b - a - 1)].push_back(at);
      // field 3 (other) follows subject_components
      const auto c = line.find(',', b + 1),
                 d = c == std::string::npos ? std::string::npos : line.find(',', c + 1);
      if (d != std::string::npos)
        byOther_[base(line.substr(c + 1, d - c - 1))].push_back(at);
    }
    at = in.tellg();
  }
}

std::vector<StoredCobordism> DatabaseIndex::ofSubject(const std::string &subject) const {
  auto it = offsets_.find(subject);
  return it == offsets_.end() ? std::vector<StoredCobordism>{} : read(it->second, 1u << 30);
}

std::vector<StoredCobordism> DatabaseIndex::byOutgoing(const std::string &b, size_t cap) const {
  auto it = byOther_.find(b);
  return it == byOther_.end() ? std::vector<StoredCobordism>{} : read(it->second, cap);
}

std::vector<StoredCobordism> DatabaseIndex::read(const std::vector<std::streamoff> &offsets,
                                                 size_t cap) const {
  std::vector<StoredCobordism> out;
  std::ifstream in(path_, std::ios::binary);
  std::string line;
  for (std::streamoff off : offsets) {
    if (out.size() >= cap)
      break;
    in.seekg(off);
    if (!std::getline(in, line))
      continue;
    StoredCobordism s;
    if (!parseCobordismLine(line, s.cobordism, /*keepPairSig=*/true, /*wantPairSigKey=*/false,
                          path_))
      continue;
    s.cobordism.fileOffset = static_cast<long long>(off);
    if (!incomingPD_.empty())
      if (auto r = incomingPD_.find(cobordisms::cobordismKey(s.cobordism.pairSig)); r != incomingPD_.end())
        s.incomingPD = r->second;
    out.push_back(std::move(s));
  }
  return out;
}

std::unordered_set<std::string>
cobordismIdentities(const std::filesystem::path &path, unsigned threads,
                  std::streamoff minRangeBytes) {
  std::unordered_set<std::string> result;
  std::ifstream in(path, std::ios::binary);
  if (!in)
    return result;
  std::string line;
  std::getline(in, line); // header
  if (line != COBORDISMS_HEADER && line != COBORDISMS_HEADER_12)
    std::cerr << "[!] " << path.string()
              << ": unexpected witness-file header; reading it anyway\n";
  const std::streamoff first = in.tellg();
  if (first < 0)
    return result; // the header alone, or not even that
  const std::streamoff size =
      static_cast<std::streamoff>(std::filesystem::file_size(path));
  // A range per thread, each starting at a line start: a cut falls just
  // after the first newline at or past its even share. Small files: one.
  const size_t n = std::clamp<size_t>(
      static_cast<size_t>((size - first) / std::max<std::streamoff>(minRangeBytes, 1)), 1,
      std::max(threads, 1u));
  std::vector<std::streamoff> cut(n + 1, size);
  cut[0] = first;
  for (size_t k = 1; k < n; ++k) {
    std::streamoff pos = first + (size - first) * static_cast<std::streamoff>(k) /
                                     static_cast<std::streamoff>(n);
    pos = std::max(pos, cut[k - 1]);
    in.clear();
    in.seekg(pos - 1);
    char c = 0;
    while (in.get(c) && c != '\n') {
    }
    cut[k] = in ? static_cast<std::streamoff>(in.tellg()) : size;
  }
  std::vector<std::unordered_set<std::string>> sets(n);
  std::vector<size_t> malformed(n, 0);
  auto work = [&](size_t k) {
    std::ifstream r(path, std::ios::binary);
    r.seekg(cut[k]);
    std::string l;
    std::streamoff at = cut[k];
    while (at < cut[k + 1]) {
      if (!std::getline(r, l))
        break;
      at += static_cast<std::streamoff>(l.size()) + 1;
      if (r.eof()) {
        // As loadWitnesses(): a last line with no newline is torn.
        if (!l.empty())
          std::cerr << "[!] " << path.string() << ": ignoring a torn last line ("
                    << l.size() << " bytes, no newline)\n";
        break;
      }
      if (l.empty())
        continue;
      cobordisms::Cobordism w;
      if (!parseCobordismLine(l, w, /*keepPairSig=*/false, /*wantPairSigKey=*/false, path)) {
        ++malformed[k];
        continue;
      }
      sets[k].insert(cobordisms::cobordismIdentity(w));
    }
  };
  parallelFor(n, static_cast<unsigned>(n), work);
  size_t bad = 0;
  for (size_t k = 0; k < n; ++k) {
    bad += malformed[k];
    if (result.empty())
      result.swap(sets[k]);
    else
      result.merge(sets[k]);
  }
  if (bad > 0)
    std::cerr << "[!] " << path.string() << ": skipped " << bad
              << " malformed witness lines\n";
  return result;
}

namespace {
using appendonly::writeAll;

// The first line of the open file `fd` (without its newline).
std::string readHeader(int fd) {
  std::string header;
  char buf[4096];
  off_t pos = 0;
  while (true) {
    ssize_t n = ::pread(fd, buf, sizeof buf, pos);
    if (n <= 0)
      break;
    const char *nl = static_cast<const char *>(std::memchr(buf, '\n', n));
    header.append(buf, nl ? static_cast<size_t>(nl - buf)
                          : static_cast<size_t>(n));
    if (nl)
      break;
    pos += n;
  }
  return header;
}
} // namespace

// Appends witnesses[from..] to `path`, then fsyncs. The file is never
// rewritten: everything already in it stays byte-for-byte, so a merge, a
// solve or a crash can never lose what an earlier write put there.
//
// Creates the file (with the header) if absent. Refuses a file with the old
// 12-column header -- the 12->13 migration was one explicit run of the
// retired verifyslicegenus --rewrite-witnesses, never something an append
// does implicitly -- and truncates a torn last line (no newline) before appending.
//
// On success each appended witness gets its fileOffset and drops its pair
// signature from memory. Throws on any failure, leaving those witnesses
// untouched in memory for the next attempt.
void appendCobordisms(const std::filesystem::path &path,
                     std::vector<cobordisms::Cobordism> &cobordisms,
                     size_t from) {
  if (from >= cobordisms.size())
    return;
  const std::string what = path.string();
  int fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
  if (fd < 0)
    throw std::runtime_error("cannot open " + what + ": " +
                             std::strerror(errno));
  struct FdCloser {
    int fd;
    ~FdCloser() { ::close(fd); }
  } closer{fd};

  off_t size = ::lseek(fd, 0, SEEK_END);
  if (size < 0)
    throw std::runtime_error("cannot seek " + what);
  if (size == 0) {
    writeAll(fd, std::string(COBORDISMS_HEADER) + "\n", what);
  } else {
    const std::string header = readHeader(fd);
    if (header != COBORDISMS_HEADER)
      throw std::runtime_error(
          what + " does not have the current witness-file header" +
          (header == COBORDISMS_HEADER_12
               ? std::string(" (it is the 12-column one: migrate it first, as "
                             "the retired verifyslicegenus --rewrite-witnesses did)")
               : std::string()) +
          "; refusing to append to it");
    appendonly::cutTornLine(fd, what);
  }

  off_t pos = ::lseek(fd, 0, SEEK_END);
  std::string buffer;
  std::vector<long long> offsets;
  offsets.reserve(cobordisms.size() - from);
  for (size_t i = from; i < cobordisms.size(); ++i) {
    offsets.push_back(static_cast<long long>(pos) +
                      static_cast<long long>(buffer.size()));
    buffer += formatCobordism(cobordisms[i]);
    buffer += '\n';
  }
  writeAll(fd, buffer, what);
  appendonly::sync(fd, what);

  for (size_t i = from; i < cobordisms.size(); ++i) {
    cobordisms::Cobordism &w = cobordisms[i];
    w.fileOffset = offsets[i - from];
    if (w.pairSigKey.empty() && !w.pairSig.empty())
      w.pairSigKey = cobordisms::cobordismKey(w.pairSig);
    std::string().swap(w.pairSig);
  }
}

} // namespace cobordisms
