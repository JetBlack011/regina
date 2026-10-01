//
//  witnessstore.cpp
//
//  The append-only witness file (cobordisms.csv) and the literature tables'
//  CSV rows, shared by verifyslicegenus and cascadesearch. Moved verbatim
//  from verifyslicegenus.cpp (2026-09-28).
//

#include "cobound/cobordisms/witnessstore.h"

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
#include "cobound/cobordisms/witnesskey.h"

namespace witnessstore {
// One row of the input knot table (Name,PD Notation,Genus-4D). No RFC-4180
// quoting appears in that file (PD Notation uses ';' internally, never a
// literal comma), so a naive two-comma split suffices -- see the input
// table's own format, confirmed during design.
bool splitInputLine(const std::string &line, std::string &name,
                    std::string &pd, std::string &genusField) {
  size_t c1 = line.find(',');
  if (c1 == std::string::npos)
    return false;
  size_t c2 = line.find(',', c1 + 1);
  if (c2 == std::string::npos)
    return false;
  name = line.substr(0, c1);
  pd = line.substr(c1 + 1, c2 - c1 - 1);
  genusField = line.substr(c2 + 1);
  if (!genusField.empty() && genusField.back() == '\r')
    genusField.pop_back();
  return true;
}

// Parses "N" or "[lo;hi]" into lo/hi (lo == hi in the plain-integer case).
void parseGenusField(const std::string &field, int &lo, int &hi) {
  if (!field.empty() && field.front() == '[') {
    size_t semi = field.find(';');
    lo = std::stoi(field.substr(1, semi - 1));
    hi = std::stoi(field.substr(semi + 1, field.size() - semi - 2));
  } else {
    lo = hi = std::stoi(field);
  }
}

// Loads a literature table for its names and bounds only, skipping the PD
// code entirely. Used for tables that aren't this run's --input: we need
// their names (to expand orientation-blind identifications into candidate
// sets) and their bounds, but never build a triangulation from them, so
// there is no reason to pay parsePDCode()'s cost across 12k+ rows.
size_t loadNameTable(const std::filesystem::path &path,
                     cobordismgraph::NameTable &names) {
  std::ifstream in(path);
  if (!in)
    throw std::runtime_error("Cannot open name table: " + path.string());

  size_t loaded = 0;
  std::string line;
  std::getline(in, line); // header
  while (std::getline(in, line)) {
    if (line.empty())
      continue;
    std::string name, pd, genusField;
    if (!splitInputLine(line, name, pd, genusField))
      continue;
    int lo = 0, hi = 0;
    try {
      parseGenusField(genusField, lo, hi);
    } catch (const std::exception &) {
      continue;
    }
    names.addLiterature(name, lo, hi);
    ++loaded;
  }
  return loaded;
}

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
std::string formatWitness(const cobordismgraph::Witness &w) {
  std::ostringstream candidates;
  for (size_t i = 0; i < w.otherCandidates.size(); ++i) {
    if (i)
      candidates << ';';
    candidates << w.otherCandidates[i];
  }
  std::ostringstream out;
  out << (w.kind == cobordismgraph::WitnessKind::direct ? "direct"
                                                        : "cobordism")
      << ',' << csvField(w.subject) << ',' << w.subjectComponents << ','
      << csvField(w.other) << ',' << csvField(candidates.str()) << ','
      << w.otherComponents << ',' << w.genus << ','
      << (w.tubed ? "true" : "false") << ',' << csvField(w.pairSig) << ','
      << csvField(w.sourceRow) << ',' << w.thickenLayers << ',' << w.maxFaces
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
bool witnessFromFields(std::vector<std::string> f, cobordismgraph::Witness &w,
                       bool keepPairSig, bool wantPairSigKey,
                       const std::filesystem::path &path) {
  if (f.size() < 12)
    return false;
  w.kind = f[0] == "direct" ? cobordismgraph::WitnessKind::direct
                            : cobordismgraph::WitnessKind::cobordism;
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
  if (!f[4].empty()) {
    std::istringstream candidates(f[4]);
    std::string one;
    while (std::getline(candidates, one, ';'))
      if (!one.empty())
        w.otherCandidates.push_back(one);
  }
  w.tubed = f[7] == "true";
  if (wantPairSigKey && !f[8].empty())
    w.pairSigKey = witnesskey::witnessKey(f[8]);
  if (keepPairSig)
    w.pairSig = std::move(f[8]);
  w.sourceRow = f[9];
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

namespace {
// Parses one witness line (12 or 13 fields) into `w`, keeping the pair
// signature only if `keepPairSig`. Returns false for a malformed line.
bool parseWitnessLine(const std::string &line, cobordismgraph::Witness &w,
                      bool keepPairSig, bool wantPairSigKey,
                      const std::filesystem::path &path) {
  return witnessFromFields(parseCsvLine(line), w, keepPairSig, wantPairSigKey,
                           path);
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
std::vector<cobordismgraph::Witness>
loadWitnesses(const std::filesystem::path &path, bool wantPairSigKeys) {
  std::vector<cobordismgraph::Witness> result;
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
    cobordismgraph::Witness w;
    if (!parseWitnessLine(line, w, /*keepPairSig=*/false, wantPairSigKeys,
                          path)) {
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

std::unordered_set<std::string>
witnessIdentities(const std::filesystem::path &path, unsigned threads,
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
      cobordismgraph::Witness w;
      if (!parseWitnessLine(l, w, /*keepPairSig=*/false, /*wantPairSigKey=*/false, path)) {
        ++malformed[k];
        continue;
      }
      sets[k].insert(cobordismgraph::witnessIdentity(w));
    }
  };
  std::vector<std::thread> pool;
  for (size_t k = 1; k < n; ++k)
    pool.emplace_back(work, k);
  work(0);
  for (auto &t : pool)
    t.join();
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
// write(2) until done, or throw.
void writeAll(int fd, const std::string &data, const std::string &what) {
  const char *p = data.data();
  size_t left = data.size();
  while (left > 0) {
    ssize_t n = ::write(fd, p, left);
    if (n < 0) {
      if (errno == EINTR)
        continue;
      throw std::runtime_error("write to " + what + " failed: " +
                               std::strerror(errno));
    }
    p += n;
    left -= static_cast<size_t>(n);
  }
}

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
// 12-column header -- the 12->13 migration is one explicit
// --rewrite-witnesses run, never something an append does implicitly --
// and truncates a torn last line (no newline) before appending.
//
// On success each appended witness gets its fileOffset and drops its pair
// signature from memory. Throws on any failure, leaving those witnesses
// untouched in memory for the next attempt.
void appendWitnesses(const std::filesystem::path &path,
                     std::vector<cobordismgraph::Witness> &witnesses,
                     size_t from) {
  if (from >= witnesses.size())
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
               ? std::string(" (it is the 12-column one: run "
                             "--rewrite-witnesses once to migrate it)")
               : std::string()) +
          "; refusing to append to it");
    char last = 0;
    if (::pread(fd, &last, 1, size - 1) != 1)
      throw std::runtime_error("cannot read " + what);
    if (last != '\n') {
      // A torn append: find the last complete line and cut back to it.
      off_t cut = size;
      char c = 0;
      while (cut > 0 && ::pread(fd, &c, 1, cut - 1) == 1 && c != '\n')
        --cut;
      std::cerr << "[!] " << what << ": truncating a torn last line ("
                << (size - cut) << " bytes)\n";
      if (::ftruncate(fd, cut) != 0)
        throw std::runtime_error("cannot truncate " + what);
    }
  }

  off_t pos = ::lseek(fd, 0, SEEK_END);
  std::string buffer;
  std::vector<long long> offsets;
  offsets.reserve(witnesses.size() - from);
  for (size_t i = from; i < witnesses.size(); ++i) {
    offsets.push_back(static_cast<long long>(pos) +
                      static_cast<long long>(buffer.size()));
    buffer += formatWitness(witnesses[i]);
    buffer += '\n';
  }
  writeAll(fd, buffer, what);
  if (::fsync(fd) != 0)
    throw std::runtime_error("fsync of " + what + " failed: " +
                             std::strerror(errno));

  for (size_t i = from; i < witnesses.size(); ++i) {
    cobordismgraph::Witness &w = witnesses[i];
    w.fileOffset = offsets[i - from];
    if (w.pairSigKey.empty() && !w.pairSig.empty())
      w.pairSigKey = witnesskey::witnessKey(w.pairSig);
    std::string().swap(w.pairSig);
  }
}

// --rewrite-witnesses: the one full rewrite, used for the 12->13 column
// migration. Streams the file line by line (never holding it in memory),
// re-emitting each witness through formatWitness(), and checks the round
// trip as it goes: a 13-field line must come back byte-identical, a
// 12-field line must come back as itself plus the trailing ',' of an empty
// resolved_vertices. Anything else means formatWitness() is not a faithful
// inverse of the loader, and the rewrite is abandoned before it replaces
// anything. Returns the number of witnesses written.
size_t rewriteWitnessFile(const std::filesystem::path &path) {
  std::ifstream in(path, std::ios::binary);
  if (!in)
    throw std::runtime_error("cannot open " + path.string());
  std::filesystem::path tmp = path;
  tmp += ".rewrite.tmp";
  std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
  if (!out)
    throw std::runtime_error("cannot open " + tmp.string() + " for writing");

  std::string line;
  std::getline(in, line); // old header
  out << COBORDISMS_HEADER << "\n";
  size_t written = 0, migrated = 0;
  while (std::getline(in, line)) {
    if (in.eof()) {
      if (!line.empty())
        std::cerr << "[!] " << path.string()
                  << ": dropping a torn last line (" << line.size()
                  << " bytes)\n";
      break;
    }
    if (line.empty())
      continue;
    cobordismgraph::Witness w;
    if (!parseWitnessLine(line, w, /*keepPairSig=*/true,
                          /*wantPairSigKey=*/false, path))
      throw std::runtime_error("malformed witness line " +
                               std::to_string(written + 2) + " of " +
                               path.string() + "; nothing rewritten");
    std::string again = formatWitness(w);
    const size_t fields = parseCsvLine(line).size();
    const bool faithful =
        fields == 12 ? again == line + "," : again == line;
    if (!faithful)
      throw std::runtime_error(
          "line " + std::to_string(written + 2) + " of " + path.string() +
          " does not round-trip through formatWitness(); nothing rewritten");
    migrated += fields == 12;
    out << again << '\n';
    ++written;
  }
  out.flush();
  if (!out)
    throw std::runtime_error("writing " + tmp.string() + " failed");
  out.close();
  std::filesystem::rename(tmp, path);
  std::cout << "[+] --rewrite-witnesses: " << written << " witnesses, "
            << migrated << " migrated from 12 to 13 columns, every line "
            << "round-tripped\n";
  return written;
}

} // namespace witnessstore
