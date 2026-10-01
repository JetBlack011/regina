// readbackcache.cpp

#include "cobound/outgoing/fromdatabase.h"

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <sys/file.h>
#include <unistd.h>

#include "cobound/cobordisms/appendonly.h"
#include "cobound/cobordisms/cobordismkey.h"
#include "surfer/report/atomicwrite.h"

namespace fs = std::filesystem;

namespace cascade {

namespace {

void joinSizes(std::ostringstream &o, const std::vector<size_t> &v) {
  for (size_t i = 0; i < v.size(); ++i) o << (i ? "," : "") << v[i];
}

bool splitSizes(const std::string &s, std::vector<size_t> &out) {
  out.clear();
  if (s.empty()) return true;
  std::istringstream in(s);
  for (std::string t; std::getline(in, t, ',');) {
    if (t.empty()) return false;
    try {
      size_t used = 0;
      out.push_back(std::stoul(t, &used));
      if (used != t.size()) return false;
    } catch (const std::exception &) {
      return false;
    }
  }
  return true;
}

} // namespace

// "e+,e-,...;e+,...|sc,sc|firstEdge,...|surfaceComponent,..."
std::string serialiseLink(const farside::OutgoingLink &link) {
  std::ostringstream o;
  for (size_t c = 0; c < link.curves.size(); ++c) {
    o << (c ? ";" : "");
    for (size_t i = 0; i < link.curves[c].size(); ++i)
      o << (i ? "," : "") << link.curves[c][i].edge << (link.curves[c][i].reversed ? '-' : '+');
  }
  o << '|';
  joinSizes(o, link.surfaceComponent);
  o << '|';
  joinSizes(o, link.incomingFirstEdge);
  o << '|';
  joinSizes(o, link.incomingSurfaceComponent);
  return o.str();
}

std::optional<farside::OutgoingLink> parseLink(const std::string &text) {
  std::vector<std::string> parts;
  {
    std::istringstream in(text);
    for (std::string p; std::getline(in, p, '|');) parts.push_back(p);
    if (!text.empty() && text.back() == '|') parts.emplace_back();
  }
  if (parts.size() != 4) return std::nullopt;
  farside::OutgoingLink link;
  if (!parts[0].empty()) {
    std::istringstream curves(parts[0]);
    for (std::string c; std::getline(curves, c, ';');) {
      knotbuilder::EdgeCycle cycle;
      std::istringstream edges(c);
      for (std::string e; std::getline(edges, e, ',');) {
        if (e.size() < 2 || (e.back() != '+' && e.back() != '-')) return std::nullopt;
        try {
          size_t used = 0;
          const size_t edge = std::stoul(e.substr(0, e.size() - 1), &used);
          if (used != e.size() - 1) return std::nullopt;
          cycle.push_back({edge, e.back() == '-'});
        } catch (const std::exception &) {
          return std::nullopt;
        }
      }
      link.curves.push_back(std::move(cycle));
    }
  }
  if (!splitSizes(parts[1], link.surfaceComponent) ||
      !splitSizes(parts[2], link.incomingFirstEdge) ||
      !splitSizes(parts[3], link.incomingSurfaceComponent))
    return std::nullopt;
  if (link.surfaceComponent.size() != link.curves.size() ||
      link.incomingFirstEdge.size() != link.incomingSurfaceComponent.size())
    return std::nullopt;
  return link;
}

RowReadBacks::RowReadBacks(const std::string &dir, const std::string &rowPD, int layers,
                           const std::string &buildDigest)
    : digest_(buildDigest) {
  if (dir.empty()) return;
  fs::create_directories(dir);
  path_ = dir + "/" + witnesskey::witnessKey(rowPD + "|" + std::to_string(layers)) + ".readback";
  std::ifstream in(path_, std::ios::binary);
  std::string line;
  if (!in || !std::getline(in, line) || in.eof()) {
    rewrite_ = true; // absent, empty or a torn header
    return;
  }
  if (line != "readback 1 " + digest_) {
    rewrite_ = true; // another build of the row: its edge numbers mean nothing here
    return;
  }
  while (std::getline(in, line)) {
    if (in.eof()) break; // no newline: a torn last line
    const size_t a = line.find('\t'), b = a == std::string::npos ? a : line.find('\t', a + 1);
    if (b == std::string::npos) continue;
    const std::string key = line.substr(0, a), kind = line.substr(a + 1, b - a - 1),
                      rest = line.substr(b + 1);
    CachedReadBack r;
    if (kind == "ok") {
      r.link = parseLink(rest);
      if (!r.link) continue; // unreadable: recomputed
    } else if (kind == "fail") {
      r.why = rest;
    } else {
      continue;
    }
    entries_[key] = std::move(r);
  }
}

const CachedReadBack *RowReadBacks::get(const std::string &witnessKey) const {
  auto it = entries_.find(witnessKey);
  if (it == entries_.end()) return nullptr;
  ++hits_;
  return &it->second;
}

void RowReadBacks::put(const std::string &witnessKey, const CachedReadBack &r) {
  if (path_.empty()) return;
  if (!entries_.emplace(witnessKey, r).second) return;
  std::string why = r.why;
  for (char &ch : why)
    if (ch == '\n' || ch == '\t') ch = ' ';
  pending_ += witnessKey + (r.link ? "\tok\t" + serialiseLink(*r.link) : "\tfail\t" + why) + '\n';
}

void RowReadBacks::flush() {
  if (path_.empty() || (pending_.empty() && !rewrite_)) return;
  appendonly::FileLock lock(path_ + ".lock");
  std::string out;
  if (rewrite_) {
    // Start the file afresh with this build's digest and everything known now.
    out = "readback 1 " + digest_ + "\n";
    for (const auto &[key, r] : entries_) {
      std::string why = r.why;
      for (char &ch : why)
        if (ch == '\n' || ch == '\t') ch = ' ';
      out += key + (r.link ? "\tok\t" + serialiseLink(*r.link) : "\tfail\t" + why) + '\n';
    }
    report::atomicWrite(path_, [&](std::ostream &f) { f << out; });
    rewrite_ = false;
  } else {
    // Another run may have appended meanwhile; a duplicate key is harmless
    // (the loader keeps one). A torn last line from a killed run is cut first.
    // A cache: not fsynced.
    appendonly::append(path_, pending_, appendonly::Sync::no);
  }
  pending_.clear();
}

} // namespace cascade
