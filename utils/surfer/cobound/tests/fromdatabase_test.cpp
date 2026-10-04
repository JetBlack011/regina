// readbackcache_test.cpp
//
// Master cobordisms' read-backs kept across runs (outgoing/fromdatabase.h), on the
// real cobordisms of 10_3's Phase 0 search:
//
//   1. serialiseLink()/parseLink() round-trip every read-back exactly, and
//      malformed text is refused.
//   2. A cache written by one ReadBacks is read back whole by the next:
//      every link equal to a fresh outgoingLinkFast(), failures kept too.
//   3. A file with another build digest is ignored and replaced (its edge
//      numbers belong to another thickening).
//   4. A torn last line (a run killed mid-append) is skipped, then cut
//      before the next append, and the file stays readable.
//   5. An empty directory means no cache: nothing found, nothing written.

#include <filesystem>
#include <fstream>
#include <string>
#include <vector>

#include <unistd.h>

#include "cobound/cobordisms/cobordismkey.h"
#include "cobound/outgoing/fromdatabase.h"
#include "linknaming/tests/check.h"

using namespace outgoing;
namespace fs = std::filesystem;

namespace {

const char *PD_10_3 =
    "[[2;10;3;9];[4;17;5;18];[6;15;7;16];[8;4;9;3];[10;2;11;1];[12;20;13;19];"
    "[14;7;15;8];[16;5;17;6];[18;14;19;13];[20;12;1;11]]";

std::vector<std::string> pairsigs(const std::string &path) {
  std::ifstream in(path);
  std::string line;
  std::getline(in, line); // header
  std::vector<std::string> out;
  while (std::getline(in, line)) {
    std::vector<std::string> f;
    std::string cur;
    bool q = false;
    for (char c : line) {
      if (c == '"') q = !q;
      else if (c == ',' && !q) { f.push_back(cur); cur.clear(); }
      else cur += c;
    }
    f.push_back(cur);
    out.push_back(f[8]);
  }
  return out;
}

bool sameLink(const outgoing::OutgoingLink &a, const outgoing::OutgoingLink &b) {
  if (a.curves.size() != b.curves.size()) return false;
  for (size_t c = 0; c < a.curves.size(); ++c) {
    if (a.curves[c].size() != b.curves[c].size()) return false;
    for (size_t i = 0; i < a.curves[c].size(); ++i)
      if (a.curves[c][i].edge != b.curves[c][i].edge ||
          a.curves[c][i].reversed != b.curves[c][i].reversed)
        return false;
  }
  return a.surfaceComponent == b.surfaceComponent &&
         a.incomingFirstEdge == b.incomingFirstEdge &&
         a.incomingSurfaceComponent == b.incomingSurfaceComponent;
}

std::string tempDir() {
  std::string d = (fs::temp_directory_path() / ("readbackcache_test_" + std::to_string(::getpid())))
                      .string();
  fs::remove_all(d);
  return d;
}

} // namespace

int main() {
  const std::vector<std::string> sigs =
      pairsigs(std::string(COBOUND_TEST_DATA) + "/search_10_3_cobordisms.csv");
  CHECK_EQ(static_cast<int>(sigs.size()), 9, "nine 10_3 witnesses");
  outgoing::OutgoingReader redraw(PD_10_3, 2);
  const std::string digest = redraw.buildChecksum();
  std::vector<std::optional<outgoing::OutgoingLink>> fresh;
  for (const std::string &s : sigs) {
    std::string why;
    fresh.push_back(redraw.outgoingLinkFast(s, why));
    CHECK(fresh.back().has_value(), "a 10_3 witness reads back: " + why);
  }

  // 1. Round trip, and refusals.
  for (const auto &l : fresh) {
    if (!l) continue;
    auto back = parseLink(serialiseLink(*l));
    CHECK(back && sameLink(*back, *l), "serialise/parse round-trips a read-back");
  }
  CHECK(!parseLink("1+,2x|0|3|0").has_value(), "a bad direction is refused");
  CHECK(!parseLink("1+,2-|0,1|3|0").has_value(), "a component list of the wrong length is refused");
  CHECK(!parseLink("1+|0|3").has_value(), "three fields are refused");

  // 2. Written by one, read whole by the next; a failure is kept too.
  const std::string dir = tempDir();
  {
    ReadBacks c(dir, PD_10_3, 2, digest);
    CHECK_EQ(static_cast<int>(c.loaded()), 0, "a new cache is empty");
    for (size_t i = 0; i < sigs.size(); ++i)
      c.put(cobordisms::cobordismKey(sigs[i]), {fresh[i], ""});
    c.put("feedfacecafe", {std::nullopt, "no isomorphism carries its incoming curve"});
    c.flush();
  }
  {
    ReadBacks c(dir, PD_10_3, 2, digest);
    CHECK_EQ(static_cast<int>(c.loaded()), 10, "every entry comes back");
    bool all = true;
    for (size_t i = 0; i < sigs.size(); ++i) {
      const CachedReadBack *r = c.get(cobordisms::cobordismKey(sigs[i]));
      all = all && r && r->link && fresh[i] && sameLink(*r->link, *fresh[i]);
    }
    CHECK(all, "every cached read-back equals a fresh one");
    const CachedReadBack *f = c.get("feedfacecafe");
    CHECK(f && !f->link && f->why == "no isomorphism carries its incoming curve",
          "a failure is kept with its reason");
    CHECK_EQ(static_cast<int>(c.hits()), 10, "hits counted");
    // The same diagram at other layers is another file.
    ReadBacks other(dir, PD_10_3, 3, digest);
    CHECK_EQ(static_cast<int>(other.loaded()), 0, "layers are part of the row");
  }

  // 3. Another build digest: ignored, then replaced.
  {
    ReadBacks c(dir, PD_10_3, 2, "not-this-build");
    CHECK_EQ(static_cast<int>(c.loaded()), 0, "a stale digest's entries are not used");
    c.put(cobordisms::cobordismKey(sigs[0]), {fresh[0], ""});
    c.flush();
    ReadBacks again(dir, PD_10_3, 2, "not-this-build");
    CHECK_EQ(static_cast<int>(again.loaded()), 1, "the file was started afresh for that build");
    ReadBacks orig(dir, PD_10_3, 2, digest);
    CHECK_EQ(static_cast<int>(orig.loaded()), 0, "and the old build's entries are gone");
  }

  // 4. A torn last line.
  {
    fs::remove_all(dir);
    {
      ReadBacks c(dir, PD_10_3, 2, digest);
      c.put(cobordisms::cobordismKey(sigs[0]), {fresh[0], ""});
      c.flush();
    }
    std::string path;
    for (const auto &e : fs::directory_iterator(dir))
      if (e.path().extension() == ".readback") path = e.path().string();
    {
      std::ofstream f(path, std::ios::app | std::ios::binary);
      f << cobordisms::cobordismKey(sigs[1]) << "\tok\t" << serialiseLink(*fresh[1]).substr(0, 5);
    }
    ReadBacks c(dir, PD_10_3, 2, digest);
    CHECK_EQ(static_cast<int>(c.loaded()), 1, "the torn line is skipped");
    c.put(cobordisms::cobordismKey(sigs[2]), {fresh[2], ""});
    c.flush();
    ReadBacks after(dir, PD_10_3, 2, digest);
    CHECK_EQ(static_cast<int>(after.loaded()), 2, "the torn line was cut, the new one kept");
    const CachedReadBack *r = after.get(cobordisms::cobordismKey(sigs[2]));
    CHECK(r && r->link && sameLink(*r->link, *fresh[2]), "the appended read-back is intact");
  }

  // 5. No cache.
  {
    ReadBacks none("", PD_10_3, 2, digest);
    none.put(cobordisms::cobordismKey(sigs[0]), {fresh[0], ""});
    none.flush();
    CHECK(none.get(cobordisms::cobordismKey(sigs[0])) == nullptr, "an empty dir keeps nothing");
  }
  fs::remove_all(dir);
  return checks::finish("readbackcache_test");
}
