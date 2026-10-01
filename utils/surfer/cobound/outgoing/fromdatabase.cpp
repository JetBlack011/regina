//
//  fromdatabase.cpp
//
//  WitnessRedrawer (was farsideredraw.cpp) and RowReadBacks (was
//  readbackcache.cpp).
//

#include "cobound/outgoing/fromdatabase.h"

#include <algorithm>
#include <cerrno>
#include <chrono>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <map>
#include <sstream>
#include <stdexcept>
#include <sys/file.h>
#include <unistd.h>

#include <triangulation/dim2.h>
#include <utilities/sigutils.h>

#include "cobound/cobordisms/appendonly.h"
#include "cobound/cobordisms/cobordismkey.h"
#include "surfer/pairsig/pairsig.h"
#include "surfer/pairsig/sha1.h"
#include "surfer/report/atomicwrite.h"
#include "surfer/submanifold/submanifold.h"
#include "surfer/submanifold/vertexlinks.h"

namespace farside {

namespace {

// Face index of `t` (a triangle of `from`) under `iso`, in `to`.
int carryTriangle(const regina::Triangle<4> *t, const regina::Isomorphism<4> &iso,
                  const regina::Triangulation<4> &to) {
    auto emb = t->front();
    size_t p = emb.pentachoron()->index();
    regina::Perm<5> v = iso.facetPerm(p) * emb.vertices();
    return static_cast<int>(
        to.pentachoron(iso.simpImage(p))->triangle(regina::Face<4, 2>::faceNumber(v))->index());
}

} // namespace

std::optional<std::vector<int>> WitnessRedrawer::pinned_(const regina::Triangulation<4> &ambient,
                                                         const std::vector<int> &faces,
                                                         const IsoSource &isos) const {
    const regina::Triangulation<4> &W = thickening();
    std::optional<std::vector<int>> carried;
    isos([&](const regina::Isomorphism<4> &iso) {
        std::vector<int> image;
        image.reserve(faces.size());
        for (int f : faces) image.push_back(carryTriangle(ambient.triangle(f), iso, W));
        if (boundaryEdgesOf(W, image, rb_.searchSideBC) != rb_.searchEdges) return false;
        carried = std::move(image);
        return true;
    });
    return carried;
}

WitnessRedrawer::WitnessRedrawer(const std::string &rowPD, int layers) {
    if (layers < 1) throw regina::InvalidArgument("WitnessRedrawer: layers must be >= 1");
    // The row's thickening, exactly as verifyslicegenus builds it: collared
    // through every layer, no cone.
    rowsearch::buildRow(rowPD, layers, layers, /*useCone=*/false, rb_);
    const regina::Triangulation<4> &W = rb_.tri;
    outgoing_ = std::make_unique<OutgoingMap>(rb_.link.tri, *rb_.cob);
    drawer_ = std::make_unique<knotbuilder::DiagramDrawer>(rb_.link.tri, rb_.pdcode.size());
    skeleton_ = std::make_unique<Skeleton<4, 2>>(W);
    rowCycles_ = knotbuilder::DiagramDrawer::cyclesOf(rb_.link.edges, rb_.link.reversed);
    {
        // A search-side edge -> its row edge (the row map) -> its T edge ->
        // the cycle holding it.
        std::unordered_map<size_t, size_t> cycleOfT;
        for (size_t c = 0; c < rowCycles_.size(); ++c)
            for (const auto &de : rowCycles_[c]) cycleOfT[de.edge] = c;
        for (const auto &[edge, rowEdge] : rb_.orientation->rowIndexOf)
            rowComponentOf_[edge] = cycleOfT.at(rb_.link.edges[rowEdge]->index());
    }
    for (size_t c = 0; c < W.countBoundaryComponents(); ++c) {
        const regina::BoundaryComponent<4> *bc = W.boundaryComponent(c);
        // Same-indexed edges of bc and its build() are numbered alike
        // (BoundaryComponent<4>::build(); KnottedSurface relies on it too).
        boundaries_.push_back(bc->build());
        for (size_t k = 0; k < bc->countEdges(); ++k) boundaryEdge_[bc->edge(k)] = {c, k};
    }
}

std::optional<std::vector<int>> WitnessRedrawer::carry(const std::string &pairsig,
                                                       std::string &why) const {
    const regina::Triangulation<4> &W = thickening();
    auto t0 = std::chrono::steady_clock::now();
    DecodedKnottedSurfaceSig dec = fromKnottedSurfaceSig(pairsig);
    const std::vector<int> faces = dec.surface->markedFaces();
    auto t1 = std::chrono::steady_clock::now();
    msDecode_ += std::chrono::duration<double, std::milli>(t1 - t0).count();

    // The isomorphism onto W sending the witness's incoming curve onto L x {0},
    // found while the ambient's isomorphisms are enumerated.
    std::optional<std::vector<int>> carried =
        pinned_(*dec.ambient, faces, [&](const IsoVisitor &visit) {
            dec.ambient->findAllIsomorphisms(W, visit);
        });
    msIso_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    if (carried) return carried;

    // Say why: no isomorphism at all (another thickening), or isomorphisms
    // whose incoming curve is not L x {0}.
    size_t isos = 0;
    std::string sizes;
    dec.ambient->findAllIsomorphisms(W, [&](const regina::Isomorphism<4> &iso) {
        std::vector<int> image;
        for (int f : faces) image.push_back(carryTriangle(dec.ambient->triangle(f), iso, W));
        if (isos < 4) {
            auto in = boundaryEdgesOf(W, image, rb_.searchSideBC);
            auto out = boundaryEdgesOf(W, image, outgoing_->boundaryComponent());
            std::vector<size_t> common;
            std::ranges::set_intersection(in, rb_.searchEdges, std::back_inserter(common));
            sizes += " [in " + std::to_string(in.size()) + " edges, " +
                     std::to_string(common.size()) + " on L; out " + std::to_string(out.size()) + "]";
        }
        ++isos;
        return false;
    });
    why = "no isomorphism carries its incoming curve onto L (" + std::to_string(isos) +
          " isomorphisms onto the thickening; L has " + std::to_string(rb_.searchEdges.size()) +
          " edges;" + sizes + ")";
    return std::nullopt;
}

namespace {

// Directed boundary edges chained head to tail into closed curves (as
// KnottedSurface::orientedBoundaryLinks() chains them), refusing a dead end.
std::vector<OrientedCurve> chain(const std::vector<OrientedEdge> &directed) {
    auto curves = chainIntoCurves(directed, edgecycles::OpenChain::refuse);
    if (!curves)
        throw regina::InvalidArgument("the surface's boundary is not closed curves");
    return std::move(*curves);
}

} // namespace

std::optional<OutgoingLink> WitnessRedrawer::outgoingLinkFast(const std::string &pairsig,
                                                              std::string &why) const {
    const regina::Triangulation<4> &W = thickening();
    auto t0 = std::chrono::steady_clock::now();
    const size_t pos = pairsig.find(regina::Base64Encoder::spare[0]); // pairsig.cpp's delimiter
    if (pos == std::string::npos) throw regina::InvalidArgument("not a pair signature");
    const std::string ambientSig = pairsig.substr(0, pos);
    if (ambientSig != ambientSig_) {
        ambient_ = std::make_unique<regina::Triangulation<4>>(
            regina::Triangulation<4>::fromSig(ambientSig));
        const size_t m = ambient_->countTriangles();
        faceWidth_ = regina::Base64Encoder::integerWidth(m == 0 ? 0 : m - 1);
        auto ti = std::chrono::steady_clock::now();
        isos_.clear();
        ambient_->findAllIsomorphisms(W, [&](const regina::Isomorphism<4> &iso) {
            isos_.push_back(iso);
            return false;
        });
        msIso_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - ti).count();
        ambientSig_ = ambientSig;
    }
    std::vector<int> faces;
    {
        const std::string suffix = pairsig.substr(pos + 1);
        if (suffix.size() % static_cast<size_t>(faceWidth_) != 0)
            throw regina::InvalidArgument("malformed pair signature suffix");
        regina::Base64Decoder decoder(suffix.begin(), suffix.end());
        for (size_t i = 0; i < suffix.size() / static_cast<size_t>(faceWidth_); ++i)
            faces.push_back(decoder.template decodeInt<int>(faceWidth_));
    }
    auto t1 = std::chrono::steady_clock::now();
    msDecode_ += std::chrono::duration<double, std::milli>(t1 - t0).count();

    // The isomorphism carrying the incoming curve onto L x {0}, from the
    // ambient's isomorphisms found once per row.
    const std::optional<std::vector<int>> pinned =
        pinned_(*ambient_, faces, [&](const IsoVisitor &visit) {
            for (const regina::Isomorphism<4> &iso : isos_)
                if (visit(iso)) break;
        });
    auto t2 = std::chrono::steady_clock::now();
    msIso_ += std::chrono::duration<double, std::milli>(t2 - t1).count();
    if (!pinned) {
        why = "no isomorphism carries its incoming curve onto L (" +
              std::to_string(isos_.size()) + " isomorphisms onto the thickening)";
        return std::nullopt;
    }
    const std::vector<int> &carried = *pinned;

    // The surface as a plain 2-triangulation: one triangle per face, glued
    // along the edges of W two faces share (an embedded surface has no edge
    // in three).
    // Keyed by edge index, so the gluings are made in an order that depends
    // on the surface alone, never on addresses.
    regina::Triangulation<2> surf;
    std::map<size_t, std::vector<std::pair<size_t, int>>> byEdge;
    for (size_t k = 0; k < carried.size(); ++k) {
        surf.newSimplex();
        const regina::Triangle<4> *t = W.triangle(carried[k]);
        for (int i = 0; i < 3; ++i) byEdge[t->edge(i)->index()].push_back({k, i});
    }
    for (const auto &[edge, uses] : byEdge) {
        if (uses.size() == 1) continue;
        if (uses.size() != 2) throw regina::InvalidArgument("the witness is not a surface");
        const auto [ka, ia] = uses[0];
        const auto [kb, ib] = uses[1];
        const regina::Perm<5> p = W.triangle(carried[ka])->edgeMapping(ia);
        const regina::Perm<5> q = W.triangle(carried[kb])->edgeMapping(ib);
        int image[3];
        image[ia] = ib;
        image[p[0]] = q[0];
        image[p[1]] = q[1];
        surf.simplex(ka)->join(ia, surf.simplex(kb), regina::Perm<3>(image[0], image[1], image[2]));
    }

    // Oriented boundary, per surface component, exactly as
    // KnottedSurface::orientedBoundaryLinks()/boundaryEdgeSurfaceComponent().
    std::vector<std::vector<OrientedEdge>> directed(boundaries_.size());
    std::map<const regina::Edge<3> *, size_t> surfaceOf;
    for (size_t k = 0; k < carried.size(); ++k) {
        const regina::Simplex<2> *simplex = surf.simplex(k);
        const int sign = simplex->orientation();
        const regina::Triangle<4> *t = W.triangle(carried[k]);
        for (int i = 0; i < 3; ++i) {
            if (simplex->adjacentSimplex(i) != nullptr) continue;
            auto be = boundaryEdge_.find(t->edge(i));
            if (be == boundaryEdge_.end()) continue;
            const auto [c, local] = be->second;
            const int headLocal = sign > 0 ? (i + 2) % 3 : (i + 1) % 3;
            const regina::Perm<5> p = t->edgeMapping(i);
            const regina::Edge<3> *edge = boundaries_[c].edge(local);
            directed[c].push_back({edge, p[0] == headLocal});
            surfaceOf[edge] = simplex->component()->index();
        }
    }
    // Then oriented against the row, as for a rebuilt surface
    // (orientedOutgoingLink()).
    const std::vector<std::pair<size_t, std::vector<OrientedCurve>>> oriented = {
        {rb_.searchSideBC, chain(directed[rb_.searchSideBC])},
        {outgoing_->boundaryComponent(), chain(directed[outgoing_->boundaryComponent()])}};
    std::optional<OutgoingLink> out = orientedOutgoingLink(
        oriented, surfaceOf, *outgoing_, *rb_.orientation, rb_.searchSideBC, &why);
    msRead_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t2).count();
    return out;
}

std::optional<OutgoingLink> WitnessRedrawer::outgoingLink(const std::string &pairsig,
                                                          std::string &why) const {
    std::optional<std::vector<int>> carried = carry(pairsig, why);
    if (!carried) return std::nullopt;
    auto t0 = std::chrono::steady_clock::now();
    KnottedSurface surface(*skeleton_);
    auto tb = std::chrono::steady_clock::now();
    if (!surface.addFaces(*carried)) {
        why = "the witness's faces do not embed";
        return std::nullopt;
    }
    auto t1 = std::chrono::steady_clock::now();
    msBoundaryBuild_ += std::chrono::duration<double, std::milli>(tb - t0).count();
    auto link = orientedOutgoingLink(surface, *outgoing_, *rb_.orientation, rb_.searchSideBC);
    msSurface_ += std::chrono::duration<double, std::milli>(t1 - t0).count();
    msRead_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    if (!link) why = "incoming orientation is inconsistent";
    return link;
}

bool WitnessRedrawer::rebuild(const std::vector<int> &faces, KnottedSurface &surface,
                              std::string &why) const {
    const regina::Triangulation<4> &W = rb_.tri;
    for (int f : faces)
        if (f < 0 || static_cast<size_t>(f) >= W.countTriangles()) {
            why = "face " + std::to_string(f) + " is not a triangle of the thickening";
            return false;
        }
    // Before anything is built: the search side must be exactly L x {0}.
    if (boundaryEdgesOf(W, faces, rb_.searchSideBC) != rb_.searchEdges) {
        why = "its incoming boundary is not L x {0}";
        return false;
    }
    if (!surface.addFaces(faces)) {
        why = "its faces fail the search's embedding checks";
        return false;
    }
    if (!surface.satisfies(BoundaryCondition::proper)) {
        why = "it is not properly embedded";
        return false;
    }
    if (!surface.isAcceptable()) {
        why = "it is not acceptable (a self-intersection that does not resolve, "
              "or not smooth at the boundary)";
        return false;
    }
    return true;
}

std::optional<OutgoingLink> WitnessRedrawer::outgoingLinkFromFaces(const std::vector<int> &faces,
                                                                   std::string &why) const {
    PetalCache cache;
    KnottedSurface::SelfIntersectionOptions options;
    options.resolveUnlinked = true;
    KnottedSurface surface(options, *skeleton_, cache);
    if (!rebuild(faces, surface, why)) return std::nullopt;
    auto link = orientedOutgoingLink(surface, *outgoing_, *rb_.orientation, rb_.searchSideBC);
    if (!link) why = "incoming orientation is inconsistent";
    return link;
}

std::string WitnessRedrawer::buildChecksum() const {
    const regina::Triangulation<4> &W = rb_.tri;
    std::string data = std::to_string(W.size()) + '|';
    for (size_t i = 0; i < W.size(); ++i) {
        const regina::Pentachoron<4> *p = W.pentachoron(i);
        for (int f = 0; f < 5; ++f) {
            if (const regina::Pentachoron<4> *adj = p->adjacentPentachoron(f))
                data += std::to_string(adj->index()) + ':' +
                        std::to_string(p->adjacentGluing(f).S5Index());
            else
                data += '-';
            data += ',';
        }
    }
    data += '|';
    for (size_t k = 0; k < W.countTriangles(); ++k) {
        const auto &emb = W.triangle(k)->front();
        data += std::to_string(emb.simplex()->index()) + ':' + std::to_string(emb.face()) + ',';
    }
    return pairsig::sha1Hex(data).substr(0, 16);
}

} // namespace farside

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
