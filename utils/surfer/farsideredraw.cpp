//
//  farsideredraw.cpp
//

#include "farsideredraw.h"

#include <algorithm>
#include <chrono>
#include <iterator>
#include <map>

#include <triangulation/dim2.h>
#include <utilities/sigutils.h>

#include "collar.h"
#include "embeddedsubmanifold.h"
#include "pairsig.h"

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

WitnessRedrawer::WitnessRedrawer(const std::string &rowPD, int layers)
    : pd_(knotbuilder::parsePDCode(rowPD)), built_(knotbuilder::buildLink(pd_)) {
    if (layers < 1) throw regina::InvalidArgument("WitnessRedrawer: layers must be >= 1");
    // The row's thickening, exactly as verifyslicegenus builds it.
    std::vector<int> edgeIndices;
    for (const regina::Edge<3> *e : built_.edges) edgeIndices.push_back(static_cast<int>(e->index()));
    cob_ = std::make_unique<CobordismBuilder<3>>(built_.tri);
    CollarBuilder collar(edgeIndices);
    for (int i = 0; i < layers; ++i) {
        cob_->thicken();
        collar.addLayer(*cob_);
    }
    const regina::Triangulation<4> &W = cob_->getCobordism();
    std::vector<int> seed;
    for (regina::Triangle<4> *t : collar.resolve()) seed.push_back(static_cast<int>(t->index()));
    incomingBC_ = cob_->baseBoundaryComponent()->index();
    rowEdges_ = boundaryEdgesOf(W, seed, incomingBC_);
    row_ = cobordismgraph::buildRowOrientation(built_.edges, built_.reversed,
                                               W.boundaryComponent(incomingBC_)->build(), &rowEdges_);
    outgoing_ = std::make_unique<OutgoingMap>(built_.tri, *cob_);
    drawer_ = std::make_unique<knotbuilder::DiagramDrawer>(built_.tri, pd_.size());
    skeleton_ = std::make_unique<Skeleton<4, 2>>(W);
    rowCycles_ = knotbuilder::DiagramDrawer::cyclesOf(built_.edges, built_.reversed);
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

    // The isomorphism onto W sending the witness's incoming curve onto L x {0}.
    std::vector<int> carried;
    bool found = false;
    dec.ambient->findAllIsomorphisms(W, [&](const regina::Isomorphism<4> &iso) {
        std::vector<int> image;
        image.reserve(faces.size());
        for (int f : faces) image.push_back(carryTriangle(dec.ambient->triangle(f), iso, W));
        if (boundaryEdgesOf(W, image, incomingBC_) != rowEdges_) return false;
        carried = std::move(image);
        found = true;
        return true;
    });
    msIso_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    if (found) return carried;

    // Say why: no isomorphism at all (another thickening), or isomorphisms
    // whose incoming curve is not L x {0}.
    size_t isos = 0;
    std::string sizes;
    dec.ambient->findAllIsomorphisms(W, [&](const regina::Isomorphism<4> &iso) {
        std::vector<int> image;
        for (int f : faces) image.push_back(carryTriangle(dec.ambient->triangle(f), iso, W));
        if (isos < 4) {
            auto in = boundaryEdgesOf(W, image, incomingBC_);
            auto out = boundaryEdgesOf(W, image, outgoing_->boundaryComponent());
            std::vector<size_t> common;
            std::ranges::set_intersection(in, rowEdges_, std::back_inserter(common));
            sizes += " [in " + std::to_string(in.size()) + " edges, " +
                     std::to_string(common.size()) + " on L; out " + std::to_string(out.size()) + "]";
        }
        ++isos;
        return false;
    });
    why = "no isomorphism carries its incoming curve onto L (" + std::to_string(isos) +
          " isomorphisms onto the thickening; L has " + std::to_string(rowEdges_.size()) +
          " edges;" + sizes + ")";
    return std::nullopt;
}

namespace {

// Directed boundary edges chained head to tail into closed curves (as
// embeddedsubmanifold.cpp's chainIntoCurves()).
std::vector<OrientedCurve> chain(const std::vector<OrientedEdge> &directed) {
    std::unordered_map<const regina::Vertex<3> *, OrientedEdge> outFrom;
    for (const OrientedEdge &oe : directed)
        outFrom[oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0)] = oe;
    std::vector<OrientedCurve> curves;
    std::unordered_map<const regina::Edge<3> *, bool> visited;
    for (const OrientedEdge &oe : directed) {
        if (visited.contains(oe.edge)) continue;
        OrientedCurve curve;
        const regina::Vertex<3> *start = oe.reversed ? oe.edge->vertex(1) : oe.edge->vertex(0);
        const regina::Vertex<3> *curr = start;
        for (size_t step = 0; step <= directed.size(); ++step) {
            auto it = outFrom.find(curr);
            if (it == outFrom.end())
                throw regina::InvalidArgument("the surface's boundary is not closed curves");
            const OrientedEdge &next = it->second;
            visited[next.edge] = true;
            curve.push_back(next);
            curr = next.reversed ? next.edge->vertex(0) : next.edge->vertex(1);
            if (curr == start) break;
        }
        curves.push_back(std::move(curve));
    }
    return curves;
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

    // The isomorphism carrying the incoming curve onto L x {0}.
    std::vector<int> carried;
    for (const regina::Isomorphism<4> &iso : isos_) {
        std::vector<int> image;
        image.reserve(faces.size());
        for (int f : faces) image.push_back(carryTriangle(ambient_->triangle(f), iso, W));
        if (boundaryEdgesOf(W, image, incomingBC_) == rowEdges_) {
            carried = std::move(image);
            break;
        }
    }
    auto t2 = std::chrono::steady_clock::now();
    msIso_ += std::chrono::duration<double, std::milli>(t2 - t1).count();
    if (carried.empty()) {
        why = "no isomorphism carries its incoming curve onto L (" +
              std::to_string(isos_.size()) + " isomorphisms onto the thickening)";
        return std::nullopt;
    }

    // The surface as a plain 2-triangulation: one triangle per face, glued
    // along the edges of W two faces share (an embedded surface has no edge
    // in three).
    regina::Triangulation<2> surf;
    std::unordered_map<const regina::Edge<4> *, std::vector<std::pair<size_t, int>>> byEdge;
    for (size_t k = 0; k < carried.size(); ++k) {
        surf.newSimplex();
        const regina::Triangle<4> *t = W.triangle(carried[k]);
        for (int i = 0; i < 3; ++i) byEdge[t->edge(i)].push_back({k, i});
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
    std::vector<OrientedCurve> incoming = chain(directed[incomingBC_]);
    std::optional<std::map<size_t, int>> flips = incomingFlips(row_, incoming, surfaceOf);
    if (!flips) {
        why = "incoming orientation is inconsistent";
        return std::nullopt;
    }
    OutgoingLink out;
    for (const OrientedCurve &curve : incoming) {
        if (curve.empty()) continue;
        out.incomingFirstEdge.push_back(curve.front().edge->index());
        out.incomingSurfaceComponent.push_back(surfaceOf.at(curve.front().edge));
    }
    for (const OrientedCurve &curve : chain(directed[outgoing_->boundaryComponent()])) {
        if (curve.empty()) continue;
        const size_t comp = surfaceOf.at(curve.front().edge);
        auto f = flips->find(comp);
        if (f == flips->end()) {
            why = "a surface component misses the row";
            return std::nullopt;
        }
        knotbuilder::EdgeCycle cyc = outgoing_->carry(curve);
        if (f->second < 0) {
            std::reverse(cyc.begin(), cyc.end());
            for (auto &de : cyc) de.reversed = !de.reversed;
        }
        out.curves.push_back(std::move(cyc));
        out.surfaceComponent.push_back(comp);
    }
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
    auto link = orientedOutgoingLink(surface, *outgoing_, row_, incomingBC_);
    msSurface_ += std::chrono::duration<double, std::milli>(t1 - t0).count();
    msRead_ += std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t1).count();
    if (!link) why = "incoming orientation is inconsistent";
    return link;
}

} // namespace farside
