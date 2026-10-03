//
//  rowsearch_test.cpp
//
//  The row pipeline verifyslicegenus and cascadesearch share
//  (rowsearch.h), each piece pinned on its own:
//
//    1. conditionFor(): `connected` only ever for a knot.
//    2. gateReason(): the names --rejection-sample-log has always written.
//    3. RowAccounting: every bucket, every failure message, and the exact
//       `accounting:` body tools/orchestrate/dispatch.py parses.
//    4. RowWatchdog: the surface target before the clocks, each limit's
//       reason, at most one call, no thread without a limit.
//    5. gateSurface() on real exhaustive searches at face cap 3 (the
//       canaries' shape): 3_1 accepts all 1,752 surfaces, and L2a1{0}
//       accepts 795 of 945 and turns 150 away as the other orientation --
//       the numbers tools/orchestrate/canaries.expected pins through
//       verifyslicegenus, reached here without it and without any naming.
//    6. SurfaceBoundaryInfo::captureFaces: a surface's pair signature, taken
//       later from (ambient, faces), is the one the search would have
//       captured -- what lets the cascade defer signatures to the few
//       witnesses a proof uses.
//

#include <atomic>
#include <chrono>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#include "surfer/pairsig/pairsig.h"
#include "cobound/search/incoming.h"
#include "cobound/search/preconditions.h"
#include "cobound/search/search.h"
#include "surfer/enumeration/surfacesearch.h"
#include "linknaming/census/censusnaming.h"

using namespace rowsearch;

namespace {

int passed = 0;
int failedCount = 0;

void expect(bool ok, const std::string &desc) {
    if (ok) {
        ++passed;
    } else {
        ++failedCount;
        std::cout << "  FAIL: " << desc << "\n";
    }
}

template <typename A, typename E>
void expectEq(const A &actual, const E &expected, const std::string &desc) {
    if (actual == expected) {
        ++passed;
    } else {
        ++failedCount;
        std::cout << "  FAIL: " << desc << "\n        expected " << expected
                  << ", got " << actual << "\n";
    }
}

void test_condition_for() {
    using M = BoundaryConditionMode;
    expect(conditionFor(M::automatic, 1) == BoundaryCondition::connected,
           "automatic: a knot is searched under connected");
    expect(conditionFor(M::automatic, 2) == BoundaryCondition::proper,
           "automatic: a link under proper");
    expect(conditionFor(M::connected, 1) == BoundaryCondition::connected,
           "connected: honoured for a knot");
    expect(conditionFor(M::connected, 3) == BoundaryCondition::proper,
           "connected: never for a link, which cannot meet it on its own "
           "search side");
    expect(conditionFor(M::proper, 1) == BoundaryCondition::proper,
           "proper: for a knot too (the only way to find knot-to-link "
           "cobordisms)");
    expect(conditionFor(M::proper, 2) == BoundaryCondition::proper,
           "proper: for a link");
}

void test_gate_reasons() {
    expectEq(std::string(gateReason(Gate::nonOrientable)),
             std::string("non-orientable"), "non-orientable");
    expectEq(std::string(gateReason(Gate::unnamedSide)),
             std::string("unnamed-side"), "unnamed-side");
    expectEq(std::string(gateReason(Gate::searchSideElsewhere)),
             std::string("search-side-elsewhere"), "search-side-elsewhere");
    expectEq(std::string(gateReason(Gate::searchSideBroken)),
             std::string("search-side-broken"), "search-side-broken");
    expectEq(std::string(gateReason(Gate::orientation)),
             std::string("orientation"), "orientation");
    expectEq(std::string(gateReason(Gate::orientationBroken)),
             std::string("orientation-broken"), "orientation-broken");
    expectEq(std::string(gateReason(Gate::multiFarSide)),
             std::string("multi-far-side"), "multi-far-side");
}

void test_accounting() {
    {
        RowAccounting a;
        a.described = 10;
        a.recorded = 3;
        a.duplicate = 4;
        a.reject(Gate::orientation);
        a.reject(Gate::orientation);
        a.reject(Gate::searchSideElsewhere);
        expectEq(a.bucketed(), 10LL, "every gate lands in its own bucket");
        expectEq(a.impossible(), 0LL, "no impossible bucket touched");
        expectEq(a.failure(10, 0, false), std::string(),
                 "balanced: no failure");
        expect(!a.nothingExamined(), "some surfaces examined");
        // Exactly what dispatch.py's RE_ACCOUNTING reads.
        expectEq(a.summary(10, false),
                 std::string("accepted 10, described 10, recorded 3, "
                             "duplicate 4, other-orientation 2, "
                             "search-side-elsewhere 1, impossible 0, drain "
                             "complete, ok"),
                 "the accounting line's body, byte for byte");
        expectEq(a.failure(11, 0, false),
                 std::string("11 surfaces accepted but 10 described"),
                 "an accepted surface the drain never described");
        expectEq(a.failure(11, 0, true), std::string(),
                 "unless the drain was deliberately cut short");
        expectEq(a.failure(10, 2, false),
                 std::string("2 accepted surfaces failed to rebuild in the "
                             "drain"),
                 "a rebuild failure");
    }
    {
        RowAccounting a;
        a.described = 5;
        a.recorded = 4;
        expectEq(a.failure(5, 0, false),
                 std::string("5 surfaces described but 4 accounted for"),
                 "a described surface in no bucket");
    }
    {
        RowAccounting a;
        a.described = 6;
        a.recorded = 1;
        for (Gate g : {Gate::nonOrientable, Gate::searchSideBroken,
                       Gate::orientationBroken, Gate::multiFarSide,
                       Gate::unnamedSide})
            a.reject(g);
        expectEq(a.impossible(), 5LL, "each impossible gate counts");
        expectEq(a.failure(6, 0, false),
                 std::string("5 surfaces hit a state that cannot occur "
                             "(non-orientable 1, search side 1, orientation 1, "
                             "multi far side 1, unnamed side 1)"),
                 "impossible states are a failure");
    }
    {
        RowAccounting a;
        a.described = 3;
        a.reject(Gate::orientation);
        a.reject(Gate::orientation);
        a.reject(Gate::orientation);
        expect(a.nothingExamined(),
               "every surface turned away: nothing examined");
        expectEq(a.summary(3, true),
                 std::string("accepted 3, described 3, recorded 0, "
                             "duplicate 0, other-orientation 3, "
                             "search-side-elsewhere 0, impossible 0, drain "
                             "skipped, WARNING"),
                 "a skipped drain and a WARNING");
    }
}

// Runs a watchdog until it fires or `wait` passes; returns every reason given.
std::vector<std::string> watch(WatchdogLimits limits, long long satisfying,
                               std::chrono::milliseconds wait) {
    std::mutex m;
    std::vector<std::string> reasons;
    RowWatchdog w(std::move(limits), [&](const char *why) {
        std::lock_guard<std::mutex> lock(m);
        reasons.emplace_back(why);
    });
    w.publishSatisfying(satisfying);
    std::this_thread::sleep_for(wait);
    w.stop();
    w.stop(); // idempotent
    return reasons;
}

void test_watchdog() {
    using namespace std::chrono_literals;
    auto r = watch({}, 0, 450ms);
    expect(r.empty(), "no limit: never fires");

    r = watch({.surfaceTarget = 100, .rowSeconds = 0.0}, 100, 450ms);
    expect(r.size() == 1 && r[0] == "surface-target",
           "target and deadline in one tick: the target is recorded");

    r = watch({.surfaceTarget = 100, .rowSeconds = 0.0}, 99, 450ms);
    expect(r.size() == 1 && r[0] == "timeout",
           "short of the target: the deadline");

    r = watch({.rowSeconds = 10.0}, 0, 450ms);
    expect(r.empty(), "before the deadline: nothing");
}

// Names every far side "far", so a search needs no tables and no complement
// identification: the gates never consult a name beyond its presence.
class FarNamer : public BoundaryNamer {
  public:
    explicit FarNamer(size_t searchSide) : searchSide_(searchSide) {}
    std::string nameLink(size_t bc, const Link &curves) const override {
        if (bc != searchSide_) return "far";
        return curves.comps_.size() == 1 ? identify::identify(curves.comps_.front())
                                         : identify::identify(curves);
    }
    std::string nameCurve(size_t, const Knot &curve) const override {
        return identify::identify(curve);
    }

  private:
    size_t searchSide_;
};

struct GateRun {
    long long accepted = 0;
    RowAccounting acct;
    std::string failure;
    long long drainSkipped = 0;
    int signaturesChecked = 0;
    int signaturesAgreeing = 0;
    int facesMatchingCount = 0;
};

// An exhaustive search at face cap 3 on the canaries' shape: two layers,
// collared through both, `proper`. Every accepted surface is gated;
// the first few also have their deferred pair signature compared.
void gateRun(const std::string &pd, const std::string &name, GateRun &out,
             int signatures) {
    RowBuild rb;
    buildRow(pd, 2, 2, rb);
    SurfaceSearchLimits limits;
    limits.capturePairSig = true;
    SurfaceSearch e(rb.tri, rb.seedFaces, rb.searchSideBC);
    e.configureLimits(limits);
    FarNamer namer(rb.searchSideBC);
    e.setBoundaryNamer(namer);
    e.primeBoundaryName(rb.searchSideBC, rb.searchEdges, name);

    std::mutex m;
    SurfaceSearchCallbacks callbacks;
    callbacks.onSurfaceBoundaryProcessed = [&](const SurfaceBoundaryInfo &info) {
        out.acct.described.fetch_add(1);
        const GatedSurface g = gateSurface(info, rb);
        if (!g.accepted()) {
            out.acct.reject(g.gate);
            return;
        }
        out.acct.recorded.fetch_add(1);
        std::vector<int> faces = info.captureFaces();
        {
            std::lock_guard<std::mutex> lock(m);
            if (static_cast<long long>(faces.size()) == info.triangleCount)
                ++out.facesMatchingCount;
            if (out.signaturesChecked >= signatures)
                return;
            ++out.signaturesChecked;
        }
        const bool same =
            pairSig<4, 2>(rb.tri, faces) == info.capturePairSig();
        std::lock_guard<std::mutex> lock(m);
        if (same)
            ++out.signaturesAgreeing;
    };
    const SearchStats stats = e.search(
        4, conditionFor(BoundaryConditionMode::proper, rb.componentCount),
        callbacks, 0, 0, std::nullopt, std::nullopt,
        /*orientableOnly=*/true, 3, 0, 2);
    out.accepted = stats.satisfyingCount;
    out.failure = out.acct.failure(out.accepted, e.rebuildFailures(),
                                   e.boundaryProcessingSkipped());
}

void test_gates_on_real_searches() {
    {
        GateRun r;
        gateRun("[[1;5;2;4];[3;1;4;6];[5;3;6;2]]", "3_1", r, 3);
        expectEq(r.accepted, 1752LL, "3_1 at cap 3: accepted (canaries)");
        expectEq(r.acct.recorded.load(), 1752LL,
                 "3_1: every surface passes the gates");
        expectEq(r.failure, std::string(), "3_1: the accounting balances");
        expectEq(r.facesMatchingCount, 1752,
                 "3_1: captureFaces() gives every triangle of the surface");
        expectEq(r.signaturesAgreeing, 3,
                 "3_1: pairSig(ambient, captureFaces()) is the captured "
                 "signature");
    }
    {
        GateRun r;
        gateRun("PD[X[4; 1; 3; 2]; X[2; 3; 1; 4]]", "L2a1{0}", r, 2);
        expectEq(r.accepted, 945LL, "L2a1{0} at cap 3: accepted (canaries)");
        expectEq(r.acct.recorded.load(), 795LL,
                 "L2a1{0}: the row's own orientation (canaries)");
        expectEq(r.acct.orientation.load(), 150LL,
                 "L2a1{0}: the other orientation, turned away (canaries)");
        expectEq(r.acct.impossible(), 0LL, "L2a1{0}: nothing impossible");
        expectEq(r.failure, std::string(),
                 "L2a1{0}: the accounting balances");
        expectEq(r.signaturesAgreeing, 2,
                 "L2a1{0}: deferred signatures agree");
    }
}

} // namespace

int main() {
    test_condition_for();
    test_gate_reasons();
    test_accounting();
    test_watchdog();
    test_gates_on_real_searches();
    std::cout << passed << " passed, " << failedCount << " failed\n";
    return failedCount == 0 ? 0 : 1;
}
