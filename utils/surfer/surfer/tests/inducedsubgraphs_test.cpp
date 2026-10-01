//
//  enumerator_test.cpp
//
//  ConnectedInducedSubgraphEnumerator against brute force.
//
//  The search can only ever be as complete as this enumerator, and a
//  surface it never visits is lost without a trace. So on many random small
//  graphs this checks that the enumerator visits EXACTLY the connected
//  induced subgraphs a powerset sweep finds -- each one once -- in every
//  configuration the search uses: unseeded and seeded (via contractSeed()),
//  under an anti-monotonic predicate, under a depth cap (the enumerator's
//  own, setMaxSize()), and under per-root work budgets, each pass carrying on
//  from where the last stopped the way EmbeddingSearch::runSearch_() drives
//  it. With a budget it also checks that the passes together visit exactly
//  what one unbudgeted pass does, in the same order, down to a ration of one
//  attempt per pass.
//

#include <algorithm>
#include <functional>
#include <iostream>
#include <map>
#include <optional>
#include <random>
#include <set>
#include <string>
#include <vector>

#include "../enumerate_cis.h"

namespace {

int passed = 0;
int failed_count = 0;

void expect(bool ok, const std::string &desc) {
    if (ok) {
        ++passed;
    } else {
        ++failed_count;
        std::cout << "  FAIL: " << desc << "\n";
    }
}

struct Graph {
    int n;
    std::vector<std::vector<int>> adj; // 1-indexed
};

Graph randomGraph(std::mt19937 &rng, int n, double p) {
    Graph g{n, std::vector<std::vector<int>>(n + 1)};
    std::bernoulli_distribution edge(p);
    for (int u = 1; u <= n; ++u)
        for (int v = u + 1; v <= n; ++v)
            if (edge(rng)) {
                g.adj[u].push_back(v);
                g.adj[v].push_back(u);
            }
    return g;
}

bool connected(const Graph &g, const std::vector<int> &set) {
    if (set.empty())
        return false;
    std::set<int> in(set.begin(), set.end()), seen{set.front()};
    std::vector<int> stack{set.front()};
    while (!stack.empty()) {
        int u = stack.back();
        stack.pop_back();
        for (int v : g.adj[u])
            if (in.contains(v) && seen.insert(v).second)
                stack.push_back(v);
    }
    return seen.size() == in.size();
}

// A hereditary predicate: no two vertices of a forbidden pair together.
// Stateful, as the search's own predicates are (ConditionalPredicate).
struct ForbiddenPairs : ConditionalPredicate {
    std::set<std::pair<int, int>> forbidden; // original ids, u < v
    std::function<std::vector<int>(int)> expand; // graph id -> original ids
    std::vector<int> members;                    // original ids in the set
    std::vector<size_t> marks;

    bool clashes(int a, int b) const {
        return forbidden.contains({std::min(a, b), std::max(a, b)});
    }
    bool admits(const std::vector<int> &set) const {
        for (size_t i = 0; i < set.size(); ++i)
            for (size_t j = i + 1; j < set.size(); ++j)
                if (clashes(set[i], set[j]))
                    return false;
        return true;
    }
    bool tryAdd(int v) override {
        std::vector<int> add = expand(v);
        for (int a : add)
            for (int m : members)
                if (clashes(a, m))
                    return false;
        if (!admits(add))
            return false;
        marks.push_back(members.size());
        members.insert(members.end(), add.begin(), add.end());
        return true;
    }
    void undo(int) override {
        members.resize(marks.back());
        marks.pop_back();
    }
};

struct Config {
    bool seeded = false;
    bool predicate = false;
    std::optional<int> cap;       // faces added beyond the seed
    std::optional<long long> budgetStart;
};

// Every connected induced subgraph brute force admits: containing the seed
// and strictly larger than it when seeded, within the cap, and passing the
// predicate.
std::set<std::vector<int>> bruteForce(const Graph &g,
                                      const std::vector<int> &seed,
                                      const ForbiddenPairs *pred,
                                      std::optional<int> cap) {
    std::set<std::vector<int>> out;
    for (unsigned mask = 1; mask < (1u << g.n); ++mask) {
        std::vector<int> set;
        for (int v = 1; v <= g.n; ++v)
            if (mask & (1u << (v - 1)))
                set.push_back(v);
        if (!seed.empty()) {
            if (!std::includes(set.begin(), set.end(), seed.begin(),
                               seed.end()) ||
                set.size() == seed.size())
                continue;
        }
        const int added = static_cast<int>(set.size() - seed.size());
        if (cap && added > *cap)
            continue;
        if (!connected(g, set))
            continue;
        if (pred && !pred->admits(set))
            continue;
        out.insert(set);
    }
    return out;
}

// Runs the enumerator the way runSearch_() drives it, returning every set it
// reports (in original ids), in the order reported.
std::vector<std::vector<int>> enumerateLikeSearch(
    const Graph &g, const std::vector<int> &seed, ForbiddenPairs *pred,
    const Config &config) {
    std::optional<ConnectedInducedSubgraphEnumerator::SeededGraph> sg;
    const Graph *graph = &g;
    Graph contracted;
    if (config.seeded) {
        sg = ConnectedInducedSubgraphEnumerator::contractSeed(g.n, g.adj,
                                                              seed);
        contracted = Graph{sg->n, sg->adj};
        graph = &contracted;
    }
    auto expand = [&](int v) -> std::vector<int> {
        if (!config.seeded)
            return {v};
        if (v == 1)
            return seed;
        return {sg->originalOf[v]};
    };

    struct AlwaysYes : ConditionalPredicate {
        bool tryAdd(int) override { return true; }
        void undo(int) override {}
    } yes;
    ConditionalPredicate *inner = &yes;
    if (pred) {
        pred->expand = expand;
        pred->members.clear();
        pred->marks.clear();
        inner = pred;
    }
    BudgetedPredicate budgeted(*inner, -1);

    std::optional<ConnectedInducedSubgraphEnumerator> e;
    if (config.seeded)
        e.emplace(graph->n, graph->adj, true, budgeted);
    else
        e.emplace(graph->n, graph->adj);
    // The cap as the search sets it: the enumerator's own, in vertices, the
    // seed's contracted vertex counting as one.
    if (config.cap)
        e->setMaxSize(static_cast<size_t>(
            std::max(1, config.seeded ? *config.cap + 1 : *config.cap)));

    std::vector<int> roots = e->getRoots();
    std::sort(roots.begin(), roots.end(), std::greater<>());

    std::vector<std::vector<int>> seen;
    auto record = [&](const std::vector<int> &U) {
        std::vector<int> set;
        for (int v : U)
            for (int o : expand(v))
                set.push_back(o);
        std::sort(set.begin(), set.end());
        seen.push_back(std::move(set));
    };

    // Each root's passes carry on from where the last stopped, the budget
    // growing as in the search; the re-adds back down are uncharged.
    for (int s : roots) {
        ConnectedInducedSubgraphEnumerator::Position position;
        for (unsigned level = 0;; ++level) {
            long long budget = -1;
            if (config.budgetStart)
                budget = *config.budgetStart << level;
            e->resetCandidateOrder();
            budgeted.reset(budget);
            budgeted.freeAttempts(position.rebuildAttempts());
            const auto outcome =
                e->enumerateFromRootFiltered(s, record, budgeted, &position);
            if (outcome != ConnectedInducedSubgraphEnumerator::Outcome::suspended)
                break;
        }
    }
    return seen;
}

void runOne(std::mt19937 &rng, const Config &config, const std::string &tag) {
    std::uniform_int_distribution<int> sizeDist(4, 11);
    const int n = sizeDist(rng);
    Graph g = randomGraph(rng, n, 0.35);

    std::vector<int> seed;
    if (config.seeded) {
        // A connected seed of 1-3 vertices, grown from a random vertex.
        std::uniform_int_distribution<int> vDist(1, n);
        seed.push_back(vDist(rng));
        std::uniform_int_distribution<int> want(1, 3);
        const int target = want(rng);
        while (static_cast<int>(seed.size()) < target) {
            std::vector<int> frontier;
            for (int u : seed)
                for (int v : g.adj[u])
                    if (std::find(seed.begin(), seed.end(), v) == seed.end())
                        frontier.push_back(v);
            if (frontier.empty())
                break;
            std::uniform_int_distribution<size_t> pick(0, frontier.size() - 1);
            seed.push_back(frontier[pick(rng)]);
        }
        std::sort(seed.begin(), seed.end());
    }

    std::optional<ForbiddenPairs> pred;
    if (config.predicate) {
        pred.emplace();
        std::bernoulli_distribution forbid(0.15);
        std::set<int> seedSet(seed.begin(), seed.end());
        for (int u = 1; u <= n; ++u)
            for (int v = u + 1; v <= n; ++v)
                if (!(seedSet.contains(u) && seedSet.contains(v)) &&
                    forbid(rng))
                    pred->forbidden.insert({u, v});
    }

    auto expected =
        bruteForce(g, seed, pred ? &*pred : nullptr, config.cap);
    const auto sequence =
        enumerateLikeSearch(g, seed, pred ? &*pred : nullptr, config);
    std::map<std::vector<int>, int> got;
    for (const auto &set : sequence)
        ++got[set];

    size_t duplicates = 0, extra = 0, missing = 0;
    for (const auto &[set, count] : got) {
        duplicates += count > 1;
        extra += !expected.contains(set);
    }
    for (const auto &set : expected)
        missing += !got.contains(set);
    expect(duplicates == 0 && extra == 0 && missing == 0,
           tag + " (n=" + std::to_string(n) + "): " +
               std::to_string(missing) + " missing, " +
               std::to_string(extra) + " extra, " +
               std::to_string(duplicates) + " reported more than once");

    // Budgeted passes that carry on from one another must visit exactly
    // what one unbudgeted pass does, in the same order -- at any ration,
    // however small (a ration of 1 suspends after every attempt).
    if (config.budgetStart) {
        Config unbudgeted = config;
        unbudgeted.budgetStart.reset();
        const auto reference =
            enumerateLikeSearch(g, seed, pred ? &*pred : nullptr, unbudgeted);
        expect(sequence == reference,
               tag + " (n=" + std::to_string(n) +
                   "): budgeted passes visit in the unbudgeted order (" +
                   std::to_string(sequence.size()) + " vs " +
                   std::to_string(reference.size()) + " visits)");
    }
}

} // namespace

int main() {
    std::mt19937 rng(20260926);
    const std::vector<std::pair<std::string, Config>> configs = {
        {"unseeded", {}},
        {"seeded", {.seeded = true}},
        {"unseeded + predicate", {.predicate = true}},
        {"seeded + predicate", {.seeded = true, .predicate = true}},
        {"seeded + cap 2", {.seeded = true, .cap = 2}},
        {"unseeded + cap 3", {.cap = 3}},
        {"seeded + predicate + budget", {.seeded = true, .predicate = true,
                                         .budgetStart = 2}},
        {"unseeded + budget", {.budgetStart = 1}},
        {"seeded + cap 3 + budget",
         {.seeded = true, .cap = 3, .budgetStart = 3}},
        {"seeded + predicate + cap 4 + budget 1",
         {.seeded = true, .predicate = true, .cap = 4, .budgetStart = 1}},
        {"unseeded + predicate + budget 7",
         {.predicate = true, .budgetStart = 7}},
    };
    for (const auto &[tag, config] : configs)
        for (int trial = 0; trial < 60; ++trial)
            runOne(rng, config, tag);
    std::cout << passed << " passed, " << failed_count << " failed\n";
    return failed_count == 0 ? 0 : 1;
}
