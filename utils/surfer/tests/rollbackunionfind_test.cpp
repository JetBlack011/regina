// rollbackunionfind_test.cpp
//
// Direct unit tests of RollbackUnionFind, independent of EmbeddedSubmanifold:
// unite()/find() behave like an ordinary union-find, and rollbackTo() undoes
// unite() calls back to a checkpoint exactly as if they had never happened
// -- including that further unite() calls after a rollback behave correctly,
// not merely that find() happens to report the old roots.

#include <iostream>
#include <unistd.h>

#include "../rollbackunionfind.h"

static int passed = 0, failed_count = 0;

namespace {
bool colorEnabled() {
    static bool enabled = isatty(fileno(stdout));
    return enabled;
}
std::ostream &green(std::ostream &os) {
    return colorEnabled() ? os << "\033[32m" : os;
}
std::ostream &red(std::ostream &os) {
    return colorEnabled() ? os << "\033[31m" : os;
}
std::ostream &bold(std::ostream &os) {
    return colorEnabled() ? os << "\033[1m" : os;
}
std::ostream &resetColor(std::ostream &os) {
    return colorEnabled() ? os << "\033[0m" : os;
}
} // namespace

#define EXPECT_EQ(actual, expected, desc)                                      \
    do {                                                                       \
        auto _a = (actual);                                                    \
        auto _e = (expected);                                                  \
        if (_a == _e) {                                                        \
            std::cout << green << "  PASS: " << resetColor << (desc) << "\n";  \
            ++passed;                                                          \
        } else {                                                               \
            std::cout << red << "  FAIL: " << (desc) << "\n"                   \
                      << "        expected " << _e << ", got " << _a           \
                      << resetColor << "\n";                                   \
            ++failed_count;                                                    \
        }                                                                      \
    } while (0)

void test_basic_unite_and_find() {
    std::cout << "\n--- basic unite()/find() ---\n";

    RollbackUnionFind dsu(5);
    for (int i = 0; i < 5; ++i)
        EXPECT_EQ(dsu.find(i), i, "element " + std::to_string(i) +
                                       " starts as its own root");

    EXPECT_EQ(dsu.unite(0, 1), true, "uniting 0 and 1 succeeds");
    EXPECT_EQ(dsu.find(0), dsu.find(1), "0 and 1 share a root after uniting");
    EXPECT_EQ(dsu.unite(0, 1), false,
              "uniting already-same-set elements is a no-op");

    dsu.unite(2, 3);
    EXPECT_EQ(dsu.find(2), dsu.find(3), "2 and 3 share a root");
    EXPECT_EQ(dsu.find(0) == dsu.find(2), false,
              "{0,1} and {2,3} remain separate classes");

    dsu.unite(1, 2);
    EXPECT_EQ(dsu.find(0), dsu.find(3),
              "uniting a member of each class merges {0,1} and {2,3}");
    EXPECT_EQ(dsu.find(4) == dsu.find(0), false,
              "untouched element 4 stays in its own class");
}

void test_rollback_single_union() {
    std::cout << "\n--- rollback a single unite() ---\n";

    RollbackUnionFind dsu(4);
    size_t mark = dsu.checkpoint();
    dsu.unite(0, 1);
    EXPECT_EQ(dsu.find(0), dsu.find(1), "0 and 1 united");

    dsu.rollbackTo(mark);
    EXPECT_EQ(dsu.find(0), 0, "0 is its own root again after rollback");
    EXPECT_EQ(dsu.find(1), 1, "1 is its own root again after rollback");
}

void test_rollback_multiple_unions_lifo() {
    std::cout << "\n--- rollback several unite() calls, LIFO ---\n";

    RollbackUnionFind dsu(6);
    size_t mark0 = dsu.checkpoint();
    dsu.unite(0, 1);
    size_t mark1 = dsu.checkpoint();
    dsu.unite(1, 2);
    dsu.unite(3, 4);
    EXPECT_EQ(dsu.find(0), dsu.find(2), "{0,1,2} fully merged");
    EXPECT_EQ(dsu.find(3), dsu.find(4), "{3,4} merged");

    // Roll back only to mark1: undoes unite(3,4) and unite(1,2), but NOT
    // unite(0,1) (which happened before mark1 was taken).
    dsu.rollbackTo(mark1);
    EXPECT_EQ(dsu.find(0), dsu.find(1), "{0,1} survives rollback to mark1");
    EXPECT_EQ(dsu.find(1) == dsu.find(2), false,
              "1 and 2 are separate again after rollback to mark1");
    EXPECT_EQ(dsu.find(3) == dsu.find(4), false,
              "3 and 4 are separate again after rollback to mark1");
    EXPECT_EQ(dsu.find(5), 5, "untouched element 5 unaffected");

    // Roll back the rest, to mark0: undoes unite(0,1) too.
    dsu.rollbackTo(mark0);
    EXPECT_EQ(dsu.find(0), 0, "0 is its own root again after full rollback");
    EXPECT_EQ(dsu.find(1), 1, "1 is its own root again after full rollback");
}

void test_reunite_after_rollback_behaves_fresh() {
    std::cout << "\n--- unite() after rollback behaves as if history never "
                 "happened ---\n";

    RollbackUnionFind dsu(4);
    size_t mark = dsu.checkpoint();
    dsu.unite(0, 1);
    dsu.unite(2, 3);
    dsu.unite(1, 2); // merges everything into one class
    dsu.rollbackTo(mark);

    // Re-derive an entirely different partition after rollback: {0,2} and
    // {1,3} instead of the {0,1,2,3} single class from before.
    dsu.unite(0, 2);
    dsu.unite(1, 3);
    EXPECT_EQ(dsu.find(0), dsu.find(2), "{0,2} merged in the new history");
    EXPECT_EQ(dsu.find(1), dsu.find(3), "{1,3} merged in the new history");
    EXPECT_EQ(dsu.find(0) == dsu.find(1), false,
              "{0,2} and {1,3} remain separate -- no leftover state from the "
              "rolled-back {0,1,2,3} merge");
}

void run(const std::string &name, void (*fn)()) {
    std::cout << bold << "\n=== " << name << " ===" << resetColor << "\n";
    fn();
}

int main() {
    run("basic_unite_and_find", test_basic_unite_and_find);
    run("rollback_single_union", test_rollback_single_union);
    run("rollback_multiple_unions_lifo", test_rollback_multiple_unions_lifo);
    run("reunite_after_rollback_behaves_fresh",
        test_reunite_after_rollback_behaves_fresh);

    std::cout << bold << "\n=== Summary: " << passed << " passed, "
              << failed_count << " failed ===" << resetColor << "\n";
    return failed_count > 0 ? 1 : 0;
}
