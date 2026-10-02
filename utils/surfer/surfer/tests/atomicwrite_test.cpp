// atomicwrite_test.cpp: report::atomicWrite() (../report/atomicwrite.h)
// replaces a file whole -- new contents, nothing left beside it -- and on
// any failure leaves the old file as it was.

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>

#include <stdlib.h>

#include "surfer/report/atomicwrite.h"

namespace fs = std::filesystem;

static int passed = 0, failed = 0;

static void check(bool ok, const std::string &what) {
    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "FAIL: " << what << "\n";
    }
}

static std::string contents(const fs::path &p) {
    std::ifstream in(p, std::ios::binary);
    return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

static size_t entries(const fs::path &dir) {
    size_t n = 0;
    for ([[maybe_unused]] const auto &e : fs::directory_iterator(dir))
        ++n;
    return n;
}

int main() {
    char tmpl[] = "/tmp/atomicwrite_test.XXXXXX";
    const fs::path dir = mkdtemp(tmpl);
    const fs::path file = dir / "out.txt";

    report::atomicWrite(file, [](std::ostream &out) { out << "first\n"; });
    check(contents(file) == "first\n", "a new file holds what was written");
    check(entries(dir) == 1, "and nothing is left beside it");

    report::atomicWrite(file, [](std::ostream &out) { out << "second\nline\n"; });
    check(contents(file) == "second\nline\n", "a second write replaces the first whole");
    check(entries(dir) == 1, "still nothing beside it");

    bool threw = false;
    try {
        report::atomicWrite(file, [](std::ostream &out) {
            out << "half";
            throw std::runtime_error("writer failed");
        });
    } catch (const std::runtime_error &) {
        threw = true;
    }
    check(threw, "a failing writer's exception propagates");
    check(contents(file) == "second\nline\n", "and the old file is untouched");
    check(entries(dir) == 1, "and its temporary file is removed");

    threw = false;
    try {
        report::atomicWrite(dir / "no-such-dir" / "x", [](std::ostream &out) { out << "x"; });
    } catch (const std::runtime_error &) {
        threw = true;
    }
    check(threw, "an unwritable path throws");

    report::atomicWrite(file, [](std::ostream &) {});
    check(contents(file).empty(), "an empty write leaves an empty file");

    report::atomicWrite(file, [](std::ostream &out) { out << "cache\n"; },
                        report::Durability::cache);
    check(contents(file) == "cache\n", "a cache's write replaces the file whole too");
    check(entries(dir) == 1, "and leaves nothing beside it");

    fs::remove_all(dir);
    std::cout << "atomicwrite_test: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
