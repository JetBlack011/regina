// progress_test.cpp: report::RollingReport (../report/progress.h) redraws
// its block in place -- rewinding exactly the previous block's lines -- and
// commitLine()/forget() leave what is on screen alone.

#include <iostream>
#include <sstream>
#include <string>

#include "surfer/report/progress.h"

static int passed = 0, failed = 0;

static void check(bool ok, const std::string &what) {
    if (ok) {
        ++passed;
    } else {
        ++failed;
        std::cout << "FAIL: " << what << "\n";
    }
}

// What `steps` writes to std::cerr.
template <class Steps> static std::string captured(Steps steps) {
    std::ostringstream out;
    std::streambuf *old = std::cerr.rdbuf(out.rdbuf());
    steps();
    std::cerr.rdbuf(old);
    return out.str();
}

int main() {
    {
        report::RollingReport r;
        const std::string got = captured([&] {
            r.draw("a\nb\n");
            r.draw("c\n");
            r.draw("d\ne\nf\n");
        });
        check(got == "a\nb\n" "\x1b[2F\x1b[0J" "c\n" "\x1b[1F\x1b[0J" "d\ne\nf\n",
              "each draw rewinds the previous block's lines, then prints");
    }
    {
        report::RollingReport r;
        const std::string got = captured([&] {
            r.draw("a\nb\n");
            r.commitLine("kept\n");
            r.draw("c\n");
        });
        check(got == "a\nb\n" "kept\n" "c\n",
              "a committed line is permanent: the next draw erases nothing");
    }
    {
        report::RollingReport r;
        const std::string got = captured([&] {
            r.draw("a\n");
            r.forget();
            r.draw("b\n");
            r.draw("c\n");
        });
        check(got == "a\n" "b\n" "\x1b[1F\x1b[0J" "c\n",
              "forget(): the next draw starts fresh, the one after rewinds again");
    }
    {
        report::RollingReport r;
        const std::string got = captured([&] { r.draw(""); r.draw("x\n"); });
        check(got == "x\n", "an empty block has no lines to rewind");
    }
    std::cout << "progress_test: " << passed << " passed, " << failed << " failed\n";
    return failed == 0 ? 0 : 1;
}
