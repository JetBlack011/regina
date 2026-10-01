// appendonly_test.cpp: appendonly (../cobordisms/appendonly.h) appends whole
// lines, cuts a torn last line before appending, and holds a lock.

#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>

#include <fcntl.h>
#include <stdlib.h>
#include <unistd.h>

#include "cobound/cobordisms/appendonly.h"
#include "linknaming/tests/check.h"

namespace fs = std::filesystem;

static std::string contents(const fs::path &p) {
  std::ifstream in(p, std::ios::binary);
  return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

int main() {
  char tmpl[] = "/tmp/appendonly_test.XXXXXX";
  const fs::path dir = mkdtemp(tmpl);
  const std::string file = (dir / "log.csv").string();

  appendonly::append(file, "a,1\n", appendonly::Sync::yes);
  CHECK_EQ(contents(file), std::string("a,1\n"), "append creates the file");
  appendonly::append(file, "b,2\nc,3\n", appendonly::Sync::no);
  CHECK_EQ(contents(file), std::string("a,1\nb,2\nc,3\n"), "and appends after it");

  // A crash mid-append leaves a torn last line: the next append cuts it.
  std::ofstream(file, std::ios::app | std::ios::binary) << "d,tor";
  appendonly::append(file, "e,5\n", appendonly::Sync::yes);
  CHECK_EQ(contents(file), std::string("a,1\nb,2\nc,3\ne,5\n"),
           "a torn last line is cut before appending");

  // cutTornLine itself: the bytes it cut, and nothing on a whole file.
  {
    std::ofstream(file, std::ios::app | std::ios::binary) << "xyz";
    const int fd = ::open(file.c_str(), O_RDWR);
    CHECK_EQ(appendonly::cutTornLine(fd, file), 3L, "three torn bytes cut");
    CHECK_EQ(appendonly::cutTornLine(fd, file), 0L, "nothing on a whole file");
    ::close(fd);
  }
  {
    const std::string one = (dir / "one.csv").string();
    std::ofstream(one, std::ios::binary) << "no newline at all";
    const int fd = ::open(one.c_str(), O_RDWR);
    CHECK_EQ(appendonly::cutTornLine(fd, one), 17L, "a file of one torn line is emptied");
    ::close(fd);
    CHECK_EQ(contents(one), std::string(""), "empty");
  }

  // The lock: taken and released (a second taker after the first is gone).
  {
    appendonly::FileLock lock(file + ".lock");
  }
  {
    appendonly::FileLock again(file + ".lock");
  }
  CHECK(fs::exists(file + ".lock"), "the lock file is created");

  bool threw = false;
  try {
    appendonly::append((dir / "no-such-dir" / "x").string(), "x\n", appendonly::Sync::yes);
  } catch (const std::runtime_error &) {
    threw = true;
  }
  CHECK(threw, "an unopenable file throws");

  fs::remove_all(dir);
  return cascadetest::finish("appendonly_test");
}
