//
//  atomicwrite.cpp
//

#include "surfer/report/atomicwrite.h"

#include <atomic>
#include <cerrno>
#include <cstring>
#include <fstream>
#include <stdexcept>
#include <string>
#include <system_error>

#include <fcntl.h>
#include <unistd.h>

namespace report {

namespace {

std::atomic<unsigned long> counter{0};

// fsync()s `path` (a file, or with O_DIRECTORY a directory) by name.
void syncPath(const std::filesystem::path &path, int flags) {
    const int fd = ::open(path.c_str(), flags | O_CLOEXEC);
    if (fd < 0)
        throw std::runtime_error("cannot open " + path.string() + " to sync it: " +
                                 std::strerror(errno));
    const int synced = ::fsync(fd);
    const int err = errno;
    ::close(fd);
    if (synced != 0)
        throw std::runtime_error("fsync of " + path.string() + " failed: " +
                                 std::strerror(err));
}

} // namespace

void atomicWrite(const std::filesystem::path &path,
                 const std::function<void(std::ostream &)> &write, Durability durability) {
    const bool durable = durability == Durability::durable;
    std::filesystem::path tmp = path;
    tmp += ".tmp." + std::to_string(::getpid()) + "." + std::to_string(counter++);
    try {
        {
            std::ofstream out(tmp, std::ios::binary | std::ios::trunc);
            if (!out)
                throw std::runtime_error("cannot write " + tmp.string());
            write(out);
            out.flush();
            if (!out)
                throw std::runtime_error("error writing " + tmp.string());
        }
        if (durable)
            syncPath(tmp, O_RDONLY);
        std::error_code ec;
        std::filesystem::rename(tmp, path, ec);
        if (ec)
            throw std::runtime_error("cannot rename " + tmp.string() + " to " +
                                     path.string() + ": " + ec.message());
    } catch (...) {
        std::error_code ignored;
        std::filesystem::remove(tmp, ignored);
        throw;
    }
    if (!durable)
        return;
    const std::filesystem::path dir =
        path.has_parent_path() ? path.parent_path() : std::filesystem::path(".");
    syncPath(dir, O_RDONLY | O_DIRECTORY);
}

} // namespace report
