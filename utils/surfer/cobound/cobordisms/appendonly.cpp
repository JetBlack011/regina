//
//  appendonly.cpp
//

#include "cobound/cobordisms/appendonly.h"

#include <cerrno>
#include <cstring>
#include <iostream>
#include <stdexcept>

#include <fcntl.h>
#include <sys/file.h>
#include <unistd.h>

namespace appendonly {

FileLock::FileLock(const std::string &path) {
    fd_ = ::open(path.c_str(), O_RDWR | O_CREAT | O_CLOEXEC, 0644);
    if (fd_ < 0)
        throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
    while (::flock(fd_, LOCK_EX) != 0)
        if (errno != EINTR) {
            const int err = errno;
            ::close(fd_);
            throw std::runtime_error("cannot lock " + path + ": " + std::strerror(err));
        }
}

FileLock::~FileLock() {
    if (fd_ >= 0) {
        ::flock(fd_, LOCK_UN);
        ::close(fd_);
    }
}

void writeAll(int fd, const std::string &bytes, const std::string &what) {
    const char *p = bytes.data();
    size_t left = bytes.size();
    while (left > 0) {
        const ssize_t n = ::write(fd, p, left);
        if (n < 0) {
            if (errno == EINTR)
                continue;
            throw std::runtime_error("write to " + what + " failed: " + std::strerror(errno));
        }
        p += n;
        left -= static_cast<size_t>(n);
    }
}

void sync(int fd, const std::string &what) {
    if (::fsync(fd) != 0)
        throw std::runtime_error("fsync of " + what + " failed: " + std::strerror(errno));
}

long cutTornLine(int fd, const std::string &what) {
    const off_t size = ::lseek(fd, 0, SEEK_END);
    if (size < 0)
        throw std::runtime_error("cannot seek " + what);
    if (size == 0)
        return 0;
    char last = 0;
    if (::pread(fd, &last, 1, size - 1) != 1)
        throw std::runtime_error("cannot read " + what);
    if (last == '\n')
        return 0;
    // A torn append: find the last complete line and cut back to it.
    off_t cut = size;
    char c = 0;
    while (cut > 0 && ::pread(fd, &c, 1, cut - 1) == 1 && c != '\n')
        --cut;
    std::cerr << "[!] " << what << ": truncating a torn last line (" << (size - cut)
              << " bytes)\n";
    if (::ftruncate(fd, cut) != 0)
        throw std::runtime_error("cannot truncate " + what);
    return static_cast<long>(size - cut);
}

void append(const std::string &path, const std::string &bytes, Sync syncIt) {
    const int fd = ::open(path.c_str(), O_RDWR | O_CREAT | O_APPEND | O_CLOEXEC, 0644);
    if (fd < 0)
        throw std::runtime_error("cannot open " + path + ": " + std::strerror(errno));
    try {
        cutTornLine(fd, path);
        writeAll(fd, bytes, path);
        if (syncIt == Sync::yes)
            sync(fd, path);
    } catch (...) {
        ::close(fd);
        throw;
    }
    ::close(fd);
}

} // namespace appendonly
