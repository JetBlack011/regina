//
//  progress.cpp
//

#include "surfer/report/progress.h"

#include <algorithm>
#include <iostream>

namespace report {

void RollingReport::draw(const std::string &text) {
    std::lock_guard<std::mutex> lock(mutex_);
    if (prevLines_ > 0)
        std::cerr << "\x1b[" << prevLines_ << "F\x1b[0J";
    std::cerr << text;
    prevLines_ = static_cast<std::size_t>(std::count(text.begin(), text.end(), '\n'));
}

void RollingReport::commitLine(const std::string &text) {
    std::lock_guard<std::mutex> lock(mutex_);
    std::cerr << text;
    prevLines_ = 0;
}

void RollingReport::forget() {
    std::lock_guard<std::mutex> lock(mutex_);
    prevLines_ = 0;
}

} // namespace report
