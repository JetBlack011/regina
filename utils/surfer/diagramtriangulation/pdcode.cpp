//
//  pdcode.cpp
//
//  parsePDCode() moved here from knotbuilder.cpp (now fromdiagram.cpp).
//

#include "diagramtriangulation/pdcode.h"

#include <algorithm>
#include <cctype>
#include <sstream>

diagramtriangulation::PDCode diagramtriangulation::parsePDCode(std::string pdcode_str) {
    std::vector<std::array<int, 4>> pdcode;

    for (char &c : pdcode_str) {
        if (!std::isdigit(c)) {
            c = ' ';
        }
    }

    std::stringstream ss(pdcode_str);
    std::vector<int> pdlist;
    int token;
    while (ss >> token) {
        pdlist.push_back(token);
    }

    // PD codes conventionally 1-index strand labels; detect 0-indexed
    // input (a literal 0 can only appear in that case) and normalize to
    // 0-indexed internally either way.
    bool isZeroIndexed = false;
    if (std::ranges::find(pdlist, 0) != pdlist.end()) {
        isZeroIndexed = true;
    }

    if (!isZeroIndexed) {
        for (int &i : pdlist) {
            --i;
        }
    }

    for (int i = 0; i < pdlist.size(); i += 4) {
        std::array<int, 4> crossing;
        for (int j = 0; j < 4; j++) {
            crossing[j] = pdlist[i + j];
        }
        pdcode.push_back(crossing);
    }

    return pdcode;
}
