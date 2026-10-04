//
//  tableclasses.cpp
//
//  The table's link classes: every knot and link table entry whose
//  canonical name is another entry's, as linknaming::LinkNamer names it.
//
//  Usage:
//    tableclasses --knots <table.csv> --links <table.csv> [--symmetry <csv>]
//                 > table_link_classes.csv
//
//  stdout, CSV: name,canonical,proof -- one line per entry whose canonical
//  name differs from its own. `proof` is `diagram` when a version of one
//  table diagram is the other (Tables::canonical()), or `isometry` when
//  an isometry of the complements carries meridians to meridians with a
//  uniform orientation sign (isometry/isometry.h): either way the two entries
//  are one oriented link up to mirror and global reversal, so one link of
//  the graph. The solvers read the file opt-in (C++ --link-classes, frontier.py
//  --link-classes). A class whose members' literature 4-genera differ is a
//  table error or a bug: it is reported on stderr and the run fails.
//

#include <iostream>
#include <string>

#include "linknaming/linknamer.h"
#include "linknaming/tables.h"

int main(int argc, char **argv) {
    std::string knots, links, symmetry;
    for (int i = 1; i + 1 < argc; i += 2) {
        const std::string a = argv[i], v = argv[i + 1];
        if (a == "--knots") knots = v;
        else if (a == "--links") links = v;
        else if (a == "--symmetry") symmetry = v;
        else { std::cerr << "unknown option " << a << "\n"; return 2; }
    }
    if (knots.empty() || links.empty()) {
        std::cerr << "usage: tableclasses --knots <csv> --links <csv> [--symmetry <csv>]\n";
        return 2;
    }
    const linknaming::Tables tables = linknaming::Tables::load(knots, links, symmetry);
    const linknaming::LinkNamer namer(tables);
    std::cout << "name,canonical,proof\n";
    size_t merged = 0;
    for (const linknaming::TableEntry &e : tables.entries()) {
        const std::string &c = namer.canonicalName(e);
        if (c == e.name)
            continue;
        ++merged;
        std::cout << e.name << ',' << c << ','
                  << (tables.canonical(e.name) == c ? "diagram" : "isometry") << '\n';
    }
    const auto conflicts = namer.classConflicts();
    for (const std::string &conflict : conflicts)
        std::cerr << "[!] class with two literature values: " << conflict << "\n";
    std::cerr << merged << " of " << tables.size() << " entries named by another entry's name\n";
    return conflicts.empty() ? 0 : 1;
}
