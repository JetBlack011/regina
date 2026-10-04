//
//  name.cpp
//
//  cobound name (was farsidename): names for stored cobordisms'
//  outgoing links, from their pair signatures.
//
//  Usage:
//    cobound name --knots <table.csv> --links <table.csv> [--symmetry <csv>]
//                 [--search-height H] [--search-visits N] [--simplify-tries N]
//                 [--exhaustive-height H] [--max-search-crossings N]
//                 [--deep-height H] [--deep-visits N] [--max-deep-crossings N]
//                 [--profile] [--reference] < cobordisms
//
//  The flags are farsidename's, which the atlas's run_farsidename.py passes;
//  each is a config key (knot_table, link_table, knot_symmetry, namer_*,
//  name_profile, name_reference), so --config and --set work too. The limits
//  default to linknaming::NamerLimits.
//
//  stdin: one "<id>\t<incoming>\t<thicken_layers>\t<pair signature>" per
//  line, grouped by incoming diagram (each one's thickening is built once).
//  <incoming> is a table name, or a PD code itself ("[...]" or "PD[...]":
//  a goal run's search is on a graph link's own diagram). stdout, one line
//  per cobordism, tab-separated:
//
//    <id> <incoming> ok <name> <exact> <pinned> <components> <split_unknots>
//        <pieces> <proof> <drawn_crossings> <ms>
//    <id> <incoming> FAILED <reason>
//
//  The outgoing link is recovered exactly as the search saw it
//  (outgoing::OutgoingReader: the pair carried onto the incoming
//  diagram's own thickening by an isomorphism pinned by L x {0}), oriented as a cobordism
//  from the oriented incoming link, drawn (knotbuilder::DiagramDrawer, which
//  refuses any non-planar drawing) and named by linknaming::LinkNamer.
//  <exact> is 1 when the name is an identity, <pinned> when every piece's
//  oriented variant is proved; <pieces> lists "display/by/crossings" per
//  piece, joined by " & " (a link tag contains ';'). Names are cached per
//  drawn diagram (its exact signature), so an outgoing link seen again costs
//  one drawing.
//

#include <chrono>
#include <fstream>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>

#include "linknaming/linknamer.h"
#include "linknaming/tables.h"
#include "cobound/outgoing/fromdatabase.h"
#include "cobound/driver/commands.h"
#include "cobound/driver/config.h"

namespace {

std::unordered_map<std::string, std::string> pdCodes(const std::string &path) {
    std::unordered_map<std::string, std::string> out;
    for (const linknaming::TableRow &row : linknaming::readTableRows(path))
        out.emplace(row.name, row.pd);
    return out;
}

const char *byName(linknaming::PieceName::By by) {
    switch (by) {
        case linknaming::PieceName::By::exactDiagram: return "diagram";
        case linknaming::PieceName::By::isometry: return "isometry";
        case linknaming::PieceName::By::searchAndInvariants: return "search";
        default: return "untabulated";
    }
}

} // namespace

int commands::name(const std::vector<std::string> &args) {
    std::string knots, links, symmetry;
    linknaming::NamerLimits limits;
    bool profile = false, reference = false;
    try {
        const config::Config cfg = config::forCommand(
            "name", config::Context::name, args,
            {{"--knots", "knot_table"},
             {"--links", "link_table"},
             {"--symmetry", "knot_symmetry"},
             {"--search-height", "namer_search_height"},
             {"--search-visits", "namer_search_visits"},
             {"--simplify-tries", "namer_simplify_tries"},
             {"--exhaustive-height", "namer_exhaustive_height"},
             {"--max-search-crossings", "namer_max_search_crossings"},
             {"--deep-height", "namer_deep_height"},
             {"--deep-visits", "namer_deep_visits"},
             {"--max-deep-crossings", "namer_max_deep_crossings"},
             {"--profile", "name_profile", false, "1"},
             {"--reference", "name_reference", false, "1"}});
        knots = cfg.text("knot_table");
        links = cfg.text("link_table");
        symmetry = cfg.text("knot_symmetry");
        limits.searchHeight = static_cast<int>(cfg.integer("namer_search_height"));
        limits.searchVisits = static_cast<size_t>(cfg.integer("namer_search_visits"));
        limits.simplifyTries = static_cast<int>(cfg.integer("namer_simplify_tries"));
        limits.exhaustiveHeight = static_cast<int>(cfg.integer("namer_exhaustive_height"));
        limits.maxSearchCrossings = static_cast<size_t>(cfg.integer("namer_max_search_crossings"));
        limits.deepHeight = static_cast<int>(cfg.integer("namer_deep_height"));
        limits.deepVisits = static_cast<size_t>(cfg.integer("namer_deep_visits"));
        limits.maxDeepCrossings = static_cast<size_t>(cfg.integer("namer_max_deep_crossings"));
        profile = cfg.flag("name_profile");
        reference = cfg.flag("name_reference");
    } catch (const config::Error &e) {
        std::cerr << e.what() << "\n"
                  << "usage: cobound name --knots <csv> --links <csv> [--symmetry <csv>]\n"
                     "       [--search-height H] [--search-visits N] [--simplify-tries N]\n"
                     "       [--exhaustive-height H] [--max-search-crossings N]\n"
                     "       [--deep-height H] [--deep-visits N] [--max-deep-crossings N]\n"
                     "       [--profile] [--reference] < witnesses\n"
                     "(defaults: exactnaming::NamerLimits)\n";
        return 2;
    }
    const linknaming::Tables tables = linknaming::Tables::load(knots, links, symmetry);
    for (const std::string &c : tables.inconsistentClasses())
        std::cerr << "[!] table classes with two literature values: " << c << "\n";
    const linknaming::LinkNamer namer(tables, limits);
    std::unordered_map<std::string, std::string> pd = pdCodes(knots);
    for (auto &[k, v] : pdCodes(links)) pd.emplace(k, v);

    // --profile: cumulative ms per step, printed to stderr at the end.
    double msIncoming = 0, msOrient = 0, msDraw = 0, msName = 0, msDecode = 0, msIso = 0;
    double msSurface = 0, msRead = 0, msBoundaryBuild = 0;
    long cobordisms = 0, diagrams = 0, named = 0;
    auto ms = [](std::chrono::steady_clock::time_point t) {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t).count();
    };
    auto flushIncoming = [&](const outgoing::OutgoingReader *r) {
        if (r) {
            msDecode += r->msDecode(); msIso += r->msIsomorphism();
            msSurface += r->msSurfaceBuild(); msRead += r->msBoundaryRead();
            msBoundaryBuild += r->msBoundaryBuild();
        }
    };
    std::unique_ptr<outgoing::OutgoingReader> redraw;
    std::string redrawKey;
    std::unordered_map<std::string, linknaming::LinkName> cache; // per incoming diagram
    std::string line;
    while (std::getline(std::cin, line)) {
        std::istringstream in(line);
        std::string id, incoming, layers, sig;
        if (!std::getline(in, id, '\t') || !std::getline(in, incoming, '\t') ||
            !std::getline(in, layers, '\t') || !std::getline(in, sig))
            continue;
        const auto start = std::chrono::steady_clock::now();
        try {
            if (redrawKey != incoming + "\t" + layers) {
                // The incoming diagram is a table name, or the PD code itself: a
                // goal run's search is on a graph link's own diagram, which no
                // table holds (its database's <cobordisms>.rows.csv gives it;
                // cobordisms/pending.h).
                std::string code;
                if (!incoming.empty() && (incoming.front() == '[' || incoming.rfind("PD[", 0) == 0)) {
                    code = incoming;
                } else {
                    auto it = pd.find(incoming);
                    if (it == pd.end()) throw regina::InvalidArgument("no PD code for row " + incoming);
                    code = it->second;
                }
                for (char &ch : code)
                    if (ch == ';') ch = ',';
                flushIncoming(redraw.get());
                redraw.reset();
                const auto tr = std::chrono::steady_clock::now();
                redraw = std::make_unique<outgoing::OutgoingReader>(code, std::stoi(layers));
                msIncoming += ms(tr);
                ++diagrams;
                redrawKey = incoming + "\t" + layers;
                cache.clear();
            }
            ++cobordisms;
            std::string why;
            const double before = redraw->msDecode() + redraw->msIsomorphism();
            const auto to = std::chrono::steady_clock::now();
            auto link = reference ? redraw->outgoingLink(sig, why)
                                  : redraw->outgoingLinkFast(sig, why);
            msOrient += ms(to) - (redraw->msDecode() + redraw->msIsomorphism() - before);
            if (!link) {
                std::cout << id << '\t' << incoming << "\tFAILED\t" << why << '\n';
                continue;
            }
            const auto td = std::chrono::steady_clock::now();
            const knotbuilder::Diagram d = redraw->drawer().draw(link->curves);
            const regina::Link drawn = d.link();
            const std::string key = drawn.sig<2>(false, false, true);
            msDraw += ms(td);
            auto hit = cache.find(key);
            if (hit == cache.end()) {
                const auto tn = std::chrono::steady_clock::now();
                hit = cache.emplace(key, namer.name(drawn)).first;
                msName += ms(tn);
                ++named;
            }
            const linknaming::LinkName &n = hit->second;
            std::ostringstream pieces;
            for (size_t i = 0; i < n.pieces.size(); ++i)
                pieces << (i ? " & " : "") << n.pieces[i].display() << '/' << byName(n.pieces[i].by)
                       << '/' << n.pieces[i].crossings;
            const long ms = std::chrono::duration_cast<std::chrono::milliseconds>(
                                std::chrono::steady_clock::now() - start).count();
            std::cout << id << '\t' << incoming << "\tok\t" << n.name << '\t' << n.isName << '\t'
                      << n.pinned << '\t' << d.components << '\t' << n.splitUnknots << '\t'
                      << pieces.str() << '\t' << n.proof() << '\t' << drawn.size() << '\t'
                      << ms << '\n';
        } catch (const std::exception &e) {
            std::cout << id << '\t' << incoming << "\tFAILED\t" << e.what() << '\n';
        }
    }
    flushIncoming(redraw.get());
    if (profile)
        std::cerr << "profile: " << cobordisms << " witnesses, " << diagrams << " rows, "
                  << named << " distinct diagrams named; ms: build row " << msIncoming
                  << ", decode pairsig " << msDecode << ", isomorphism search " << msIso
                  << ", surface+orient " << msOrient << " (KnottedSurface build " << msSurface
                  << " = constructor/boundary builds " << msBoundaryBuild << " + addFaces "
                  << (msSurface - msBoundaryBuild)
                  << ", oriented read-out " << msRead << "), draw " << msDraw
                  << ", exact naming " << msName << "\n";
    return 0;
}
