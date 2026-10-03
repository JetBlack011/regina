#include "surfer/enumeration/searchstrategy.h"

#include "surfer/report/atomicwrite.h"

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <istream>
#include <sstream>
#include <stdexcept>

// The file format, one record per line:
//
//   surfer-search-frontier <1|2>
//   fingerprint <hex>
//   round <r> <R> cap <c|-> suppress_below <s> deepest_exhausted <e|-> complete <0|1>
//   runs <n> satisfying <S> found <F> attempts <A>
//   pending <bytes> <path>                      (format 2 only; the path relative
//                                               to the frontier file's directory)
//   roots <N>
//   r <idx> <done> <level> <deepest> <levels> {<childIndex> <child> <npruned> <pruned>...}
//   ...
//   end
//
// Only roots with something to say (done, a pass had, or a position) get an
// `r` line; the rest start from scratch. A complete frontier lists none.

namespace {

const char *kMagic = "surfer-search-frontier";
// 1: no pending line; 2: a pending line. A frontier without one is written
// as format 1, byte for byte as before.
constexpr int kFormatPlain = 1, kFormatPending = 2;

std::string optionalField(const std::optional<long long> &v) {
    return v ? std::to_string(*v) : std::string("-");
}

std::optional<long long> parseOptional(const std::string &s) {
    if (s == "-")
        return std::nullopt;
    return std::stoll(s);
}

[[noreturn]] void bad(const std::string &what) {
    throw std::runtime_error("SearchFrontier: " + what);
}

// Reads the next line and checks its leading keyword.
std::istringstream expectLine(std::istream &in, const std::string &keyword) {
    std::string line;
    if (!std::getline(in, line))
        bad("expected '" + keyword + "', found end of file");
    std::istringstream fields(line);
    std::string word;
    fields >> word;
    if (word != keyword)
        bad("expected '" + keyword + "', found '" + line + "'");
    return fields;
}

} // namespace

size_t SearchFrontier::rootsDone() const {
    size_t n = 0;
    for (const Root &r : roots)
        n += r.done;
    return n;
}

size_t SearchFrontier::rootsStarted() const {
    size_t n = 0;
    for (const Root &r : roots)
        n += !r.done && (r.level > 0 || !r.position.empty());
    return n;
}

unsigned SearchFrontier::maxLevel() const {
    unsigned m = 0;
    for (const Root &r : roots)
        m = std::max(m, r.level);
    return m;
}

std::string SearchFrontier::summary() const {
    std::ostringstream s;
    s << "round " << round << "/" << rounds << " cap "
      << (cap ? std::to_string(*cap) : std::string("none")) << ", roots "
      << rootsDone() << " done + " << rootsStarted() << " part-walked of "
      << roots.size() << ", max level " << maxLevel() << "; runs " << runs
      << ", satisfying " << satisfying << ", found " << found
      << ", attempts " << attempts;
    if (complete)
        s << ", complete";
    return s.str();
}

void SearchFrontier::write(std::ostream &out) const {
    out << kMagic << ' ' << (pending ? kFormatPending : kFormatPlain) << '\n'
        << "fingerprint " << fingerprint << '\n'
        << "round " << round << ' ' << rounds << " cap " << optionalField(cap)
        << " suppress_below " << suppressBelow << " deepest_exhausted "
        << optionalField(deepestExhausted) << " complete " << (complete ? 1 : 0)
        << '\n'
        << "runs " << runs << " satisfying " << satisfying << " found " << found
        << " attempts " << attempts << '\n';
    if (pending)
        out << "pending " << pending->bytes << ' ' << pending->path << '\n';
    out << "roots " << roots.size() << '\n';
    if (!complete)
        for (size_t i = 0; i < roots.size(); ++i) {
            const Root &r = roots[i];
            if (!r.done && r.level == 0 && r.position.empty())
                continue;
            out << "r " << i << ' ' << (r.done ? 1 : 0) << ' ' << r.level << ' '
                << r.position.deepest << ' ' << r.position.levels.size();
            for (const auto &l : r.position.levels) {
                out << ' ' << l.childIndex << ' ' << l.child << ' '
                    << l.pruned.size();
                for (int p : l.pruned)
                    out << ' ' << p;
            }
            out << '\n';
        }
    out << "end\n";
}

SearchFrontier SearchFrontier::read(std::istream &in) {
    SearchFrontier f;
    int format = 0;
    {
        auto fields = expectLine(in, kMagic);
        if (!(fields >> format) || (format != kFormatPlain && format != kFormatPending))
            bad("unsupported format");
    }
    {
        auto fields = expectLine(in, "fingerprint");
        if (!(fields >> f.fingerprint))
            bad("no fingerprint");
    }
    {
        auto fields = expectLine(in, "round");
        std::string k1, capField, k2, k3, deepest, k4;
        int complete = 0;
        if (!(fields >> f.round >> f.rounds >> k1 >> capField >> k2 >>
              f.suppressBelow >> k3 >> deepest >> k4 >> complete) ||
            k1 != "cap" || k2 != "suppress_below" || k3 != "deepest_exhausted" ||
            k4 != "complete")
            bad("malformed round line");
        f.cap = parseOptional(capField);
        f.deepestExhausted = parseOptional(deepest);
        f.complete = complete != 0;
    }
    {
        auto fields = expectLine(in, "runs");
        std::string k1, k2, k3;
        if (!(fields >> f.runs >> k1 >> f.satisfying >> k2 >> f.found >> k3 >>
              f.attempts) ||
            k1 != "satisfying" || k2 != "found" || k3 != "attempts")
            bad("malformed runs line");
    }
    if (format == kFormatPending) {
        auto fields = expectLine(in, "pending");
        Pending p;
        if (!(fields >> p.bytes))
            bad("malformed pending line");
        std::getline(fields >> std::ws, p.path);
        if (p.path.empty())
            bad("malformed pending line");
        f.pending = std::move(p);
    }
    {
        auto fields = expectLine(in, "roots");
        size_t n = 0;
        if (!(fields >> n))
            bad("malformed roots line");
        f.roots.assign(n, Root{});
        if (f.complete)
            for (Root &r : f.roots)
                r.done = true;
    }
    std::string line;
    while (std::getline(in, line)) {
        if (line == "end")
            return f;
        std::istringstream fields(line);
        std::string tag;
        size_t idx = 0, levels = 0;
        int done = 0;
        Root r;
        if (!(fields >> tag >> idx >> done >> r.level >> r.position.deepest >>
              levels) ||
            tag != "r" || idx >= f.roots.size())
            bad("malformed root line '" + line + "'");
        r.done = done != 0;
        r.position.levels.resize(levels);
        for (auto &l : r.position.levels) {
            size_t pruned = 0;
            if (!(fields >> l.childIndex >> l.child >> pruned))
                bad("malformed position in '" + line + "'");
            l.pruned.resize(pruned);
            for (int &p : l.pruned)
                if (!(fields >> p))
                    bad("malformed pruned list in '" + line + "'");
        }
        f.roots[idx] = std::move(r);
    }
    bad("no 'end' line: the file is truncated");
}

namespace {

namespace fs = std::filesystem;

// The directory a frontier file at `path` lives in, absolute and lexically
// normal.
fs::path frontierDirectory(const std::string &path) {
    return fs::absolute(fs::path(path)).lexically_normal().parent_path();
}

// The pending path as save() records it: relative to the frontier's
// directory, so a work tree that is packed, synced or copied keeps its
// resume; the file's own name when no relative path exists.
std::string recordedPending(const std::string &pending, const std::string &frontierPath) {
    const fs::path file = fs::absolute(fs::path(pending)).lexically_normal();
    const fs::path rel = file.lexically_relative(frontierDirectory(frontierPath));
    return rel.empty() ? file.filename().string() : rel.string();
}

// The pending file a loaded frontier names: a relative record resolved
// against the frontier's directory (an absolute one, as 4(b)'s format 2
// wrote it, as it stands); when that file does not exist but one of its
// name lies beside the frontier, that one.
std::string resolvedPending(const std::string &recorded, const std::string &frontierPath) {
    const fs::path p(recorded);
    const fs::path dir = frontierDirectory(frontierPath);
    fs::path file = p.is_absolute() ? p : (dir / p).lexically_normal();
    std::error_code ec;
    if (!fs::exists(file, ec)) {
        const fs::path beside = dir / p.filename();
        if (fs::exists(beside, ec))
            file = beside;
    }
    return file.string();
}

} // namespace

void SearchFrontier::save(const std::string &path) const {
    try {
        if (pending) {
            SearchFrontier recorded = *this;
            recorded.pending->path = recordedPending(pending->path, path);
            report::atomicWrite(path, [&recorded](std::ostream &out) { recorded.write(out); },
                                report::Durability::cache);
        } else {
            report::atomicWrite(path, [this](std::ostream &out) { write(out); },
                                report::Durability::cache);
        }
    } catch (const std::runtime_error &e) {
        bad(e.what());
    }
}

std::optional<SearchFrontier> SearchFrontier::load(const std::string &path) {
    std::ifstream in(path);
    if (!in)
        return std::nullopt;
    SearchFrontier f = read(in);
    if (f.pending)
        f.pending->path = resolvedPending(f.pending->path, path);
    return f;
}
