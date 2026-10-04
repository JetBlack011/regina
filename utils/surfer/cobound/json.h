//
//  json.h
//
//  The one JSON writer: escaping, and the arrays the run records hold.
//

#ifndef SURFER_COBOUND_JSON_H
#define SURFER_COBOUND_JSON_H

#include <sstream>
#include <string>

/*! \file utils/surfer/cobound/json.h
 *  \brief What every JSON line cobound writes is built from (cascade.jsonl,
 *  node_bounds.jsonl, the certificates, the cobordism graph's partition
 *  genera (profiles.jsonl), `cobound draw`'s `incoming=` field). The run directory's formats are
 *  frozen, so each helper writes exactly what its callers wrote before.
 */

namespace json {

/**
 * `s` escaped for a JSON string: `"` and `\` get a backslash and a newline
 * becomes `\n`. Nothing else is escaped: no frozen file holds another control
 * character, and changing the escaping would change those files.
 */
inline std::string escape(const std::string &s) {
    std::string o;
    o.reserve(s.size());
    for (char c : s) {
        if (c == '"' || c == '\\') o += '\\';
        if (c == '\n') {
            o += "\\n";
            continue;
        }
        o += c;
    }
    return o;
}

/** `"` + escape(s) + `"`. */
inline std::string quote(const std::string &s) { return '"' + escape(s) + '"'; }

/** A sequence of numbers as `[a,b,c]` (no spaces). */
template <typename Seq>
std::string array(const Seq &v) {
    std::ostringstream o;
    o << '[';
    bool first = true;
    for (const auto &x : v) {
        if (!first) o << ',';
        o << x;
        first = false;
    }
    o << ']';
    return o.str();
}

/** A sequence of sequences of numbers as `[[a,b],[c,d]]`. */
template <typename Rows>
std::string matrix(const Rows &m) {
    std::string o = "[";
    bool first = true;
    for (const auto &row : m) {
        if (!first) o += ',';
        o += array(row);
        first = false;
    }
    return o + ']';
}

} // namespace json

#endif // SURFER_COBOUND_JSON_H
