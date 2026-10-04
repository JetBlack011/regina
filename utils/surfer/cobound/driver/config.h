//
//  config.h
//
//  cobound's configuration: one schema for every command.
//

#ifndef SURFER_COBOUND_CONFIG_H
#define SURFER_COBOUND_CONFIG_H

#include <map>
#include <optional>
#include <ostream>
#include <stdexcept>
#include <string>
#include <vector>

/*! \file utils/surfer/cobound/driver/config.h
 *  \brief The one schema of cobound's configuration, its file format, and
 *  the one parser of every command's arguments.
 *
 *  A configuration is a set of `key = value` assignments, read from config
 *  files (`--config FILE`, any number, in order) and then `--set key=value`
 *  (any number, last wins). Every key has a type, the commands it applies
 *  to, and in each of those a default, or no default (unset), or none at
 *  all (required). Defaults reproduce the reference behaviour: where the
 *  two retired drivers differed, a key's default without a goal is
 *  verifyslicegenus's and with a goal cascadesearch's (schema() lists both;
 *  config_test checks every one).
 *
 *  The file format: one assignment per line, `key = value`, the value
 *  running to the end of the line with surrounding spaces trimmed (a PD
 *  code may hold spaces); blank lines and lines starting with `#` are
 *  skipped. `none` unsets a key that may be unset; a list (`paths`) is
 *  comma-separated; a flag is 0/1 (or true/false, yes/no, on/off).
 *
 *  `run` writes the configuration it ran with to <work>/cobound.conf
 *  (writeEffective()): every key it reads, its value, and where the value
 *  came from.
 */

namespace config {

/** What a command runs as, as far as keys go. `run` is goal-directed when
 *  goal_genus or goal_lower is set (runContext()). */
enum class Context { run, goal, solve, sign, name, draw };

/** "run", "run with a goal", "solve", ... */
const char *contextName(Context c);

enum class Type {
    flag,    ///< 0 or 1
    integer, ///< a whole number
    real,    ///< a number
    threads, ///< a whole number > 0, or `auto` (the machine's hardware threads)
    text,    ///< anything
    path,    ///< a file or directory
    paths,   ///< comma-separated paths
    choice,  ///< one of Key::choices
};

/** A key's rule in one context. */
struct Rule {
    enum class Kind {
        unset,    ///< no default: unset unless assigned (`none` unsets it)
        value,    ///< `value` unless assigned
        required, ///< must be assigned
    };
    Kind kind = Kind::unset;
    std::string value;
};

struct Key {
    std::string name;
    /** Earlier spellings of the key, still accepted in config files and by
     *  --set (phase 6 renamed the keys that held a retired term: phase-5
     *  cobound.conf files spell them so). `cobound help keys` lists them; the
     *  effective configuration writes `name`. */
    std::vector<std::string> oldNames;
    Type type = Type::text;
    std::vector<std::string> choices; ///< Type::choice
    /** The contexts the key applies to, with its rule in each. */
    std::vector<std::pair<Context, Rule>> rules;
    /** The retired options it replaces, as "--x (driver)". */
    std::string replaces;
    std::string doc;

    const Rule *rule(Context c) const;
};

/** Every key, in the order the effective configuration lists them. */
const std::vector<Key> &schema();
/** The key named `name` (or spelled so formerly: Key::oldNames), or nullptr. */
const Key *findKey(const std::string &name);

/** A configuration or command line that cannot be used. */
struct Error : std::runtime_error {
    using std::runtime_error::runtime_error;
};

/** One assignment, and where it came from ("--set", "<file>:<line>", or a
 *  command's own flag). */
struct Assignment {
    std::string key, value, source;
};

/**
 * A command line: `--config FILE` and `--set key=value` read into
 * assignments (files first, in order, then every --set in order, so --set
 * wins), each of `aliases`' flags into its key, and the rest positional.
 * `draw` and `name` keep their reference tools' flags as aliases: they are
 * stdin tools the atlas and its checker (cascade_check.py) drive.
 */
struct Alias {
    std::string flag;  ///< "--knots"
    std::string key;   ///< "knot_table"
    bool takesValue = true;
    std::string value; ///< the value assigned when !takesValue ("1")
};
struct CommandLine {
    std::vector<Assignment> assignments;
    std::vector<std::string> positional;
};
CommandLine parseCommandLine(const std::vector<std::string> &args,
                             const std::vector<Alias> &aliases = {});

/** `run`'s context: `goal` when goal_genus or goal_lower is assigned (and
 *  not `none`), else `run`. */
Context runContext(const std::vector<Assignment> &assignments);

class Config {
  public:
    /**
     * Every key applying to `context`: its last assignment, or its default.
     * \exception Error an unknown key, a value of the wrong type, or a
     *            required key unassigned. A known key that does not apply to
     *            `context` is reported in ignored() and otherwise ignored.
     */
    Config(Context context, const std::vector<Assignment> &assignments);

    Context context() const { return context_; }
    /** Whether `key` has a value (a default, or assigned and not `none`). */
    bool has(const std::string &key) const;
    bool flag(const std::string &key) const;
    long long integer(const std::string &key) const;
    std::optional<long long> optionalInteger(const std::string &key) const;
    double real(const std::string &key) const;
    std::optional<double> optionalReal(const std::string &key) const;
    /** A text or path value; "" when unset. */
    std::string text(const std::string &key) const;
    std::optional<std::string> optionalText(const std::string &key) const;
    std::vector<std::string> paths(const std::string &key) const;
    /** threads: the value, or the machine's hardware threads for `auto`. */
    unsigned threads(const std::string &key = "threads") const;

    /** Keys assigned that do not apply to this context (with their source). */
    const std::vector<Assignment> &ignored() const { return ignored_; }

    /** Every key of this context, `key = value  # source`, in schema
     *  order, under a header naming the context. */
    void writeEffective(std::ostream &out, const std::string &command) const;

  private:
    struct Value {
        std::optional<std::string> text; ///< nullopt: unset
        std::string source;              ///< "default", "--set", "<file>:<line>"
        std::string spelling;            ///< the old name it was given as, if one
    };
    const Value &value_(const std::string &key, Type expected) const;
    Context context_;
    std::map<std::string, Value> values_;
    std::vector<Assignment> ignored_;
};

/** Parses a flag's text (0/1, true/false, yes/no, on/off); nullopt if none. */
std::optional<bool> parseFlag(const std::string &text);

/**
 * A command's configuration from its arguments: parseCommandLine(`args`,
 * `aliases`), then Config in `context` (unset: runContext() of what was
 * assigned). Positional arguments go to `positional`, or are refused when it
 * is null. Each assigned key that does not apply is reported on stderr.
 * \exception Error as parseCommandLine() and Config().
 */
Config forCommand(const std::string &command, std::optional<Context> context,
                  const std::vector<std::string> &args, const std::vector<Alias> &aliases = {},
                  std::vector<std::string> *positional = nullptr);

} // namespace config

#endif // SURFER_COBOUND_CONFIG_H
