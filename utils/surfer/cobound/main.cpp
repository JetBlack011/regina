//
//  main.cpp
//
//  cobound: bounding the smooth 4-genus by cobordisms.
//
//    cobound run       [--config FILE]... [--set key=value]...
//    cobound solve     [--config FILE]... [--set key=value]...
//    cobound sign      [--config FILE]... [--set key=value]...
//    cobound draw      [--layers N] [--gauss] [--faces [--pairsig [--sig-cache DIR]]] '<incoming PD>'
//    cobound name      --knots CSV --links CSV [--symmetry CSV] [namer limits] [--profile]
//                      [--reference]
//    cobound meridians sig|dump|dump-link|dump-subset|slope
//    cobound help [keys]
//
//  See driver/commands.h for what each command does and replaces, and
//  driver/config.h for the configuration (`cobound help keys` lists every
//  key, its type, its default in each command, and the option it replaces).
//

#include <iostream>
#include <string>
#include <vector>

#include "cobound/driver/commands.h"
#include "cobound/driver/config.h"

namespace {

void help(std::ostream &out) {
  out << "Usage: cobound <command> [arguments]\n\n"
         "  run        search from each target: without a goal each once (was\n"
         "             verifyslicegenus), with one (goal_genus or goal_lower) outwards\n"
         "             from the target until the goal has a proof (was cascadesearch)\n"
         "  solve      re-derive every verdict from the database (was verifyslicegenus\n"
         "             --solve-only)\n"
         "  sign       sign a work directory's pending cobordisms into the database (was\n"
         "             cascadesearch --sign-only)\n"
         "  draw       stored cobordisms' outgoing links, drawn (was farsidediagram)\n"
         "  name       stored cobordisms' outgoing links, named (was farsidename)\n"
         "  meridians  sig | dump | dump-link | dump-subset | slope (was peripheral_slopes)\n"
         "  help keys  every configuration key\n\n"
         "run, solve and sign read their configuration from --config FILE (any number;\n"
         "`key = value` lines) and --set key=value (any number; the last wins). A run\n"
         "writes the configuration it ran with to <work>/cobound.conf.\n";
}

void helpKeys(std::ostream &out) {
  static const char *const types[] = {"flag",  "integer", "real",  "threads",
                                      "text",  "path",    "paths", "choice"};
  for (const config::Key &k : config::schema()) {
    out << k.name << " (" << types[static_cast<int>(k.type)];
    if (k.type == config::Type::choice) {
      out << ":";
      for (const std::string &c : k.choices) out << " " << c;
    }
    out << ")\n";
    for (const std::string &old : k.oldNames)
      out << "    formerly " << old << " (still accepted)\n";
    for (const auto &[ctx, rule] : k.rules) {
      out << "    " << config::contextName(ctx) << ": ";
      switch (rule.kind) {
      case config::Rule::Kind::required: out << "required"; break;
      case config::Rule::Kind::unset: out << "unset by default"; break;
      case config::Rule::Kind::value: out << "default " << rule.value; break;
      }
      out << "\n";
    }
    out << "    replaces " << k.replaces << "\n    " << k.doc << "\n";
  }
}

} // namespace

int main(int argc, char **argv) {
  if (argc < 2) {
    help(std::cerr);
    return 2;
  }
  const std::string command = argv[1];
  const std::vector<std::string> args(argv + 2, argv + argc);
  if (command == "run") return commands::run(args);
  if (command == "solve") return commands::solve(args);
  if (command == "sign") return commands::sign(args);
  if (command == "draw") return commands::draw(args);
  if (command == "name") return commands::name(args);
  if (command == "meridians") return commands::meridians(args);
  if (command == "help" || command == "-h" || command == "--help") {
    if (!args.empty() && args.front() == "keys") helpKeys(std::cout);
    else help(std::cout);
    return 0;
  }
  std::cerr << "cobound: unknown command '" << command << "'\n\n";
  help(std::cerr);
  return 2;
}
