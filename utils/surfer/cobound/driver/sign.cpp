//
//  sign.cpp
//
//  cobound sign: a work directory's pending cobordisms, signed into the
//  database.
//

#include <iostream>

#include "cobound/cobordisms/pending.h"
#include "cobound/driver/commands.h"
#include "cobound/driver/config.h"
#include "cobound/solver/literature.h"

// `sign` (was cascadesearch --sign-only): every pending file under `work`
// (hop_*/kept.csv), from where it was last signed, signed and appended to
// `cobordisms` -- each cobordism whose identity neither it nor any of
// `dedupe_against` holds -- and each file's <file>.signed record set (plan
// divergence 1). What a run's own end does, for a run that was killed.
int commands::sign(const std::vector<std::string> &args) {
  try {
    const config::Config cfg = config::forCommand("sign", config::Context::sign, args);
    const std::string store = cfg.text("cobordisms");
    const solver::NameTable names =
        solver::loadTableNames(cfg.text("knot_table"), cfg.text("link_table"), "");
    const cobordisms::StoreResult s =
        cobordisms::signPending(cfg.text("work"), store, cfg.paths("dedupe_against"), names,
                             cfg.threads(), cfg.text("pair_sig_cache"));
    std::cout << "[+] witness store: " << s.kept << " kept, " << s.fresh << " new, "
              << s.appended << " appended to " << store << "\n";
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "cobound sign: " << e.what() << "\n";
    return 2;
  }
}
