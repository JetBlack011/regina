// searchjudge.h
//
// A depth-0 search's own cobordism graph, which judges its finds.

#pragma once

#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "cobound/bounds/axioms.h"
#include "cobound/bounds/cobordismgraph.h"
#include "cobound/bounds/links.h"
#include "cobound/bounds/searchcobordisms.h"
#include "cobound/outgoing/fromdatabase.h"
#include "linknaming/linknamer.h"
#include "linknaming/tables.h"

namespace bounds {

/**
 * The cobordism graph of one search without a goal (plan divergences 6 and
 * 10): it holds the searched link, with its literature lower bound, and each
 * of the search's finds as it is kept -- its outgoing link interned and given
 * its outside facts (NodeAxioms) exactly as a goal run's graph does. From
 * those and the literature alone it judges what the finds prove.
 *
 *   - A contradiction: some derived bound below a proved lower bound (the
 *     graph's gates, ProofGraph::contradictions()). Something the search or
 *     the naming computed is wrong, and the run halts.
 *   - Whether the searched link now has a constructive bound at its
 *     literature lower bound: a proof resting on no literature value.
 *
 * No database is read: a link constructive only through another search's
 * cobordism is not judged so here (plan divergence 10, intended).
 *
 * The search runs in row().rowBuild(), the row this judge reads its finds'
 * outgoing links from (search::buildRow(pd, layers, layers)).
 */
class SearchJudge {
public:
  /// \param name the searched link's name; \param pd its row's PD, as
  /// searched; \param literatureLo its literature lower bound.
  /// \throws std::runtime_error when the row's own link does not redraw as
  /// its diagram (HopAssembler's certification).
  SearchJudge(const std::string &name, const std::string &pd, int layers, int literatureLo,
              const linknaming::ExactTables &tables, const linknaming::ExactNamer &namer,
              const linknaming::SymmetryTable &symmetries, unsigned threads);
  SearchJudge(const SearchJudge &) = delete;
  SearchJudge &operator=(const SearchJudge &) = delete;

  /// The row the search runs in, and reads its finds' outgoing links from.
  const outgoing::OutgoingReader &reader() const { return assembler_->redrawer(); }

  struct Verdict {
    /// The graph's contradictions so far (each a sentence); none, normally.
    std::vector<std::string> contradictions;
    /// The searched link's least connected genus, when a proof resting on
    /// no literature value reaches its literature lower bound.
    std::optional<int> constructive;
  };

  /// One find, as the search kept it: its oriented outgoing link on row(),
  /// its genus (tubed), and a key naming it. Not thread-safe: the search
  /// calls it under its own lock. A find the graph cannot take (it breaks
  /// an invariant of its reading) is counted in failures() and judged as
  /// absent.
  Verdict add(const outgoing::OutgoingLink &link, int genus, const std::string &key);

  long long finds() const { return finds_; }
  long long failures() const { return failures_; }
  const ProofGraph &graph() const { return g_; }

private:
  ProofGraph g_;
  NodeRegistry reg_;
  NodeAxioms axioms_;
  std::unique_ptr<CobordismAssembler> assembler_;
  NodeId target_ = -1;
  int components_ = 1;
  int literatureLo_ = 0;
  long long finds_ = 0, failures_ = 0;
};

} // namespace bounds
