// diagramiso.h
//
// Isomorphisms between oriented link diagrams given as signed Gauss data,
// with the component map they induce. See README.md, "Component maps".

#pragma once

#include <optional>
#include <vector>

#include "exactnaming/gaussdiagram.h"

namespace cascade {

/**
 * A combinatorial isomorphism between two diagrams, as a map on components.
 * Diagram `a`'s component i goes to `b`'s component `componentMap[i]`.
 * `mirrored` means b is a's mirror image: crossing signs negated and over and
 * under swapped. `reversed` means every component of b runs against a's.
 */
struct DiagramIsomorphism {
  std::vector<int> componentMap;
  bool mirrored = false;
  bool reversed = false;
};

/**
 * Whether the diagrams are the same up to relabelling crossings, permuting
 * components, and choosing each component's starting point, optionally also
 * up to mirror image and/or reversing EVERY component. An isomorphic pair of
 * diagrams is the same oriented link up to the allowed symmetries, and the
 * component map is the one that isomorphism realises.
 *
 * Reversing SOME components is never allowed: that is a different oriented
 * link in general.
 *
 * Exact and exhaustive (backtracking over rotations), for the small diagrams
 * the cascade meets. Returns the first isomorphism found; which one, when
 * there are several (symmetric diagrams), is unspecified.
 */
std::optional<DiagramIsomorphism>
findDiagramIsomorphism(const exactnaming::GaussDiagram &a,
                       const exactnaming::GaussDiagram &b,
                       bool allowMirror, bool allowReverse);

/// The diagram with every component reversed.
exactnaming::GaussDiagram reverseAll(const exactnaming::GaussDiagram &d);
/// The mirror image: signs negated, over and under swapped.
exactnaming::GaussDiagram mirrorImage(const exactnaming::GaussDiagram &d);
/// The diagram with ONE component reversed (a different oriented link in
/// general). The signs of that component's crossings with others flip.
exactnaming::GaussDiagram reverseComponent(const exactnaming::GaussDiagram &d,
                                           size_t c);

} // namespace cascade
