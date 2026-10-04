// diagramiso.h
//
// Isomorphisms between oriented link diagrams given as signed Gauss data,
// with the component map they induce. See README.md, "Component maps".

#pragma once

#include <optional>
#include <vector>

#include "linknaming/diagrams/gaussdiagram.h"

namespace linknaming {

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
 * a goal run meets. Returns the first isomorphism found; which one, when
 * there are several (symmetric diagrams), is unspecified -- unless
 * `componentMap` is given, when only isomorphisms taking a's component i to
 * b's componentMap[i] are considered. With a == b, that asks whether a
 * permutation of components is a symmetry of the link.
 */
std::optional<DiagramIsomorphism>
findDiagramIsomorphism(const linknaming::GaussDiagram &a,
                       const linknaming::GaussDiagram &b,
                       bool allowMirror, bool allowReverse,
                       const std::vector<int> *componentMap = nullptr);

/// The diagram with every component reversed.
linknaming::GaussDiagram reverseAll(const linknaming::GaussDiagram &d);
/// The mirror image: signs negated, over and under swapped.
linknaming::GaussDiagram mirrorImage(const linknaming::GaussDiagram &d);
/// The diagram with ONE component reversed (a different oriented link in
/// general). The signs of that component's crossings with others flip.
linknaming::GaussDiagram reverseComponent(const linknaming::GaussDiagram &d,
                                           size_t c);

} // namespace linknaming
