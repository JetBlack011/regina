"""
enumerate_surfaces.py
======================

Enumerate all connected, closed surfaces (2-manifolds) embedded in the
2-skeleton of a triangulated compact 4-manifold.

TOPOLOGICAL BACKGROUND
-----------------------
Given a triangulation T of a 4-manifold, its 2-skeleton is the set of all
triangles (2-simplices) together with the edges and vertices they touch.
We want to find every subset S of triangles such that S, as a simplicial
complex in its own right, is a closed 2-manifold (i.e. a disjoint union of
closed surfaces -- spheres, tori, projective planes, Klein bottles, genus-g
surfaces, etc., possibly non-orientable).

A subset S of triangles is the triangle-set of a closed 2-manifold iff BOTH
of the following local conditions hold:

  (1) EDGE CONDITION. Every edge of T lies in exactly 0 or exactly 2 of the
      selected triangles. (In a closed surface, every edge is shared by
      exactly two faces; an edge with 1 or with >=3 incident selected
      triangles is a boundary edge or a non-manifold "book" edge, both of
      which are forbidden.)

  (2) VERTEX CONDITION. At every vertex v of T, look at the "link" of v
      restricted to the selected triangles, and ask whether it forms a
      valid 2-manifold neighborhood there (no two or more separate
      sheets of surface passing through the same point -- a "pinch
      point" / non-manifold singularity, which we must forbid). See
      ENFORCING THE VERTEX CONDITION below for how this is actually
      checked -- it turns out to be subtler than a simple graph-theoretic
      condition on v's link whenever the triangulation reuses vertices
      heavily (see below).

This is encoded as a Boolean satisfiability problem and solved with an
All-SAT loop in Z3, then the (possibly disconnected) solutions are filtered
down to genuinely connected surfaces with a BFS over triangle-adjacency.

ENFORCING THE VERTEX CONDITION: LAZY VS. EAGER
-----------------------------------------------
Two ways to enforce condition (2) are implemented, selected via
`mode="lazy"` (default) or `mode="eager"`:

  * EAGER (`build_vertex_constraints`): for every vertex v, precompute
    every simple cycle in v's link graph (nodes = edges of T incident to
    v, arcs = triangles incident to v), and add a Z3 clause up front
    forbidding each node-disjoint pair of cycles from being
    simultaneously selected. Sound *only* when no triangle is
    self-identified (find_self_identified_triangles(tri) is empty -- see
    below); the enumeration is also worst-case super-exponential in the
    number of distinct edges touching a vertex, so a single
    densely-linked vertex can hang or exhaust memory even in an
    otherwise tiny triangulation.

  * LAZY (`is_valid_closed_surface`, the default): add only the edge
    constraints up front, then let Z3 propose candidates freely. For
    each candidate, build the actual standalone abstract
    regina.Triangulation2 it would produce (see
    build_surface_triangulation2) and ask REGINA DIRECTLY whether it's a
    valid closed 2-manifold. A candidate that fails is rejected (a plain
    "don't propose this exact assignment again" cut -- deliberately not
    a cleverer generalizing cut, see below) and the solver is re-queried.
    This never enumerates a vertex link's cycle set, so it sidesteps
    eager mode's performance cliff entirely, AND it has no known
    soundness caveats, for reasons explained next.

  WHY NOT JUST CHECK CONNECTIVITY OF EACH VERTEX'S ACTIVE LINK GRAPH?
  An earlier version of lazy mode did exactly that: BFS the selected
  edges of each vertex's link graph and reject if it split into more
  than one connected component (see find_vertex_violations, kept in this
  file for reference). That is WRONG in both directions once a
  triangulation has self-identified triangles (two or more of one
  triangle's own vertices coinciding at the same ambient vertex -- common
  in Regina's minimal/efficient example triangulations, which
  deliberately reuse few vertices):
    - FALSE NEGATIVE (eager mode's failure, same root cause): a single
      self-identified triangle can make one genuinely-connected piece
      look, under a naive decomposition, like two separate pieces that a
      third triangle later bridges -- confirmed on
      regina.Example4.ballBundle(), where eager mode finds only 6 of the
      12 genuinely valid connected surfaces.
    - FALSE POSITIVE (the ambient-vertex BFS approach's own failure):
      two components that look disjoint in one ambient vertex's combined
      link graph do NOT always indicate a pinch point -- they can be two
      entirely separate abstract vertices that simply happen to coincide
      at the same ambient point, which is perfectly valid. Confirmed on
      regina.Example4.s2xs2(), triangles {2, 3}: the ambient-vertex BFS
      flags a "violation", yet Regina confirms this is a valid sphere
      with 3 abstract vertices -- the two flagged pieces are each their
      own separate, valid, degree-1 vertex. Cutting on this false
      positive would have permanently excluded a valid surface.
  The fundamental issue: whether two triangle-corners belong to the same
  abstract vertex, or to two different abstract vertices that happen to
  coincide at the same ambient vertex, depends on the FULL transitive
  gluing structure -- not on simple connectivity of one ambient vertex's
  arc graph. Rather than re-deriving that transitive structure ourselves
  (and risking a third subtle bug), is_valid_closed_surface just asks
  Regina, whose own vertex-link machinery already gets this right.

  Because Regina's check can occasionally say "invalid" for reasons that
  don't decompose neatly into "these specific triangles are the
  problem", all_sat_enumerate uses only the provably-safe "block this
  exact assignment" cut for rejected candidates, not a cleverer
  generalizing one -- see find_vertex_violations's docstring for the two
  counterexamples above showing why a generalizing cut based on ambient
  vertex structure cannot be trusted here. This trades away some
  convergence speed for guaranteed correctness; the correctness test
  suite (test_enumerate_surfaces.py) cross-checks lazy mode against
  brute-force ground truth on many triangulations (including both of the
  counterexamples above) to confirm it.

DEPENDENCIES
------------
    pip install regina z3-solver networkx

Author: (generated)
"""

import argparse
import sys
from collections import Counter, defaultdict, deque
from itertools import combinations

import networkx as nx
import z3

try:
    import regina
except ImportError:  # pragma: no cover
    regina = None


# ---------------------------------------------------------------------------
# 1. Variable initialization
# ---------------------------------------------------------------------------

def build_variables(tri):
    """
    Create one Z3 Boolean variable x_i per triangle (2-simplex) of the
    4-manifold triangulation `tri`. x_i == True means triangle i is part
    of the candidate surface.

    Returns
    -------
    list[z3.BoolRef]
        xs[i] is the indicator variable for triangle i (indexed the same
        way as tri.triangle(i).index()).
    """
    n = tri.countTriangles()
    return [z3.Bool(f"x_{i}") for i in range(n)]


def forbid_empty_surface(xs, solver):
    """At least one triangle must be selected -- forbid the trivial/empty
    surface, which trivially satisfies every local constraint below."""
    solver.add(z3.Sum([z3.If(v, 1, 0) for v in xs]) >= 1)


# ---------------------------------------------------------------------------
# 2. The Edge Condition (closedness)
# ---------------------------------------------------------------------------

def compute_edge_to_triangles(tri):
    """
    Build a map: edge_index -> list of triangle indices incident to it.

    Regina's face-numbering convention for a Triangle<4> `t` is that
    `t.edge(i)` is the edge OPPOSITE `t.vertex(i)`, i.e. the edge that
    does *not* touch local vertex i. Iterating i in {0,1,2} therefore
    visits all three edges of the triangle exactly once.
    """
    edge_to_tris = defaultdict(list)
    for t in tri.triangles():
        ti = t.index()
        for local_e in range(3):
            e = t.edge(local_e)
            edge_to_tris[e.index()].append(ti)
    return edge_to_tris


def build_edge_constraints(tri, xs, solver):
    """
    For every edge e of the triangulation, the number of *selected*
    triangles containing e must be exactly 0 (the edge is not touched by
    the surface at all) or exactly 2 (the edge is a normal interior edge
    of the surface, shared by precisely two faces, as in any closed
    2-manifold). Any other count (1, 3, 4, ...) would make e a boundary
    edge or a branch/"book" edge, neither of which is allowed in a closed
    manifold surface.
    """
    edge_to_tris = compute_edge_to_triangles(tri)
    for _eidx, tri_list in edge_to_tris.items():
        incident_vars = [xs[i] for i in tri_list]
        count = z3.Sum([z3.If(v, 1, 0) for v in incident_vars])
        solver.add(z3.Or(count == 0, count == 2))
    return edge_to_tris


# ---------------------------------------------------------------------------
# 3. The Vertex Condition (manifold verification at each vertex)
# ---------------------------------------------------------------------------

def build_vertex_link_graphs(tri):
    """
    For every vertex v of the 4-manifold, build a multigraph L_v (the
    "link graph restricted to candidate surface triangles"):

        * NODES of L_v are the edges of the 4-manifold incident to v.
        * Each TRIANGLE t incident to v contributes exactly one edge to
          L_v, connecting the two edges-of-t that touch v, labelled with
          t's own index (so we can look up its Z3 variable later).

    Why: a triangle t = {v, a, b} incident to v is "seen" locally at v as
    a connecting arc between the two 4-manifold-edges {v,a} and {v,b}.
    If triangle t is selected, that arc is present in the local picture of
    the surface at v. The edge condition guarantees each node of L_v has
    degree 0 or 2 among the *selected* arcs, so the selected sub-graph of
    L_v is always a disjoint union of simple cycles: exactly what a
    2-manifold's vertex link should look like (a disjoint union of
    circles), and we still need to make sure there's at most 1 of them.

    Returns
    -------
    dict[int, networkx.MultiGraph]
        vertex_index -> L_v
    """
    link_graphs = {v.index(): nx.MultiGraph() for v in tri.vertices()}

    for t in tri.triangles():
        ti = t.index()
        for j in range(3):
            v = t.vertex(j)
            others = [k for k in range(3) if k != j]
            e1 = t.edge(others[0]).index()
            e2 = t.edge(others[1]).index()
            # The multigraph key must be unique per (u, v) node pair, not
            # just globally-ish unique. Using key=ti alone is NOT safe: in
            # extremely degenerate/minimal triangulations (very few global
            # vertices, heavy self-identification -- e.g. a triangle all
            # three of whose vertices coincide at one global vertex), two
            # DIFFERENT corners (j values) of the SAME triangle can land
            # on the exact same (e1, e2) node pair. Reusing key=ti for
            # both would make networkx silently overwrite the first arc
            # instead of keeping both -- so we key on (ti, j) instead,
            # which is always distinct per corner, and keep the actual
            # triangle index in the 'var' attribute where it's needed.
            link_graphs[v.index()].add_edge(e1, e2, key=(ti, j), var=ti)

    return link_graphs


def enumerate_link_cycles(graph):
    """
    Enumerate all simple cycles in a (small) multigraph via exhaustive
    DFS backtracking. We roll our own instead of relying purely on
    networkx.simple_cycles for two reasons:

      * we must correctly handle *parallel edges* (two distinct triangles
        connecting the same pair of link-graph nodes), which show up as
        valid length-2 cycles ("bigons" -- e.g. two triangles glued along
        both shared edges, as in a pillow/lens-shaped piece of surface);
      * we must correctly handle *self-loops* (a single triangle whose
        own two edges at this vertex coincide -- possible in highly
        degenerate/minimal triangulations with very few global vertices),
        which are a complete, valid 1-node cycle all by themselves and
        are registered directly rather than via the general DFS (which
        never revisits its own start node at path-length 1).

    Parameters
    ----------
    graph : networkx.MultiGraph
        Nodes are 4-manifold edge indices, multi-edges are keyed by
        triangle index with a 'var' attribute equal to that same index.

    Returns
    -------
    list[tuple[frozenset[int], frozenset[int]]]
        A list of (node_set, triangle_var_set) pairs, one per distinct
        simple cycle found (deduplicated by node/variable-set identity,
        independent of traversal direction or start point).
    """
    adj = defaultdict(list)  # node -> list of (neighbor, edge_key, var)
    for u, v, key, data in graph.edges(keys=True, data=True):
        adj[u].append((v, key, data["var"]))
        adj[v].append((u, key, data["var"]))

    found = set()  # dedup by (node_set, var_set)

    # Self-loops (u == v -- a triangle whose own two edges at this vertex
    # coincide, which happens in highly degenerate/minimal triangulations)
    # are a valid, complete 1-node cycle all by themselves: that single
    # triangle already saturates the node's degree to 2. The general DFS
    # below never revisits `start` at path-length 1, so it can never
    # "close" a pure self-loop into a cycle -- we register these directly.
    for u, v, key, data in graph.edges(keys=True, data=True):
        if u == v:
            found.add((frozenset({u}), frozenset({data["var"]})))

    def dfs(start, current, path_nodes, path_node_set, path_vars, used_keys):
        for nbr, key, var in adj[current]:
            if key in used_keys:
                continue
            if nbr == start and len(path_nodes) >= 2:
                cycle_nodes = frozenset(path_node_set)
                cycle_vars = frozenset(path_vars + [var])
                found.add((cycle_nodes, cycle_vars))
            elif nbr not in path_node_set:
                dfs(
                    start,
                    nbr,
                    path_nodes + [nbr],
                    path_node_set | {nbr},
                    path_vars + [var],
                    used_keys | {key},
                )

    for start in graph.nodes():
        dfs(start, start, [start], {start}, [], set())

    return list(found)


def find_self_identified_triangles(tri):
    """
    Return the indices of all triangles that have a repeated vertex label
    among their own 3 vertices (i.e. two or more of a single triangle's
    corners map to the *same* ambient vertex).

    This matters because it is exactly the condition under which EAGER
    mode's pairwise disjoint-cycle exclusion (build_vertex_constraints)
    can become UNSOUND. A self-identified triangle contributes more than
    one corner-arc to the *same* vertex's link graph, which can make a
    single link-graph node have selected-degree greater than 2 -- and a
    graph with degree->2 nodes cannot be soundly analyzed by decomposing
    it into simple cycles and forbidding node-disjoint pairs: two "cycles"
    found in isolation can be bridged into one genuinely connected
    component by a *third* triangle's arcs, which the pairwise-exclusion
    approach has no way to account for (it reasons about each pair of
    cycles independently of everything else that might be selected).
    Concretely, this can make EAGER mode incorrectly reject a valid
    closed surface -- confirmed by brute-force cross-check (see the test
    suite) on regina.Example4.ballBundle(), where EAGER finds only 6 of
    the 12 genuinely valid connected surfaces.

    LAZY mode has no such issue: it checks the *actual* connectivity of
    the concrete active subgraph in each candidate model directly (via
    BFS), which remains correct regardless of node degree.

    This is common in Regina's minimal/efficient example triangulations
    (which deliberately reuse few vertices via heavy identification) --
    in a survey of Regina's built-in Example4 constructions, 13 of 18 had
    at least one self-identified triangle.

    Returns
    -------
    list[int]
        Indices of every triangle with a repeated own-vertex label.
        Empty list means EAGER mode is safe to use on this triangulation.
    """
    flagged = []
    for t in tri.triangles():
        labels = [t.vertex(i).index() for i in range(3)]
        if len(set(labels)) < 3:
            flagged.append(t.index())
    return flagged


def build_vertex_constraints(tri, xs, solver, link_graphs=None):
    """
    EAGER mode. For every vertex v, enumerate the simple cycles of its
    link graph L_v and forbid any pair of *node-disjoint* cycles from
    being fully selected simultaneously.

    The reasoning behind this constraint: IF every node of L_v has
    selected-degree at most 2 (guaranteed by the edge condition whenever
    no triangle is self-identified -- see find_self_identified_triangles),
    fully selecting one simple cycle C saturates every one of its nodes at
    degree 2, so no other triangle touching a node of C can also be
    selected. Then "all of C's variables are true" exactly captures "C is
    a whole connected component of the vertex link", and forbidding any
    two *node-disjoint* cycles from both being fully true forbids v from
    ever having 2+ separate circle-components in its link (a pinch
    point), with pairwise exclusion sufficing for 3+ as well.

    CORRECTNESS CAVEAT (confirmed by brute-force counterexample): this
    argument relies on every link-graph node having selected-degree at
    most 2, which can FAIL when a triangle is self-identified (two or
    more of its own vertices coincide -- see
    find_self_identified_triangles). Such a triangle can push a node's
    degree above 2, and a graph with degree->2 nodes can't be soundly
    analyzed via simple-cycle decomposition + pairwise exclusion: two
    "cycles" that look node-disjoint in isolation can be bridged into one
    genuinely-connected piece by a third triangle's arcs, which this
    pairwise reasoning has no way to detect. Concretely, on
    regina.Example4.ballBundle() (which has self-identified triangles),
    this function causes EAGER mode to find only 6 of the 12 genuinely
    valid connected surfaces (verified by brute force over all 2^12
    triangle subsets) -- it incorrectly excludes valid surfaces rather
    than accepting invalid ones. LAZY mode (find_vertex_violations) has
    no such issue, since it checks the actual connectivity of the
    concrete active subgraph directly rather than reasoning about
    precomputed cycles in isolation.

    Cost warning: this enumerates EVERY simple cycle in every vertex's
    link graph up front. That enumeration is worst-case super-exponential
    in the number of distinct edges touching a vertex (i.e. that vertex's
    degree in the 1-skeleton) -- a vertex whose link graph looks like a
    dense/complete graph on more than about a dozen nodes can make this
    hang or exhaust memory, *regardless* of how small the rest of the
    triangulation is. See build_solver(mode="lazy") for a version that
    entirely avoids this enumeration.
    """
    if link_graphs is None:
        link_graphs = build_vertex_link_graphs(tri)

    n_constraints = 0
    for _vidx, g in link_graphs.items():
        cycles = enumerate_link_cycles(g)
        for (nodes_a, vars_a), (nodes_b, vars_b) in combinations(cycles, 2):
            if nodes_a.isdisjoint(nodes_b):
                cycle_a_selected = z3.And([xs[i] for i in vars_a])
                cycle_b_selected = z3.And([xs[i] for i in vars_b])
                solver.add(z3.Not(z3.And(cycle_a_selected, cycle_b_selected)))
                n_constraints += 1
    return n_constraints


def find_vertex_violations(model, xs, link_graphs):
    """
    NOT SOUND IN GENERAL -- kept for reference and diagnostic use only.
    See is_valid_closed_surface for the robust replacement actually used
    by all_sat_enumerate.

    This checks, per ambient vertex, whether the selected link-graph
    edges form more than one connected component, on the theory that
    >1 component means a pinch point. That theory has TWO failure modes,
    both confirmed by brute-force counterexample on triangulations with
    self-identified triangles:

      * FALSE NEGATIVES (eager mode's failure): a single triangle with
        multiple corners at the same ambient vertex can make one
        genuinely-connected piece look, under a naive cycle
        decomposition, like two disjoint pieces that a third triangle
        then bridges together -- see build_vertex_constraints's
        docstring (regina.Example4.ballBundle() example).

      * FALSE POSITIVES (this function's own failure): conversely, two
        components that are genuinely disjoint in the ambient link graph
        do NOT always mean a pinch point -- they can be two entirely
        separate abstract vertices that simply happen to coincide at the
        same ambient vertex, which is perfectly valid (no different from
        two points of a manifold happening to land at the same
        coordinates under an embedding). Confirmed on
        regina.Example4.s2xs2(), solution (triangles 2, 3): this function
        reports a "violation" (two disjoint self-loop components at
        ambient vertex 1), yet Regina confirms the reconstructed surface
        is a perfectly valid sphere -- it has 3 abstract vertices, and
        the two flagged components are each their own separate, entirely
        valid, degree-1 abstract vertex. Treating this as a violation and
        cutting on it would permanently and incorrectly forbid triangles
        2 and 3 from ever being selected together.

    The fundamental problem: connectivity of the *ambient*-vertex-grouped
    ARC graph is not the same question as "how many *abstract* vertices
    does this correspond to, and does each have a valid circular link" --
    those coincide only when no triangle is self-identified. There is no
    way to fix this by adjusting the ambient-vertex grouping; the
    correct check has to operate on the actual abstract reconstruction,
    which is exactly what is_valid_closed_surface does by asking Regina
    directly.

    Returns
    -------
    list[tuple[int, list[frozenset[int]]]]
        One entry per *flagged* vertex (see caveats above):
        (vertex_index, components), each component a frozenset of
        triangle (Z3 variable) indices.
    """
    violations = []
    for vidx, g in link_graphs.items():
        active_adj = defaultdict(list)
        for u, v, _key, data in g.edges(keys=True, data=True):
            var = data["var"]
            if z3.is_true(model.eval(xs[var], model_completion=True)):
                active_adj[u].append((v, var))
                active_adj[v].append((u, var))
        if not active_adj:
            continue  # surface doesn't touch v at all -- nothing to check

        visited = set()
        components = []
        for start in active_adj:
            if start in visited:
                continue
            comp_vars = set()
            queue = deque([start])
            visited.add(start)
            while queue:
                cur = queue.popleft()
                for nbr, var in active_adj[cur]:
                    comp_vars.add(var)
                    if nbr not in visited:
                        visited.add(nbr)
                        queue.append(nbr)
            components.append(frozenset(comp_vars))

        if len(components) > 1:
            violations.append((vidx, components))

    return violations


def add_refinement_cuts(solver, xs, violations):
    """
    NOT SOUND IN GENERAL -- kept for reference only alongside
    find_vertex_violations. See that function's docstring: the
    "violations" it detects are not always genuine, so cutting on them
    can permanently and incorrectly exclude valid surfaces. Not used by
    all_sat_enumerate's default path; see is_valid_closed_surface.
    """
    n_cuts = 0
    for _vidx, components in violations:
        for comp_a, comp_b in combinations(components, 2):
            cond_a = z3.And([xs[i] for i in comp_a])
            cond_b = z3.And([xs[i] for i in comp_b])
            solver.add(z3.Not(z3.And(cond_a, cond_b)))
            n_cuts += 1
    return n_cuts


def is_valid_closed_surface(tri, triangle_indices):
    """
    THE ROBUST CHECK actually used by LAZY mode. Build the standalone
    abstract regina.Triangulation2 for `triangle_indices` (see
    build_surface_triangulation2, defined later in this file -- the
    forward reference is fine since Python resolves it at call time) and
    ask Regina directly whether every vertex link is valid.

    Why this is correct where the ambient-vertex-grouping approach isn't:
    Regina's isValid() operates on the actual abstract vertex classes
    formed by the real gluing structure (transitively closing over every
    edge identification), not on an ad-hoc grouping by which ambient
    vertex each corner happens to map to. Whether two corners end up as
    the same abstract vertex, or as two different abstract vertices that
    simply coincide at the same ambient point, is exactly what Regina's
    own combinatorial machinery is designed to get right -- so we lean on
    it rather than re-deriving an equivalent (and, as demonstrated twice
    above, easy-to-get-subtly-wrong) criterion ourselves.

    Returns
    -------
    bool
        True iff `triangle_indices` forms a valid closed 2-manifold
        (every edge shared by exactly 2 triangles -- already guaranteed
        by the Z3 edge constraints -- and every vertex link a circle).
    """
    t2 = build_surface_triangulation2(tri, triangle_indices)
    return t2.isValid() and t2.isClosed()


# ---------------------------------------------------------------------------
# 4. All-SAT enumeration loop
# ---------------------------------------------------------------------------

def build_solver(tri, mode="lazy"):
    """
    Assemble the Z3 solver for the given 4-manifold triangulation.

    Parameters
    ----------
    mode : {"lazy", "eager"}
        "lazy" (default, recommended): only the (cheap) edge constraints
            are added up front. The vertex condition is instead enforced
            during the All-SAT loop itself via a CEGAR-style loop: each
            candidate model is checked with is_valid_closed_surface,
            which builds the actual abstract Triangulation2 and asks
            Regina directly whether it's a valid closed surface; a
            candidate that fails is rejected (plain blocking, not
            counted as a solution) and the solver is re-queried. This
            never enumerates a vertex link's full cycle set, so it
            scales to triangulations where that enumeration would be
            intractable, AND it has no known soundness caveats: unlike
            the ambient-vertex-grouping heuristic this project started
            with (see find_vertex_violations's docstring for two
            confirmed counterexamples), asking Regina directly about the
            actual abstract reconstruction is correct regardless of
            self-identified triangles.
        "eager": the original approach -- every simple cycle in every
            vertex's link graph is enumerated up front and all pairwise
            disjoint-cycle exclusions are added before search begins.
            Kept for comparison/debugging on triangulations where it is
            known to be sound. Two caveats:
              (1) Cost: can hang or exhaust memory on triangulations with
                  a vertex whose link graph is large and dense,
                  regardless of the triangulation's overall size.
              (2) Correctness: provably sound only when NO triangle is
                  self-identified (find_self_identified_triangles(tri) is
                  empty). Otherwise it can UNDER-count -- incorrectly
                  reject valid closed surfaces -- as confirmed by a
                  brute-force counterexample on
                  regina.Example4.ballBundle() (see
                  build_vertex_constraints's docstring for the full
                  explanation). A warning is printed if this mode is
                  requested on a triangulation with self-identified
                  triangles.

    Returns
    -------
    tuple[z3.Solver, list[z3.BoolRef], dict[int, networkx.MultiGraph]]
        The solver, the triangle indicator variables, and the per-vertex
        link graphs (kept in the return signature for backward
        compatibility and for anyone wanting to inspect link structure
        directly; not used by all_sat_enumerate in lazy mode, since that
        now checks each candidate via Regina rather than these graphs).
    """
    if mode not in ("lazy", "eager"):
        raise ValueError(f"mode must be 'lazy' or 'eager', got {mode!r}")

    xs = build_variables(tri)
    solver = z3.Solver()
    forbid_empty_surface(xs, solver)
    build_edge_constraints(tri, xs, solver)
    link_graphs = build_vertex_link_graphs(tri)

    if mode == "eager":
        self_identified = find_self_identified_triangles(tri)
        if self_identified:
            print(
                f"WARNING: mode='eager' requested, but {len(self_identified)} "
                f"triangle(s) are self-identified (e.g. triangle "
                f"{self_identified[0]} has a repeated vertex label among its "
                "own 3 corners). Eager mode's pairwise disjoint-cycle "
                "exclusion is not sound in this case and may incorrectly "
                "reject valid closed surfaces (under-count). Consider "
                "mode='lazy' (the default) instead -- see "
                "build_vertex_constraints's docstring for details."
            )
        build_vertex_constraints(tri, xs, solver, link_graphs)

    return solver, xs, link_graphs


def all_sat_enumerate(tri, xs, solver, use_lazy_check=False, max_solutions=None):
    """
    All-SAT loop: repeatedly ask Z3 for a satisfying model, record which
    triangles are selected, then add a "blocking clause" (the negation of
    that exact truth assignment) so the next solver.check() is forced to
    find a genuinely different configuration.

    If `use_lazy_check` is True, the vertex condition is checked and
    enforced LAZILY: every candidate model is passed through
    is_valid_closed_surface (which builds the actual abstract
    Triangulation2 and asks Regina directly). A model that fails this
    check is rejected -- blocked with a plain blocking clause and NOT
    counted as a solution -- and the solver is re-queried. Plain blocking
    (rather than a "smarter" generalizing cut) is used deliberately: see
    find_vertex_violations's docstring for two confirmed counterexamples
    showing that ambient-vertex-based generalizing cuts can be WRONG
    (permanently excluding valid surfaces), so only the provably-safe
    "reject exactly this assignment" cut is used here.

    If `use_lazy_check` is False, the vertex condition is assumed to
    already be fully encoded in `solver` (i.e. eager mode via
    build_vertex_constraints), and every satisfying model is accepted
    directly.

    Returns
    -------
    list[list[int]]
        Each entry is the sorted list of selected triangle indices for one
        satisfying assignment that is a genuine closed 2-manifold.
    """
    solutions = []
    while solver.check() == z3.sat:
        model = solver.model()
        selected = [
            i for i, v in enumerate(xs)
            if z3.is_true(model.eval(v, model_completion=True))
        ]
        is_genuine = (not use_lazy_check) or is_valid_closed_surface(tri, selected)

        if is_genuine:
            solutions.append(selected)

        # Either way, this exact assignment must not reappear: if genuine,
        # we've recorded it and want a different one next; if not, it was
        # never a valid surface and must not be re-proposed. Blocking the
        # literal assignment (rather than a broader cut) is always safe.
        selected_set = set(selected)
        blocking_clause = z3.Or([
            z3.Not(v) if i in selected_set else v
            for i, v in enumerate(xs)
        ])
        solver.add(blocking_clause)

        if is_genuine and max_solutions is not None and len(solutions) >= max_solutions:
            break

    return solutions


# ---------------------------------------------------------------------------
# 5. Post-processing: connectivity filter
# ---------------------------------------------------------------------------

def build_triangle_adjacency(tri, edge_to_tris=None):
    """
    Two triangles are adjacent iff they share a common edge of the
    4-manifold. Returns dict: triangle_index -> set of adjacent triangle
    indices (over the *whole* triangulation; the BFS below only follows
    edges between triangles that are both selected in a given solution).
    """
    if edge_to_tris is None:
        edge_to_tris = compute_edge_to_triangles(tri)

    adjacency = defaultdict(set)
    for _eidx, tri_list in edge_to_tris.items():
        for a, b in combinations(tri_list, 2):
            adjacency[a].add(b)
            adjacency[b].add(a)
    return adjacency


def connected_components(selected, adjacency):
    """BFS over the induced subgraph on `selected` triangles."""
    selected_set = set(selected)
    visited = set()
    components = []
    for start in selected:
        if start in visited:
            continue
        comp = set()
        queue = deque([start])
        visited.add(start)
        while queue:
            cur = queue.popleft()
            comp.add(cur)
            for nb in adjacency.get(cur, ()):
                if nb in selected_set and nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        components.append(comp)
    return components


def filter_connected_surfaces(solutions, adjacency):
    """Keep only those solutions whose selected triangles form a single
    connected component (i.e. discard disjoint unions of surfaces such
    as "two separate spheres")."""
    connected = []
    for sol in solutions:
        comps = connected_components(sol, adjacency)
        if len(comps) == 1:
            connected.append(sol)
    return connected


# ---------------------------------------------------------------------------
# 6. Invariants: Euler characteristic, orientability, homeomorphism type
# ---------------------------------------------------------------------------
#
# Every closed connected surface is homeomorphic to exactly one of:
#   - the sphere S^2                                 (chi = 2, orientable)
#   - the orientable genus-g surface Sigma_g          (chi = 2 - 2g, orientable)
#   - the non-orientable surface with k crosscaps N_k (chi = 2 - k, non-or.)
# (the classification theorem of compact surfaces). So computing just two
# invariants -- chi and orientability -- is enough to name the exact
# homeomorphism type of every surface our SAT search finds.
#
# Rather than re-deriving Euler characteristic (V - E + F) and
# orientability (sign-propagation across shared edges) by hand, we lean on
# Regina itself: build a genuine, independent regina.Triangulation2 that
# mirrors the selected sub-complex's own gluings, and let Regina's own
# 2-manifold engine compute eulerChar()/isOrientable() for us. This is
# both less code and a stronger guarantee of correctness (it's the same
# code Regina uses everywhere else), and it comes with a free bonus: we
# also get isValid()/isClosed()/isConnected() as independent sanity checks
# on the whole SAT encoding above.

def build_surface_triangulation2(tri, triangle_indices):
    """
    Construct a standalone regina.Triangulation2 combinatorially identical
    to the closed surface spanned by `triangle_indices` inside the ambient
    4-manifold triangulation `tri`.

    One Triangle2 is created per selected triangle. Two of these new
    triangles are glued along a facet exactly where the ambient
    4-manifold already glues the corresponding shared edge.

    Correctness note: the gluing permutation is derived from each ambient
    triangle's own `edgeMapping()` (Regina's intrinsic embedding of the
    shared edge's two vertices into that triangle's local vertex slots),
    NOT by naively matching global vertex labels. Naive label-matching
    silently breaks in highly degenerate/minimal triangulations, where a
    triangle's own vertices can be partly or fully self-identified (e.g.
    a triangulation with very few global vertices can have a triangle all
    three of whose corners map to the *same* global vertex) -- in that
    situation "the local vertex whose global label matches X" is
    ambiguous, but edgeMapping() still gives the unique correct answer
    because it comes from Regina's own combinatorial data, not from
    re-deriving it via labels that may coincide.

    This also correctly handles *self-gluing*: a single triangle whose
    own two edges get identified to each other by the ambient
    triangulation (so the same triangle index appears twice in one
    edge's incidence list, once per local facet). We track incidences at
    the (triangle, local-facet) level rather than by triangle index alone
    so this case is represented faithfully instead of being collapsed.
    """
    local_set = set(triangle_indices)
    t2 = regina.Triangulation2()
    local = {ti: t2.newTriangle() for ti in triangle_indices}

    # Map each global edge index to the (triangle_index, local_facet)
    # pairs, among the selected triangles, that touch it. Built directly
    # from (triangle, local facet) -- NOT from a pre-aggregated
    # edge_to_tris -- so that a triangle touching the same global edge
    # twice (self-gluing) yields two separate entries rather than one.
    edge_incidences = defaultdict(list)
    for ti in triangle_indices:
        amb = tri.triangle(ti)
        for local_e in range(3):
            eidx = amb.edge(local_e).index()
            edge_incidences[eidx].append((ti, local_e))

    for eidx, incidences in edge_incidences.items():
        # The edge condition guarantees this is always exactly 2 for a
        # valid closed-surface solution (whether that's two distinct
        # triangles, or the same triangle touching this edge via two of
        # its own facets).
        assert len(incidences) == 2, (
            f"edge {eidx} touched by {len(incidences)} (triangle, facet) "
            f"incidences among selected triangles -- not a closed surface"
        )
        (t1_idx, facet1), (t2_idx, facet2) = incidences
        amb1, amb2 = tri.triangle(t1_idx), tri.triangle(t2_idx)

        # edgeMapping(local_facet) is a Perm5 (Regina uses Perm<dim+1>
        # uniformly); its first two images tell us which local
        # triangle-vertex corresponds to the shared edge's own intrinsic
        # vertex 0 and vertex 1 respectively. Composing the two mappings
        # gives the correct on-edge vertex correspondence regardless of
        # any vertex/edge self-identification.
        map1 = amb1.edgeMapping(facet1)
        map2 = amb2.edgeMapping(facet2)

        images = [None, None, None]
        images[map1[0]] = map2[0]
        images[map1[1]] = map2[1]
        # The one local vertex not on the shared edge (facet1 itself)
        # maps to whichever local index of t2 is left over -- this choice
        # is free (it doesn't affect the resulting surface's combinatorics)
        # since it's off the glued facet.
        used = {images[map1[0]], images[map1[1]]}
        images[facet1] = ({0, 1, 2} - used).pop()

        local[t1_idx].join(facet1, local[t2_idx], regina.Perm3(*images))

    return t2


def classify_surface(euler_characteristic, orientable):
    """
    Apply the classification theorem of compact surfaces to name the
    homeomorphism type corresponding to a given (chi, orientable) pair.
    (Regina doesn't expose a single "name this surface" call beyond
    isSphere(), so this small lookup is the one place we still apply the
    classification theorem ourselves.)
    """
    chi = euler_characteristic
    if orientable:
        if chi == 2:
            return "Sphere (S^2)"
        if chi == 0:
            return "Torus (T^2)"
        genus = (2 - chi) // 2
        return f"Orientable genus-{genus} surface (Sigma_{genus})"
    else:
        if chi == 1:
            return "Projective plane (RP^2)"
        if chi == 0:
            return "Klein bottle (K)"
        crosscaps = 2 - chi
        return f"Non-orientable surface with {crosscaps} crosscaps (N_{crosscaps})"


def compute_surface_invariants(tri, triangle_indices):
    """
    Compute the full invariant summary for one connected closed surface by
    building an independent regina.Triangulation2 (see
    build_surface_triangulation2) and reading its invariants straight from
    Regina's own engine.
    """
    t2 = build_surface_triangulation2(tri, triangle_indices)

    # The SAT constraints should already guarantee all three of these, but
    # asking Regina to independently confirm validity/closedness/
    # connectivity is a free, strong sanity check on the whole pipeline.
    assert t2.isValid(), "Regina flagged this surface as an invalid triangulation"
    assert t2.isClosed(), "Regina flagged this surface as having boundary"
    assert t2.isConnected(), "Regina flagged this surface as disconnected"

    chi = t2.eulerChar()
    orientable = t2.isOrientable()
    surface_type = classify_surface(chi, orientable)

    return {
        "triangles": sorted(triangle_indices),
        "num_vertices": t2.countVertices(),
        "num_edges": t2.countEdges(),
        "num_faces": t2.countTriangles(),
        "euler_characteristic": chi,
        "orientable": orientable,
        "type": surface_type,
    }


def compute_all_invariants(tri, connected_surfaces):
    """Run compute_surface_invariants over a whole list of connected
    surfaces, returning one invariant-summary dict per surface (same
    order as the input list)."""
    return [
        compute_surface_invariants(tri, sol)
        for sol in connected_surfaces
    ]


# ---------------------------------------------------------------------------
# Top-level convenience wrapper
# ---------------------------------------------------------------------------

def enumerate_closed_surfaces(tri, max_solutions=None, mode="lazy"):
    """
    Full pipeline: build the SAT model for `tri`, run All-SAT, filter down
    to connected closed surfaces, and compute each one's invariants and
    homeomorphism type.

    Parameters
    ----------
    mode : {"lazy", "eager"}
        See build_solver. "lazy" (default) enforces the vertex condition
        via on-demand CEGAR refinement and never enumerates a vertex
        link's full cycle set; "eager" enumerates all cycles up front
        (kept for comparison/debugging, can be intractable on dense
        vertex links).

    Returns
    -------
    dict
        {
          "all_solutions": list of all valid (possibly disconnected) closed
                            2-manifold triangle-subsets found,
          "connected_surfaces": the subset of those that are connected,
          "invariants": list of invariant-summary dicts, one per entry of
                        "connected_surfaces" (same order), each with keys
                        triangles, num_vertices, num_edges, num_faces,
                        euler_characteristic, orientable, type,
          "mode": the mode string that was used,
        }
    """
    solver, xs, link_graphs = build_solver(tri, mode=mode)
    all_solutions = all_sat_enumerate(
        tri, xs, solver,
        use_lazy_check=(mode == "lazy"),
        max_solutions=max_solutions,
    )

    edge_to_tris = compute_edge_to_triangles(tri)
    adjacency = build_triangle_adjacency(tri, edge_to_tris)
    connected_surfaces = filter_connected_surfaces(all_solutions, adjacency)
    invariants = compute_all_invariants(tri, connected_surfaces)

    return {
        "all_solutions": all_solutions,
        "connected_surfaces": connected_surfaces,
        "invariants": invariants,
        "mode": mode,
    }


# ---------------------------------------------------------------------------
# Demo / smoke test
# ---------------------------------------------------------------------------

def _load_triangulation_from_args(args):
    """
    Resolve the --isosig / --example CLI arguments (see main()) into a
    regina.Triangulation4, or fall back to the default demo triangulation
    if neither was given. Exits with a clear error message on bad input
    rather than letting a raw Regina exception/traceback surface.
    """
    if args.isosig is not None:
        try:
            tri = regina.Triangulation4.fromIsoSig(args.isosig)
        except regina.InvalidArgument as e:
            sys.exit(
                f"error: '{args.isosig}' is not a valid isomorphism "
                f"signature for a 4-manifold triangulation ({e}). Note "
                "isomorphism signatures are dimension-specific -- a "
                "signature produced for a 3-manifold or 5-manifold "
                "triangulation will not decode here."
            )
        if tri is None:
            sys.exit(
                f"error: '{args.isosig}' could not be decoded into a "
                "triangulation"
            )
        print(f"Loaded triangulation from isosig: {args.isosig}")
        return tri

    if args.example is not None:
        try:
            constructor = getattr(regina.Example4, args.example)
        except AttributeError:
            sys.exit(
                f"error: '{args.example}' is not a known regina.Example4 "
                "constructor. Run with --list-examples to see options."
            )
        try:
            tri = constructor()
        except TypeError as e:
            sys.exit(
                f"error: regina.Example4.{args.example}() needs "
                f"arguments this script doesn't provide ({e}); try a "
                "different --example, or use --isosig instead"
            )
        print(f"Loaded triangulation from regina.Example4.{args.example}()")
        return tri

    # Default: the standard simplicial 4-sphere used throughout this
    # file's docstrings and the correctness test suite -- boundary of
    # the 5-simplex, 6 vertices / 15 edges / 20 triangles / 15
    # tetrahedra / 6 pentachora. Its 2-skeleton is the complete graph's
    # triangle complex on 6 vertices, rich enough to contain several
    # genuinely different closed surfaces.
    print(
        "No --isosig or --example given; using the default demo "
        "triangulation (regina.Example4.simplicialFourSphere())."
    )
    return regina.Example4.simplicialFourSphere()


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=(
            "Enumerate all connected, closed surfaces (2-manifolds) "
            "embedded in the 2-skeleton of a triangulated compact "
            "4-manifold."
        )
    )
    source = parser.add_mutually_exclusive_group()
    source.add_argument(
        "--isosig", metavar="SIG",
        help=(
            "Isomorphism signature of a 4-manifold triangulation to "
            "load, as produced by Regina's Triangulation4.isoSig() "
            "(e.g. copied from Regina's census data or another Regina "
            "session). Takes priority over --example."
        ),
    )
    source.add_argument(
        "--example", metavar="NAME",
        help=(
            "Name of one of Regina's built-in no-argument "
            "regina.Example4 constructors to use instead, e.g. "
            "'k3', 'rp4', 's2xs2'. See --list-examples for the full list."
        ),
    )
    parser.add_argument(
        "--list-examples", action="store_true",
        help="List regina.Example4 constructor names and exit "
             "(without running anything). Not all of them take zero "
             "arguments; this script will report a clear error if the "
             "one you pick needs parameters.",
    )
    parser.add_argument(
        "--mode", choices=["lazy", "eager"], default="lazy",
        help="Vertex-condition enforcement strategy (default: lazy, "
             "recommended). See build_solver's docstring for the "
             "correctness/performance tradeoffs of 'eager'.",
    )
    parser.add_argument(
        "--max-solutions", type=int, default=500, metavar="N",
        help="Stop after finding this many valid closed-surface "
             "triangle-subsets (default: 500). Use a value <= 0 for no "
             "limit -- may not terminate on large or combinatorially "
             "rich triangulations.",
    )
    parser.add_argument(
        "--top", type=int, default=10, metavar="N",
        help="Number of connected surfaces to print full invariants for "
             "(default: 10).",
    )
    args = parser.parse_args(argv)

    if args.list_examples:
        names = sorted(
            name for name in dir(regina.Example4)
            if not name.startswith("_") and name != "equalityType"
        )
        print("regina.Example4 constructors (pass the name after "
              "--example):")
        for name in names:
            print(f"  {name}")
        return

    tri = _load_triangulation_from_args(args)

    if not tri.isValid():
        print(
            "WARNING: this triangulation is not itself Regina-valid "
            "(tri.isValid() is False). The search below still runs over "
            "its 2-skeleton, and every reported surface is independently "
            "Regina-validated on its own terms, but it's worth knowing "
            "the ambient space itself is unusual."
        )

    print(f"  vertices={tri.countVertices()} edges={tri.countEdges()} "
          f"triangles={tri.countTriangles()} tetrahedra={tri.countTetrahedra()} "
          f"pentachora={tri.countPentachora()}")
    print(f"  mode={args.mode}")

    max_solutions = args.max_solutions if args.max_solutions > 0 else None
    result = enumerate_closed_surfaces(tri, max_solutions=max_solutions, mode=args.mode)

    hit_cap = (
        max_solutions is not None
        and len(result["all_solutions"]) >= max_solutions
    )
    print(
        f"\nTotal valid closed-surface triangle-subsets found: "
        f"{len(result['all_solutions'])}"
        f"{'  (hit --max-solutions cap; more may exist)' if hit_cap else ''}"
    )
    print(f"Of which connected (genuine single surfaces): "
          f"{len(result['connected_surfaces'])}")

    invariants = result["invariants"]

    print("\nSample connected surfaces with computed invariants "
          "(V, E, F, chi=Euler characteristic, orientable?, homeomorphism type):")
    for inv in invariants[:args.top]:
        orient_str = "orientable" if inv["orientable"] else "non-orientable"
        print(
            f"  triangles={inv['triangles']}\n"
            f"    V={inv['num_vertices']} E={inv['num_edges']} "
            f"F={inv['num_faces']}  chi={inv['euler_characteristic']:+d}  "
            f"{orient_str}  ->  {inv['type']}"
        )

    print("\nSummary: homeomorphism type distribution across all "
          f"{len(invariants)} connected surfaces found:")
    type_counts = Counter(inv["type"] for inv in invariants)
    for surface_type, count in sorted(type_counts.items(), key=lambda kv: -kv[1]):
        print(f"  {count:4d}  x  {surface_type}")


if __name__ == "__main__":
    if regina is None:
        raise SystemExit(
            "This demo requires the `regina` package (pip install regina)."
        )
    main()
