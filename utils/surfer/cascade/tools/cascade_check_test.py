"""Tests of cascade_check.py's own diagram code (independent of the C++
tests, as the checker is of the C++).

usage: ~/.venvs/atlas/bin/python cascade_check_test.py"""
import os, sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import cascade_check as cc  # noqa: E402

G = cc.Gauss
failures = 0


def check(ok, what):
    global failures
    if not ok:
        failures += 1
        print('FAIL:', what)


def same_link(a, b):
    return a.link().jones() == b.link().jones() and a.lk() == b.lk()


def relabel(g, perm):
    """g with crossing k renamed perm[k] (a diagram isomorphism)."""
    signs = [0] * len(g.signs)
    for k, s in enumerate(g.signs):
        signs[perm[k]] = s
    return G(signs, [[(1 if x > 0 else -1) * (perm[abs(x) - 1] + 1) for x in w] for w in g.comps])


def test_kinks():
    for s in (1, -1):
        for word in ([1, -1], [-1, 1]):
            r = cc.remove_nugatory(G([s], [word]))
            check(r.signs == [] and r.comps == [[]], f'kink {s} {word} -> {r.comps}')


def test_sums_through_a_twist():
    tA, tB = [1, -2, 3, -1, 2, -3], [4, -5, 6, -4, 5, -6]
    for sc in (1, -1):
        for ha, hb in ((1, 1), (1, -1), (-1, -1)):
            for word in ([7] + tA + [-7] + tB, [-7] + tA + [7] + tB):
                g = G([ha] * 3 + [hb] * 3 + [sc], [word])
                r = cc.remove_nugatory(g)
                check(len(r.signs) == 6, f'sum {sc} {ha} {hb}: {len(r.signs)} crossings')
                check(same_link(g, r), f'sum {sc} {ha} {hb}: link changed')
                check(cc.remove_nugatory(r).comps == r.comps, 'a reduced sum changed again')
    # Deleting the crossing WITHOUT turning a side over also gives a diagram
    # of the same link (the sum of the two sides along the same components),
    # so the invariants above cannot see the turn. The turn fixes WHICH
    # diagram, the one the C++ writes for a node: pin it. The side between
    # c's two visits swaps over and under, and keeps its signs.
    r = cc.remove_nugatory(G([1] * 6 + [-1], [[7] + tA + [-7] + tB]))
    check(r.comps == [[-x for x in tA] + tB] and r.signs == [1] * 6,
          f'the side between the visits was not turned over: {r.comps}')


def test_side_holding_a_component():
    # a Hopf clasp on one loop of a nugatory crossing: the whole second
    # component is on the side turned over
    g = G([1, 1, 1], [[3, 1, -2, -3], [-1, 2]])
    r = cc.remove_nugatory(g)
    check(len(r.signs) == 2, f'clasp on a loop: {len(r.signs)} crossings')
    check(same_link(g, r), 'clasp on a loop: link changed')
    check(r.lk() == [[0, 1], [1, 0]], f'clasp on a loop: lk {r.lk()}')
    check(r.comps == [[-1, 2], [1, -2]], f'clasp on a loop: the component was not turned over: {r.comps}')


def test_crossing_joined_through_another_component():
    # component 1 crosses both loops of crossing 5: not nugatory
    g = G([1] * 5, [[5, 1, -2, -5, 3, -4], [-1, 2, 4, -3]])
    check(cc.nugatory_side(g, 4) is None, 'crossing 5 called nugatory')
    check(len(cc.remove_nugatory(g).signs) == 5, 'a crossing was removed')


def test_reduced_untouched():
    for g in (G([1, 1, -1, -1], [[1, -2, 3, -4, 2, -1, 4, -3]]),   # 4_1
              G([1, 1], [[1, -2], [-1, 2]]),                        # Hopf
              G([1, 1, 1], [[1, -2, 3, -1, 2, -3]])):               # 3_1
        check(cc.remove_nugatory(g).comps == g.comps, f'reduced diagram changed: {g.comps}')


def test_same_as_node_through_a_twist():
    tA, tB = [1, -2, 3, -1, 2, -3], [4, -5, 6, -4, 5, -6]
    checker = cc.Checker({'records': [], 'nodes': []}, None, {})
    for sc in (1, -1):
        g = G([1] * 6 + [sc], [[7] + tA + [-7] + tB])
        node = relabel(cc.remove_nugatory(g), [3, 0, 5, 1, 4, 2])
        ok, why = checker.same_as_node(g, node, [0], False, False)
        check(ok and why == 'diagram', f'sum through a {sc} twist not matched: {why}')
        # the node the other side's choice gives: the reduced diagram turned over
        ok, why = checker.same_as_node(g, node.turned_over(), [0], False, False)
        check(ok and why == 'diagram', f'turned-over node not matched: {why}')


for t in (test_kinks, test_sums_through_a_twist, test_side_holding_a_component,
          test_crossing_joined_through_another_component, test_reduced_untouched,
          test_same_as_node_through_a_twist):
    t()
print(f'{"FAILED" if failures else "passed"}: {failures} failures')
sys.exit(1 if failures else 0)
