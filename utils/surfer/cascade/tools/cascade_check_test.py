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


def test_lift_split():
    # a far side a hop refused (11a_239's run): component 1 passes over both
    # its crossings, so it is a split unknot lying above the rest
    g = G([1, 1, 1, -1, -1, -1], [[5, -6], [4, 2], [1, -4, -5, 6, -2, -3], [-1, 3]])
    r = cc.lift_split(g)
    check(r.comps[1] == [] and len(r.signs) == 4, f'lift: {r.comps}')
    check(same_link(g, r), 'lift: link changed')
    check(len(cc.split_pieces(r)) >= 2, 'lift: the unknot is not split off')
    check(cc.lift_split(r).comps == r.comps, 'lift: lifted twice')
    h = G([1, 1], [[1, -2], [-1, 2]])
    check(cc.lift_split(h).comps == h.comps, 'lift: the Hopf link changed')


def test_literature_leaf_of_the_target_under_another_name():
    # L11a397{1;0} and L11a397{0;0} are one link (a meridian-carrying
    # isometry; data/table_link_classes.csv), and {1;1} is another. A
    # literature leaf on the class-mate would prove the target by its own
    # literature value: refused with the class table, and without it by the
    # isometry. The other variant is not the target, so is not refused.
    atlas = cc.ATLAS + '/data/'
    knots = atlas + '4d_smooth_slice_genus_13_crossings_pd_codes.csv'
    links = atlas + 'links_4d_smooth_slice_genus_11_crossings_pd_codes.csv'
    target = 'L11a397{1;0}'
    for label, classes in (('class table', None), ('isometry', '/nonexistent')):
        tables = cc.load_tables(knots, links, classes)
        cert = {'target': target, 'records': [],
                'nodes': [{'id': 0, 'label': 'target ' + target, 'components': 3,
                           'pd': tables['pd'][target]}]}
        ch = cc.Checker(cert, None, tables)
        own, why = ch.is_target_link('L11a397{0;0}')
        check(own, f'{label}: the class-mate was not found to be the target ({why})')
        other, why = ch.is_target_link('L11a397{1;1}')
        check(not other, f'{label}: another variant was taken for the target ({why})')
        leaf = {'id': 1, 'node': 0, 'kind': 'leaf', 'genus': 1, 'partition': '{0,1,2}',
                'source': 'literature L11a397{0;0} ' + tables['g4']['L11a397{0;0}']}
        check(not ch.check_leaf(leaf), f'{label}: a class-mate literature leaf was accepted')
        check(any("the target's own link" in p for p in ch.problems),
              f'{label}: refused for another reason: {ch.problems}')


def test_connected_summands():
    # The square knot as the cascade drew it (10_99's node 68, 2026-09-28):
    # a right trefoil and a left one joined at a visible sphere. Cut into two
    # trefoils of opposite sign; each is a trefoil (Jones); a prime diagram
    # is left whole; a granny (same signs) cuts the same way.
    square = G([1, 1, 1, -1, -1, -1], [[6, -5, 4, -6, -3, 2, -1, 3, -2, 1, 5, -4]])
    parts = cc.connected_summands(square)
    check(len(parts) == 2, f'square knot: {len(parts)} summands')
    trefoil = G([1, 1, 1], [[1, -2, 3, -1, 2, -3]])
    for p in parts:
        check(len(p.signs) == 3 and abs(sum(p.signs)) == 3, f'summand signs {p.signs}')
        check(p.link().jones() in (trefoil.link().jones(), trefoil.mirror().link().jones()),
              'a summand is not a trefoil')
    check(sorted(sum(p.signs) for p in parts) == [-3, 3], 'the two trefoils are mirrors')
    check(len(cc.connected_summands(trefoil)) == 1, 'a prime diagram is left whole')
    granny = G([1, 1, 1, 1, 1, 1], [[6, -5, 4, -6, -3, 2, -1, 3, -2, 1, 5, -4]])
    check(sorted(sum(p.signs) for p in cc.connected_summands(granny)) == [3, 3], 'granny')


def test_elementary_slice():
    sym = {'3_1': 'rev', '4_1': 'full', '8_17': 'chiral', '9_32': 'chiral', '5_1': 'rev'}
    check(cc.composite_summands('3_1#m3_1') == [('3_1', False, False), ('3_1', True, False)],
          'summand parsing')
    check(cc.composite_summands('mr8_17#8_17')[0] == ('8_17', True, True), 'mr parsing')
    check(cc.elementary_slice('3_1#m3_1', sym)[0], 'square knot is slice')
    check(not cc.elementary_slice('3_1#3_1', sym)[0], 'granny is not')
    check(cc.elementary_slice('4_1#4_1', sym)[0], '4_1#4_1 is slice')
    check(cc.elementary_slice('m3_1#3_1#5_1#m5_1', sym)[0], 'two cancelling pairs')
    check(not cc.elementary_slice('m3_1#3_1#5_1', sym)[0], 'an unpaired summand')
    check(not cc.elementary_slice('8_17#m8_17', sym)[0],
          'a non-invertible knot with its plain mirror is NOT the pattern')
    check(cc.elementary_slice('8_17#mr8_17', sym)[0], 'with its reversed mirror it is')
    check(not cc.elementary_slice('3_1#m3_1#9_32#m9_32', {'3_1': 'rev'})[0],
          'an unknown symmetry type refuses')
    check(not cc.elementary_slice('3_1', sym)[0], 'a prime knot is no composite')


for t in (test_lift_split, test_kinks, test_sums_through_a_twist, test_side_holding_a_component,
          test_crossing_joined_through_another_component, test_reduced_untouched,
          test_same_as_node_through_a_twist,
          test_literature_leaf_of_the_target_under_another_name,
          test_connected_summands, test_elementary_slice):
    t()
print(f'{"FAILED" if failures else "passed"}: {failures} failures')
sys.exit(1 if failures else 0)
