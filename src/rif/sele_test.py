from rif.sele import parse_atom_names, parse_ray_atom_names


def test_atom_sele():
    anames = parse_atom_names('foo bar or foo')
    assert anames == [('FOO',), ('BAR', 'FOO')]
    anames = parse_atom_names('foo (bar or foo)')
    assert anames == [('FOO',), ('BAR', 'FOO')]
    anames = parse_atom_names('(foo) bar or foo')
    assert anames == [('FOO',), ('BAR', 'FOO')]


def test_ray_sele():
    rays = parse_ray_atom_names('a->b')
    assert rays == [(('A',), ('B',))]
    rays = parse_ray_atom_names('a->b c->a')
    assert rays == [(('A',), ('B',)), (('C',), ('A',))]
    # rays = parse_ray_atom_names('foo->(bar or baz) (a or C)->b')
    # assert rays == [(('FOO',), ('BAR', 'BAZ')), (('A', 'C'), ('B',))]
