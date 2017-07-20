from rif.sele import parse_atom_names


def test_AtomSele():
    anames = parse_atom_names('foo (bar or foo)')
    assert anames == [('FOO',), ('BAR', 'FOO')]
