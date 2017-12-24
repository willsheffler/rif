import rif
from rif.io import pdbio


def test_pdbio(pdbfname):
    df = pdbio.readpdb(pdbfname)
    assert df.shape == (1698, 12)
    assert list(df) == 'het ai an rn ch ri x y z occ bfac elem'.split()
