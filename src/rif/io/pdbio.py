import pandas as pd
import numpy as np


def readpdb(fname):
    n = 'het ai an rn ch ri x y z occ bfac elem'.split()
    w = (6, 5, 5, 5, 2, 6, 9, 9, 9, 5, 5, 99)
    assert len(n) is len(w)
    df = pd.read_fwf(fname, widths=w, names=n)
    df = df[np.logical_or(df.het == 'ATOM', df.het == 'HETATM')]
    df.het = df.het == 'HETATM'
    df.ai = df.ai.astype('i4')
    # df.an = df.an.astype('S4')  # f*ck you, pandas!
    # df.rn = df.rn.astype('S3')  # f*ck you, pandas!
    # df.ch = df.ch.astype('S1')  # f*ck you, pandas!
    df.ri = df.ri.astype('i4')
    df.x = df.x.astype('f4')
    df.y = df.y.astype('f4')
    df.z = df.z.astype('f4')
    df.occ = df.occ.astype('f4')
    df.bfac = df.bfac.astype('f4')
    # df.elem = df.elem.astype('S4')  # f*ck you, pandas!
    return df
