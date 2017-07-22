from _rif._test import _pandas_interop as pi

import numpy as np
import pandas as pd


def test_can_treat_pd_as_dict():
    s = pd.Series([1, 3, 5, np.nan, 6, 8])
    dates = pd.date_range('20130101', periods=6)
    df = pd.DataFrame(np.random.randn(6, 4), index=dates, columns=list('ABCD'))
    sp = pi.print_dict(df)
    assert 28 == len(sp.splitlines())
    kn, totlen = pi.key_names_and_len(df)
    assert totlen == 24
