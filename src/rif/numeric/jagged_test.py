import numpy as np
import pandas as pd
import xarray as xr

from rif.dtypes import RifOperatorsDisabled


def test_indexing():
    with RifOperatorsDisabled():
        midx = pd.MultiIndex.from_product([list('abc'), [0, 1]],
                                          names=('one', 'two'))
        mda = xr.DataArray(np.random.rand(6, 3),
                           [('x', midx), ('y', range(3))])

        assert mda.shape == (6, 3)
        assert mda.sel(x=(list('ab'), [0])).shape == (2, 3)
        assert mda.sel(x=[('a', 0), ('b', 1)]).shape == (2, 3)
        assert mda.sel(x={'one': 'a', 'two': 0}).shape == (3,)
        assert mda.sel(one='a', two=0).shape == (3,)
        assert mda.loc[{'one': 'a'}, ...].shape == (2, 3)
