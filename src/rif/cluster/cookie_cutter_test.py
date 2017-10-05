from rif.cluster.cookie_cutter import *
import pandas as pd

# store = pd.HDFStore(
# '/home/sheffler/Dropbox/project/hbsat/test_data/scott_427_0001_sol.h5')
# x = store['s'].copy()
# x = x.drop(x.columns[np.min(x, axis=0) == np.max(x, axis=0)], axis=1)
# x = x.head(10000)


def test_cookie_cutter():
    x = np.arange(100).reshape((10, 10))
    x = np.concatenate((x, x), axis=0)
    clust = cookie_cutter(x, 3, 'ndiff')
    clust = x[clust, :]
    assert clust.shape == (10, 10)


def test_cookie_cutter_cpp():
    x = np.arange(10000).reshape((100, 100))
    x = np.concatenate((x, x, x, x, x), axis=0)
    clust = cookie_cutter_i2(x, 3, 'ndiff')
    clust = x[clust, :]
    assert clust.shape == (100, 100)
