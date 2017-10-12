from rif.cluster.cookie_cutter import *
import pandas as pd


def test_cookie_cutter():
    x = np.arange(100).reshape((10, 10))
    x = np.concatenate((x, x), axis=0)
    clusters = cookie_cutter(x, 3, 'ndiff')
    assert len(clusters) == 10


def test_cookie_cutter_cpp():
    x = np.arange(10000).reshape((100, 100))
    x = np.concatenate((x, x, x, x, x), axis=0)
    clusters = cookie_cutter_i2(x, 3, 'ndiff')
    assert len(clusters) == 100


def test_cookie_cutter_update_cpp():
    x = np.arange(100).reshape((10, 10))
    x = np.concatenate((x, x, x, x, x), axis=0)
    clusters = cookie_cutter_i2(x, 3, 'ndiff')
    assert len(clusters) == 10

    x = np.concatenate((x[clusters, :],
                        np.arange(60).reshape((6, 10))), axis=0)
    clusters = range(len(clusters))  # renumber clusters
    clusters = cookie_cutter_update_i2(x, 3, clusters, 'ndiff')
    assert len(clusters) == 10

    x = np.concatenate((x[clusters, :],
                        np.arange(300).reshape((30, 10))), axis=0)
    clusters = range(len(clusters))
    clusters = cookie_cutter_update_i2(x, 3, clusters, 'ndiff')
    assert len(clusters) == 30

    x = np.concatenate((x[clusters, :],
                        np.arange(10).reshape((1, 10))), axis=0)
    clusters = range(len(clusters))
    print(x)
    clusters = cookie_cutter_update_i2(x, 3, clusters, 'ndiff')
    print(clusters)
    assert len(clusters) == 30

    x = np.concatenate((x[clusters, :],
                        np.arange(300, 1000).reshape((70, 10))), axis=0)
    clusters = range(len(clusters))
    clusters = cookie_cutter_update_i2(x, 3, clusters, 'ndiff')
    assert len(clusters) == 100
