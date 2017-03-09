import os
from rif_cpp.sampling.orientations import read_karney_orientation_file
import pandas as pd
# import numpy as np


import sys
if sys.version_info[0] < 3:
    from io import BytesIO
    StrIO = BytesIO
else:
    from io import StringIO
    StrIO = StringIO


DATA_PATH = 'data/orientations/karney/'

with open(DATA_PATH + 'index.dat') as fin:
    karney_index_str = fin.read()
karney_index = pd.read_csv(StrIO(karney_index_str), sep='\s+')


def quats_from_karney_file(fname):
    q, w = read_karney_orientation_file(fname)
    return q, w


def karney_name_by_radius(cr):
    i = sum(karney_index.radius > cr)
    if i == karney_index.shape[0]:
        i -= 1
    return karney_index.iloc[i, 0]


def quaternion_set_with_covering_radius_degrees(cr=63):
    print(os.getcwd())
    fname = DATA_PATH + karney_name_by_radius(cr) + '.grid.gz'
    return quats_from_karney_file(fname)


def quaternion_set_by_name(name):
    fname = DATA_PATH + name + '.grid.gz'
    assert name in karney_index.name.values
    return quats_from_karney_file(fname)


def filter_quaternion_set_axis_within(quats, axis, angle):
    raise NotImplemented
    return quats
