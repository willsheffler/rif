
"""docstring for rif.util
"""

from rif_cpp.util import *
from numpy import dtype

sa_dtype = dict()
for i in range(1, 11):
    for dt in "f4 f8 u4 u8 i4 i8".split():
        sa_dtype[(i, dtype(dt))] = vars()['sa' + str(i) + dt + '_t']
        sa_dtype[(i, dt)] = vars()['sa' + str(i) + dt + '_t']
