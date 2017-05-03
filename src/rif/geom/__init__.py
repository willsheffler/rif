
"""docstring for rif.geom
"""

from rif_cpp.geom import *


def Ray_init(self, orig=[0, 0, 0], dirn=[1, 0, 0]):
    for i in (0, 1, 2):
        self.orig[i] = orig[i]
        self.dirn[i] = dirn[i]
        self.dirn.normalize()


Ray.__init__ = Ray_init
