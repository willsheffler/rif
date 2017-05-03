from rif_cpp.eigen_types import *


def V3_init(self, arg):
    self[0] = arg[0]
    self[1] = arg[1]
    self[2] = arg[2]


def V3_repr(self):
    return "V3([%f, %f, %f])" % (self[0], self[1], self[2])


V3.__init__ = V3_init
V3.__repr__ = V3_repr
