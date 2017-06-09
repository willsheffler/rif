import rif_cpp
from rif_cpp.eigen_types import *

__all__ = 'V3 M3 X3'.split()


for cls in vars(rif_cpp.eigen_types):
    prefix = "V3".split()
    if not any(cls.startswith(p) for p in prefix):
        continue

    def cls_init(self, arg=[0, 0, 0], b=None, c=None):
        if c is not None:
            assert b is not None
            assert not hasattr(arg, '__iter__')
            assert not hasattr(b, '__iter__')
            arg = [arg, b, c]
        self[0] = arg[0]
        self[1] = arg[1]
        self[2] = arg[2]

    def cls_repr(self):
        return "V3([%f, %f, %f])" % (self[0], self[1], self[2])

    vars()[cls].__init__ = cls_init
    vars()[cls].__repr__ = cls_repr
