from rif_cpp.numeric.lattice import *
import rif

BCC3 = BCC3f8i8
BCC4 = BCC4f8i8
BCC5 = BCC5f8i8
BCC6 = BCC6f8i8
BCC7 = BCC7f8i8
BCC8 = BCC8f8i8
BCC9 = BCC9f8i8
BCC10 = BCC10f8i8


for name in dir(rif.numeric.lattice):
    if name.startswith("BCC"):
        # print(name)
        cls = getattr(rif.numeric.lattice, name)

        def bcc_center(self, idx):
            r = self._center_impl(idx)
            return r
        cls.center = bcc_center

        def bcc_index(self, arg):
            if not arg.dtype.fields:
                arg = arg.view(rif.util.sa_dtype[self.dim, arg.dtype])
            r = self._index_impl(arg)
            assert r.shape[-1] == 1
            return r.reshape(r.shape[:-1])
        cls.index = bcc_index
