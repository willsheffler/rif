import numpy as np
import pandas as pd
import xarray as xr


def numpy_array_dtype():
    x = np.zeros(2, dtype='(4,4)f4')
    p = np.zeros(2, dtype='(4,)f4')
    print(x[1])
    print('x.shape', x.shape)
    print('xp')


def make_ray():
    r = np.random.randn(2, 4).astype('f4') * 3
    r[0, 3] = 1
    r[1, 3] = 0
    r[1, :] /= np.linalg.norm(r[1, :])
    # print(np.linalg.norm(r[0, :3]))
    # print(np.linalg.norm(r[1, :3]))
    return r


def xr_attempt2():
    N = 5
    da_1 = xr.DataArray([float(i) for i in range(N)], [
        ('index', range(N)), ])
    da_r = xr.DataArray([make_ray() for i in range(N)], [
        ('index', range(N)),
        ('r_odr', ['orig', 'dirn']),
        ('r_hmo', list('xyzw'))])
    da_x = xr.DataArray(
        data=[np.random.randn(4, 4).astype('f4') for i in range(N)],
        # coords=[('index', range(N)), ('x_row', list('xyzw')), ('x_col', list('xyzw')), ],
        dims='index x_row x_col'.split())
    ds = xr.Dataset()
    ds['scalar'] = da_1
    ds['ray'] = da_r
    ds['xform'] = da_x
    # ds = xr.Dataset({'scalar': da_1, 'ray': da_r, 'xform': da_x})

    print('================================')
    print(ds)
    print('================================')
    print(ds.isel(index=4))  # row
    print(ds[dict(index=2)])
    print('================================')
    print(np.sqrt((ds.ray**2).sum('r_hmo')))
    print(type(ds.ray), ds.ray.shape)
    print('================================')
    ds['scalar2'] = xr.DataArray(np.arange(N, dtype='f4'))
    print(ds)


def xr_attempt1():
    x1 = np.arange(16).reshape(4, 4)
    x2 = np.arange(16, 32).reshape(4, 4)
    r1 = np.arange(6).reshape(2, 3)
    r2 = np.arange(6, 12).reshape(2, 3)
    pdat = pd.DataFrame({'index': [42, 51], 'ray': [r1, r2]})
    print(pdat)
    print(pdat.dtypes)
    print()
    xdat1 = xr.DataArray([37, 45], coords=[('x', [1, 2])])
    print(xdat1)

    xdat2 = xr.DataArray(
        [r1, r2], [('x', [1, 2]), ('y0', [1, 2, 3]), ('y1', [1, 2])])
    print(xdat2)


def main():
    # xr_attempt2()
    numpy_array_dtype()


if __name__ == '__main__':
    main()
