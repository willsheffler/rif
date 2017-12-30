import numpy as np


def guess_is_degrees(angle):
    return np.max(np.abs(angle)) > 2 * np.pi


def fast_axis_of(xforms):
    return np.stack((
        xforms[..., 2, 1] - xforms[..., 1, 2],
        xforms[..., 0, 2] - xforms[..., 2, 0],
        xforms[..., 1, 0] - xforms[..., 0, 1]), axis=-1)


def axis_angle_of(xforms):
    axis = fast_axis_of(xforms)
    four_sin2 = np.sum(axis**2, axis=-1)
    axis = axis / np.linalg.norm(axis, axis=-1)[..., np.newaxis]
    sin_angl = np.sqrt(four_sin2 / 4)
    cos_angl = np.trace(xforms, axis1=-1, axis2=-2) / 2 - 1
    angl = np.arctan2(sin_angl, cos_angl)
    return axis, angl


def angle_of(xforms):
    axis = fast_axis_of(xforms)
    four_sin2 = np.sum(axis**2, axis=-1)
    sin_angl = np.sqrt(four_sin2 / 4)
    cos_angl = np.trace(xforms, axis1=-1, axis2=-2) / 2 - 1
    angl = np.arctan2(sin_angl, cos_angl)
    return angl


def center_of(xforms):
    # (X-I) * x == 0
    c = np.identity(4)
    raise NotImplementedError
    return c


def rot(axis, angle, degrees='auto', dtype='f4', shape=(3, 3)):
    axis = np.array(axis, dtype=dtype)
    angle = np.array(angle, dtype=dtype)
    if degrees is 'auto': degrees = guess_is_degrees(angle)
    angle = angle * np.pi / 180.0 if degrees else angle
    if axis.shape and angle.shape and (axis.shape[:-1] != angle.shape):
        raise ValueError('axis and angle not compatible')
    axis /= np.linalg.norm(axis, axis=-1)[..., np.newaxis]
    a = np.cos(angle / 2.0)
    tmp = axis * -np.sin(angle / 2)[..., np.newaxis]
    b, c, d = tmp[..., 0], tmp[..., 1], tmp[..., 2]
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    outshape = angle.shape if angle.shape else axis.shape[:-1]
    rot3 = np.zeros(outshape + shape, dtype=dtype)
    rot3[..., 0, 0] = aa + bb - cc - dd
    rot3[..., 0, 1] = 2 * (bc + ad)
    rot3[..., 0, 2] = 2 * (bd - ac)
    rot3[..., 1, 0] = 2 * (bc - ad)
    rot3[..., 1, 1] = aa + cc - bb - dd
    rot3[..., 1, 2] = 2 * (cd + ab)
    rot3[..., 2, 0] = 2 * (bd + ac)
    rot3[..., 2, 1] = 2 * (cd - ab)
    rot3[..., 2, 2] = aa + dd - bb - cc
    return rot3


def hrot(axis, angle, center=None, dtype='f4', **args):
    axis = np.array(axis, dtype=dtype)
    angle = np.array(angle, dtype=dtype)
    center = (np.array([0, 0, 0], dtype=dtype) if center is None
              else np.array(center, dtype=dtype))
    r = rot(axis, angle, dtype=dtype, shape=(4, 4), **args)
    x, y, z = center[..., 0], center[..., 1], center[..., 2]
    r[..., 0, 3] = x - r[..., 0, 0] * x - r[..., 0, 1] * y - r[..., 0, 2] * z
    r[..., 1, 3] = y - r[..., 1, 0] * x - r[..., 1, 1] * y - r[..., 1, 2] * z
    r[..., 2, 3] = z - r[..., 2, 0] * x - r[..., 2, 1] * y - r[..., 2, 2] * z
    r[..., 3, 3] = 1
    return r


def hray(origin, direction):
    s = np.broadcast(origin, direction).shape
    r = np.empty(s[:-1] + (2, 4))
    r[..., 0, :3] = origin
    r[..., 0, 3] = 1
    r[..., 1, :3] = hnormalized(direction)
    r[..., 1, 3] = 0
    return r


def htrans(trans, dtype='f4'):
    trans = np.asanyarray(trans)
    if trans.shape[-1] != 3:
        raise ValueError('trans should be shape (..., 3)')
    tileshape = trans.shape[:-1] + (1, 1)
    t = np.tile(np.identity(4, dtype), tileshape)
    t[..., :3, 3] = trans
    return t


def hdot(a, b):
    a = np.asanyarray(a)
    b = np.asanyarray(b)
    return np.sum(a * b, axis=-1)


def hnorm(a):
    a = np.asanyarray(a)
    return np.sqrt(np.sum(a * a, axis=-1))


def hnorm2(a):
    a = np.asanyarray(a)
    return np.sum(a * a, axis=-1)


def hnormalized(a):
    a = np.asanyarray(a)
    return a / hnorm(a)[..., None]


def is_valid_rays(r):
    r = np.asanyarray(r)
    if r.shape[-2:] != (2, 4): return False
    if np.any(r[..., :, 3] != (1, 0)): return False
    if np.any(abs(np.linalg.norm(r[..., 1, :3], axis=-1) - 1) > 0.000001):
        return False
    return True


def random_rays(shape=(), cen=(0, 0, 0), sdev=1):
    cen = np.asanyarray(cen)
    cen = cen + np.random.randn(*(shape + (3,))) * sdev
    norm = np.random.randn(*(shape + (3,)))
    norm /= np.linalg.norm(norm, axis=-1)[..., np.newaxis]
    a = np.stack([cen, norm], axis=-2)
    b = np.tile(np.array([1, 0]).reshape(2, 1), a.shape[:-2] + (1, 1))
    return np.concatenate([a, b], axis=-1)


def proj_perp(u, v):
    u = np.asanyarray(u)
    v = np.asanyarray(v)
    return v - hdot(u, v)[..., None] / hnorm2(u)[..., None] * u


def point_in_plane(plane, pt):
    return np.abs(hdot(plane[..., 1, :3], pt - plane[..., 0, :3])) < 0.000001


def ray_in_plane(plane, ray):
    return (point_in_plane(plane, ray[..., 0, :3]) *
            point_in_plane(plane, ray[..., 0, :3] + ray[..., 1, :3]))


def intersect_planes(plane1, plane2):
    """intersect_Planes: find the 3D intersection of two planes
       Input:  two planes represented by rays shape=(..., 2, 4)
       Output: *L = the intersection line (when it exists)
       Return: rays shape=(...,2,4), status
               0 = intersection returned
               1 = disjoint (no intersection)
               2 = the two  planes coincide

    """
    if not is_valid_rays(plane1): raise ValueError('invalid plane1')
    if not is_valid_rays(plane2): raise ValueError('invalid plane2')
    shape1, shape2 = np.array(plane1.shape), np.array(plane2.shape)
    if np.any((shape1 != shape2) * (shape1 != 1) * (shape2 != 1)):
        raise ValueError('incompatible shapes for plane1, plane2:')

    p1, n1 = plane1[..., 0, :3], plane1[..., 1, :3]
    p2, n2 = plane2[..., 0, :3], plane2[..., 1, :3]
    shape = tuple(np.maximum(plane1.shape, plane2.shape))

    # print('n1', n1.shape)
    # print('n2', n2.shape)
    u = np.cross(n1, n2)
    abs_u = np.abs(u)
    # print('u', u.shape)
    planes_parallel = np.sum(abs_u, axis=-1) < 0.000001
    p2_in_plane1 = point_in_plane(plane1, p2)
    # print('pp', planes_parallel.shape)
    # print('pinp', p2_in_plane1.shape)
    status = np.zeros(shape[:-2])
    status[planes_parallel] = 1
    status[planes_parallel * p2_in_plane1] = 2

    d1 = -hdot(n1, p1)
    d2 = -hdot(n2, p2)
    amax = np.argmax(abs_u, axis=-1)
    sel0, sel1, sel2 = amax == 0, amax == 1, amax == 2
    # print('u', u, 'sels', np.sum(sel0), np.sum(sel1), np.sum(sel2))
    n1a, n2a, d1a, d2a, ua = (x[sel0] for x in (n1, n2, d1, d2, u))
    n1b, n2b, d1b, d2b, ub = (x[sel1] for x in (n1, n2, d1, d2, u))
    n1c, n2c, d1c, d2c, uc = (x[sel2] for x in (n1, n2, d1, d2, u))
    ay = (d2a * n1a[..., 2] - d1a * n2a[..., 2]) / ua[..., 0]
    az = (d1a * n2a[..., 1] - d2a * n1a[..., 1]) / ua[..., 0]
    bz = (d2b * n1b[..., 0] - d1b * n2b[..., 0]) / ub[..., 1]
    bx = (d1b * n2b[..., 2] - d2b * n1b[..., 2]) / ub[..., 1]
    cx = (d2c * n1c[..., 1] - d1c * n2c[..., 1]) / uc[..., 2]
    cy = (d1c * n2c[..., 0] - d2c * n1c[..., 0]) / uc[..., 2]
    # print(shape[:-2] + (3,))
    isect_pt = np.empty(shape[:-2] + (3,), dtype=plane1.dtype)
    # print('sel0', sel0.shape)
    # print('isect_pt', isect_pt.shape)
    isect_pt[sel0, 0] = 0
    isect_pt[sel0, 1] = ay
    isect_pt[sel0, 2] = az
    isect_pt[sel1, 0] = bx
    isect_pt[sel1, 1] = 0
    isect_pt[sel1, 2] = bz
    isect_pt[sel2, 0] = cx
    isect_pt[sel2, 1] = cy
    isect_pt[sel2, 2] = 0
    isect = hray(isect_pt, u)
    return isect, status


#    switch (maxc) {  # select max coordinate
#      case 1:        # intersect with x=0
#        iP.x() = 0;
#        iP.y() = (d2 * norm1.z() - d1 * norm2.z()) / u.x();
#        iP.z() = (d1 * norm2.y() - d2 * norm1.y()) / u.x();
#        break;
#      case 2:  # intersect with y=0
#        iP.x() = (d1 * norm2.z() - d2 * norm1.z()) / u.y();
#        iP.y() = 0;
#        iP.z() = (d2 * norm1.x() - d1 * norm2.x()) / u.y();
#        break;
#      case 3:  # intersect with z=0
#        iP.x() = (d2 * norm1.y() - d1 * norm2.y()) / u.z();
#        iP.y() = (d1 * norm2.x() - d2 * norm1.x()) / u.z();
#        iP.z() = 0;
#    }
#    L->P0 = iP;
#    L->P1 = iP + u;
#    return 2;
#  }
#
#  void rotation_axis(Vector &axis, Vector &cen, T &angle) const {
#    axis = numeric::rotation_axis<T>(R, angle);
#    Vector const p1((T)-32.09501046777237, (T)03.36227004372687,
#                    (T)35.34672781477340);  # random...
#    Vector const p2((T)21.15113978202345, (T)12.55664537217840,
#                    (T)-37.48294301885574);  # random...
#    Vector const q1((*this) * p1);
#    Vector const q2((*this) * p2);
#    Vector const n1 = (q1 - p1).normalized();
#    Vector const n2 = (q2 - p2).normalized();
#    Vector const c1 = (p1 + q1) / T(2.0);
#    Vector const c2 = (p2 + q2) / T(2.0);
#    Plane Pn1 = {n1, c1};
#    Plane Pn2 = {n2, c2};
#    Line Linter;
#    int inter_case = intersect3D_2Planes(Pn1, Pn2, &Linter);
#    switch (inter_case) {
#      case 0:
#        cen =
#            Vector(std::numeric_limits<T>::max(), std::numeric_limits<T>::max(),
#                   std::numeric_limits<T>::max());
#        break;
#      case 1:
#        cen =
#            Vector(std::numeric_limits<T>::min(), std::numeric_limits<T>::min(),
#                   std::numeric_limits<T>::min());
#        break;
#      case 2:
#        Vector Laxis = (Linter.P1 - Linter.P0).normalized();
#        if (-0.9999 < Laxis.dot(axis) && Laxis.dot(axis) < 0.9999) {
#          # std::cout << Laxis << std::endl;
#          # std::cout << axis  << std::endl;
#          # std::cout << angle << std::endl;
#          utility_exit_with_message("bad axis");
#        }
#        cen = Linter.P0;
#        break;
#    }
#  }
