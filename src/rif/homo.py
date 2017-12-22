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


def htrans(trans, dtype='f4'):
    trans = np.asanyarray(trans)
    if trans.shape[-1] != 3:
        raise ValueError('trans should be shape (..., 3)')
    tileshape = trans.shape[:-1] + (1, 1)
    t = np.tile(np.identity(4, dtype), tileshape)
    t[..., :3, 3] = trans
    return t
