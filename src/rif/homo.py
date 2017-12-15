import numpy as np


def guess_is_degrees(angle):
    return np.max(np.abs(angle)) > 2 * np.pi


def axis_of(xforms):
    return np.stack((
        xforms[..., 2, 1] - xforms[..., 1, 2],
        xforms[..., 0, 2] - xforms[..., 2, 0],
        xforms[..., 1, 0] - xforms[..., 0, 1]), axis=-1)


def rotation33(axis, angle, degrees='auto', dtype='f4', shape=(3, 3)):
    axis = np.array(axis, dtype=dtype)
    angle = np.array(angle, dtype=dtype)
    if degrees is 'auto': degrees = guess_is_degrees(angle)
    angle = angle * np.pi / 180.0 if degrees else angle
    if (axis.shape[:-1] != angle.shape):
        raise ValueError('axis and angle not compatible')
    axis /= np.linalg.norm(axis, axis=-1)[..., np.newaxis]
    a = np.cos(angle / 2.0)
    axis *= -np.sin(angle / 2)[..., np.newaxis]
    b, c, d = axis[..., 0], axis[..., 1], axis[..., 2]
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    rot33 = np.zeros(angle.shape + shape, dtype=dtype)
    rot33[..., 0, 0] = aa + bb - cc - dd
    rot33[..., 0, 1] = 2 * (bc + ad)
    rot33[..., 0, 2] = 2 * (bd - ac)
    rot33[..., 1, 0] = 2 * (bc - ad)
    rot33[..., 1, 1] = aa + cc - bb - dd
    rot33[..., 1, 2] = 2 * (cd + ab)
    rot33[..., 2, 0] = 2 * (bd + ac)
    rot33[..., 2, 1] = 2 * (cd - ab)
    rot33[..., 2, 2] = aa + dd - bb - cc
    return rot33


def homo_rotation(axis, angle, center=None, dtype='f4', **args):
    center = (np.array([0, 0, 0], dtype=dtype) if center is None
              else np.array(center, dtype=dtype))
    axis = np.array(axis, dtype=dtype)
    angle = np.array(angle, dtype=dtype)
    r = rotation33(axis, angle, dtype=dtype, shape=(4, 4), **args)
    x, y, z = center[..., 0], center[..., 1], center[..., 2]
    r[..., 0, 3] = x - r[..., 0, 0] * x - r[..., 0, 1] * y - r[..., 0, 2] * z
    r[..., 1, 3] = y - r[..., 1, 0] * x - r[..., 1, 1] * y - r[..., 1, 2] * z
    r[..., 2, 3] = z - r[..., 2, 0] * x - r[..., 2, 1] * y - r[..., 2, 2] * z
    r[..., 3, 3] = 1
    return r
