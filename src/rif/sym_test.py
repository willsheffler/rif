from .sym import *
from homog import angle_degrees


def test_symm():
    assert tetrahedral_frames.shape == (12, 4, 4)
    assert octahedral_frames.shape == (24, 4, 4)
    assert icosahedral_frames.shape == (60, 4, 4)
    x = np.concatenate([tetrahedral_frames,
                        octahedral_frames,
                        icosahedral_frames])
    assert np.all(x[..., 3, 3] == 1)
    assert np.all(x[..., 3, :3] == 0)
    assert np.all(x[..., :3, 3] == 0)
