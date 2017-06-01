import pytest
from rif.rotamer.richardson import get_rotamer_space


# @pytest.mark.xfail
def test_richardson_space():
    rotspace = get_rotamer_space()
    assert rotspace.shape == (163, 29)
    rotspace2 = get_rotamer_space()
    assert rotspace is rotspace2
