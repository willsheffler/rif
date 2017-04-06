import rif
import rif_cpp


def test_import_rif():
    print((rif.__file__))
    print((rif_cpp.__file__))
    assert hasattr(rif, 'kinematics')
    assert hasattr(rif, 'sampling')
    # assert hasattr(rif, 'numeric')
    assert hasattr(rif, 'test')
