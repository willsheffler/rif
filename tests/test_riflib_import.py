import riflib


def test_import_riflib():
    print(riflib.__file__)
    assert riflib.__version__ == '0.0.1'


def test_math():
    import riflib, sys
    print('------------------------')
    for p in sys.path:
        print('"{}"'.format(p))
    print('------------------------')
    print(riflib.__file__)
    assert riflib.add(1, 2) == 4
    assert riflib.subtract(1, 2) == -1
    assert riflib.mult(3, 3) == 9
