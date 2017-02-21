import riflib


def test_import_riflib():
    print(riflib.__file__)
    assert hasattr(riflib, '__version__')


def test_math():
    # print(riflib.__file__)
    assert riflib.test.example.add(1, 2) == 3
