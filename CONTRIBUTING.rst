Contributing
================

.. inclusion-marker-do-not-remove

Write good code
----------------

clang-format and pep8 everything please. please adhere to the following:

- better to ask forgiveness than permission
- use numpy or eigen over custom data structures as much as possible
- use hdf5 for all data

Build System
-------------

This project uses a mix of python setuptools, cmake, and custom python scripts. It's a bit janky, but there seems to be no good solution.

Adding wrapped c++ code
--------------------------

Do not put code in __init__.py files, these are auto-generated.

Create a file called path/to/my_module.pybind.cpp which defiles a function

.. code-block:: C++

    void RIFLIB_PYBIND_something_unique_path_is_good(py::module &m) {
        // pybind11 code here, adding to py::module m
    }

http://pybind11.readthedocs.io/en/master/

This function will be auto-detected by the build system and called with a pybind module rif.path.to.my_module, i.e. module namespaces based on the file path. If your filename starts in caps: path/to/MyClass.pybind.cpp, the module m will instead be rif.path.to. That is, build system will assume this is a class that you will bind  and should not be a module namespace.
