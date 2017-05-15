Contributing
================

.. inclusion-marker-do-not-remove



conding conventions?
--------------------
- clang-format and pep8 everything please (hook it to your editor)
- python style as laid out in Fluent Python
- use libraries hdf5, numpy, pandas

Build System
-------------

This project *should* support
- python setup.py
- pip
- cmake (needs help finding a numpy library)
- tools/build_and_test.py
 - should be refactored to use click (or other)
 - should dispatch on file pattern (use a library?) *.pybind.cpp *_test.py

Adding wrapped c++ code
--------------------------

This section needs an audit.

Create a file called path/to/my_module.pybind.cpp which defines a function

.. code-block:: C++

    void RIFLIB_PYBIND_something_unique_path_is_good(py::module &m) {
        // pybind11 code here, adding to py::module m
    }

http://pybind11.readthedocs.io/en/master/

This function will be auto-detected by the build system and called with a pybind module rif.path.to.my_module, i.e. module namespaces based on the file path. If your filename starts in caps: path/to/MyClass.pybind.cpp, the module m will instead be rif.path.to. That is, build system will assume this is a class that you will bind  and should not be a module namespace.

contributors
-------------
- will sheffler
- alex ford (name and advice)
- david baker (funding)
