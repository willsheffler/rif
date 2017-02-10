#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>

namespace py = pybind11;

template<class C>
std::string str(C const & c){
    return boost::lexical_cast<std::string>(c);
}

PYBIND11_PLUGIN(pigen) {

    using namespace Eigen;

    py::module m("pigen", R"pbdoc(
        pigen docs
        -----------------------

        .. currentmodule:: pigen

        .. autosummary::
           :toctree: _generate

    )pbdoc");

    typedef Vector3f v3f;
    typedef Matrix<float,3,3> m33f;
    typedef Transform<float,3,Affine> x3f;

    py::class_<Vector3f>( m, "Vector3f" )
        .def( py::init<float,float,float>() )
        .def( py::self == py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( "__getitem__", [](const v3f & v, int i) { return v[i]; } )
        // .def( "isApprox", &v3f::isApprox, py::arg("value"), py::arg("prec")=NumTraits<float>::dummy_precision() )
        .def( "__repr__", [](const v3f & v) {
            return "Vector3f("+str(v[0])+", "+str(v[1])+", "+str(v[2])+")"; } );
    ;

    py::class_<m33f>( m, "Matrix33f" )
        .def( "__init__",
            [](m33f &instance ) { instance = m33f::Identity(); } )
        .def( py::self == py::self )
        .def( py::self + py::self )
        .def( py::self - py::self )
        .def( py::self * py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self *= py::self )
        .def( "__repr__", &str<m33f> )
        .def( "isApprox", [](m33f const & m, m33f const & n, float prec){
            return m.isApprox(n,prec)        } )
    ;

    py::class_<x3f>( m, "Transform3f" )
        .def( "__init__",
            [](x3f &instance ) { instance = x3f::Identity(); } )
        // .def( py::self == py::self )
        .def( py::self * py::self )
        // .def( py::self *= py::self )
        .def( "inverse", [](x3f const & x){ return x.inverse(); } )
        .def( "isApprox", &x3f::isApprox, py::arg("value"), py::arg("prec")=NumTraits<float>::dummy_precision() )
        .def( "__repr__", [](x3f const & x){
            return str(x.rotation()) + '\n' + str(x.translation().transpose());
        } )
    ;

    return m.ptr();

}
