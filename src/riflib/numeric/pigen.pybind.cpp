#include <pybind11/pybind11.h>
#include <pybind11/operators.h>
#include <boost/lexical_cast.hpp>
#include <Eigen/Dense>

namespace py = pybind11;

template<class C>
std::string str(C const & c){
    return boost::lexical_cast<std::string>(c);
}


void RIFLIB_PYBIND_numeric_pigen(py::module & pigen) {

    using namespace Eigen;

    typedef Vector3f v3f;
    typedef Matrix<float,3,3> m33f;
    typedef Transform<float,3,Affine> x3f;

    float const fprec = NumTraits<float>::dummy_precision();

    py::class_<Vector3f>( pigen, "Vector3f" )
        .def( py::init<float,float,float>() )
        .def( py::self == py::self )
        .def( "__add__", [](v3f u, v3f v)->v3f{ return u+v; } )
        .def( "__sub__", [](v3f u, v3f v)->v3f{ return u-v; } )
        .def( "__mul__", [](float u, v3f v)->v3f{ return u*v; } )
        .def( "__div__", [](float u, v3f v)->v3f{ return u*v; } )
        .def( "__mul__", [](v3f u, float v)->v3f{ return u*v; } )
        .def( "__div__", [](v3f u, float v)->v3f{ return u*v; } )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        // .def( py::self + float() )
        // .def( py::self - float() )
        .def( "norm", [](v3f const & v){ return v.norm(); } )
        .def( "__getitem__", [](const v3f & v, int i) { return v[i]; } )
        .def( "isApprox", [](v3f const & m, v3f const & n, float prec){
            return m.isApprox(n,prec); }, py::arg("value"), py::arg("prec")=fprec )
        .def( "__repr__", [](const v3f & v) {
            return "Vector3f("+str(v[0])+", "+str(v[1])+", "+str(v[2])+")"; } );
    ;

    py::class_<m33f>( pigen, "Matrix33f" )
        .def( "__init__",
            [](m33f &instance ) { instance = m33f::Identity(); } )
        .def( "__add__", [](m33f u, m33f v)->m33f{ return u+v; } )
        .def( "__sub__", [](m33f u, m33f v)->m33f{ return u-v; } )
        .def( "__mul__", [](m33f u, m33f v)->m33f{ return u*v; } )
        .def( "__mul__", [](m33f u, v3f v)->v3f{ return u*v; } )
        .def( py::self == py::self )
        .def( py::self += py::self )
        .def( py::self -= py::self )
        .def( py::self *= py::self )
        .def( "__repr__", &str<m33f> )
        .def( "inverse", [](m33f const & m){ return m33f(m.inverse()); } )
        .def( "isApprox", [](m33f const & m, m33f const & n, float prec){
            return m.isApprox(n,prec); }, py::arg("value"), py::arg("prec")=fprec )
    ;

    py::class_<x3f>( pigen, "Transform3f" )
        .def( "__init__",
            [](x3f &instance ) { instance = x3f::Identity(); } )
        // .def( py::self == py::self )
        .def( py::self * py::self )
        // .def( py::self *= py::self )
        .def( "inverse", [](x3f const & x){ return x.inverse(); } )
        .def( "isApprox", &x3f::isApprox, py::arg("value"), py::arg("prec")=fprec )
        .def( "__repr__", [](x3f const & x){
            return str(x.rotation()) + '\n' + str(x.translation().transpose());
        } )
    ;

}
