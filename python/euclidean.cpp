#include "../cpp/include/tbgal/using_Eigen.hpp"
#include "../cpp/include/tbgal/assuming_EuclideanD.hpp"
#include "wrapper_functions.hpp"

namespace python = boost::python;

//TODO preciso setar isso antes?
auto py_vector(python::tuple args, python::dict kwargs) {
    tbgal::EuclideanD::SPACE = EuclideanMetricSpace<Dynamic, Dynamic>(len(args));
    std::vector<double> args_as_container = py_list_to_std_vector<double>(args);
    auto a = tbgal::EuclideanD::vector(args_as_container.begin(), args_as_container.end());
    return a;
}

BOOST_PYTHON_MODULE(euclidean) {
    python::class_< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >("FactoredMultivector")
        // =========== just for printing =========== 
        .def("__repr__", &tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >::repr)
        // =========== OP =========== 
        // MV_GP ^ MV_OP
        .def(python::self ^ python::other<tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >>())
        // MV_GP ^ MV_GP
        .def(python::self ^ python::self)
        // MV_GP ^ int
        .def(python::self ^ python::other<int>())
        // int ^ MV_GP
        .def(python::other<int>() ^ python::self)
        // MV_GP ^ double
        .def(python::self ^ python::other<double>())
        // double ^ MV_GP
        .def(python::other<double>() ^ python::self)
        // =========== GP =========== 
        // MV_GP * MV_OP
        .def(python::self * python::other<tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >>())
        // MV_GP * MV_GP
        .def(python::self * python::other<tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >>())
        // MV_GP * int
        .def(python::self * python::other<int>())
        // int * MV_GP
        .def(python::other<int>() * python::self)
        // MV_GP * double
        .def(python::self * python::other<double>())
        // double * MV_GP
        .def(python::other<double>() * python::self)
        // =========== UNARY =========== 
        // + MV_GP
		.def(+python::self)
        // - MV_GP
		.def(-python::self)
        // MV_GP +  MV_GP
		.def(python::self + python::self)
        // MV_GP +  MV_OP
        .def(python::self + python::other<tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>>>())
        // ~ MV_GP
        .def(~python::self)
        ;
        


    python::class_< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >("FactoredMultivector")
        // =========== just for printing =========== 
        .def("__repr__", &tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >::repr)
        // =========== OP =========== 
        // MV_OP ^ MV_GP
        .def(python::self ^ python::other<tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >>())
        // MV_OP ^ MV_OP
        .def(python::self ^ python::self)
        // MV_OP ^ int
        .def(python::self ^ python::other<int>())
        // int ^ MV_OP
        .def(python::other<int>() ^ python::self)
        // MV_OP ^ double
        .def(python::self ^ python::other<double>())
        // double ^ MV_OP
        .def(python::other<double>() ^ python::self)
        // =========== GP =========== 
        // MV_OP * MV_GP
        .def(python::self * python::other<tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic> > >>())
        // MV_OP * MV_OP
        .def(python::self * python::self)
        // MV_OP * int
        .def(python::self * python::other<int>())
        // int * MV_OP
        .def(python::other<int>() * python::self)
        // MV_OP * double
        .def(python::self * python::other<double>())
        // double * MV_OP
        .def(python::other<double>() * python::self)
        // =========== UNARY =========== 
        // + MV_OP
		.def(+python::self)
        // - MV_OP
		.def(-python::self)
        // MV_OP +  MV_OP
		.def(python::self + python::self)
        // MV_OP +  MV_GP
        .def(python::self + python::other<tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>>>())
        // ~ MV_OP
        .def(~python::self)
        ;

    // vector(scalar...)
    python::def("vector", python::raw_function(py_vector) );

    // hip(MV_OP, MV_OP)
    python::def("hip", &py_hip< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // dot(MV_OP, MV_OP)
    python::def("dot", &py_dot< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // lcont(MV_OP, MV_OP)
    python::def("lcont", &py_lcont< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // rcont(MV_OP, MV_OP)
    python::def("rcont", &py_rcont< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // reverse(MV_OP)
    python::def("reverse", &py_reverse< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);
    // reverse(MV_GP)
    python::def("reverse", &py_reverse< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // inverse(MV_OP)
    python::def("inverse", &py_inverse< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);
    // inverse(MV_GP)
    python::def("inverse", &py_inverse< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // rnorm(MV_OP)
    python::def("rnorm", &py_rnorm< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);
    // rnorm(MV_GP)
    python::def("rnorm", &py_rnorm< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // rnorm_sqr(MV_OP)
    python::def("rnorm_sqr", &py_rnorm_sqr< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);
    // rnorm_sqr(MV_GP)
    python::def("rnorm_sqr", &py_rnorm_sqr< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // dual(MV_OP)
    python::def("dual", &py_dual< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // dual(MV_GP)
    python::def("undual", &py_undual< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);

    // sp(MV_OP, MV_OP)
    python::def("sp", &py_sp< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>> >);


}