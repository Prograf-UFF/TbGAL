#include <boost/python/module.hpp> 
#include <boost/python/def.hpp>  

#include "utils.cpp"  
#include "../cpp/include/tbgal/using_Eigen.hpp"  
#include "../cpp/include/tbgal/assuming_EuclideanD.hpp"  
#include <boost/python.hpp>  
#include <boost/python/object.hpp>   
#include <boost/python/args.hpp>   
#include <boost/python/list.hpp>   
#include <boost/python/tuple.hpp>  
#include <boost/python/dict.hpp>  
#include "boost/python/stl_iterator.hpp"   
#include <boost/python/raw_function.hpp>  

// PRODUCT = tbgal::OuterProduct or tbgal::GeometricProduct  \ 
// OTHER = tbgal::OuterProduct or tbgal::GeometricProduct  \ 
// METRIC = tbgal::EuclideanMetricSpace<Dynamic, Dynamic>  \ 
// SCALAR_MAIN = SCALAR_MAIN  \ 
// SCALAR_EXTRA = SCALAR_EXTRA  \ 

#define BUILD_FACTORED_MULTIVECTOR_PYTHON(PRODUCT, OTHER_PRODUCT, METRIC, SCALAR_MAIN, SCALAR_EXTRA) \
    python::class_< tbgal::FactoredMultivector<SCALAR_MAIN, PRODUCT<METRIC<Dynamic, Dynamic>>> >("FactoredMultivector")  \ 
        // =========== just for prSCALAR_EXTRAing ===========   \ 
        .def("__repr__", &tbgal::FactoredMultivector<SCALAR_MAIN, PRODUCT<METRIC<Dynamic, Dynamic>> >::repr)  \ 
        // =========== OP ===========   \ 
        // MV_OP ^ MV_GP  \ 
        .def(python::self ^ python::other<tbgal::FactoredMultivector<SCALAR_MAIN, OTHER_PRODUCT<METRIC<Dynamic, Dynamic>> >>())  \ 
        // MV_OP ^ MV_OP  \ 
        .def(python::self ^ python::self)  \ 
        // MV_OP ^ SCALAR_EXTRA  \ 
        .def(python::self ^ python::other<SCALAR_EXTRA>())  \ 
        // SCALAR_EXTRA ^ MV_OP  \ 
        .def(python::other<SCALAR_EXTRA>() ^ python::self)  \ 
        // MV_OP ^ SCALAR_MAIN  \ 
        .def(python::self ^ python::other<SCALAR_MAIN>())  \ 
        // SCALAR_MAIN ^ MV_OP  \ 
        .def(python::other<SCALAR_MAIN>() ^ python::self)  \ 
        // =========== GP ===========   \ 
        // MV_OP * MV_GP  \ 
        .def(python::self * python::other<tbgal::FactoredMultivector<SCALAR_MAIN, OTHER_PRODUCT<METRIC<Dynamic, Dynamic>> >>())  \ 
        // MV_OP * MV_OP  \ 
        .def(python::self * python::self)  \ 
        // MV_OP * SCALAR_EXTRA  \ 
        .def(python::self * python::other<SCALAR_EXTRA>())  \ 
        // SCALAR_EXTRA * MV_OP  \ 
        .def(python::other<SCALAR_EXTRA>() * python::self)  \ 
        // MV_OP * SCALAR_MAIN  \ 
        .def(python::self * python::other<SCALAR_MAIN>())  \ 
        // SCALAR_MAIN * MV_OP  \ 
        .def(python::other<SCALAR_MAIN>() * python::self)  \ 
        // =========== UNARY ===========   \ 
        // + MV_OP  \ 
		.def(+python::self)  \ 
        // - MV_OP  \ 
		.def(-python::self)  \ 
        // MV_OP +  MV_OP  \ 
		.def(python::self + python::self)  \ 
        // MV_OP +  MV_GP  \ 
        .def(python::self + python::other<tbgal::FactoredMultivector<SCALAR_MAIN, OTHER_PRODUCT<METRIC<Dynamic, Dynamic>> > >())  \ 
        // ~ MV_OP  \ 
        .def(~python::self);  \ 
