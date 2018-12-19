#ifndef BUILDER_CU
#define BUILDER_CU

#include "../../common.cu"
#include "../metric/metric.cu"
#include "multivector.cu"
#include "operations.cu"
#include <boost/python.hpp>
#include <vector>
#include <boost/python/object.hpp>
#include <boost/python/stl_iterator.hpp>


BOOST_PYTHON_MODULE(multivector) {
	namespace python = boost::python;
	namespace operations = MultivectorOperations;


	python::class_<Multivector>("Multivector", python::no_init)
		.def("__repr__", &Multivector::to_string)
		.def("set_N", &Multivector::set_N)
			.staticmethod("set_N")
		.def("get_N", &Multivector::get_N)
			.staticmethod("get_N")
		.def("get_N_FULL", &Multivector::get_N_FULL)
			.staticmethod("get_N_FULL")
		.def("getComponent", &Multivector::getComponent)
		.def("getGrade", &Multivector::get_grade_blade)
		.def("get_component_max_projection", &Multivector::get_component_max_projection)

		.def("__eq__", &operations::is_equals)

		.def("__pos__", &operations::UNARY_PLUS)
		.def("__add__", &operations::ADD_SCALAR<int>)
		.def("__radd__", &operations::ADD_SCALAR<int>)
		.def("__add__", &operations::ADD_SCALAR<float>)
		.def("__radd__", &operations::ADD_SCALAR<float>)
		.def("__add__", &operations::ADD_SCALAR<double>)
		.def("__radd__", &operations::ADD_SCALAR<double>)
		.def("__add__", &operations::ADD)
		.def("__radd__", &operations::ADD)

		.def("__neg__", &operations::UNARY_MINUS)

		.def("__sub__", &operations::SUB_SCALAR<int>)
		.def("__rsub__", &operations::R_SUB_SCALAR<int>)
		.def("__sub__", &operations::SUB_SCALAR<float>)
		.def("__rsub__", &operations::R_SUB_SCALAR<float>)
		.def("__sub__", &operations::SUB_SCALAR<double>)
		.def("__rsub__", &operations::R_SUB_SCALAR<double>)
		.def("__sub__", &operations::SUB)
		.def("__rsub__", &operations::SUB)

		.def("__mul__", &operations::PROD<int>)
		.def("__rmul__", &operations::PROD<int>)
		.def("__mul__", &operations::PROD<float>)
		.def("__rmul__", &operations::PROD<float>)
		.def("__mul__", &operations::PROD<double>)
		.def("__rmul__", &operations::PROD<double>)

		.def("__xor__", &operations::OP)

		// .def("__invert__", &operations::REVERSE, python::return_value_policy<python::manage_new_object>())

		;


	python::def("e", &e);
	python::def("generate_T", &generate_T<EuclideanMetric>);
	python::def("build_tensor", &build_tensor<EuclideanMetric>, python::return_value_policy<python::manage_new_object>());
	python::def("extract_tensor", &extract_tensor, python::return_value_policy<python::manage_new_object>());
	python::def("get_GP_T", &Multivector::get_GP_T, python::return_value_policy<python::manage_new_object>());

	python::def("REVERSE", &operations::REVERSE);
	python::def("INVOLUTION", &operations::INVOLUTION);
	python::def("CONJUGATE", &operations::CONJUGATE);
	python::def("INVERSE", &operations::INVERSE);

	python::def("GP", &operations::GP_tensor);
	python::def("GP", &operations::GP);

	python::def("LCONT", &operations::LCONT);
	python::def("RCONT", &operations::RCONT);

	python::def("DOT", &operations::dot);
	python::def("NATIVE", &operations::native);

	python::def("SQR_NORM", &operations::SQR_NORM);//, python::return_value_policy<python::manage_new_object>());
	python::def("NORM", &operations::NORM);//, python::return_value_policy<python::manage_new_object>());
	python::def("IGP", &operations::IGP);

	python::def("SCP", &operations::SCP_tensor);
	python::def("SCP", &operations::SCP);

	python::def("take_grade", &operations::take_grade);
	python::def("fact_blade", &operations::FACT_BLADE<python::list>);
	python::def("fact_versor", &operations::FACT_VERSOR<python::list>);
	python::def("set_device", &setDevice);

}
#endif
