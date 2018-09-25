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


	python::class_<Multivector, Multivector*>("Multivector", python::no_init)
		.def("__repr__", &Multivector::to_string)
		.def("set_N", &Multivector::set_N)
			.staticmethod("set_N")
		.def("get_N", &Multivector::get_N)
			.staticmethod("get_N")
		.def("get_N_FULL", &Multivector::get_N_FULL)
			.staticmethod("get_N_FULL")
		.def("getComponent", &Multivector::getComponent, python::return_value_policy<python::manage_new_object>())

		.def("get_component_max_projection", &Multivector::get_component_max_projection, python::return_value_policy<python::manage_new_object>())

		.def("__eq__", &operations::is_equals)

		.def("__pos__", &operations::UNARY_PLUS, python::return_value_policy<python::manage_new_object>())
		.def("__add__", &operations::ADD_SCALAR<int>, python::return_value_policy<python::manage_new_object>())
		.def("__radd__", &operations::ADD_SCALAR<int>, python::return_value_policy<python::manage_new_object>())
		.def("__add__", &operations::ADD_SCALAR<float>, python::return_value_policy<python::manage_new_object>())
		.def("__radd__", &operations::ADD_SCALAR<float>, python::return_value_policy<python::manage_new_object>())
		.def("__add__", &operations::ADD_SCALAR<double>, python::return_value_policy<python::manage_new_object>())
		.def("__radd__", &operations::ADD_SCALAR<double>, python::return_value_policy<python::manage_new_object>())
		.def("__add__", &operations::ADD, python::return_value_policy<python::manage_new_object>())
		.def("__radd__", &operations::ADD, python::return_value_policy<python::manage_new_object>())

		.def("__neg__", &operations::UNARY_MINUS, python::return_value_policy<python::manage_new_object>())

		.def("__sub__", &operations::SUB_SCALAR<int>, python::return_value_policy<python::manage_new_object>())
		.def("__rsub__", &operations::R_SUB_SCALAR<int>, python::return_value_policy<python::manage_new_object>())
		.def("__sub__", &operations::SUB_SCALAR<float>, python::return_value_policy<python::manage_new_object>())
		.def("__rsub__", &operations::R_SUB_SCALAR<float>, python::return_value_policy<python::manage_new_object>())
		.def("__sub__", &operations::SUB_SCALAR<double>, python::return_value_policy<python::manage_new_object>())
		.def("__rsub__", &operations::R_SUB_SCALAR<double>, python::return_value_policy<python::manage_new_object>())
		.def("__sub__", &operations::SUB, python::return_value_policy<python::manage_new_object>())
		.def("__rsub__", &operations::SUB, python::return_value_policy<python::manage_new_object>())

		.def("__mul__", &operations::PROD<int>, python::return_value_policy<python::manage_new_object>())
		.def("__rmul__", &operations::PROD<int>, python::return_value_policy<python::manage_new_object>())
		.def("__mul__", &operations::PROD<float>, python::return_value_policy<python::manage_new_object>())
		.def("__rmul__", &operations::PROD<float>, python::return_value_policy<python::manage_new_object>())
		.def("__mul__", &operations::PROD<double>, python::return_value_policy<python::manage_new_object>())
		.def("__rmul__", &operations::PROD<double>, python::return_value_policy<python::manage_new_object>())

		.def("__xor__", &operations::OP, python::return_value_policy<python::manage_new_object>())

		// .def("__invert__", &operations::REVERSE, python::return_value_policy<python::manage_new_object>())

		;


	python::def("e", &e, python::return_value_policy<python::manage_new_object>());
	python::def("generate_T", &generate_T<EuclideanMetric>);
	python::def("build_tensor", &build_tensor<EuclideanMetric>, python::return_value_policy<python::manage_new_object>());
	python::def("extract_tensor", &extract_tensor, python::return_value_policy<python::manage_new_object>());
	python::def("get_GP_T", &Multivector::get_GP_T, python::return_value_policy<python::manage_new_object>());

	python::def("REVERSE", &operations::REVERSE, python::return_value_policy<python::manage_new_object>());
	python::def("INVOLUTION", &operations::INVOLUTION, python::return_value_policy<python::manage_new_object>());
	python::def("CONJUGATE", &operations::CONJUGATE, python::return_value_policy<python::manage_new_object>());
	python::def("INVERSE", &operations::INVERSE, python::return_value_policy<python::manage_new_object>());

	python::def("GP", &operations::GP_tensor, python::return_value_policy<python::manage_new_object>());
	python::def("GP", &operations::GP, python::return_value_policy<python::manage_new_object>());

	python::def("LCONT", &operations::LCONT, python::return_value_policy<python::manage_new_object>());
	python::def("RCONT", &operations::RCONT, python::return_value_policy<python::manage_new_object>());

	python::def("DOT", &operations::dot, python::return_value_policy<python::manage_new_object>());

	python::def("SQR_NORM", &operations::SQR_NORM);//, python::return_value_policy<python::manage_new_object>());
	python::def("NORM", &operations::NORM);//, python::return_value_policy<python::manage_new_object>());
	python::def("IGP", &operations::IGP, python::return_value_policy<python::manage_new_object>());

	python::def("SCP", &operations::SCP_tensor, python::return_value_policy<python::manage_new_object>());
	python::def("SCP", &operations::SCP, python::return_value_policy<python::manage_new_object>());

	python::def("take_grade", &operations::take_grade, python::return_value_policy<python::manage_new_object>());
	python::def("fact_blade", &operations::FACT_BLADE);
	python::def("fact_versor", &operations::FACT_VERSOR);


}
#endif
