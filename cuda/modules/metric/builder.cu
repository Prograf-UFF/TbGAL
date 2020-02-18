#ifndef BUILDER_CU
#define BUILDER_CU

#include "../../common.cu"
#include "metric.cu"
#include <boost/python.hpp>

namespace python = boost::python;

BOOST_PYTHON_MODULE(metric) {

	python::class_<Metric, boost::noncopyable>("Metric", python::no_init);

	python::def("teste", teste<EuclideanMetric>);
	python::def("teste", teste<SignedMetric>);

	python::class_<OrthogonalMetric, boost::noncopyable>("OrthogonalMetric", python::no_init);

	python::class_<EuclideanMetric, python::bases<OrthogonalMetric> >("EuclideanMetric")
  	.def("diagonal_entry", &Metric::diagonal_entry)
		.def("metric_factor", &Metric::metric_factor)
		.def("matrix_entry", &Metric::matrix_entry)
    ;

	python::class_<SignedMetric, python::bases<OrthogonalMetric> >("SignedMetric", python::init<IndexType, IndexType>())
  	.def("diagonal_entry", &Metric::diagonal_entry)
		.def("metric_factor", &Metric::metric_factor)
		.def("matrix_entry", &Metric::matrix_entry)
    ;


}
#endif
