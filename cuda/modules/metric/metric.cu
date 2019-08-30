#ifndef METRIC_CUH
#define METRIC_CUH

#include "../../common.cu"

#include <iostream>
#include <limits>

// typedef unsigned long long int IndexType;
// typedef float CoeffType;

class Metric {

	public:
		__host__ __device__ virtual CoeffType matrix_entry(const IndexType &row, const IndexType &col) = 0;
		__host__ __device__ virtual CoeffType diagonal_entry(const IndexType &index) = 0;
		__host__ __device__ virtual CoeffType metric_factor(const IndexType &basis_A, const IndexType &basis_B) = 0;
};

class OrthogonalMetric : public Metric {

	public:

		__host__ __device__ virtual CoeffType diagonal_entry(const IndexType &index) = 0;

		__host__ __device__ virtual CoeffType metric_factor(const IndexType & basis_A, const IndexType & basis_B) = 0;

		__host__ __device__ virtual CoeffType matrix_entry(const IndexType & row, const IndexType & col) {
			if (row == col) {
				return diagonal_entry(row);
			} else {
				return 0;
			}
		}
};

class EuclideanMetric : public OrthogonalMetric {

	public:

		__host__ __device__ CoeffType diagonal_entry(const IndexType &index) {
			return 1;
		}

		__host__ __device__ CoeffType metric_factor(const IndexType &basis_A, const IndexType &basis_B) {
			return 1;
		}
};

class SignedMetric : public OrthogonalMetric {

	public:

		SignedMetric(const IndexType p, const IndexType q) : p(p), q(q) {
		}

		__host__ __device__ CoeffType diagonal_entry(const IndexType &index) {
			if (index <= p) {
				return +1;
			} else if (index <= p + q) {
				return -1;
			}
			return 0;
		}

		__host__ __device__ CoeffType metric_factor(const IndexType &basis_A, const IndexType &basis_B) {
			return this->diagonal_entry(basis_A) * this->diagonal_entry(basis_B);
		}

	private:
		IndexType p;
		IndexType q;

};
template<typename T>
CoeffType teste(T m) {
	return m.matrix_entry(2, 2);
}

// BOOST_PYTHON_MODULE(metric)
// {
// 	namespace python = boost::python;
//
// 	python::class_<Metric, boost::noncopyable>("Metric", python::no_init);
//
// 	python::def("teste", teste<EuclideanMetric>);
// 	python::def("teste", teste<SignedMetric>);
//
// 	python::class_<OrthogonalMetric, boost::noncopyable>("OrthogonalMetric", python::no_init);
//
// 	python::class_<EuclideanMetric, python::bases<OrthogonalMetric> >("EuclideanMetric")
//   	.def("diagonal_entry", &Metric::diagonal_entry)
// 		.def("metric_factor", &Metric::metric_factor)
// 		.def("matrix_entry", &Metric::matrix_entry)
//     ;
//
// 	python::class_<SignedMetric, python::bases<OrthogonalMetric> >("SignedMetric", python::init<IndexType, IndexType>())
//   	.def("diagonal_entry", &Metric::diagonal_entry)
// 		.def("metric_factor", &Metric::metric_factor)
// 		.def("matrix_entry", &Metric::matrix_entry)
//     ;
//
//
// }
#endif
