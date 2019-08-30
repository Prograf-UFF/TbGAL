#pragma once
#include <SparseTensor.h>
#include <cusp\print.h>
#include <cusp\array1d.h>
#include <thrust\tuple.h>
#include "utils.cuh"
#include <vector>
#include <chrono>
#include <type_traits>

#ifndef N
#error "Please define N before include multivector header"
#endif

#define N_FULL (N * (N+1)) >> 1

template <typename IndexType, typename DataType, typename ReturnType, typename MetricType>
struct GeometricProductTensorFunctor {
	GeometricProductTensorFunctor(const MetricType &metric, const short &axis) : metric(metric), axis(axis) {}

	short axis;
	MetricType metric;

	__host__ __device__ ReturnType operator() (const IndexType &j, const IndexType &i) {
		//		uint64_t a = dynamic_cast(metric)->metric_factor(1, 3);
		thrust::tuple<IndexType, IndexType, IndexType, DataType> to_return;
		IndexType K = 0;
		DataType val = 0;

		auto my_functor = idx2ji<IndexType>();
		thrust::tuple<IndexType, IndexType> basis_i = my_functor(i);
		thrust::tuple<IndexType, IndexType> basis_j = my_functor(j);

		IndexType u = thrust::get<0>(basis_j);
		IndexType v = thrust::get<1>(basis_j);
		IndexType x = thrust::get<0>(basis_i);
		IndexType y = thrust::get<1>(basis_i);

		if (i == 0) {
			K = j;
			val = 1;
		}
		else if (j == 0) {
			K = i;
			val = 1;
		}
		else if (u == x) {
			K = 0;
			if (u == 0) {
				if (v < y) {
					val = +1 * metric.diagonal_entry(u);
					K = ji2idx<IndexType>(v, y);
				}
				else if (v > y) {
					val = -metric.diagonal_entry(u);
					K = ji2idx<IndexType>(y, v);
				}
				else if (v == y) {
					val = metric.metric_factor(u, v);
				}
			}
			else {
				if (v < y) {
					val = -1 * metric.diagonal_entry(u);
					K = ji2idx<IndexType>(v, y);
				}
				else if (v > y) {
					val = +1 * metric.diagonal_entry(u);
					K = ji2idx<IndexType>(y, v);
				}
				else if (v == y) {
					val = -metric.metric_factor(u, v);
				}
			}
		}
		else if (v == y) {
			if (u < x) {
				val = -1 * metric.diagonal_entry(v);
				K = ji2idx<IndexType>(u, x);

			}
			else {
				val = +1 * metric.diagonal_entry(v);
				K = ji2idx<IndexType>(x, u);

			}
		}
		else if (u == y) {
			if (v > x) {
				val = -1 * metric.diagonal_entry(u);
				K = ji2idx<IndexType>(x, v);
			}
			else {
				val = +1 * metric.diagonal_entry(v);
				K = ji2idx<IndexType>(v, x);
			}
		}
		else if (v == x) {
			if (u > y) {
				val = -1 * metric.diagonal_entry(v);
				K = ji2idx<IndexType>(y, u);
			}
			else {
				val = +1 * metric.diagonal_entry(v);
				K = ji2idx<IndexType>(u, y);
			}
		}

		if (axis == 0) {
			return K;
		}
		else if (axis == 1) {
			return j;
		}
		else if (axis == 2) {
			return i;
		}
		else if (axis == 3) {
			return val;
		}
		return 0;
	}
};

template<typename IndexType, typename CoeffType, class MetricType, typename MemoryType> // TODO: assert MetricType is a implementation of Metric
SparseTensor<IndexType, CoeffType, cusp::device_memory> build_geometric_product_tensor(const MetricType &metric) {
	IndexType full_size = (N * (N + 1) >> 1) + 1;

	cusp::array1d<IndexType, cusp::device_memory> all_basis = cusp::counting_array<IndexType>(full_size);

	typedef typename thrust::device_vector<IndexType>::iterator Iterator;
	repeated_range<Iterator> all_basis_repeated(all_basis.begin(), all_basis.end(), full_size);
	tiled_range<Iterator>    all_basis_tiled(all_basis.begin(), all_basis.end(), full_size);

	cusp::array1d<IndexType, cusp::device_memory> K(full_size * full_size);

	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
		all_basis_tiled.begin(),
		K.begin(),
		GeometricProductTensorFunctor<IndexType, CoeffType, IndexType, MetricType>(metric, 0));

	cusp::array1d<IndexType, cusp::device_memory> J(full_size * full_size);

	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
		all_basis_tiled.begin(),
		J.begin(),
		GeometricProductTensorFunctor<IndexType, CoeffType, IndexType, MetricType>(metric, 1));

	cusp::array1d<IndexType, cusp::device_memory> I(full_size * full_size);

	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
		all_basis_tiled.begin(),
		I.begin(),
		GeometricProductTensorFunctor<IndexType, CoeffType, IndexType, MetricType>(metric, 2));


	cusp::array1d<CoeffType, cusp::device_memory> values(full_size * full_size);

	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
		all_basis_tiled.begin(),
		values.begin(),
		GeometricProductTensorFunctor<IndexType, CoeffType, CoeffType, MetricType>(metric, 3));

	std::vector<IndexType> shape = { full_size, full_size, full_size };
	//print(I);
	//print(J);
	//print(K);
	//print(values);

	// remover de I, J, K e values sempre que values == 0

	//typename cusp::array1d<IndexType, MemoryType>::iterator new_end_I = thrust::remove_if(&I[0], &I[0] + I.size(), &values[0], is_zero<CoeffType>());
	//cusp::array1d<IndexType, MemoryType> new_I(I.begin(), new_end_I);
	//typename cusp::array1d<IndexType, MemoryType>::iterator new_end_J = thrust::remove_if(&J[0], &J[0] + J.size(), &values[0], is_zero<CoeffType>());
	//cusp::array1d<IndexType, MemoryType> new_J(J.begin(), new_end_J);
	//typename cusp::array1d<IndexType, MemoryType>::iterator new_end_K = thrust::remove_if(&K[0], &K[0] + K.size(), &values[0], is_zero<CoeffType>());
	//cusp::array1d<IndexType, MemoryType> new_K(K.begin(), new_end_K);
	//typename cusp::array1d<CoeffType, MemoryType>::iterator new_end_values = thrust::remove_if(&values[0], &values[0] + values.size(), &values[0], is_zero<CoeffType>());
	//cusp::array1d<CoeffType, MemoryType> new_values(values.begin(), new_end_values);

	//for (IndexType i = 0; i < new_I.size(); i++) {
	//	std::cout << new_K[i] << " " << new_J[i] << " " << new_I[i] << " " << new_values[i] << std::endl;
	//}

	return SparseTensor<IndexType, CoeffType, cusp::device_memory>(I, J, K, values, shape);
}

template<typename IndexType, typename CoeffType, typename MemoryType>
SparseTensor<IndexType, CoeffType, MemoryType> build_convolution_tensor(IndexType N_DIMS) {
	std::vector<IndexType> K;
	std::vector<IndexType> I;
	std::vector<IndexType> J;

	SparseTensor<IndexType, CoeffType, MemoryType> convolution_tensor;

	for (IndexType k = 0; k < N_DIMS; k++) {
		for (IndexType idx = 0; idx <= k; idx++) {
			K.push_back(k);
			J.push_back(idx);
			I.push_back(k - idx);
		}
	}

	cusp::array1d<IndexType, MemoryType> r(J);
	cusp::array1d<IndexType, MemoryType> c(I);
	cusp::array1d<IndexType, MemoryType> l(K);
	cusp::array1d<CoeffType, MemoryType> d = cusp::constant_array<CoeffType>(K.size(), 1.0);

	std::vector<IndexType> shape = { N_DIMS, N_DIMS, N_DIMS };
	return SparseTensor<IndexType, CoeffType, MemoryType>(c, r, l, d, shape);
}

