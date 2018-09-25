#pragma once
#include <SparseTensor.h>
#include <cusp\print.h>
#include <cusp\array1d.h>
#include <thrust\tuple.h>
#include "utils.cuh"
#include <vector>
#include <chrono>
#include <type_traits>
#include "metric.cuh"
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>

#ifndef N
#error "Please define N before include multivector header"
#endif

#define N_FULL ((N * (N+1)) >> 1) + 1

typedef uint64_t IndexType ;
typedef float CoeffType;
typedef cusp::device_memory MemorySpace;

class Multivector {

	public:
		Multivector() {
		}

		Multivector(const SparseTensor<IndexType, CoeffType, MemorySpace> &core) {
			this->core = core;
		}


		Multivector(IndexType index) {
			cusp::array1d<IndexType, MemorySpace> indices(1);
			cusp::array1d<CoeffType, MemorySpace> data(1);
			std::vector<IndexType> shape = { 1, N_FULL };

			data[0] = 1;
			indices[0] = index;

			this->core = SparseTensor<IndexType, CoeffType, MemorySpace>(indices, data, shape);

		}

		SparseTensor<IndexType, CoeffType, MemorySpace> getCore() const {
			return this->core;
		}

		template<typename IndexType, typename CoeffType, typename MemorySpace> friend std::ostream& operator << (std::ostream &, const Multivector<IndexType, CoeffType, MemorySpace> &);
		template<typename IndexType, typename CoeffType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> operator + (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs);
		template<typename IndexType, typename CoeffType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> operator - (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs);
		template<typename IndexType, typename CoeffType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> operator ^ (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs);
		template<typename IndexType, typename CoeffType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> take_grade (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, short grade);
		template<typename IndexType, typename CoeffType, typename MetricType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> GP (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs, const MetricType &metric);

		template<typename IndexType, typename CoeffType, typename MetricType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> SCP(const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs, const MetricType &metric);

		template<typename IndexType, typename CoeffType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> operator * (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const CoeffType &rhs);
		template<typename IndexType, typename CoeffType, typename MemorySpace> friend Multivector<IndexType, CoeffType, MemorySpace> operator * (const CoeffType &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs);


private:
		SparseTensor<IndexType, CoeffType, MemorySpace> core;
};

Multivector<IndexType, CoeffType, MemorySpace> operator * (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const CoeffType &rhs) {
	return Multivector<IndexType, CoeffType, MemorySpace>(rhs * lhs.getCore());
}

Multivector<IndexType, CoeffType, MemorySpace> operator * (const CoeffType &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs) {
	SparseTensor<IndexType, CoeffType, MemorySpace> core = lhs * rhs.getCore();
	return Multivector<IndexType, CoeffType, MemorySpace>(core);
}

std::ostream& operator << (std::ostream &os, const Multivector<IndexType, CoeffType, MemorySpace> &m) {
	if (m.core.getCore().values.size() == 0) {
		os << "0";
	}
	else {
		idx2ji<IndexType> func;
		for (uint64_t i = 0; i < m.core.getCore().values.size(); i++) {
			if (m.core.getCore().column_indices[i] == 0) {
				os << (m.core.getCore().values[i] > 0.0 ? "+" : "") << m.core.getCore().values[i];
			}
			else if (m.core.getCore().column_indices[i] <= N){
				os << (m.core.getCore().values[i] > 0.0 ? "+" : "") << m.core.getCore().values[i] << "e(" << m.core.getCore().column_indices[i] << ") ";
			}
			else {
				thrust::tuple<IndexType, IndexType> basis = func(m.core.getCore().column_indices[i]);
				os << (m.core.getCore().values[i] > 0.0 ? "+" : "") << m.core.getCore().values[i] << "e(" << thrust::get<0>(basis) << ")^e(" << thrust::get<1>(basis) << ")";
			}
		}
	}
	return os;
}

Multivector<IndexType, CoeffType, MemorySpace> operator + (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs) {
	Multivector<IndexType, CoeffType, MemorySpace> ret;
	ret.core = lhs.getCore() + rhs.getCore();
	return ret;
}

Multivector<IndexType, CoeffType, MemorySpace> operator - (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs) {
	Multivector<IndexType, CoeffType, MemorySpace> ret;
	ret.core = lhs.getCore() - rhs.getCore();
	return ret;
}


template<typename IndexType, typename CoeffType, typename MetricType, typename MemorySpace>
Multivector<IndexType, CoeffType, MemorySpace> GP(const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs, const MetricType &metric) {
	auto tensor = build_geometric_product_tensor<IndexType, CoeffType, MetricType, MemorySpace>(metric);
	Multivector<IndexType, CoeffType, MemorySpace> ret(SparseTensor<IndexType, CoeffType, MemorySpace>::tensor_dot(lhs.getCore(), rhs.getCore().t(), tensor).t());
	//std::cout << lhs.getCore().getRank() << std::endl;
	//std::cout << rhs.getCore().getRank() << std::endl;
	//std::cout << tensor.getRank() << std::endl;
	//std::cout << ret.getCore().getRank() << std::endl;
	return ret;
}

Multivector<IndexType, CoeffType, MemorySpace> operator ^ (const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs) {
	Multivector<IndexType, CoeffType, MemorySpace> ret = GP(lhs, rhs, EuclideanMetric<IndexType, CoeffType>());
	// TODO: Take_Grade
	return take_grade(ret, 2);
}

template<typename IndexType, typename CoeffType, typename MetricType, typename MemorySpace>
Multivector<IndexType, CoeffType, MemorySpace> SCP(const Multivector<IndexType, CoeffType, MemorySpace> &lhs, const Multivector<IndexType, CoeffType, MemorySpace> &rhs, const MetricType &metric) {
	Multivector<IndexType, CoeffType, MemorySpace> ret = GP(lhs, rhs, metric);
	return take_grade(ret, 0);
}


Multivector<IndexType, CoeffType, MemorySpace> take_grade(const Multivector<IndexType, CoeffType, MemorySpace> &m, short grade) {

	if (m.getCore().getRank() == 1) {
		//thrust::device_vector<IndexType> indices(m.getCore().getCore().column_indices.size());
		//thrust::device_vector<CoeffType> values(m.getCore().getCore().column_indices.size());

		//cusp::array1d<IndexType, MemorySpace> idx = m.getCore().getCore().column_indices;
		//cusp::array1d<CoeffType, MemorySpace> val = m.getCore().getCore().values;

		//thrust::device_vector<IndexType> thrust_column_indices(idx.begin(), idx.end());
		//thrust::device_vector<CoeffType> thrust_values(val.begin(), val.end());

		//thrust::device_ptr<IndexType> ptr_column_indices = &thrust_column_indices[0];//m.getCore().getCore().column_indices.data();
		//thrust::device_ptr<CoeffType> ptr_values = &thrust_values[0];//m.getCore().getCore().values.data();

		//thrust::device_ptr<IndexType> ptr_new_column_indices = &indices[0];
		//thrust::device_ptr<CoeffType> ptr_new_values = &values[0];

		//thrust::copy_if(ptr_column_indices, ptr_column_indices + m.getCore().getCore().column_indices.size(), ptr_column_indices, ptr_new_column_indices, is_grade<IndexType>(grade));
		//thrust::copy_if(ptr_values, ptr_values + m.getCore().getCore().values.size(), ptr_column_indices, ptr_new_values, is_grade<CoeffType>(grade));


		cusp::array1d<IndexType, MemorySpace> indices = m.getCore().getCore().column_indices;
		cusp::array1d<CoeffType, MemorySpace> values = m.getCore().getCore().values;

		cusp::array1d<IndexType, MemorySpace> new_indices(indices.size());
		cusp::array1d<CoeffType, MemorySpace> new_values(values.size());


		thrust::copy_if(indices.begin(), indices.end(), &indices[0], &new_indices[0], is_grade<IndexType>(grade));
		thrust::copy_if(values.begin(), values.end(), &indices[0], &new_values[0], is_grade<IndexType>(grade));
		//		thrust::copy_if(ptr_values, ptr_values + m.getCore().getCore().values.size(), ptr_column_indices, ptr_new_values, is_grade<CoeffType>(grade));


		std::vector<IndexType> dense_shape = m.getCore().getDenseShape();
//		cusp::array1d<IndexType, MemorySpace> new_indices = indices;// (ptr_new_column_indices, ptr_new_column_indices + m.getCore().getCore().column_indices.size());
//		cusp::array1d<CoeffType, MemorySpace> new_values = values;// (ptr_new_values, ptr_new_values + m.getCore().getCore().values.size());
		SparseTensor<IndexType, CoeffType, MemorySpace> new_core(new_indices, new_values, dense_shape);
		return Multivector<IndexType, CoeffType, MemorySpace>(new_core);


//		return m;

	}
}


//template<typename IndexType, typename CoeffType, typename MemorySpace>
Multivector<IndexType, CoeffType, MemorySpace> e(uint64_t index) {
	Multivector<IndexType, CoeffType, MemorySpace> m(index);
	return m;
}



//
//template<class MetricType> // TODO: assert MetricType is a implementation of Metric
//SparseTensor<uint64_t, float, cusp::device_memory> build_geometric_product_tensor(const uint64_t &N, const MetricType &metric) {
//	uint64_t full_size = (N * (N + 1) >> 1) + 1;
//
//	cusp::array1d<uint64_t, cusp::device_memory> all_basis = cusp::counting_array<uint64_t>(full_size);
//
//	typedef thrust::device_vector<uint64_t>::iterator Iterator;
//	repeated_range<Iterator> all_basis_repeated(all_basis.begin(), all_basis.end(), full_size);
//	tiled_range<Iterator>    all_basis_tiled(all_basis.begin(), all_basis.end(), full_size);
//
//	cusp::array1d<uint64_t, cusp::device_memory> K(full_size * full_size);
//
//	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
//		all_basis_tiled.begin(),
//		K.begin(),
//		EuclideanGeometricProductTensorA<uint64_t, float, uint64_t, MetricType>(N, metric, 0));
//
//	cusp::array1d<uint64_t, cusp::device_memory> J(full_size * full_size);
//
//	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
//		all_basis_tiled.begin(),
//		J.begin(),
//		EuclideanGeometricProductTensorA<uint64_t, float, uint64_t, MetricType>(N, metric, 1));
//
//	cusp::array1d<uint64_t, cusp::device_memory> I(full_size * full_size);
//
//	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
//		all_basis_tiled.begin(),
//		I.begin(),
//		EuclideanGeometricProductTensorA<uint64_t, float, uint64_t, MetricType>(N, metric, 2));
//
//
//	cusp::array1d<float, cusp::device_memory> values(full_size * full_size);
//
//	thrust::transform(all_basis_repeated.begin(), all_basis_repeated.end(),
//		all_basis_tiled.begin(),
//		values.begin(),
//		EuclideanGeometricProductTensorA<uint64_t, float, uint64_t, MetricType>(N, metric, 3));
//
//	std::vector<uint64_t> shape = { full_size, full_size, full_size };
//	return SparseTensor<uint64_t, float, cusp::device_memory>(I, J, K, values, shape);
//}
//
//template<typename IndexType, typename CoeffType, typename MemoryType>
//SparseTensor<IndexType, CoeffType, MemoryType> build_convolution_tensor(IndexType N_DIMS) {
//	std::vector<IndexType> K;
//	std::vector<IndexType> I;
//	std::vector<IndexType> J;
//
//	SparseTensor<IndexType, CoeffType, MemoryType> convolution_tensor;
//
//	for (IndexType k = 0; k < N_DIMS; k++) {
//		for (IndexType idx = 0; idx <= k; idx++) {
//			K.push_back(k);
//			J.push_back(idx);
//			I.push_back(k - idx);
//		}
//	}
//
//	cusp::array1d<IndexType, MemoryType> r(J);
//	cusp::array1d<IndexType, MemoryType> c(I);
//	cusp::array1d<IndexType, MemoryType> l(K);
//	cusp::array1d<CoeffType, MemoryType> d = cusp::constant_array<CoeffType>(K.size(), 1.0);
//
//	std::vector<uint64_t> shape = { N_DIMS, N_DIMS, N_DIMS };
//	return SparseTensor<IndexType, CoeffType, MemoryType>(c, r, l, d, shape);
//}
