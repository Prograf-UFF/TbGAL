#ifndef SPARSE_TENSOR_CU
#define SPARSE_TENSOR_CU

#include <cusp/coo_matrix.h>
#include <cusp/csr_matrix.h>
#include <cusp/print.h>
#include <cusp/multiply.h>
#include <iostream>
#include <vector>
#include <cusp/functional.h>
#include <cusp/print.h>
#include <cusp/array1d.h>
#include <thrust/tuple.h>
#include <thrust/find.h>
#include <thrust/device_vector.h>
#include <thrust/execution_policy.h>
#include <cusp/transpose.h>
#include <functional>
#include <numeric>
#include <cusp/elementwise.h>
#include <cusp/array2d.h>
#include "util_functors.cu"
#include <fstream>

template <typename T>
struct is_zero {
	__host__ __device__ bool operator()(const T &i) {
		return abs(i) <= 1e-4;
	}
};

//typedef long long unsigned int IndexType;

template<typename IndexType, typename DataType, class MemorySpace>
class SparseTensor {

private:
	typedef cusp::array1d<DataType, MemorySpace> dataArray_t;
	typedef cusp::array1d_view<typename dataArray_t::iterator> dataArray_view_t;
	typedef cusp::array1d<IndexType, MemorySpace> indexArray_t;
	typedef cusp::array1d_view<typename indexArray_t::iterator> indexArray_view_t;
	typedef cusp::coo_matrix<IndexType, DataType, MemorySpace> spmcore_t;


	short rank;
	bool initialized;
	std::vector<IndexType> dense_shape;

public:
	indexArray_t I;
	indexArray_t J;
	indexArray_t K;
	dataArray_t data;
	spmcore_t core;
	/* DEFAULT CONSTRUCTOR */
	SparseTensor() {
		this->initialized = false;
	}

	SparseTensor(const SparseTensor& s) {
		rank = s.getRank();
		core = s.getCore();
		dense_shape = s.getDenseShape();
		initialized = s.isInitialized();
	}
	// TODO change to copy_if
	template<typename T1, typename T2>
	T1 erase_zeros(T1 indices, T2 data, bool erase=true) {
		if (!erase) {
			return T1(indices.begin(), indices.end());
		}
		typename T1::iterator new_end_indices = thrust::remove_if(indices.begin(), indices.end(), &data[0], is_zero<DataType>());
		return T1(indices.begin(), new_end_indices);
	}

	/* RANK 1 CONSTRUCTOR */
	SparseTensor(indexArray_t &indices, dataArray_t &data, std::vector<IndexType> dense_shape, bool erase=true) {

		this->dense_shape = dense_shape;
		this->rank = 1;
		this->core.resize(dense_shape[0], dense_shape[1], dataArray_view_t(data).size());
		if (dense_shape[0] == 1) {
			this->core.column_indices = erase_zeros<indexArray_t, dataArray_t>(indices, data, erase);
			this->I = this->core.column_indices;
		}
		else if (dense_shape[1] == 1) {
			this->core.row_indices = erase_zeros<indexArray_t, dataArray_t>(indices, data, erase);
			this->J = this->core.row_indices;
		}
		else {
			throw std::invalid_argument("Wrong arguments for constructor");
		}
		this->core.values = erase_zeros<dataArray_t, dataArray_t>(data, data, erase);
		this->initialized = true;
	}

	/* RANK 2 CONSTRUCTOR */
	SparseTensor(indexArray_t &cols, indexArray_t &rows, dataArray_t &data, std::vector<IndexType> dense_shape) {

		// TODO assert values and N
		this->dense_shape = dense_shape;
		this->rank = 2;
		this->core.resize(dense_shape[0], dense_shape[1], data.size());
		this->core.column_indices = erase_zeros<indexArray_t, dataArray_t>(cols, data);
		this->core.row_indices = erase_zeros<indexArray_t, dataArray_t>(rows, data);

		this->I = this->core.column_indices;
		this->J = this->core.row_indices;

		this->core.values = erase_zeros<dataArray_t, dataArray_t>(data, data);
		this->initialized = true;
	}

	/* RANK 3 CONSTRUCTOR */
	SparseTensor(indexArray_t &cols, indexArray_t &rows, indexArray_t &layers, dataArray_t &data, std::vector<IndexType> dense_shape) {

		indexArray_t c = indexArray_t(erase_zeros<indexArray_t, dataArray_t>(cols, data));
		indexArray_t r = indexArray_t(erase_zeros<indexArray_t, dataArray_t>(rows, data));
		indexArray_t l = indexArray_t(erase_zeros<indexArray_t, dataArray_t>(layers, data));
		dataArray_t d = dataArray_t(erase_zeros<dataArray_t, dataArray_t>(data, data));

		this->I = c;
		this->J = r;
		this->K = l;
		this->data = d;

		// std::cout << "BEFORE: " << d.size() << std::endl;
//		std::ofstream myfile;
//		myfile.open("tensor.dat");
//		myfile << dense_shape[0] << std::endl;
//		for (IndexType i = 0; i < d.size(); i++) {
//			myfile << l[i] << ";" << r[i] << ";" << c[i] << ";" << d[i] << std::endl;
//		}
//		myfile.close();

		// TODO assert values and N
		this->dense_shape = dense_shape;
		this->rank = 3;
		this->core.resize(dense_shape[0], dense_shape[1] * dense_shape[2], d.size());
		this->core.row_indices = r;
		thrust::transform(l.begin(), l.end(), this->core.column_indices.begin(), cusp::multiplies_value<IndexType>(dense_shape[1]));
		thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(c.begin(), this->core.column_indices.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(c.begin(), this->core.column_indices.begin())) + this->core.column_indices.size(),
			this->core.column_indices.begin(),
			cusp::sum_pair_functor<IndexType>());
		this->core.values = d;
		// std::cout << "BEFORE: " << d.size() << std::endl;
		this->core.sort_by_row_and_column();
		this->initialized = true;

	}

	/* CORE CONSTRUCTOR */
	SparseTensor(const std::vector<IndexType>& dense_shape, short rank, const spmcore_t &core) {
		spmcore_t new_core(core);
		new_core.column_indices = erase_zeros<indexArray_t, dataArray_t>(core.column_indices, core.values);
		new_core.row_indices = erase_zeros<indexArray_t, dataArray_t>(core.row_indices, core.values);
		new_core.values = erase_zeros<dataArray_t, dataArray_t>(core.values, core.values);

		this->dense_shape = dense_shape;
		this->rank = rank;
		this->core = new_core;
		this->initialized = true;
	}

	/* DEFAULT DESTRUCTOR */
	~SparseTensor() {

	}

	template<typename IndexT, typename DataT, typename MemoryS>
	friend std::ostream& operator << (std::ostream &, const SparseTensor<IndexT, DataT, MemoryS> &);

	template<typename IndexT, typename DataT, typename MemoryS>
	friend SparseTensor<IndexT, DataT, MemoryS> operator* (const DataT &, const SparseTensor<IndexT, DataT, MemoryS> &);
//	template<typename IndexType, typename DataType, typename MemorySpace> friend SparseTensor<IndexType, DataType, MemorySpace> operator* (const SparseTensor<IndexType, DataType, MemorySpace> &, const DataType &);

	void print_shape() {
		std::cout << "( ";
		for (short i = 0; i < this->dense_shape.size(); i++) {
			std::cout << this->dense_shape[i] << " ";
		}
		std::cout << ")" << std::endl;
	}

	bool isInitialized() const {
		return this->initialized;
	}

	const std::vector<IndexType>& getDenseShape() const {
		return this->dense_shape;
	}

	const spmcore_t& getCore() const {
		return this->core;
	}

	const short& getRank() const {
		return this->rank;
	}

	DataType operator() (IndexType i) {
		if (this->getRank() == 1) {
			indexArray_t *array;
			if (this->dense_shape[0] == 1) {
				array = &core.column_indices;
			}
			else {
				array = &core.row_indices;
			}
			auto iter = thrust::find(thrust::device, array->begin(), array->end(), i);
			if (iter == array->end()) {
				return 0;
			}
			return core.values[iter - array->begin()];
		}
		else {
			throw std::logic_error("Not implemented for rank > 1 tensors");
		}
	}

	SparseTensor operator+ (const SparseTensor &lhs) const {
		spmcore_t new_core;

		cusp::add(this->getCore(), lhs.getCore(), new_core);
		return SparseTensor(this->getDenseShape(), this->getRank(), new_core);
	}
	SparseTensor operator- (const SparseTensor &lhs) const {

		spmcore_t new_core;
		cusp::subtract(this->getCore(), lhs.getCore(), new_core);
		return SparseTensor(this->getDenseShape(), this->getRank(), new_core);
	}

	bool operator==(const SparseTensor &t) {
		return (this->getCore().values.size() == t.getCore().values.size()
		&& thrust::equal(this->getCore().values.begin(), this->getCore().values.end(), t.getCore().values.begin(), compare_given_threshold<CoeffType>(epsilon))
		&& thrust::equal(this->getCore().column_indices.begin(), this->getCore().column_indices.end(), t.getCore().column_indices.begin())
		&& thrust::equal(this->getCore().row_indices.begin(), this->getCore().row_indices.end(), t.getCore().row_indices.begin()));
	}

	SparseTensor operator* (const SparseTensor &rhs) const {

		IndexType inner_dim_lhs = this->getDenseShape()[this->getDenseShape().size() - 1];
		IndexType inner_dim_rhs = rhs.getDenseShape()[0];

		if (inner_dim_lhs == inner_dim_rhs) {

			std::vector<IndexType> output_shape;
			output_shape.reserve(this->getDenseShape().size() + rhs.getDenseShape().size() - 2);
			output_shape.insert(output_shape.end(), this->getDenseShape().begin(), this->getDenseShape().begin() + this->getDenseShape().size() - 1);
			output_shape.insert(output_shape.end(), rhs.getDenseShape().begin() + 1, rhs.getDenseShape().end());

			short output_rank = output_shape.size();

			cusp::coo_matrix<IndexType, DataType, MemorySpace> output_core;

			cusp::coo_matrix<IndexType, DataType, MemorySpace> lhs_core = this->getCore();
//			cusp::sort_by_row(lhs_core.row_indices, lhs_core.column_indices, lhs_core.values);

			cusp::coo_matrix<IndexType, DataType, MemorySpace> rhs_core = rhs.getCore();
//			cusp::sort_by_row(rhs_core.row_indices, rhs_core.column_indices, rhs_core.values);
			// rhs_core.sort_by_row_and_column();
			// lhs_core.sort_by_row_and_column();

			cusp::multiply(lhs_core, rhs_core, output_core);


//			if (output_rank == 3) {
//				if (output_core.row_indices.size() > output_core.column_indices.size()) {
//					cusp::transpose(output_core, output_core);
//				}
//			}
			return SparseTensor(output_shape, output_rank - std::count(output_shape.begin(), output_shape.end(), 1), output_core);
		}
		else {
			throw std::invalid_argument("Incompatible shapes");
		}
	}

	/* Reshape */
	SparseTensor<IndexType, DataType, MemorySpace> reshape(std::vector<IndexType> new_shape) {
		IndexType current_n_elements = std::accumulate(std::begin(this->dense_shape), std::end(this->dense_shape), 1, std::multiplies<IndexType>());
		IndexType new_n_elements = std::accumulate(std::begin(new_shape), std::end(new_shape), 1, std::multiplies<IndexType>());

		spmcore_t core = this->getCore();

		IndexType w_old = this->getRank() == 3 ? this->getDenseShape()[1] + this->getDenseShape()[2] : this->getDenseShape()[1];
		IndexType w_new = new_shape.size() == 3 ? new_shape[1] + new_shape[2] : new_shape[1];

		cusp::array1d<IndexType, MemorySpace> idx(current_n_elements);
		cusp::array1d<IndexType, MemorySpace> row_indices(current_n_elements);
		cusp::array1d<IndexType, MemorySpace> column_indices(current_n_elements);
		cusp::array1d<DataType, MemorySpace> data = this->getCore().values;

		if (current_n_elements == new_n_elements) {

			thrust::transform(core.row_indices.begin(), core.row_indices.end(), idx.begin(), cusp::multiplies_value<IndexType>(w_old));
			thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(idx.begin(), core.column_indices.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(idx.begin(), core.column_indices.begin())) + core.column_indices.size(),
				idx.begin(),
				cusp::sum_pair_functor<IndexType>());

			thrust::transform(idx.begin(), idx.end(), row_indices.begin(), cusp::divide_value<IndexType>(w_new));
			thrust::device_vector<IndexType> minus_jw(current_n_elements);
			thrust::transform(row_indices.begin(), row_indices.end(), minus_jw.begin(), cusp::multiplies_value<IndexType>(-w_new));
			thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(idx.begin(), minus_jw.begin())),
				thrust::make_zip_iterator(thrust::make_tuple(idx.begin(), minus_jw.begin())) + idx.size(),
				column_indices.begin(),
				cusp::sum_pair_functor<IndexType>());
		}
		else {
			throw std::invalid_argument("Incompatible shapes");
		}
		return SparseTensor<IndexType, DataType, MemorySpace>(column_indices, row_indices, data, new_shape);
	}

	/* Transpose */
	SparseTensor t() {
		if (this->getRank() < 3) {
			spmcore_t coreT;
			std::vector<IndexType> denseShapeT(this->getDenseShape().size());
			std::reverse_copy(this->dense_shape.begin(), this->dense_shape.end(), denseShapeT.begin());
			cusp::transpose(this->core, coreT);
			SparseTensor<IndexType, DataType, MemorySpace> tensorT(denseShapeT, this->getRank(), coreT);
			return tensorT;
		}
		throw std::logic_error("Not implemented");
	}

	static SparseTensor<IndexType, DataType, MemorySpace> tensor_dot(const SparseTensor &op1, const SparseTensor &op2, const SparseTensor &T) {


		IndexType inner_dim_op1 = op1.getDenseShape()[op1.getDenseShape().size() - 1];
		IndexType inner_dim_op1_T = T.getDenseShape()[0];
		IndexType inner_dim_T_op2 = T.getDenseShape()[T.getDenseShape().size() - 1];
		IndexType inner_dim_op2 = op2.getDenseShape()[0];

		// if (inner_dim_op1 == inner_dim_op1_T && inner_dim_T_op2 == inner_dim_op2) {
		// 	SparseTensor<IndexType, DataType, MemorySpace> Z = op1 * T;
		// 	std::vector<IndexType> new_shape = { op2.getDenseShape()[0], op2.getDenseShape()[0] };
		// 	SparseTensor<IndexType, DataType, MemorySpace> K = Z.reshape(new_shape);
		// 	auto a =  K * op2;
		// 	return a;
		// }


		if (inner_dim_op1 == inner_dim_op1_T && inner_dim_T_op2 == inner_dim_op2) {
			SparseTensor<IndexType, DataType, MemorySpace> Z = op1 * T;
			auto w = inner_dim_op1_T;
			cusp::array1d<IndexType, MemorySpace> c = Z.getCore().column_indices;
			cusp::array1d<IndexType, MemorySpace> I(c.size());
			cusp::array1d<IndexType, MemorySpace> J(c.size());
			cusp::array1d<DataType, MemorySpace> data(Z.getCore().values);

			thrust::transform(c.begin(), c.end(), I.begin(), cusp::modulus_value<IndexType>(w));
			thrust::transform(c.begin(), c.end(), J.begin(), cusp::divide_value<IndexType>(w));

			std::vector<IndexType> dims = {w, w};
			SparseTensor<IndexType, DataType, MemorySpace> Z_reshaped(I, J, data, dims);

			auto a =  Z_reshaped * op2;

			return a;
		}

		else {
			throw std::invalid_argument("Incompatible shapes");
		}
	}
};

template<typename IndexType, typename DataType, typename MemorySpace>
SparseTensor<IndexType, DataType, MemorySpace> operator* (const DataType &lhs, const SparseTensor<IndexType, DataType, MemorySpace> &rhs) {
	cusp::coo_matrix<IndexType, DataType, MemorySpace> new_core = rhs.getCore();
	thrust::transform(new_core.values.begin(), new_core.values.end(), new_core.values.begin(), cusp::multiplies_value<DataType>(lhs));
//	cusp::multiply(lhs, rhs.getCore(), new_core);
	return SparseTensor<IndexType, DataType, MemorySpace>(rhs.getDenseShape(), rhs.getRank(), new_core);
}

//template<typename IndexType, typename DataType, typename MemorySpace>
//SparseTensor<IndexType, DataType, MemorySpace> operator* (const SparseTensor<IndexType, DataType, MemorySpace> &lhs, const DataType &rhs) {
//	cusp::coo_matrix<IndexType, DataType, MemorySpace> new_core;
//	cusp::multiply(rhs, lhs.getCore(), new_core);
//	return SparseTensor<IndexType, DataType, MemorySpace>(lhs.getDenseShape(), lhs.getRank(), new_core);
//}


template<typename IndexType, typename DataType, typename MemorySpace>
std::ostream& operator << (std::ostream &os, const SparseTensor<IndexType, DataType, MemorySpace> &T) {

	if (T.getRank() == 3) {

		cusp::array1d<IndexType, MemorySpace> K(T.getCore().values.size());
		cusp::array1d<IndexType, MemorySpace> I(T.getCore().values.size());

		thrust::transform(T.getCore().column_indices.begin(), T.getCore().column_indices.end(), K.begin(), cusp::multiplies_value<IndexType>(1 / T.getDenseShape()[1]));

		cusp::array1d<IndexType, MemorySpace> kw(T.getCore().values.size());
		thrust::transform(K.begin(), K.end(), kw.begin(), cusp::multiplies_value<IndexType>(-T.getDenseShape()[1]));

		thrust::transform(thrust::make_zip_iterator(thrust::make_tuple(T.getCore().column_indices.begin(), kw.begin())),
			thrust::make_zip_iterator(thrust::make_tuple(T.getCore().column_indices.begin(), kw.begin())) + kw.size(),
			I.begin(),
			cusp::sum_pair_functor<IndexType>());

		os << "K" << '\t' << "J" << '\t' << "I" << '\t' << "VALUE" << std::endl;
		os << "-----------------------------" << std::endl;
		for (IndexType i = 0; i < K.size(); i++) {
			os << K[i] << '\t' << T.getCore().row_indices[i] << '\t' << I[i] << '\t' << T.getCore().values[i] << std::endl;
		}
	}
	return os;
}
#endif
