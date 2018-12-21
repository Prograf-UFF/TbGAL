#ifndef MULTIVECTOR_CU
#define MULTIVECTOR_CU

#include "../../common.cu"

#include "../sparse_tensor/SparseTensor.cu"
#include <cusp/print.h>
#include <cusp/array1d.h>
#include <thrust/tuple.h>
#include "../sparse_tensor/utils.cu"
#include <vector>
#include <chrono>
#include <type_traits>
#include "../metric/metric.cu"
#include <thrust/device_vector.h>
#include <thrust/device_ptr.h>
#include <type_traits>
#include <set>

#include <boost/python.hpp>

SparseTensor<IndexType, CoeffType, MemorySpace> *GP_T;
SparseTensor<IndexType, CoeffType, MemorySpace> *OP_T;


class Multivector {
	private:
		SparseTensor<IndexType, CoeffType, MemorySpace> core;
		static IndexType N;
		static IndexType N_FULL;

	public:
		static IndexType get_N();
		static IndexType get_N_FULL();
		static void set_N(const IndexType N);
		static SparseTensor<IndexType, CoeffType, MemorySpace>* get_OP_T();
		static SparseTensor<IndexType, CoeffType, MemorySpace>* get_GP_T();

		Multivector();
		Multivector(const Multivector& m);
		Multivector(const SparseTensor<IndexType, CoeffType, MemorySpace> &core);
	    Multivector(IndexType index, CoeffType coeff = 1);

		SparseTensor<IndexType, CoeffType, MemorySpace> getCore() const;
		std::vector<int> get_grade() const;
		int get_grade_blade();
	    cusp::array1d<IndexType, MemorySpace> getComponentIndexes() const;

	    Multivector getComponent(IndexType idx);
		Multivector get_component_max_projection();

		// REVERSE
		Multivector operator ~();
		// UNARY PLUS
		Multivector operator +();
		// UNARY MINUS
		Multivector operator -();


		// SUM
		friend Multivector operator +(const Multivector &lhs, const Multivector &rhs);
		template<typename ScalarType, typename>
		friend Multivector operator +(const Multivector &lhs, const ScalarType &rhs);
		template<typename ScalarType, typename>
		friend Multivector operator +(const ScalarType &lhs, const Multivector &rhs);

		// DIFF
		friend Multivector operator -(const Multivector &lhs, const Multivector &rhs);
		template<typename ScalarType, typename>
		friend Multivector operator -(const Multivector &lhs, const ScalarType &rhs);
		template<typename ScalarType, typename>
		friend Multivector operator -(const ScalarType &lhs, const Multivector &rhs);

		// OUTER PRODUCT
		friend Multivector operator ^(const Multivector &lhs, const Multivector &rhs);

		// PRODUCT
		template<typename ScalarType, typename>
		friend Multivector operator *(const Multivector &lhs, const ScalarType &rhs);
		template<typename ScalarType, typename>
		friend Multivector operator *(const ScalarType &lhs, const Multivector &rhs);

		// OPERATOR ==
		friend bool operator ==(const Multivector &lhs, const Multivector &rhs);

		// OPERATOR <<
		friend std::ostream& operator <<(std::ostream& os, Multivector& m);


		// Multivector* REVERSE();
		// Multivector* INVOLUTION();
	    // Multivector* take_grade(IndexType grade);
	    std::string to_string();

};

IndexType Multivector::N = 0;
IndexType Multivector::N_FULL = 0;


/************************** END OF MULTIVECTOR.H ***************************/

/************************** GENERAL OPERATIONS.H ***************************/

template<typename MetricType>
SparseTensor<IndexType, CoeffType, MemorySpace> *build_tensor(MetricType metric);

template<typename MetricType>
void generate_T(const MetricType metric);

template<typename MetricType>
SparseTensor<IndexType, CoeffType, MemorySpace> *build_tensor(MetricType metric) {
  return new SparseTensor<IndexType, CoeffType, cusp::device_memory>(build_geometric_product_tensor<MetricType>(Multivector::get_N(), metric));
}

void setDevice(int device){
  cudaSetDevice(device);
}

template<typename MetricType>
void generate_T(const MetricType metric) {
    GP_T = new SparseTensor<IndexType, CoeffType, MemorySpace>(build_geometric_product_tensor<MetricType>(Multivector::get_N(), metric));
}


/*************** END OF GENERAL OPERATIONS.H ****************/


/******************** MULTIVECTOR.CPP **********************/

IndexType Multivector::get_N() {
	return N;
}

IndexType Multivector::get_N_FULL() {
	return N_FULL;
}

void Multivector::set_N(const IndexType N) {
	Multivector::N = N;
	Multivector::N_FULL = (N * (N+1) >> 1) + 1;
}

Multivector::Multivector() {
	this->core = SparseTensor<IndexType, CoeffType, MemorySpace>();
}

Multivector::Multivector(const Multivector& m) {
	core = m.getCore();
}

Multivector::Multivector(const SparseTensor<IndexType, CoeffType, MemorySpace> &core) {
	this->core = core;
}

Multivector::Multivector(IndexType index, CoeffType coeff) {
	cusp::array1d<IndexType, MemorySpace> indices(1);
	cusp::array1d<CoeffType, MemorySpace> data(1);
	std::vector<IndexType> shape = { 1, Multivector::get_N_FULL()};

	data[0] = coeff;
	indices[0] = index;

	this->core = SparseTensor<IndexType, CoeffType, MemorySpace>(indices, data, shape, coeff == 0);

}

SparseTensor<IndexType, CoeffType, MemorySpace> Multivector::getCore() const {
	return this->core;
}

std::vector<int> Multivector::get_grade() const {
	std::vector<int> ret;
	if (core.getRank() == 1) {
		cusp::array1d<IndexType, MemorySpace> indices = core.getCore().column_indices;
		for (int i = 0; i <= 2; i++) {
			int count_grade = thrust::count_if(indices.begin(), indices.end(), is_grade<IndexType>(Multivector::get_N(), i));
			if (count_grade != 0) {
				ret.push_back(i);
			}
		}
	}
	// TODO else for handling exception
	return ret;
}

int Multivector::get_grade_blade() {
	return this->get_grade()[0];
}

cusp::array1d<IndexType, MemorySpace> Multivector::getComponentIndexes() const {
	return this->getCore().getCore().column_indices;
}

Multivector Multivector::getComponent(IndexType idx) {
	cusp::array1d<IndexType, MemorySpace> indices = this->getCore().getCore().column_indices;
	cusp::array1d<CoeffType, MemorySpace> values = this->getCore().getCore().values;

	cusp::array1d<IndexType, MemorySpace> new_indices(1);
	cusp::array1d<CoeffType, MemorySpace> new_values(1);

	thrust::copy_if(indices.begin(), indices.end(), &indices[0], &new_indices[0], is_component<IndexType>(idx));
	thrust::copy_if(values.begin(), values.end(), &indices[0], &new_values[0], is_component<IndexType>(idx));

	std::vector<IndexType> dense_shape = this->getCore().getDenseShape();

	SparseTensor<IndexType, CoeffType, MemorySpace> new_core(new_indices, new_values, dense_shape, false);
	return Multivector(new_core);
}

Multivector Multivector::get_component_max_projection() {
	cusp::array1d<IndexType, MemorySpace> indices = this->getComponentIndexes();
	cusp::array1d<CoeffType, MemorySpace> values = this->getCore().getCore().values;

	if (indices.size() == 0) {
		// throw
		// return NULL;
	}

	thrust::device_vector<CoeffType>::iterator iter = thrust::max_element(values.begin(), values.end());
	auto coeff_max = *iter;
	auto idx_max = iter - values.begin();

	iter = thrust::min_element(values.begin(), values.end());
	auto coeff_min = *iter;
	auto idx_min = iter - values.begin();

	if (abs(coeff_min) > coeff_max) {
		coeff_max = coeff_min;
		idx_max = idx_min;
	}

	// TODO fixop
	// return new Multivector(indices[idx_max]);
	return Multivector(indices[idx_max], coeff_max);
}

std::string Multivector::to_string() {
	std::string repr = "";
	bool first = true;
	if (core.getCore().values.size() == 0) {
		repr += "0";
	} else {
		idx2ji<IndexType> func(Multivector::get_N());
		for (IndexType i = 0; i < core.getCore().values.size(); i++) {
			if (core.getCore().column_indices[i] == 0) {
				repr += (core.getCore().values[i] > 0.0 ? (first ? "" : "+" ) : "") + std::to_string(core.getCore().values[i]);
				first = false;
			} else if (core.getCore().column_indices[i] <= Multivector::get_N()){
				repr += (core.getCore().values[i] > 0.0 ? (first ? "" : "+" ) : "") + (core.getCore().values[i] == 1.0 ? "" : std::to_string(core.getCore().values[i]) + "*") + "e(" + std::to_string(core.getCore().column_indices[i]) + ")";
				first = false;
			} else {
				thrust::tuple<IndexType, IndexType> basis = func(core.getCore().column_indices[i]);
				repr += (core.getCore().values[i] > 0.0 ? (first ? "" : "+" ) : "") + (core.getCore().values[i] == 1.0 ? "" : std::to_string(core.getCore().values[i]) + "*") + "e(" + std::to_string(thrust::get<0>(basis)) + ")^e(" + std::to_string(thrust::get<1>(basis)) + ")";
				first = false;
			}
		}
	}
	return repr;
}

Multivector e(IndexType index) {
  if (index > Multivector::get_N()) {
    throw std::logic_error("Can't allocate basis blade e(" + std::to_string(index) + ") in a " + std::to_string(Multivector::get_N()) + "-dimensional space");
  }
	return Multivector(index);
}

SparseTensor<IndexType, CoeffType, MemorySpace>* Multivector::get_GP_T() {
    return GP_T;
}

SparseTensor<IndexType, CoeffType, MemorySpace>* Multivector::get_OP_T() {;
    return OP_T;
}


#include "operations.cu"

// REVERSE
Multivector Multivector::operator ~() {
	return MultivectorOperations::REVERSE(*this);
}
// UNARY PLUS
Multivector Multivector::operator +() {
	return MultivectorOperations::UNARY_PLUS(*this);
}
// UNARY MINUS
Multivector Multivector::operator -() {
	return MultivectorOperations::UNARY_MINUS(*this);
}
// SUM
Multivector operator +(const Multivector &lhs, const Multivector &rhs) {
	return MultivectorOperations::ADD(lhs, rhs);
}
template<typename ScalarType>
Multivector operator +(const Multivector &lhs, const ScalarType &rhs) {
	return MultivectorOperations::ADD_SCALAR(lhs, rhs);
}
template<typename ScalarType>
Multivector operator +(const ScalarType &lhs, const Multivector &rhs) {
	return MultivectorOperations::ADD_SCALAR(rhs, lhs);
}
// DIFF
Multivector operator -(const Multivector &lhs, const Multivector &rhs) {
	return MultivectorOperations::SUB(lhs, rhs);
}
template<typename ScalarType>
Multivector operator -(const Multivector &lhs, const ScalarType &rhs) {
	return MultivectorOperations::SUB_SCALAR(lhs, rhs);
}
template<typename ScalarType>
Multivector operator -(const ScalarType &lhs, const Multivector &rhs) {
	return MultivectorOperations::SUB_SCALAR(rhs, lhs);
}
// OUTER PRODUCT
Multivector operator ^(const Multivector &lhs, const Multivector &rhs) {
	return MultivectorOperations::OP(lhs, rhs);
}
// PRODUCT
template<typename ScalarType>
Multivector operator *(const Multivector &lhs, const ScalarType &rhs) {
	return MultivectorOperations::PROD(lhs, rhs);
}
template<typename ScalarType>
Multivector operator *(const ScalarType &lhs, const Multivector &rhs) {
	return MultivectorOperations::PROD(rhs, lhs);
}
// OPERATOR ==
bool operator ==(const Multivector &lhs, const Multivector &rhs) {
	return MultivectorOperations::is_equals(lhs, rhs);
}
// OPERATOR <<
std::ostream& operator <<(std::ostream& os, Multivector& m) {
	os << m.to_string();
	return os;
}

#endif
