#ifndef OPERATIONS_CU
#define OPERATIONS_CU

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
#include "multivector.cu"

#include <type_traits>
#include <set>

#include <boost/python.hpp>

namespace MultivectorOperations {

    Multivector* GP_tensor(Multivector *lhs, Multivector *rhs, SparseTensor<IndexType, CoeffType, MemorySpace> *tensor);
    Multivector* GP(Multivector *lhs, Multivector *rhs);
    Multivector* ADD(Multivector *lhs, Multivector *rhs);
    Multivector* SUB(Multivector *lhs, Multivector *rhs);
    Multivector* SCP(Multivector *lhs, Multivector *rhs);
    Multivector* SCP_tensor(Multivector *lhs, Multivector *rhs, SparseTensor<IndexType, CoeffType, MemorySpace> *tensor);
    Multivector* OP(Multivector *lhs, Multivector *rhs);
    Multivector* LCONT(Multivector *lhs, Multivector *rhs);
    Multivector* RCONT(Multivector *lhs, Multivector *rhs);
    Multivector *take_grade(Multivector*m, IndexType grade);
    Multivector* UNARY_MINUS(Multivector *m);

    template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
    Multivector* PROD (Multivector *lhs, const ScalarType &rhs);
    template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
    Multivector* ADD_SCALAR(Multivector *lhs, const ScalarType &rhs);
    template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
    Multivector* SUB_SCALAR(Multivector *lhs, const ScalarType &rhs);
    template<typename ScalarType, typename = typename std::enable_if<std::is_arithmetic<ScalarType>::value, ScalarType>::type>
    Multivector* R_SUB_SCALAR(Multivector *lhs, const ScalarType &rhs);


    Multivector* ADD(Multivector *lhs, Multivector *rhs) {
        return new Multivector(lhs->getCore() + rhs->getCore());
    }

    template<typename ScalarType, typename>
    Multivector* ADD_SCALAR(Multivector *lhs, const ScalarType &rhs) {
        return ADD(lhs, new Multivector(0, static_cast<CoeffType>(rhs)));
    }

    template<typename ScalarType, typename>
    Multivector* PROD (Multivector *lhs, const ScalarType &rhs) {
        return new Multivector(static_cast<CoeffType>(rhs) * lhs->getCore());
    }

    Multivector* UNARY_PLUS(Multivector *m) {
        return new Multivector(m->getCore());
    }

    Multivector* UNARY_MINUS(Multivector *m) {
        return PROD<CoeffType>(m, -1);
    }

    Multivector* SUB(Multivector *lhs, Multivector *rhs) {
        return ADD(lhs, UNARY_MINUS(rhs));
    }

    template<typename ScalarType, typename>
    Multivector* SUB_SCALAR(Multivector *lhs, const ScalarType &rhs) {
        return ADD_SCALAR(lhs, -static_cast<CoeffType>(rhs));
    }

    template<typename ScalarType, typename>
    Multivector* R_SUB_SCALAR(Multivector *lhs, const ScalarType &rhs) {
        return ADD_SCALAR(UNARY_MINUS(lhs), static_cast<CoeffType>(rhs));
    }

    bool is_equals(Multivector *lhs, Multivector *rhs) {
        return (lhs->getCore() == rhs->getCore());
    }

    Multivector* REVERSE(Multivector *m) {
    	auto grade_2 = take_grade(m, 2);
    	auto new_mv = new Multivector(*m);
    	new_mv = ADD(new_mv, UNARY_MINUS(grade_2));
        new_mv = ADD(new_mv, UNARY_MINUS(grade_2));
        return new_mv;
    }

    Multivector* INVOLUTION(Multivector *m) {
    	auto grade_1 = take_grade(m, 1);
    	auto new_mv = new Multivector(*m);
    	new_mv = ADD(new_mv, UNARY_MINUS(grade_1));
        new_mv = ADD(new_mv, UNARY_MINUS(grade_1));
        return new_mv;
    }

    Multivector* GP_tensor(Multivector *lhs, Multivector *rhs, SparseTensor<IndexType, CoeffType, MemorySpace> *tensor) {
      auto rhs_t = rhs->getCore().t();
    	return new Multivector(SparseTensor<IndexType, CoeffType, MemorySpace>::tensor_dot(lhs->getCore(), rhs_t, *tensor).t());
    }

    Multivector* GP(Multivector *lhs, Multivector *rhs) {
    	return GP_tensor(lhs, rhs, Multivector::get_GP_T());
    }








    Multivector* SCP(Multivector *lhs, Multivector *rhs) {
    	Multivector *ret = GP(lhs, rhs);
    	return take_grade(ret, static_cast<IndexType>(0));
    }

    Multivector* SCP_tensor(Multivector *lhs, Multivector *rhs, SparseTensor<IndexType, CoeffType, MemorySpace> *tensor) {
    	Multivector *ret = GP_tensor(lhs, rhs, tensor);
    	return take_grade(ret, 0);
    }

    Multivector* OP(Multivector *lhs, Multivector *rhs) {
    	if (Multivector::get_OP_T() == NULL) {
    		OP_T = extract_tensor(Multivector::get_N(), Multivector::get_GP_T(), Operation::OUTER_PRODUCT); //TODO fix
    	}
    	auto r = GP_tensor(lhs, rhs, Multivector::get_OP_T());
    	return r;
    }

    Multivector* LCONT(Multivector *lhs, Multivector *rhs) {
    	auto *T = extract_tensor(Multivector::get_N(), Multivector::get_GP_T(), Operation::LEFT_CONTRACTION); //TODO fix
    	auto r = GP_tensor(lhs, rhs, T);
    	return r;
    }


    Multivector* RCONT(Multivector *lhs, Multivector *rhs) {
    	auto *T = extract_tensor(Multivector::get_N(), Multivector::get_GP_T(), Operation::RIGHT_CONTRACTION); //TODO fix
    	auto r = GP_tensor(lhs, rhs, T);
    	return r;
    }

    Multivector* dot(Multivector *lhs, Multivector *rhs) {
    	auto *T = extract_tensor(Multivector::get_N(), Multivector::get_GP_T(), Operation::DOT_PRODUCT); //TODO fix
    	auto r = GP_tensor(lhs, rhs, T);
    	return r;
    }

    CoeffType SQR_NORM(Multivector *m) {
    	auto r = GP(m, REVERSE(m));
    	return r->getCore().getCore().values[0];
    }

    Multivector *INVERSE(Multivector *m) {
    	return PROD<float>(REVERSE(m), (1.0/SQR_NORM(m)));
    }

    Multivector *CONJUGATE(Multivector *m) {
    	return REVERSE(INVOLUTION(m));
    }

    Multivector* IGP(Multivector *lhs, Multivector *rhs) {
    	auto r = GP(lhs, INVERSE(rhs));
    	return r;
    }

    CoeffType NORM(Multivector *m) {
    	return sqrt(SQR_NORM(m));
    }





    IndexType compute_max_arg_projection(Multivector *m) {
    	auto all_idx = m->getComponentIndexes();
    	auto decompose_base = idx2ji<IndexType>(Multivector::get_N());
    	std::set<IndexType> all_basis;
    	for (IndexType idx : all_idx) {
    		thrust::tuple<IndexType, IndexType> pair = decompose_base(idx);
    		if (thrust::get<0>(pair) != 0) {
    			all_basis.insert(thrust::get<0>(pair));
    		}
    		all_basis.insert(thrust::get<1>(pair));
    	}
    	CoeffType max_norm = 0;
    	IndexType max_idx = 0;
    	for (IndexType i : all_basis) {
    		auto sqr_norm = SQR_NORM(LCONT(e(i), INVERSE(m)));
    		if (sqr_norm > max_norm) {
    			max_norm = sqr_norm;
    			max_idx = i;
    		}
    	}
    	return max_idx;
    }

    Multivector *take_grade(Multivector*m, IndexType grade) {
    	if (m->getCore().getRank() == 1) {
    	    cusp::array1d<IndexType, MemorySpace> indices = m->getCore().getCore().column_indices;
    	    cusp::array1d<CoeffType, MemorySpace> values = m->getCore().getCore().values;

    		cusp::array1d<IndexType, MemorySpace> new_indices(indices.size());
    		cusp::array1d<CoeffType, MemorySpace> new_values(values.size());

    		thrust::copy_if(indices.begin(), indices.end(), &indices[0], &new_indices[0], is_grade<IndexType>(Multivector::get_N(), grade));
    		thrust::copy_if(values.begin(), values.end(), &indices[0], &new_values[0], is_grade<IndexType>(Multivector::get_N(), grade));

    	    new_indices = m->getCore().erase_zeros<cusp::array1d<IndexType, MemorySpace>, cusp::array1d<CoeffType, MemorySpace>>(new_indices, new_values);
    	    new_values = m->getCore().erase_zeros<cusp::array1d<CoeffType, MemorySpace>, cusp::array1d<CoeffType, MemorySpace>>(new_values, new_values);

    		std::vector<IndexType> dense_shape = m->getCore().getDenseShape();

    		SparseTensor<IndexType, CoeffType, MemorySpace> new_core(new_indices, new_values, dense_shape, false);
    		return new Multivector(new_core);
    	}
    	// TODO handler for exception
    	return NULL;
    }

    boost::python::list FACT_BLADE(Multivector *m) {
    	boost::python::list list;

    	// Perwass approach
    	// Multivector *A = new Multivector(*m);
    	// list.append(e(compute_max_arg_projection(A)));
    	// for (int i = 0; i < 1; i++) {
    	// 	auto *nj = LCONT(e(compute_max_arg_projection(A)), INVERSE(A));
    	// 	auto *nj_ = PROD(nj, 1.0 / NORM(nj));
    	// 	list.append(nj_);
    	// 	A = LCONT(INVERSE(nj_), A);
    	// }

    	//Fernandes et al approach
    	auto component_max_projection = m->get_component_max_projection();
    	// auto scalar = NORM(m);
    	auto temp = m;//PROD(m, 1.0 / scalar);

    	// list.append(scalar);
    	// list.append(temp);

    	auto all_idx = component_max_projection->getComponentIndexes();
    	auto decompose_base = idx2ji<IndexType>(Multivector::get_N());
    	std::set<IndexType> all_basis;
    	for (IndexType idx : all_idx) {
    		thrust::tuple<IndexType, IndexType> pair = decompose_base(idx);
    		if (thrust::get<0>(pair) != 0) {
    			all_basis.insert(thrust::get<0>(pair));
    		}
    		all_basis.insert(thrust::get<1>(pair));
    	}
    	std::vector<IndexType> basis(all_basis.begin(), all_basis.end());

    	for (IndexType i = 0; i < basis.size() - 1; i++) {
    		auto *nj = LCONT(e(basis[i]), INVERSE(temp));
    		auto *fatorj = PROD(nj, 1.0/NORM(nj));
    		list.append(fatorj);
    		temp = LCONT(INVERSE(fatorj), temp);
    	}
    	list.append(temp);
    	//
    	// for (IndexType i : idx) {
    	// }
    	// for (IndexType i = 0; i < basis.size(); i++) {
    	// }
    	//
    	//
    	// std::vector<Multivector*> v = {m1, m2};
    	return list;
    }

    boost::python::list FACT_VERSOR(Multivector *V) {
    	boost::python::list list;

    	Multivector *rev_V = REVERSE(V);
    	std::vector<int> grades_V = rev_V->get_grade();
    	// int k = 0;
    	while (!(grades_V.size() == 1 && grades_V[0] == 0)) {
    		// for (int i : grades_V) {
    			// std::cout << "GRADES_V: " << i << " IT: " << k << std::endl;
    		// }
    		// k++;
    		auto *A = take_grade(rev_V, grades_V[grades_V.size() - 1]);
    		auto vectors = FACT_BLADE(A);
    		int i = 0;
    		while (dot(boost::python::extract<Multivector*>(vectors[i]), boost::python::extract<Multivector*>(vectors[i]))->getCore().getCore().values[0] == 0) {
    			i++;
    		}
    		Multivector *n = boost::python::extract<Multivector*>(vectors[i]);
    		list.append(n);
    		rev_V = GP(rev_V, n);
    		grades_V = rev_V->get_grade();
    	}
    	list.append(rev_V);

    	return list;
    }
}

#endif
