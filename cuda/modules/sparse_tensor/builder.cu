#ifndef BUILDER_CU
#define BUILDER_CU

// typedef unsigned long long int IndexType;
// typedef float CoeffType;
#include "../../common.cu"

// #include <cusp/coo_matrix.h>

// typedef cusp::device_memory MemorySpace;
#include "SparseTensor.cu"
#include <boost/python.hpp>

namespace python = boost::python;

BOOST_PYTHON_MODULE(sparse_tensor) {
	python::class_<SparseTensor<IndexType, CoeffType, MemorySpace>, boost::noncopyable>("SparseTensor", python::no_init);
}

#endif
