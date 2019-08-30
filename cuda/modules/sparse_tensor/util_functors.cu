#ifndef UTIL_FUNCTORS
#define UTIL_FUNCTORS

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>

#include <thrust/fill.h>
#include <thrust/device_vector.h>

template <typename T>
struct is_grade {
	is_grade(IndexType N, T grade) : grade_(grade), N_(N) {}
	T grade_;
	IndexType N_;

	__host__ __device__ bool operator()(const T &i) {
		if (grade_ == 0) {
			return i == 0;
		} else if (grade_ == 1) {
			return i > 0 && i <= N_;
		} else if (grade_ == 2) {
			return i > N_;
		}
		return false;
	}
};

template <typename T>
struct is_component {
	is_component(T component) : component_(component) {}
	T component_;

	__host__ __device__ bool operator()(const T &i) {
		return component_ == i;
	}
};

template <typename T>
struct compare_given_threshold {

	compare_given_threshold(CoeffType threshold) : threshold(threshold) {}
	CoeffType threshold;

  __host__ __device__ bool operator()(T x, T y) const {
	return abs(x-y) < threshold;
  }
};
#endif
