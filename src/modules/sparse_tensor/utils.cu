#ifndef UTILS_CU
#define UTILS_CU

#include <thrust/iterator/counting_iterator.h>
#include <thrust/iterator/transform_iterator.h>
#include <thrust/iterator/permutation_iterator.h>
#include <thrust/functional.h>
#include <thrust/fill.h>
#include <thrust/device_vector.h>
#include "../metric/metric.cu"
#include "util_functors.cu"
#include <unistd.h>


enum Axis {
	I,
	J,
	K,
	EXISTS,
	VALUES
};

enum Operation {
	OUTER_PRODUCT,
	LEFT_CONTRACTION,
	RIGHT_CONTRACTION,
	DOT_PRODUCT
};


template <typename Iterator>
class repeated_range {
    public:

	    typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	    struct repeat_functor : public thrust::unary_function<difference_type,difference_type> {
	        difference_type repeats;

	        repeat_functor(difference_type repeats) : repeats(repeats) {}

	        __host__ __device__ difference_type operator()(const difference_type& i) const {
	            return i / repeats;
	        }
	    };

	    typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	    typedef typename thrust::transform_iterator<repeat_functor, CountingIterator> TransformIterator;
	    typedef typename thrust::permutation_iterator<Iterator,TransformIterator> PermutationIterator;

	    typedef PermutationIterator iterator;

	    repeated_range(Iterator first, Iterator last, difference_type repeats) : first(first), last(last), repeats(repeats) {}

	    iterator begin(void) const {
	        return PermutationIterator(first, TransformIterator(CountingIterator(0), repeat_functor(repeats)));
	    }

	    iterator end(void) const {
	        return begin() + repeats * (last - first);
	    }

    protected:
	    Iterator first;
	    Iterator last;
	    difference_type repeats;
};


template <typename Iterator>
class tiled_range {
    public:

	    typedef typename thrust::iterator_difference<Iterator>::type difference_type;

	    struct tile_functor : public thrust::unary_function<difference_type,difference_type> {
	        difference_type tile_size;

	        tile_functor(difference_type tile_size) : tile_size(tile_size) {}

	        __host__ __device__ difference_type operator()(const difference_type& i) const {
	            return i % tile_size;
	        }
	    };

	    typedef typename thrust::counting_iterator<difference_type> CountingIterator;
	    typedef typename thrust::transform_iterator<tile_functor, CountingIterator> TransformIterator;
	    typedef typename thrust::permutation_iterator<Iterator,TransformIterator> PermutationIterator;

	    typedef PermutationIterator iterator;

	    tiled_range(Iterator first, Iterator last, difference_type tiles) : first(first), last(last), tiles(tiles) {}

	    iterator begin(void) const {
	        return PermutationIterator(first, TransformIterator(CountingIterator(0), tile_functor(last - first)));
	    }

	    iterator end(void) const {
	        return begin() + tiles * (last - first);
	    }

    protected:
	    Iterator first;
	    Iterator last;
	    difference_type tiles;
};


template <typename T>
__host__ __device__ T getY(T idx, T N_) {
	return (T)(N_ - 2 - static_cast<T>(sqrt(-8 * idx + 4 * N_ *(N_ - 1) - 7) / 2.0 - 0.5));
}


template <typename T>
__host__ __device__ T getX(T idx, T N_) {
	auto Y = getY(idx, N_);
	return idx + Y + 1 - N_ * (N_ - 1) / 2 + (N_ - Y) * ((N_ - Y - 1)) / 2;
}


template <typename T>
__host__ __device__ T ji2idx(T j, T i, T N) {
	return ((N * (N - 1)) >> 1) - (((N - j) * (N - j - 1)) >> 1) + i;
}


template <typename T>
struct idx2ji {
	__host__ __device__ idx2ji(IndexType N) : N_(N + 1) {}
	T N_;
	__host__ __device__	thrust::tuple<T, T> operator()(const T& idx) const {
		if (idx == 0) {
			return thrust::make_tuple<T, T>(0, 0);
		}
		T idx_ = idx - 1;
		thrust::tuple<T, T> to_return = thrust::make_tuple<T, T>(getY(idx_, N_), getX(idx_, N_));
		return to_return;
	}
};


__host__ __device__ CoeffType canonical_sort(const IndexType &u, const IndexType &v, const IndexType &x, const IndexType &y) {
	IndexType sorter[4];
	sorter[0] = u;
	sorter[1] = v;
	sorter[2] = x;
	sorter[3] = y;

	int count_changes = 0;
	bool changed = false;
	do {
		changed = false;
		for (int i = 0; i < 3; i++) {
			if ((sorter[i] > sorter[i + 1])) {
				IndexType aux = sorter[i];
				sorter[i] = sorter[i + 1];
				sorter[i + 1] = aux;
				changed = true;
				if (sorter[i] != 0 && sorter[i + 1] != 0) {
					count_changes++;
				}
			}
		}
	} while (changed);
	return static_cast<CoeffType>(pow(-1, count_changes));
};


__host__ __device__ IndexType output_basis(const IndexType &i, const IndexType &j, const IndexType &N) {
	if (i == 0) {
		return ji2idx<IndexType>(i, j, N);
	} else if (j == 0) {
		return ji2idx<IndexType>(j, i, N);
	} else if (i < j) {
		return ji2idx<IndexType>(i, j, N);
	} else if (i > j){
		return ji2idx<IndexType>(j, i, N);
	} else {
		return static_cast<IndexType>(0);
	}
};

template < class MetricType, typename ReturnType, class = typename std::enable_if<std::is_base_of<Metric, MetricType>::value>::type>
struct GeometricProductTensorFunctor {
	GeometricProductTensorFunctor(IndexType N, const MetricType &metric, const Axis &axis) : metric(metric), axis(axis), N_(N) {}

	short axis;
	MetricType metric;
	IndexType N_;

	__host__ __device__ ReturnType operator() (const thrust::tuple<IndexType, IndexType> &it) {
		IndexType j = thrust::get<0>(it);
		IndexType i = thrust::get<1>(it);

		if (axis == Axis::I) {
			return i;
		} else if (axis == Axis::J) {
			return j;
		}

		IndexType K = 0;
		CoeffType val = 0;

		auto my_functor = idx2ji<IndexType>(N_);
		thrust::tuple<IndexType, IndexType> basis_i = my_functor(i);
		thrust::tuple<IndexType, IndexType> basis_j = my_functor(j);

		IndexType u = thrust::get<0>(basis_j);
		IndexType v = thrust::get<1>(basis_j);
		IndexType x = thrust::get<0>(basis_i);
		IndexType y = thrust::get<1>(basis_i);

		if (i == 0) {
			K = j;
			val = 1;
		} else if (j == 0) {
			K = i;
			val = 1;
		} else {
			if (i == j) {
				K = 0;
				val = canonical_sort(u, v, x, y) * metric.metric_factor(u, v);
			} else if (u == x) {
				K = output_basis(v, y, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(u); // vy
			} else if (u == y) {
				K = output_basis(v, x, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(u); // vx
			} else if (v == x) {
				K = output_basis(u, y, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(v); // uy
			} else if (v == y) {
				K = output_basis(u, x, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(v); // ux
			}
		}

		if (axis == Axis::K) {
			return K;
		} else if (axis == Axis::VALUES) {
			return val;
		} else if (axis == Axis::EXISTS) {
			return val != 0;
		}
		return 0;
	}
};

template < class MetricType, typename ReturnType, class = typename std::enable_if<std::is_base_of<Metric, MetricType>::value>::type>
struct NONGeometricProductTensorFunctor {
	NONGeometricProductTensorFunctor(IndexType N, const MetricType &metric, const Axis &axis) : metric(metric), axis(axis), N_(N) {}

	short axis;
	MetricType metric;
	IndexType N_;

	__host__ __device__ ReturnType operator() (const thrust::tuple<IndexType, IndexType> &it) {
		IndexType j = thrust::get<0>(it);
		IndexType i = thrust::get<1>(it);

		if (axis == Axis::I) {
			return i;
		} else if (axis == Axis::J) {
			return j;
		}

		IndexType K = 0;
		CoeffType val = 0;

		auto my_functor = idx2ji<IndexType>(N_);
		thrust::tuple<IndexType, IndexType> basis_i = my_functor(i);
		thrust::tuple<IndexType, IndexType> basis_j = my_functor(j);

		IndexType u = thrust::get<0>(basis_j);
		IndexType v = thrust::get<1>(basis_j);
		IndexType x = thrust::get<0>(basis_i);
		IndexType y = thrust::get<1>(basis_i);

		if (i == 0) {
			K = j;
			val = 1;
		} else if (j == 0) {
			K = i;
			val = 1;
		} else {
			if (i == j) {
				K = 0;
				val = canonical_sort(u, v, x, y) * metric.metric_factor(u, v);
			} else if (u == x) {
				K = output_basis(v, y, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(u); // vy
			} else if (u == y) {
				K = output_basis(v, x, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(u); // vx
			} else if (v == x) {
				K = output_basis(u, y, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(v); // uy
			} else if (v == y) {
				K = output_basis(u, x, N_);
				val = canonical_sort(u, v, x, y) * metric.diagonal_entry(v); // ux
			}
		}

		if (axis == Axis::K) {
			return K;
		} else if (axis == Axis::VALUES) {
			return val;
		} else if (axis == Axis::EXISTS) {
			return val == 0;
		}
		return 0;
	}
};


// template < class MetricType, class = typename std::enable_if<std::is_base_of<Metric, MetricType>::value>::type>
// SparseTensor<IndexType, CoeffType, cusp::device_memory> build_geometric_product_tensor(IndexType N, MetricType metric) {
//
// 	IndexType full_size = (N * (N + 1) >> 1) + 1;
//
// 	cusp::array1d<IndexType, cusp::host_memory> all_basis = cusp::counting_array<IndexType>(full_size);
//
// 	typedef typename thrust::host_vector<IndexType>::iterator Iterator;
// 	repeated_range<Iterator> all_basis_repeated(all_basis.begin(), all_basis.end(), full_size);
// 	tiled_range<Iterator> all_basis_tiled(all_basis.begin(), all_basis.end(), full_size);
//
// 	typedef typename repeated_range<Iterator>::iterator repeated_it;
// 	typedef typename tiled_range<Iterator>::iterator tiled_it;
//
// 	typedef thrust::tuple<repeated_it, tiled_it> IteratorTuple;
// 	typedef thrust::zip_iterator<IteratorTuple> ZipIterator;
//
//
// 	ZipIterator zbegin = thrust::make_zip_iterator(thrust::make_tuple(all_basis_repeated.begin(), all_basis_tiled.begin()));
// 	ZipIterator zend = thrust::make_zip_iterator(thrust::make_tuple(all_basis_repeated.begin(), all_basis_tiled.begin())) + (full_size * full_size);
//
// 	IndexType size = thrust::count_if<ZipIterator>(zbegin, zend, GeometricProductTensorFunctor<MetricType, bool>(N, metric, Axis::EXISTS));
//
// 	// std::cout << "FULL SIZE: " << full_size * full_size << std::endl;
// 	// std::cout << "SIZE: " << size << std::endl;
//
// 	cusp::array1d<thrust::tuple<IndexType, IndexType>, cusp::host_memory> tuples(size);
//
// 	thrust::copy_if(thrust::host, zbegin, zend, tuples.begin(), GeometricProductTensorFunctor<MetricType, bool>(N, metric, Axis::EXISTS));
//
// 	cusp::array1d<IndexType, cusp::host_memory> K(size);
// 	thrust::transform(thrust::host, tuples.begin(), tuples.end(),
// 		K.begin(),
// 		GeometricProductTensorFunctor<MetricType, IndexType>(N, metric, Axis::K));
//
// 	cusp::array1d<IndexType, cusp::host_memory> I(size);
// 	thrust::transform(thrust::host, tuples.begin(), tuples.end(),
// 		I.begin(),
// 		GeometricProductTensorFunctor<MetricType, IndexType>(N, metric, Axis::I));
//
// 	cusp::array1d<IndexType, cusp::host_memory> J(size);
// 	thrust::transform(thrust::host, tuples.begin(), tuples.end(),
// 		J.begin(),
// 		GeometricProductTensorFunctor<MetricType, IndexType>(N, metric, Axis::J));
//
//
// 	cusp::array1d<CoeffType, cusp::host_memory> values(size);
// 	thrust::transform(thrust::host, tuples.begin(), tuples.end(),
// 		values.begin(),
// 		GeometricProductTensorFunctor<MetricType, CoeffType>(N, metric, Axis::VALUES));
//
// 	// std::cout << "I: " << *(thrust::max_element(I.begin(), I.end())) << std::endl;
// 	// std::cout << "J: " << *(thrust::max_element(J.begin(), J.end())) << std::endl;
// 	// std::cout << "K: " << *(thrust::max_element(K.begin(), K.end())) << std::endl;
//
// 	std::vector<IndexType> shape = { full_size, full_size, full_size };
//
// 	FILE *pfile;
// 	pfile = fopen("/home/eduardovera/values.bin", "wb");
// 	fwrite(&(*values.begin()), sizeof(values), values.size(), pfile);
// 	fclose(pfile);
//
// 	cusp::array1d<IndexType, cusp::device_memory> I_dev(I.begin(), I.end());
// 	cusp::array1d<IndexType, cusp::device_memory> J_dev(J.begin(), J.end());
// 	cusp::array1d<IndexType, cusp::device_memory> K_dev(K.begin(), K.end());
// 	cusp::array1d<CoeffType, cusp::device_memory> values_dev(values.begin(), values.end());
//
//
// 	return SparseTensor<IndexType, CoeffType, cusp::device_memory>(I_dev, J_dev, K_dev, values_dev, shape);
// }


template < class MetricType, class = typename std::enable_if<std::is_base_of<Metric, MetricType>::value>::type>
SparseTensor<IndexType, CoeffType, cusp::device_memory> build_geometric_product_tensor(IndexType N, MetricType metric) {

	IndexType full_size = (N * (N + 1) >> 1) + 1;

	cusp::array1d<IndexType, cusp::device_memory> all_basis = cusp::counting_array<IndexType>(full_size);

	typedef typename thrust::device_vector<IndexType>::iterator Iterator;
	repeated_range<Iterator> all_basis_repeated(all_basis.begin(), all_basis.end(), full_size);
	tiled_range<Iterator> all_basis_tiled(all_basis.begin(), all_basis.end(), full_size);

	typedef typename repeated_range<Iterator>::iterator repeated_it;
	typedef typename tiled_range<Iterator>::iterator tiled_it;

	typedef thrust::tuple<repeated_it, tiled_it> IteratorTuple;
	typedef thrust::zip_iterator<IteratorTuple> ZipIterator;


	ZipIterator zbegin = thrust::make_zip_iterator(thrust::make_tuple(all_basis_repeated.begin(), all_basis_tiled.begin()));
	ZipIterator zend = thrust::make_zip_iterator(thrust::make_tuple(all_basis_repeated.begin(), all_basis_tiled.begin())) + (full_size * full_size);

	std::cout << "SIZE: " << thrust::distance(zbegin, zend) << std::endl;

	IndexType size = thrust::count_if<ZipIterator>(zbegin, zend, GeometricProductTensorFunctor<MetricType, bool>(N, metric, Axis::EXISTS));

	std::cout << "SIZE: " << size << std::endl;

	// cusp::array1d<thrust::tuple<IndexType, IndexType>, cusp::device_memory> tuples(size);
	// cusp::array1d<ZipIterator, cusp::device_memory> tuples(size);

	cusp::array1d<IndexType, cusp::device_memory> indices(size);

	thrust::counting_iterator<IndexType> first_index(0);

	thrust::counting_iterator<IndexType> last_index(thrust::distance(zbegin, zend));


	thrust::copy_if(first_index, last_index, zbegin, indices.begin(), GeometricProductTensorFunctor<MetricType, bool>(N, metric, Axis::EXISTS));
	// ZipIterator nzend = thrust::remove_copy_if(zbegin, zend, tuples.begin(), NONGeometricProductTensorFunctor<MetricType, bool>(N, metric, Axis::EXISTS));

	// std::cout << "SIZE: " << tuples.end() - tuples.begin() << std::endl;

	typedef thrust::device_vector<IndexType>::iterator IndexIterator;

	thrust::permutation_iterator<ZipIterator, IndexIterator> permiter_begin(zbegin, indices.begin());
	thrust::permutation_iterator<ZipIterator, IndexIterator> permiter_end(zbegin, indices.end());

	// delete indices?

	cusp::array1d<IndexType, cusp::device_memory> K(size);
	thrust::transform(permiter_begin, permiter_end,
		K.begin(),
		GeometricProductTensorFunctor<MetricType, IndexType>(N, metric, Axis::K));

	cusp::array1d<IndexType, cusp::device_memory> I(size);
	thrust::transform(permiter_begin, permiter_end,
		I.begin(),
		GeometricProductTensorFunctor<MetricType, IndexType>(N, metric, Axis::I));

	cusp::array1d<IndexType, cusp::device_memory> J(size);
	thrust::transform(permiter_begin, permiter_end,
		J.begin(),
		GeometricProductTensorFunctor<MetricType, IndexType>(N, metric, Axis::J));

	cusp::array1d<CoeffType, cusp::device_memory> values(size);
	thrust::transform(permiter_begin, permiter_end,
		values.begin(),
		GeometricProductTensorFunctor<MetricType, CoeffType>(N, metric, Axis::VALUES));

	// std::cout << "I: " << *(thrust::max_element(I.begin(), I.end())) << std::endl;
	// std::cout << "J: " << *(thrust::max_element(J.begin(), J.end())) << std::endl;
	// std::cout << "K: " << *(thrust::max_element(K.begin(), K.end())) << std::endl;

	std::vector<IndexType> shape = { full_size, full_size, full_size };

	return SparseTensor<IndexType, CoeffType, cusp::device_memory>(I, J, K, values, shape);
}



template <typename ReturnType>
struct get_all_grades_and_compare {
	IndexType N_;
	IndexType inv_N_FULL;
	IndexType N_FULL;
	Axis axis;
	Operation operation;
	get_all_grades_and_compare(IndexType N, Operation operation, Axis axis) {
		this->N_ = N;
		this->N_FULL = (N_ * (N_ - 1) >> 1);
		this->axis = axis;
		this->operation = operation;
	}

	__host__ __device__ ReturnType operator()(const thrust::tuple<IndexType, IndexType, IndexType, CoeffType> &t) {
		IndexType r = thrust::get<0>(t);
		IndexType c = thrust::get<1>(t);
		IndexType l = thrust::get<2>(t);

		int grade_j = (r == static_cast<IndexType>(0)) ? 0 : (r <= this->N_ ) ? 1 : 2;
		int grade_i = (c == static_cast<IndexType>(0)) ? 0 : (c <= this->N_)  ? 1 : 2;
		int grade_k = (l == static_cast<IndexType>(0)) ? 0 : (l <= this->N_)  ? 1 : 2;

		bool itFits = false;

		if ((operation == Operation::OUTER_PRODUCT) && (grade_i + grade_j == grade_k) ) {
			itFits = true;
		} else if ((operation == Operation::LEFT_CONTRACTION) && (grade_i - grade_j == grade_k) ) {
			itFits = true;
		} else if ((operation == Operation::RIGHT_CONTRACTION) && (grade_j - grade_i == grade_k) ) {
			itFits = true;
		} else if ((operation == Operation::DOT_PRODUCT) && (abs(grade_i - grade_j) == grade_k) ) {
			itFits = true;
		}

		if (itFits) {
			if (this->axis == Axis::J) {
				return r;
			} else if (this->axis == Axis::I) {
				return c;
			} else if (this->axis == Axis::K) {
				return l;
			} else if (this->axis == Axis::VALUES) {
				return thrust::get<3>(t);
			} else if (this->axis == Axis::EXISTS) {
				return true;
			}
		}
		return static_cast<ReturnType>(0);
	}
};

SparseTensor<IndexType, CoeffType, cusp::device_memory>* extract_tensor(IndexType N, SparseTensor<IndexType, CoeffType, cusp::device_memory> *tensor, const Operation &operation) {

	typedef typename cusp::array1d<IndexType, cusp::device_memory>::iterator Iterator_Index;
	typedef typename cusp::array1d<CoeffType, cusp::device_memory>::iterator Iterator_Coeff;

	typedef thrust::tuple<Iterator_Index, Iterator_Index, Iterator_Index, Iterator_Coeff> IteratorTuple;
	typedef thrust::zip_iterator<IteratorTuple> ZipIterator;

	ZipIterator zbegin = thrust::make_zip_iterator(thrust::make_tuple(tensor->J.begin(), tensor->I.begin(), tensor->K.begin(), tensor->data.begin()));
	ZipIterator zend = zbegin + tensor->data.size();

	IndexType size = thrust::count_if<ZipIterator>(zbegin, zend, get_all_grades_and_compare<bool>(N, operation, Axis::EXISTS));

	thrust::device_vector<thrust::tuple<IndexType, IndexType, IndexType, CoeffType>> tuples(size);
	thrust::copy_if(zbegin, zend, tuples.begin(), get_all_grades_and_compare<bool>(N, operation, Axis::EXISTS));

	cusp::array1d<CoeffType, cusp::device_memory> newdata(size);
	cusp::array1d<IndexType, cusp::device_memory> newc(size);
	cusp::array1d<IndexType, cusp::device_memory> newr(size);
	cusp::array1d<IndexType, cusp::device_memory> newl(size);

	thrust::transform(tuples.begin(), tuples.end(),
					&newr[0],
					get_all_grades_and_compare<IndexType>(N, operation, Axis::J));

	thrust::transform(tuples.begin(), tuples.end(),
					&newc[0],
					get_all_grades_and_compare<IndexType>(N, operation, Axis::I));

	thrust::transform(tuples.begin(), tuples.end(),
					&newl[0],
					get_all_grades_and_compare<IndexType>(N, operation, Axis::K));

	thrust::transform(tuples.begin(), tuples.end(),
					&newdata[0],
					get_all_grades_and_compare<CoeffType>(N, operation, Axis::VALUES));

	return new SparseTensor<IndexType, CoeffType, cusp::device_memory>(newc, newr, newl, newdata, tensor->getDenseShape());

}

#endif
