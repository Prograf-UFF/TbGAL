#ifndef __TBGAL_MATRIX_TRAITS_HPP__
#define __TBGAL_MATRIX_TRAITS_HPP__

namespace tbgal {
    
    namespace detail {

        template<typename MatrixType>
        struct scalar_type;

        template<typename MatrixType>
        using scalar_type_t = typename MatrixType::Scalar;

        template<typename MatrixType>
        struct index_type;

        template<typename MatrixType>
        using index_type_t = typename MatrixType::Index;

        template<typename FirstMatrixType, typename SecondMatrixType>
        struct common_matrix_type;

        template<typename FirstMatrixType, typename SecondMatrixType>
        using common_matrix_type_t = typename common_matrix_type<FirstMatrixType, SecondMatrixType>::type;

        template<typename FirstMatrixType, typename SecondScalarType>
        struct promote_to_common_scalar_type;

        template<typename FirstMatrixType, typename SecondScalarType>
        using promote_to_common_scalar_type_t = typename promote_to_common_scalar_type<FirstMatrixType, SecondScalarType>::type;

        template<typename ScalarType, typename MetricSpaceType>
        struct identity_matrix;

        template<typename ScalarType, typename MetricSpaceType>
        using identity_matrix_t = typename identity_matrix<ScalarType, MetricSpaceType>::type;
        
        template<typename ScalarType, typename MetricSpaceType>
        constexpr decltype(auto) make_identity_matrix(MetricSpaceType const &);

        template<typename MatrixType>
        constexpr decltype(auto) rows_count(MatrixType const &);

        template<typename MatrixType, typename IndexType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &, IndexType const &);

        template<typename MatrixType>
        constexpr bool is_symmetric(MatrixType const &);

}

}

#endif // __TBGAL_MATRIX_TRAITS_HPP__
