#ifndef __TBGAL_MATRIX_DECLARATIONS_HPP__
#define __TBGAL_MATRIX_DECLARATIONS_HPP__

#include <cmath>
#include <cstdint>
#include <type_traits>

#ifndef TBGAL_DEFAULT_SCALAR_TYPE
    #define TBGAL_DEFAULT_SCALAR_TYPE std::double_t
#endif // TBGAL_DEFAULT_SCALAR_TYPE

#ifndef TBGAL_DEFAULT_INDEX_TYPE
    #define TBGAL_DEFAULT_INDEX_TYPE std::int64_t
#endif // TBGAL_DEFAULT_INDEX_TYPE

namespace tbgal {
    
    using DefaultIndexType = TBGAL_DEFAULT_INDEX_TYPE;
    using DefaultScalarType = TBGAL_DEFAULT_SCALAR_TYPE;

    constexpr static DefaultIndexType Dynamic = -1;

    namespace detail {

        template<typename MatrixType>
        struct index_type;

        template<typename MatrixType>
        using index_type_t = typename index_type<MatrixType>::type;

        template<typename MatrixType>
        struct scalar_type;

        template<typename MatrixType>
        using scalar_type_t = typename scalar_type<MatrixType>::type;

        template<typename MatrixType>
        struct cols_at_compile_time;
        
        template<typename MatrixType>
        constexpr auto cols_at_compile_time_v = cols_at_compile_time<MatrixType>::value;

        template<typename MatrixType>
        struct rows_at_compile_time;
        
        template<typename MatrixType>
        constexpr auto rows_at_compile_time_v = rows_at_compile_time<MatrixType>::value;

        template<typename MatrixType>
        constexpr decltype(auto) coeff(MatrixType &, DefaultIndexType, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &) noexcept;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        struct matrix_type;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        using matrix_type_t = typename matrix_type<ScalarType, RowsAtCompileTime, ColsAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime>::type;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_matrix(DefaultIndexType, DefaultIndexType) noexcept;

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(DefaultIndexType) noexcept;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_zero_matrix(DefaultIndexType, DefaultIndexType) noexcept;

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceType, typename TargetMatrixType>
        constexpr void assign_block(SourceType const &, TargetMatrixType &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType) noexcept;

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &, DefaultIndexType, DefaultIndexType, TargetMatrixType &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr void conservative_resize(MatrixType &, DefaultIndexType, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) determinant(MatrixType const &) noexcept;
        
        template<typename MatrixType>
        constexpr decltype(auto) evaluate(MatrixType const &arg) noexcept;

        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&...) noexcept;

        template<DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime, typename IteratorType>
        constexpr decltype(auto) fill_column_matrix_using_iterators(IteratorType, IteratorType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) inverse(MatrixType const &) noexcept;

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod(FirstMatrixType const &, SecondMatrixType const &) noexcept;

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType, SecondMatrixType const &) noexcept;

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &, SecondMatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) qr_orthogonal_matrix(MatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &, DefaultIndexType) noexcept;

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) sub(FirstMatrixType const &, SecondMatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) svd_left_eigenvectors(MatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) transpose(MatrixType const &) noexcept;

    }

}

#endif // __TBGAL_MATRIX_DECLARATIONS_HPP__
