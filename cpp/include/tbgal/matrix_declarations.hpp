#ifndef __TBGAL_MATRIX_DECLARATIONS_HPP__
#define __TBGAL_MATRIX_DECLARATIONS_HPP__

namespace tbgal {
    
    namespace detail {

        template<typename FirstType, typename SecondType>
        struct common_type : std::common_type<FirstType, SecondType> {
        };

        template<typename FirstType, typename SecondType>
        using common_type_t = typename common_type<FirstType, SecondType>::type;

        template<typename MatrixType>
        struct index_type;

        template<typename MatrixType>
        using index_type_t = typename MatrixType::Index;

        template<typename MatrixType>
        struct scalar_type;

        template<typename MatrixType>
        using scalar_type_t = typename MatrixType::Scalar;

        template<typename MatrixType>
        struct cols_at_compile_time;
        
        template<typename MatrixType>
        constexpr auto cols_at_compile_time_v = cols_at_compile_time<MatrixType>::value;

        template<typename MatrixType>
        struct rows_at_compile_time;
        
        template<typename MatrixType>
        constexpr auto rows_at_compile_time_v = rows_at_compile_time<MatrixType>::value;

        template<typename MatrixType>
        constexpr decltype(auto) coeff(MatrixType const &, DefaultIndexType, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &) noexcept;

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime>
        struct identity_matrix_type;

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime>
        using identity_matrix_type_t = typename identity_matrix_type<ScalarType, SizeAtCompileTime>::type;
        
        template<typename ScalarType, DefaultIndexType SizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(DefaultIndexType) noexcept;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime>
        struct matrix_type;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime>
        using matrix_type_t = typename matrix_type<ScalarType, RowsAtCompileTime, ColsAtCompileTime>::type;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime>
        constexpr decltype(auto) make_matrix(DefaultIndexType, DefaultIndexType) noexcept;

        template<typename SourceMatrixType, typename TargetMatrixType>
        constexpr TargetMatrixType& copy_columns(SourceMatrixType const &, DefaultIndexType, TargetMatrixType &, DefaultIndexType, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) determinant(MatrixType const &) noexcept;
        
        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&... args) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) left_columns(MatrixType const &, DefaultIndexType) noexcept;
        
        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod(FirstMatrixType const &, SecondMatrixType const &) noexcept;
        
        template<typename MatrixType>
        constexpr decltype(auto) qr_decomposition(MatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &, DefaultIndexType) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) transpose(MatrixType const &) noexcept;

    }

}

#endif // __TBGAL_MATRIX_DECLARATIONS_HPP__
