#ifndef __TBGAL_MATRIX_TRAITS_HPP__
#define __TBGAL_MATRIX_TRAITS_HPP__

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
        constexpr decltype(auto) cols(MatrixType const &) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &) noexcept;

        template<typename ScalarType, std::int32_t SizeAtCompileTime>
        struct identity_matrix_type;

        template<typename ScalarType, std::int32_t SizeAtCompileTime>
        using identity_matrix_type_t = typename identity_matrix_type<ScalarType, SizeAtCompileTime>::type;
        
        template<typename ScalarType, std::int32_t SizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(std::int32_t) noexcept;

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        struct matrix_type;

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        using matrix_type_t = typename matrix_type<ScalarType, RowsAtCompileTime, ColsAtCompileTime>::type;

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        constexpr decltype(auto) make_matrix(std::int32_t, std::int32_t) noexcept;

        template<typename SourceMatrixType, typename TargetMatrixType>
        constexpr TargetMatrixType& copy_columns(SourceMatrixType const &, std::int32_t, TargetMatrixType &, std::int32_t, std::int32_t) noexcept;

        template<typename TriangularMatrixType>
        constexpr decltype(auto) determinant_triangular_matrix(TriangularMatrixType const &, std::int32_t) noexcept;

        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&... args) noexcept;

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &, std::int32_t) noexcept;

        template<typename MatrixQType_, typename MatrixRType_>
        class BaseQRDecompositionResult {
        public:

            using MatrixQType = MatrixQType_;
            using MatrixRType = MatrixRType_;
            using IndexType = index_type_t<MatrixQType>;

            constexpr BaseQRDecompositionResult() = default;
            constexpr BaseQRDecompositionResult(BaseQRDecompositionResult const &) = default;
            constexpr BaseQRDecompositionResult(BaseQRDecompositionResult &&) = default;

            virtual MatrixQType const & matrix_q() const noexcept = 0;
            virtual MatrixRType const & matrix_r() const noexcept = 0;
            virtual IndexType rank() const noexcept = 0;
        };

        template<typename MatrixType>
        constexpr decltype(auto) qr_decomposition(MatrixType const &) noexcept;

    }

}

#endif // __TBGAL_MATRIX_TRAITS_HPP__
