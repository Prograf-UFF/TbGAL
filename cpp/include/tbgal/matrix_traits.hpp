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

        template<typename ResultingMatrixType, typename InputMatrixType>
        constexpr ResultingMatrixType& copy_first_columns(ResultingMatrixType &, std::size_t, InputMatrixType const &, std::size_t);
        
        template<typename UpperTriangularMatrixType>
        constexpr decltype(auto) determinant_triangular(UpperTriangularMatrixType const &);

        template<typename MatrixType>
        constexpr MatrixType make_matrix(std::size_t, std::size_t);
        
        template<typename ScalarType, typename MetricSpaceType>
        struct identity_matrix_type;

        template<typename ScalarType, typename MetricSpaceType>
        using identity_matrix_type_t = typename identity_matrix_type<ScalarType, MetricSpaceType>::type;
        
        template<typename ScalarType, typename MetricSpaceType>
        constexpr decltype(auto) make_identity_matrix(MetricSpaceType const &);

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &);

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &);

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &, std::size_t);

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
        constexpr decltype(auto) qr_decomposition(MatrixType const &);

    }

}

#endif // __TBGAL_MATRIX_TRAITS_HPP__
