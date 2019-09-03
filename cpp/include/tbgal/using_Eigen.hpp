#ifndef __TBGAL_USING_EIGEN_HPP__
#define __TBGAL_USING_EIGEN_HPP__

#include "core.hpp"
#include <eigen3/Eigen/Dense>

namespace tbgal {

    namespace detail {

        template<typename FirstScalarType, int FirstRows, int FirstCols, int FirstOptions, int FirstMaxRows, int FirstMaxCols, typename SecondScalarType, int SecondRows, int SecondCols, int SecondOptions, int SecondMaxRows, int SecondMaxCols>
        struct common_type<Eigen::Matrix<FirstScalarType, FirstRows, FirstCols, FirstOptions, FirstMaxRows, FirstMaxCols>, Eigen::Matrix<SecondScalarType, SecondRows, SecondCols, SecondOptions, SecondMaxRows, SecondMaxCols> > {
            static_assert((FirstRows == Eigen::Dynamic) || (SecondRows == Eigen::Dynamic) || (FirstRows == SecondRows), "The given matrices have incompatible number of rows.");
            static_assert((FirstCols == Eigen::Dynamic) || (SecondCols == Eigen::Dynamic) || (FirstCols == SecondCols), "The given matrices have incompatible number of cols.");
            using type = Eigen::Matrix<
                std::common_type_t<FirstScalarType, SecondScalarType>,
                ((FirstRows == Eigen::Dynamic) || (SecondRows == Eigen::Dynamic)) ? Eigen::Dynamic : FirstRows,
                ((FirstCols == Eigen::Dynamic) || (SecondCols == Eigen::Dynamic)) ? Eigen::Dynamic : FirstCols
            >;
        };

        template<typename FirstScalarType, int FirstRows, int FirstCols, int FirstOptions, int FirstMaxRows, int FirstMaxCols, typename SecondScalarType, int SecondSizeAtCompileTime, int SecondMaxSizeAtCompileTime>
        struct common_type<Eigen::Matrix<FirstScalarType, FirstRows, FirstCols, FirstOptions, FirstMaxRows, FirstMaxCols>, Eigen::DiagonalMatrix<SecondScalarType, SecondSizeAtCompileTime, SecondMaxSizeAtCompileTime> > {
            static_assert((FirstRows == Eigen::Dynamic) || (SecondSizeAtCompileTime == Eigen::Dynamic) || (FirstRows == SecondSizeAtCompileTime), "The given matrices have incompatible number of rows.");
            static_assert((FirstCols == Eigen::Dynamic) || (SecondSizeAtCompileTime == Eigen::Dynamic) || (FirstCols == SecondSizeAtCompileTime), "The given matrices have incompatible number of cols.");
            using type = Eigen::Matrix<
                std::common_type_t<FirstScalarType, SecondScalarType>,
                ((FirstRows == Eigen::Dynamic) || (SecondSizeAtCompileTime == Eigen::Dynamic)) ? Eigen::Dynamic : FirstRows,
                ((FirstCols == Eigen::Dynamic) || (SecondSizeAtCompileTime == Eigen::Dynamic)) ? Eigen::Dynamic : FirstCols
            >;
        };

        template<typename FirstScalarType, int FirstRows, int FirstCols, int FirstOptions, int FirstMaxRows, int FirstMaxCols, typename SecondScalarType>
        struct common_type<Eigen::Matrix<FirstScalarType, FirstRows, FirstCols, FirstOptions, FirstMaxRows, FirstMaxCols>, SecondScalarType> {
            using type = Eigen::Matrix<std::common_type_t<FirstScalarType, SecondScalarType>, FirstRows, FirstCols, FirstOptions, FirstMaxRows, FirstMaxCols>;
        };

        template<typename FirstScalarType, int FirstSizeAtCompileTime, int FirstMaxSizeAtCompileTime, typename SecondScalarType, int SecondSizeAtCompileTime, int SecondMaxSizeAtCompileTime>
        struct common_type<Eigen::DiagonalMatrix<FirstScalarType, FirstSizeAtCompileTime, FirstMaxSizeAtCompileTime>, Eigen::DiagonalMatrix<SecondScalarType, SecondSizeAtCompileTime, SecondMaxSizeAtCompileTime> > {
            static_assert((FirstSizeAtCompileTime == Eigen::Dynamic) || (SecondSizeAtCompileTime == Eigen::Dynamic) || (FirstSizeAtCompileTime == SecondSizeAtCompileTime), "The given matrices have incompatible sizes.");
            using type = Eigen::DiagonalMatrix<
                std::common_type_t<FirstScalarType, SecondScalarType>,
                ((FirstSizeAtCompileTime == Eigen::Dynamic) || (SecondSizeAtCompileTime == Eigen::Dynamic)) ? Eigen::Dynamic : FirstSizeAtCompileTime
            >;
        };

        template<typename FirstScalarType, int FirstSizeAtCompileTime, int FirstMaxSizeAtCompileTime, typename SecondScalarType, int SecondRows, int SecondCols, int SecondOptions, int SecondMaxRows, int SecondMaxCols>
        struct common_type<Eigen::DiagonalMatrix<FirstScalarType, FirstSizeAtCompileTime, FirstMaxSizeAtCompileTime>, Eigen::Matrix<SecondScalarType, SecondRows, SecondCols, SecondOptions, SecondMaxRows, SecondMaxCols> > {
            static_assert((FirstSizeAtCompileTime == Eigen::Dynamic) || (SecondRows == Eigen::Dynamic) || (FirstSizeAtCompileTime == SecondRows), "The given matrices have incompatible number of rows.");
            static_assert((FirstSizeAtCompileTime == Eigen::Dynamic) || (SecondCols == Eigen::Dynamic) || (FirstSizeAtCompileTime == SecondCols), "The given matrices have incompatible number of cols.");
            using type = Eigen::Matrix<
                std::common_type_t<FirstScalarType, SecondScalarType>,
                ((FirstSizeAtCompileTime == Eigen::Dynamic) || (SecondRows == Eigen::Dynamic)) ? Eigen::Dynamic : FirstSizeAtCompileTime,
                ((FirstSizeAtCompileTime == Eigen::Dynamic) || (SecondCols == Eigen::Dynamic)) ? Eigen::Dynamic : FirstSizeAtCompileTime
            >;
        };

        template<typename FirstScalarType, int FirstSizeAtCompileTime, int FirstMaxSizeAtCompileTime, typename SecondScalarType>
        struct common_type<Eigen::DiagonalMatrix<FirstScalarType, FirstSizeAtCompileTime, FirstMaxSizeAtCompileTime>, SecondScalarType> {
            using type = Eigen::DiagonalMatrix<std::common_type_t<FirstScalarType, SecondScalarType>, FirstSizeAtCompileTime, FirstMaxSizeAtCompileTime>;
        };

        template<typename FirstScalarType, typename SecondScalarType, int SecondRows, int SecondCols, int SecondOptions, int SecondMaxRows, int SecondMaxCols>
        struct common_type<FirstScalarType, Eigen::Matrix<SecondScalarType, SecondRows, SecondCols, SecondOptions, SecondMaxRows, SecondMaxCols> > {
            using type = Eigen::Matrix<std::common_type_t<FirstScalarType, SecondScalarType>, SecondRows, SecondCols, SecondOptions, SecondMaxRows, SecondMaxCols>;
        };

        template<typename FirstScalarType, typename SecondScalarType, int SecondSizeAtCompileTime, int SecondMaxSizeAtCompileTime>
        struct common_type<FirstScalarType, Eigen::DiagonalMatrix<SecondScalarType, SecondSizeAtCompileTime, SecondMaxSizeAtCompileTime> > {
            using type = Eigen::DiagonalMatrix<std::common_type_t<FirstScalarType, SecondScalarType>, SecondSizeAtCompileTime, SecondMaxSizeAtCompileTime>;
        };

        template<typename MatrixType>
        struct index_type {
            using type = Eigen::Index;
        };

        template<typename MatrixType>
        struct scalar_type {
            using type = typename MatrixType::Scalar;
        };

        template<typename MatrixType>
        struct cols_at_compile_time {
            constexpr auto value = MatrixType::ColsAtCompileTime;
        };

        template<typename MatrixType>
        struct rows_at_compile_time {
            constexpr auto value = MatrixType::RowsAtCompileTime;
        };

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &arg) noexcept {
            return arg.cols();
        }

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &arg) noexcept {
            return arg.rows();
        }

        template<typename ScalarType, std::int32_t SizeAtCompileTime>
        struct identity_matrix_type {
            using type = Eigen::DiagonalMatrix<ScalarType, SizeAtCompileTime>;
        };

        template<typename ScalarType, std::int32_t SizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(std::int32_t size) noexcept {
            return identity_matrix_type_t<ScalarType, SizeAtCompileTime>::Identity(size);
        }

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        struct matrix_type {
            using type = Eigen::Matrix<ScalarType, RowsAtCompileTime, ColsAtCompileTime>;
        };

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        constexpr decltype(auto) make_matrix(std::int32_t rows, std::int32_t cols) noexcept {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime>(rows, cols);
        }
        
        template<typename... ScalarTypes>
        constexpr decltype(auto) make_column_matrix(ScalarTypes &&... args) noexcept {
            matrix_type_t<std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>, sizeof...(ScalarTypes), 1> result;
            result << std::move(args)...;
            return result;
        }

        template<typename InputMatrixType, typename ResultingMatrixType>
        constexpr ResultingMatrixType& copy_columns(InputMatrixType const &source, std::int32_t first_source, ResultingMatrixType &target, std::int32_t first_target, std::int32_t count) noexcept {
            for (std::int32_t offset = 0, col_source = first_source, col_target = first_target; offset < count; ++offset, ++col_source, ++col_target) {
                for (Eigen::Index row = 0; row != target.rows(); ++row) {
                    target(row, col_target) = source(row, col_source);
                }
            }
            return target;
        }

        template<typename TriangularMatrixType>
        constexpr decltype(auto) determinant_triangular_matrix(TriangularMatrixType const &arg, std::int32_t size) noexcept {
            using ScalarType = typename TriangularMatrixType::Scalar;
            ScalarType result = 1;
            for (std::int32_t ind = 0; ind != size; ++ind) {
                result *= arg(ind, ind);
            }
            return result;
        }

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &input, std::int32_t col) noexcept {
            using ResultingMatrixType = Eigen::Matrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime>;
            ResultingMatrixType result(input.rows(), input.cols());
            for (std::int32_t offset = 0, end = input.cols() - col; offset != end; ++offset) {
                for (Eigen::Index row = 0; row != input.rows(); ++row) {
                    result(row, offset) = input(row, col + offset);
                }
            }
            for (std::int32_t offset = 0, first = input.cols() - col; offset != col; ++offset) {
                for (Eigen::Index row = 0; row != input.rows(); ++row) {
                    result(row, first + offset) = input(row, offset);
                }
            }
            return result;
        }

        template<typename MatrixType_>
        class QRDecompositionResult final : public BaseQRDecompositionResult<
            typename Eigen::ColPivHouseholderQR<MatrixType_>::HouseholderSequenceType,
            typename Eigen::ColPivHouseholderQR<MatrixType_>::MatrixType
        > {
        public:

            using MatrixType = MatrixType_;

            constexpr QRDecompositionResult(QRDecompositionResult const &) = default;
            constexpr QRDecompositionResult(QRDecompositionResult &&) = default;

            constexpr QRDecompositionResult(MatrixType const &arg) :
                qr_(arg) {
            }

            constexpr MatrixQType const & matrix_q() const noexcept override {
                return qr_.householderQ();
            }

            constexpr MatrixRType const & matrix_r() const noexcept override {
                return qr_.matrixR();
            }

            constexpr IndexType rank() const noexcept override {
                return qr_.rank();
            }

        private:

            Eigen::ColPivHouseholderQR<MatrixType> qr_;
        };

        template<typename MatrixType>
        constexpr decltype(auto) qr_decomposition(MatrixType const &arg) noexcept {
             return QRDecompositionResult<MatrixType>(arg);
        }

    }

}

#endif // __TBGAL_USING_EIGEN_HPP__
