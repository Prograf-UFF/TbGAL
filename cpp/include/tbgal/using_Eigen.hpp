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
            constexpr static auto value = MatrixType::ColsAtCompileTime;
        };

        template<typename MatrixType>
        struct rows_at_compile_time {
            constexpr static auto value = MatrixType::RowsAtCompileTime;
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
            identity_matrix_type_t<ScalarType, SizeAtCompileTime> result(size);
            result.setIdentity();
            return result;
        }

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        struct matrix_type {
            using type = Eigen::Matrix<ScalarType, RowsAtCompileTime, ColsAtCompileTime>;
        };

        template<typename ScalarType, std::int32_t RowsAtCompileTime, std::int32_t ColsAtCompileTime>
        constexpr decltype(auto) make_matrix(std::int32_t rows, std::int32_t cols) noexcept {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime>(rows, cols);
        }
        
        template<typename SourceMatrixType>
        struct _copy_columns_impl {
            template<typename TargetMatrixType>
            constexpr static TargetMatrixType& eval(SourceMatrixType const &source, std::int32_t first_source, TargetMatrixType &target, std::int32_t first_target, std::int32_t count) noexcept {
                for (std::int32_t offset = 0, col_source = first_source, col_target = first_target; offset < count; ++offset, ++col_source, ++col_target) {
                    for (std::int32_t row = 0; row != target.rows(); ++row) {
                        target(row, col_target) = source(row, col_source);
                    }
                }
                return target;
            }
        };
        
        template<typename SourceScalarType, int SourceSizeAtCompileTime, int SourceMaxSizeAtCompileTime>
        struct _copy_columns_impl<Eigen::DiagonalMatrix<SourceScalarType, SourceSizeAtCompileTime, SourceMaxSizeAtCompileTime> > {
            
            template<typename TargetMatrixType>
            constexpr static TargetMatrixType& eval(Eigen::DiagonalMatrix<SourceScalarType, SourceSizeAtCompileTime, SourceMaxSizeAtCompileTime> const &source, std::int32_t first_source, TargetMatrixType &target, std::int32_t first_target, std::int32_t count) noexcept {
                for (std::int32_t offset = 0, col_source = first_source, col_target = first_target; offset < count; ++offset, ++col_source, ++col_target) {
                    for (std::int32_t row = 0; row < col_source; ++row) {
                        target(row, col_target) = 0;
                    }
                    target(col_source, col_target) = source.diagonal()[col_source];
                    for (std::int32_t row = col_source + 1; row < target.rows(); ++row) {
                        target(row, col_target) = 0;
                    }
                }
                return target;
            }
        };
        
        template<typename SourceMatrixType, typename TargetMatrixType>
        constexpr TargetMatrixType& copy_columns(SourceMatrixType const &source, std::int32_t first_source, TargetMatrixType &target, std::int32_t first_target, std::int32_t count) noexcept {
            return _copy_columns_impl<SourceMatrixType>::eval(source, first_source, target, first_target, count);
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
        constexpr MatrixType& _fill_column_matrix_impl(MatrixType &target) noexcept {
            return target;
        }
        
        template<typename MatrixType, typename FirstScalarType, typename... NextScalarTypes>
        constexpr MatrixType& _fill_column_matrix_impl(MatrixType &target, FirstScalarType &&arg1, NextScalarTypes &&... args) noexcept {
            target(target.rows() - (1 + sizeof...(NextScalarTypes)), 0) = std::move(arg1);
            return _fill_column_matrix_impl(target, std::move(args)...);
        }
        
        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&... args) noexcept {
            matrix_type_t<std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>, sizeof...(ScalarTypes), 1> result;
            return _fill_column_matrix_impl(result, std::move(args)...);
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
        private:

            using Super = BaseQRDecompositionResult<
                typename Eigen::ColPivHouseholderQR<MatrixType_>::HouseholderSequenceType,
                typename Eigen::ColPivHouseholderQR<MatrixType_>::MatrixType
            >;

        public:

            using MatrixQType = typename Super::MatrixQType;
            using MatrixRType = typename Super::MatrixRType;
            using IndexType = typename Super::IndexType;

            using MatrixType = MatrixType_;

            constexpr QRDecompositionResult(QRDecompositionResult const &) = default;
            constexpr QRDecompositionResult(QRDecompositionResult &&) = default;

            constexpr QRDecompositionResult(MatrixType const &arg) :
                qr_(arg),
                matrix_q_(qr_.householderQ()) {
            }

            constexpr MatrixQType const & matrix_q() const noexcept override {
                return matrix_q_;
            }

            constexpr MatrixRType const & matrix_r() const noexcept override {
                return qr_.matrixR();
            }

            constexpr IndexType rank() const noexcept override {
                return qr_.rank();
            }

        private:

            Eigen::ColPivHouseholderQR<MatrixType> qr_;
            MatrixQType matrix_q_;
        };

        template<typename MatrixType>
        constexpr decltype(auto) qr_decomposition(MatrixType const &arg) noexcept {
             return QRDecompositionResult<MatrixType>(arg);
        }

    }

}

#endif // __TBGAL_USING_EIGEN_HPP__
