#ifndef __TBGAL_USING_EIGEN_HPP__
#define __TBGAL_USING_EIGEN_HPP__

#include <eigen3/Eigen/Dense>

#include "core.hpp"

namespace tbgal {
    
    namespace detail {

        template<typename FirstMatrixType, typename SecondMatrixType>
        struct common_type;
        //TODO Implementar

        template<typename ScalarType, int SizeAtCompileTime, int MaxSizeAtCompileTime>
        struct index_type<Eigen::DiagonalMatrix<ScalarType, SizeAtCompileTime, MaxSizeAtCompileTime> > {
            using type = Eigen::Index;
        };

        template<typename ScalarType, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
        struct index_type<Eigen::Matrix<ScalarType, Rows, Cols, Options, MaxRows, MaxCols> > {
            using type = Eigen::Index;
        };

        template<typename ScalarType, int SizeAtCompileTime, int MaxSizeAtCompileTime>
        struct scalar_type<Eigen::DiagonalMatrix<ScalarType, SizeAtCompileTime, MaxSizeAtCompileTime> > {
            using type = typename Eigen::DiagonalMatrix<ScalarType, SizeAtCompileTime, MaxSizeAtCompileTime>::Scalar;
        };

        template<typename ScalarType, int Rows, int Cols, int Options, int MaxRows, int MaxCols>
        struct scalar_type<Eigen::Matrix<ScalarType, Rows, Cols, Options, MaxRows, MaxCols> > {
            using type = typename Eigen::Matrix<ScalarType, Rows, Cols, Options, MaxRows, MaxCols>::Scalar;
        };

        template<typename ResultingMatrixType, typename InputMatrixType>
        constexpr ResultingMatrixType& copy_first_columns(ResultingMatrixType &result, std::size_t first, InputMatrixType const &input, std::size_t cols) noexcept {
            for (std::size_t offset = 0; offset != cols; ++offset) {
                for (Eigen::Index row = 0; row != result.rows(); ++row) {
                    result(row, first + offset) = input(row, offset);
                }
            }
            return result;
        }
        
        template<typename UpperTriangularMatrixType>
        constexpr decltype(auto) determinant_triangular(UpperTriangularMatrixType const &arg) noexcept {
            using ScalarType = typename UpperTriangularMatrixType::Scalar;
            ScalarType result = arg(0, 0);
            for (Eigen::Index ind = 1; ind < arg.rows(); ++ind) {
                result *= arg(ind, ind);
            }
            return result;
        }

        template<typename MatrixType>
        constexpr MatrixType make_matrix(std::size_t rows, std::size_t cols) noexcept {
            return MatrixType(rows, cols);
        }
        
        template<typename ScalarType, typename MetricSpaceType>
        struct identity_matrix_type {
            using type = Eigen::DiagonalMatrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime>;
        };

        template<typename ScalarType, typename MetricSpaceType>
        constexpr decltype(auto) make_identity_matrix(MetricSpaceType const &space) noexcept {
            return identity_matrix_type_t<ScalarType, MetricSpaceType>::Identity(space.dimensions(), space.dimensions());
        }

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &arg) noexcept {
            return arg.cols();
        }

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &arg) noexcept {
            return arg.rows();
        }

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &input, std::size_t col) noexcept {
            using ResultingMatrixType = Eigen::Matrix<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime>;
            ResultingMatrixType result(input.rows(), input.cols());
            for (std::size_t offset = 0, end = input.cols() - col; offset != end; ++offset) {
                for (Eigen::Index row = 0; row != input.rows(); ++row) {
                    result(row, offset) = input(row, col + offset);
                }
            }
            for (std::size_t offset = 0, first = input.cols() - col; offset != col; ++offset) {
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
