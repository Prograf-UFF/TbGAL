#ifndef __TBGAL_USING_EIGEN_HPP__
#define __TBGAL_USING_EIGEN_HPP__

#include <Eigen/Dense>
#include "matrix_declarations.hpp"

#ifdef TBGAL_USING_MATRIX_DEFINITIONS
    #error "Matrix definition already included by command '#include <using_{SomeMatrixAlgebraLibrary}.hpp>.'"
#else
    #define TBGAL_USING_MATRIX_DEFINITIONS
#endif // TBGAL_USING_MATRIX_DEFINITIONS

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

        template<typename FirstScalarType, int FirstRows, int FirstCols, int FirstOptions, int FirstMaxRows, int FirstMaxCols, typename SecondScalarType>
        struct common_type<Eigen::Matrix<FirstScalarType, FirstRows, FirstCols, FirstOptions, FirstMaxRows, FirstMaxCols>, SecondScalarType> {
            using type = Eigen::Matrix<std::common_type_t<FirstScalarType, SecondScalarType>, FirstRows, FirstCols, FirstOptions, FirstMaxRows, FirstMaxCols>;
        };

        template<typename FirstScalarType, typename SecondScalarType, int SecondRows, int SecondCols, int SecondOptions, int SecondMaxRows, int SecondMaxCols>
        struct common_type<FirstScalarType, Eigen::Matrix<SecondScalarType, SecondRows, SecondCols, SecondOptions, SecondMaxRows, SecondMaxCols> > {
            using type = Eigen::Matrix<std::common_type_t<FirstScalarType, SecondScalarType>, SecondRows, SecondCols, SecondOptions, SecondMaxRows, SecondMaxCols>;
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
        constexpr decltype(auto) coeff(MatrixType &arg, DefaultIndexType row, DefaultIndexType col) noexcept {
            return arg(row, col);
        }

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &arg) noexcept {
            return arg.cols();
        }

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &arg) noexcept {
            return arg.rows();
        }

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime>
        struct matrix_type {
            using type = Eigen::Matrix<ScalarType, RowsAtCompileTime, ColsAtCompileTime>;
        };

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime>
        constexpr decltype(auto) make_matrix(DefaultIndexType rows, DefaultIndexType cols) noexcept {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime>(rows, cols);
        }
        
        template<typename ScalarType, DefaultIndexType SizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(DefaultIndexType size) noexcept {
            return matrix_type_t<ScalarType, SizeAtCompileTime, SizeAtCompileTime>::Identity(size, size);
        }

        template<DefaultIndexType ColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr TargetMatrixType& copy_columns(SourceMatrixType const &source, DefaultIndexType start_col_source, TargetMatrixType &target, DefaultIndexType start_col_target, DefaultIndexType cols) noexcept {
            target.template block<TargetMatrixType::RowsAtCompileTime, ColsAtCompileTime>(0, start_col_target, target.rows(), cols) = source.template block<SourceMatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, start_col_source, source.rows(), cols);
            return target;
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr TargetMatrixType& copy_to_block(SourceMatrixType const &source, TargetMatrixType &target, DefaultIndexType start_row, DefaultIndexType start_col, DefaultIndexType block_rows, DefaultIndexType block_cols) noexcept {
            target.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row, start_col, block_rows, block_cols) = source;
            return target;
        }

        template<typename MatrixType>
        constexpr decltype(auto) determinant(MatrixType const &arg) noexcept {
            return arg.determinant();
        }
        
        template<typename MatrixType, typename ScalarType>
        constexpr void _fill_column_matrix_impl(MatrixType &target, ScalarType &&arg) noexcept {
            target(target.rows() - 1, 0) = std::move(arg);
        }
        
        template<typename MatrixType, typename FirstScalarType, typename... NextScalarTypes>
        constexpr void _fill_column_matrix_impl(MatrixType &target, FirstScalarType &&arg1, NextScalarTypes &&... args) noexcept {
            target(target.rows() - (1 + sizeof...(NextScalarTypes)), 0) = std::move(arg1);
            _fill_column_matrix_impl(target, std::move(args)...);
        }
        
        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&... args) noexcept {
            matrix_type_t<std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>, sizeof...(ScalarTypes), 1> result;
            _fill_column_matrix_impl(result, std::move(args)...);
            return result;
        }

        template<typename MatrixType>
        constexpr decltype(auto) left_columns(MatrixType const &arg, DefaultIndexType count) noexcept {
            return arg.leftCols(count);
        }

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod(FirstMatrixType const &arg1, SecondMatrixType const &arg2) noexcept {
            return (arg1 * arg2).eval();
        }
        
        template<typename MatrixType>
        constexpr decltype(auto) qr_decomposition(MatrixType const &arg) noexcept {
            using MatrixQType = matrix_type_t<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime>;
            Eigen::ColPivHouseholderQR<MatrixType> qr(arg);
            auto rank = qr.rank();
            return std::make_tuple(MatrixQType(qr.householderQ()), qr.matrixR().topLeftCorner(rank, rank).template triangularView<Eigen::Upper>(), rank);
        }

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &input, DefaultIndexType col) noexcept {
            using ResultingMatrixType = matrix_type_t<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime>;
            ResultingMatrixType result(input.rows(), input.cols());
            result.template block<ResultingMatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, 0, result.rows(), input.cols() - col) = input.template block<MatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, col, input.rows(), input.cols() - col);
            result.template block<ResultingMatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, input.cols() - col, result.rows(), col) = input.template block<MatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, 0, input.rows(), col);
            return result;
        }

        template<typename MatrixType>
        constexpr decltype(auto) transpose(MatrixType const &arg) noexcept {
            return arg.transpose();
        }

    }

}

#endif // __TBGAL_USING_EIGEN_HPP__
