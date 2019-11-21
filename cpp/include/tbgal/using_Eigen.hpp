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

        template<typename MatrixType>
        struct index_type {
            using type = Eigen::Index;
        };

        template<typename MatrixType>
        struct scalar_type {
            using type = typename MatrixType::Scalar;
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

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        struct matrix_type {
            using type = Eigen::Matrix<
                    ScalarType,
                    RowsAtCompileTime,
                    ColsAtCompileTime,
                    (MaxRowsAtCompileTime == 1 && MaxColsAtCompileTime != 1) ? Eigen::RowMajor : Eigen::ColMajor, // This is a workaround required by Eigen library. 
                    MaxRowsAtCompileTime,
                    MaxColsAtCompileTime
                >;
        };

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_matrix(DefaultIndexType rows, DefaultIndexType cols) noexcept {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime>(rows, cols);
        }
        
        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(DefaultIndexType size) noexcept {
            return matrix_type_t<ScalarType, SizeAtCompileTime, SizeAtCompileTime, MaxSizeAtCompileTime, MaxSizeAtCompileTime>::Identity(size, size);
        }

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_zero_matrix(DefaultIndexType rows, DefaultIndexType cols) noexcept {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime>::Zero(rows, cols);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &source, DefaultIndexType start_row_source, DefaultIndexType start_col_source, TargetMatrixType &target, DefaultIndexType start_row_target, DefaultIndexType start_col_target, DefaultIndexType block_rows, DefaultIndexType block_cols) noexcept {
            target.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_target, start_col_target, block_rows, block_cols) = source.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_source, start_col_source, block_rows, block_cols);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &source, DefaultIndexType start_row_source, DefaultIndexType start_col_source, TargetMatrixType &target, DefaultIndexType block_rows, DefaultIndexType block_cols) noexcept {
            target = source.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_source, start_col_source, block_rows, block_cols);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceType, typename TargetMatrixType>
        constexpr void assign_block(SourceType const &source, TargetMatrixType &target, DefaultIndexType start_row_target, DefaultIndexType start_col_target, DefaultIndexType block_rows, DefaultIndexType block_cols) noexcept {
            target.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_target, start_col_target, block_rows, block_cols) = source;
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) block(MatrixType const &arg, DefaultIndexType start_row, DefaultIndexType start_col, DefaultIndexType block_rows, DefaultIndexType block_cols) noexcept {
            return arg.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row, start_col, block_rows, block_cols);
        }

        template<typename MatrixType>
        constexpr void conservative_resize(MatrixType &arg, DefaultIndexType rows, DefaultIndexType cols) noexcept {
            arg.conservativeResize(rows, cols);
        }

        template<typename MatrixType>
        constexpr decltype(auto) determinant(MatrixType const &arg) noexcept {
            return arg.determinant();
        }
        
        template<typename MatrixType, int BlockRowsAtCompileTime, int BlockColsAtCompileTime, bool InnerPanel>
        constexpr decltype(auto) evaluate(Eigen::Block<MatrixType, BlockRowsAtCompileTime, BlockColsAtCompileTime, InnerPanel> const &arg) noexcept {
            return arg.eval();
        }

        template<typename ScalarType, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
        constexpr decltype(auto) evaluate(Eigen::Matrix<ScalarType, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> const &arg) noexcept {
            return arg;
        }

        template<typename FirstMatrixType, typename SecondMatrixType, int Option>
        constexpr decltype(auto) evaluate(Eigen::Product<FirstMatrixType, SecondMatrixType, Option> const &arg) noexcept {
            return arg.eval();
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
            matrix_type_t<std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>, sizeof...(ScalarTypes), 1, sizeof...(ScalarTypes), 1> result;
            _fill_column_matrix_impl(result, std::move(args)...);
            return result;
        }

        template<DefaultIndexType RowsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, typename IteratorType>
        constexpr decltype(auto) fill_column_matrix_using_iterators(IteratorType begin, IteratorType end) noexcept {
            assert(RowsAtCompileTime == Eigen::Dynamic || RowsAtCompileTime == std::distance(begin, end));
            matrix_type_t<std::remove_cv_t<std::remove_reference_t<typename std::iterator_traits<IteratorType>::value_type> >, RowsAtCompileTime, 1, MaxRowsAtCompileTime, 1> result(std::distance(begin, end), 1);
            for (Eigen::Index ind = 0; begin != end; ++ind, ++begin) {
                result(ind, 0) = *begin;
            }
            return result;
        }

        template<typename MatrixType>
        constexpr decltype(auto) inverse(MatrixType const &arg) noexcept {
            return arg.inverse();
        }

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod(FirstMatrixType const &arg1, SecondMatrixType const &arg2) noexcept {
            return arg1 * arg2;
        }
        
        template<DefaultIndexType FirstBlockRowsAtCompileTime, DefaultIndexType SecondBlockRowsAtCompileTime, DefaultIndexType SecondBlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &arg1, DefaultIndexType start_row1, DefaultIndexType start_col1, DefaultIndexType block_rows1, SecondMatrixType const &arg2, DefaultIndexType start_row2, DefaultIndexType start_col2, DefaultIndexType block_rows2, DefaultIndexType block_cols2) noexcept {
            return arg1.template block<FirstBlockRowsAtCompileTime, SecondBlockRowsAtCompileTime>(start_row1, start_col1, block_rows1, block_rows2) * arg2.template block<SecondBlockRowsAtCompileTime, SecondBlockColsAtCompileTime>(start_row2, start_col2, block_rows2, block_cols2);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &arg1, DefaultIndexType start_row1, DefaultIndexType start_col1, DefaultIndexType block_rows1, DefaultIndexType block_cols1, SecondMatrixType const &arg2) noexcept {
            return arg1.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row1, start_col1, block_rows1, block_cols1) * arg2;
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &arg1, SecondMatrixType const &arg2, DefaultIndexType start_row2, DefaultIndexType start_col2, DefaultIndexType block_rows2, DefaultIndexType block_cols2) noexcept {
            return arg1 * arg2.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row2, start_col2, block_rows2, block_cols2);
        }

        template<typename MatrixType>
        constexpr decltype(auto) qr_orthogonal_matrix(MatrixType const &arg) noexcept {
            using MatrixQType = matrix_type_t<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime, MatrixType::MaxRowsAtCompileTime, MatrixType::MaxRowsAtCompileTime>;
            Eigen::ColPivHouseholderQR<MatrixType> qr(arg);
            auto rank = qr.rank();
            return std::make_tuple(MatrixQType(qr.householderQ()), rank);
        }

        template<typename MatrixType>
        constexpr decltype(auto) reverse_columns(MatrixType const &arg) noexcept {
            return Eigen::Reverse<MatrixType, Eigen::Horizontal>(arg);
        }

        //TODO Passar para dentro do código de dual e undual.
        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &input, DefaultIndexType col) noexcept {
            using ResultingMatrixType = matrix_type_t<typename MatrixType::Scalar, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime, MatrixType::MaxRowsAtCompileTime, MatrixType::MaxColsAtCompileTime>;
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
