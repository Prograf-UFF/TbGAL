/* Copyright (C) Eduardo Vera Sousa and Leandro Augusto Frata Fernandes
 * 
 * authors    : Sousa, Eduardo V.
 *              Fernandes, Leandro A. F.
 * repository : https://github.com/Prograf-UFF/TbGAL
 * 
 * This file is part of the Tensor-based Geometric Algebra Library (TbGAL).
 * 
 * TbGAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * TbGAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with TbGAL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef __TBGAL_USING_EIGEN_HPP__
#define __TBGAL_USING_EIGEN_HPP__

#ifdef TBGAL_USING_MATRIX_DEFINITIONS
    #error "Matrix definition already included by command '#include <using_{SomeMatrixAlgebraLibrary}.hpp>.'"
#else
    #define TBGAL_USING_MATRIX_DEFINITIONS
#endif // TBGAL_USING_MATRIX_DEFINITIONS

#include <Eigen/Dense>

#ifndef TBGAL_DEFAULT_INDEX_TYPE
    #define TBGAL_DEFAULT_INDEX_TYPE Eigen::Index
#endif // TBGAL_DEFAULT_INDEX_TYPE

#include "matrix_declarations.hpp"

namespace tbgal {

    namespace detail {

        template<typename MatrixType>
        struct index_type {
            using type = Eigen::Index;
        };

        template<typename MatrixType>
        struct scalar_type {
            using type = typename std::remove_cv_t<std::remove_reference_t<MatrixType> >::Scalar;
        };

        template<typename MatrixType>
        constexpr decltype(auto) coeff(MatrixType &arg, DefaultIndexType row, DefaultIndexType col) {
            return arg(row, col);
        }

        template<typename ScalarType, int SizeAtCompileTime, int MaxSizeAtCompileTime>
        constexpr decltype(auto) coeff(Eigen::DiagonalMatrix<ScalarType, SizeAtCompileTime, MaxSizeAtCompileTime> &arg, DefaultIndexType row, DefaultIndexType col) {
            assert(row == col);
            return arg.diagonal()[row];
        }

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &arg) {
            return arg.cols();
        }

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &arg) {
            return arg.rows();
        }

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        struct diagonal_matrix_type {
            using type = Eigen::DiagonalMatrix<
                    ScalarType,
                    SizeAtCompileTime,
                    MaxSizeAtCompileTime
                >;
        };

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

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        constexpr decltype(auto) make_diagonal_matrix(DefaultIndexType size) {
            return diagonal_matrix_type_t<ScalarType, SizeAtCompileTime, MaxSizeAtCompileTime>(size);
        }
        
        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(DefaultIndexType size) {
            return matrix_type_t<ScalarType, SizeAtCompileTime, SizeAtCompileTime, MaxSizeAtCompileTime, MaxSizeAtCompileTime>::Identity(size, size);
        }

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_matrix(DefaultIndexType rows, DefaultIndexType cols) {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime>(rows, cols);
        }
        
        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_zero_matrix(DefaultIndexType rows, DefaultIndexType cols) {
            return matrix_type_t<ScalarType, RowsAtCompileTime, ColsAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime>::Zero(rows, cols);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &source, DefaultIndexType start_row_source, DefaultIndexType start_col_source, TargetMatrixType &target, DefaultIndexType start_row_target, DefaultIndexType start_col_target, DefaultIndexType block_rows, DefaultIndexType block_cols) {
            target.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_target, start_col_target, block_rows, block_cols) = source.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_source, start_col_source, block_rows, block_cols);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &source, DefaultIndexType start_row_source, DefaultIndexType start_col_source, TargetMatrixType &target, DefaultIndexType block_rows, DefaultIndexType block_cols) {
            target = source.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_source, start_col_source, block_rows, block_cols);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceType, typename TargetMatrixType>
        constexpr void assign_block(SourceType const &source, TargetMatrixType &target, DefaultIndexType start_row_target, DefaultIndexType start_col_target, DefaultIndexType block_rows, DefaultIndexType block_cols) {
            target.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row_target, start_col_target, block_rows, block_cols) = source;
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) block(MatrixType const &arg, DefaultIndexType start_row, DefaultIndexType start_col, DefaultIndexType block_rows, DefaultIndexType block_cols) {
            return arg.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row, start_col, block_rows, block_cols);
        }

        template<typename MatrixType>
        constexpr void conservative_resize(MatrixType &arg, DefaultIndexType rows, DefaultIndexType cols) {
            arg.conservativeResize(rows, cols);
        }

        template<typename MatrixType>
        constexpr decltype(auto) determinant(MatrixType const &arg) {
            return arg.determinant();
        }
        
        template<typename MatrixType>
        constexpr decltype(auto) es_eigenvectors_matrix(MatrixType const &arg) {
            using EigenSolverType = Eigen::SelfAdjointEigenSolver<MatrixType>;
            using EigenVectorsMatrixType = typename EigenSolverType::EigenvectorsType;
            return EigenVectorsMatrixType(EigenSolverType(arg, Eigen::ComputeEigenvectors).eigenvectors());
        }

        template<typename MatrixType, int BlockRowsAtCompileTime, int BlockColsAtCompileTime, bool InnerPanel>
        constexpr decltype(auto) evaluate(Eigen::Block<MatrixType, BlockRowsAtCompileTime, BlockColsAtCompileTime, InnerPanel> const &arg) {
            return arg.eval();
        }

        template<typename ScalarType, int RowsAtCompileTime, int ColsAtCompileTime, int Options, int MaxRowsAtCompileTime, int MaxColsAtCompileTime>
        constexpr decltype(auto) evaluate(Eigen::Matrix<ScalarType, RowsAtCompileTime, ColsAtCompileTime, Options, MaxRowsAtCompileTime, MaxColsAtCompileTime> const &arg) {
            return arg;
        }

        template<typename FirstMatrixType, typename SecondMatrixType, int Option>
        constexpr decltype(auto) evaluate(Eigen::Product<FirstMatrixType, SecondMatrixType, Option> const &arg) {
            return arg.eval();
        }

        template<typename MatrixType>
        constexpr void _fill_column_matrix_impl(MatrixType const &, Eigen::Index const) {
        }
        
        template<typename MatrixType, typename ScalarType>
        constexpr void _fill_column_matrix_impl(MatrixType &target, Eigen::Index offset, ScalarType &&arg) {
            target(offset, 0) = std::move(arg);
        }
        
        template<typename MatrixType, typename FirstScalarType, typename... NextScalarTypes>
        constexpr void _fill_column_matrix_impl(MatrixType &target, Eigen::Index offset, FirstScalarType &&arg1, NextScalarTypes &&... args) {
            target(offset, 0) = std::move(arg1);
            _fill_column_matrix_impl(target, offset + 1, std::move(args)...);
        }
        
        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&... args) {
            matrix_type_t<std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>, sizeof...(ScalarTypes), 1, sizeof...(ScalarTypes), 1> result;
            _fill_column_matrix_impl(result, 0, std::move(args)...);
            return result;
        }

        template<DefaultIndexType RowsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, typename IteratorType, typename... ExtraScalarTypes>
        constexpr decltype(auto) fill_column_matrix_using_iterator(IteratorType begin, IteratorType end, ExtraScalarTypes &&...extra_args) {
            matrix_type_t<std::remove_cv_t<std::remove_reference_t<typename std::iterator_traits<IteratorType>::value_type> >, RowsAtCompileTime, 1, MaxRowsAtCompileTime, 1> result(std::distance(begin, end) + sizeof...(extra_args), 1);
            Eigen::Index ind = 0;
            for (; begin != end; ++ind, ++begin) {
                result(ind, 0) = *begin;
            }
            _fill_column_matrix_impl(result, ind, std::move(extra_args)...);
            return result;
        }

        template<typename MatrixType>
        constexpr decltype(auto) inverse(MatrixType const &arg) {
            return arg.inverse();
        }

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod(FirstMatrixType const &arg1, SecondMatrixType const &arg2) {
            return arg1 * arg2;
        }
        
        template<DefaultIndexType FirstBlockRowsAtCompileTime, DefaultIndexType SecondBlockRowsAtCompileTime, DefaultIndexType SecondBlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &arg1, DefaultIndexType start_row1, DefaultIndexType start_col1, DefaultIndexType block_rows1, SecondMatrixType const &arg2, DefaultIndexType start_row2, DefaultIndexType start_col2, DefaultIndexType block_rows2, DefaultIndexType block_cols2) {
            return arg1.template block<FirstBlockRowsAtCompileTime, SecondBlockRowsAtCompileTime>(start_row1, start_col1, block_rows1, block_rows2) * arg2.template block<SecondBlockRowsAtCompileTime, SecondBlockColsAtCompileTime>(start_row2, start_col2, block_rows2, block_cols2);
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &arg1, DefaultIndexType start_row1, DefaultIndexType start_col1, DefaultIndexType block_rows1, DefaultIndexType block_cols1, SecondMatrixType const &arg2) {
            return arg1.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row1, start_col1, block_rows1, block_cols1) * arg2;
        }

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &arg1, SecondMatrixType const &arg2, DefaultIndexType start_row2, DefaultIndexType start_col2, DefaultIndexType block_rows2, DefaultIndexType block_cols2) {
            return arg1 * arg2.template block<BlockRowsAtCompileTime, BlockColsAtCompileTime>(start_row2, start_col2, block_rows2, block_cols2);
        }

        template<typename MatrixType>
        constexpr decltype(auto) qr_orthogonal_matrix(MatrixType const &arg) {
            using MatrixQType = matrix_type_t<scalar_type_t<MatrixType>, MatrixType::RowsAtCompileTime, MatrixType::RowsAtCompileTime, MatrixType::MaxRowsAtCompileTime, MatrixType::MaxRowsAtCompileTime>;
            Eigen::ColPivHouseholderQR<MatrixType> qr(arg);
            auto rank = qr.rank();
            return std::make_tuple(MatrixQType(qr.householderQ()), rank);
        }

        template<typename MatrixType>
        constexpr decltype(auto) reverse_columns(MatrixType const &arg) {
            return Eigen::Reverse<MatrixType, Eigen::Horizontal>(arg);
        }

        //TODO Passar para dentro do código de dual e undual.
        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &input, DefaultIndexType col) {
            using ResultingMatrixType = matrix_type_t<scalar_type_t<MatrixType>, MatrixType::RowsAtCompileTime, MatrixType::ColsAtCompileTime, MatrixType::MaxRowsAtCompileTime, MatrixType::MaxColsAtCompileTime>;
            ResultingMatrixType result(input.rows(), input.cols());
            result.template block<ResultingMatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, 0, result.rows(), input.cols() - col) = input.template block<MatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, col, input.rows(), input.cols() - col);
            result.template block<ResultingMatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, input.cols() - col, result.rows(), col) = input.template block<MatrixType::RowsAtCompileTime, Eigen::Dynamic>(0, 0, input.rows(), col);
            return result;
        }

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) subtract(FirstMatrixType const &arg1, SecondMatrixType const &arg2) {
            return arg1 - arg2;
        }

        template<typename MatrixType>
        constexpr decltype(auto) transpose(MatrixType const &arg) {
            return arg.transpose();
        }

    }

}

#endif // __TBGAL_USING_EIGEN_HPP__
