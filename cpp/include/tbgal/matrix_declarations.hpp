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

#ifndef __TBGAL_MATRIX_DECLARATIONS_HPP__
#define __TBGAL_MATRIX_DECLARATIONS_HPP__

#include <cmath>
#include <cstdint>
#include <iostream>
#include <iterator>
#include <type_traits>
#include <tuple>
#include <utility>

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
        constexpr decltype(auto) coeff(MatrixType &, DefaultIndexType, DefaultIndexType);

        template<typename MatrixType>
        constexpr decltype(auto) cols(MatrixType const &);

        template<typename MatrixType>
        constexpr decltype(auto) rows(MatrixType const &);

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        struct diagonal_matrix_type;

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        using diagonal_matrix_type_t = typename diagonal_matrix_type<ScalarType, SizeAtCompileTime, MaxSizeAtCompileTime>::type;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        struct matrix_type;

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        using matrix_type_t = typename matrix_type<ScalarType, RowsAtCompileTime, ColsAtCompileTime, MaxRowsAtCompileTime, MaxColsAtCompileTime>::type;

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        constexpr decltype(auto) make_diagonal_matrix(DefaultIndexType);

        template<typename ScalarType, DefaultIndexType SizeAtCompileTime, DefaultIndexType MaxSizeAtCompileTime>
        constexpr decltype(auto) make_identity_matrix(DefaultIndexType);

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_matrix(DefaultIndexType, DefaultIndexType);

        template<typename ScalarType, DefaultIndexType RowsAtCompileTime, DefaultIndexType ColsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, DefaultIndexType MaxColsAtCompileTime>
        constexpr decltype(auto) make_zero_matrix(DefaultIndexType, DefaultIndexType);

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &, DefaultIndexType, DefaultIndexType, TargetMatrixType &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType);

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceMatrixType, typename TargetMatrixType>
        constexpr void assign_block(SourceMatrixType const &, DefaultIndexType, DefaultIndexType, TargetMatrixType &, DefaultIndexType, DefaultIndexType);

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename SourceType, typename TargetMatrixType>
        constexpr void assign_block(SourceType const &, TargetMatrixType &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType);

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) block(MatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType);

        template<typename MatrixType>
        constexpr void conservative_resize(MatrixType &, DefaultIndexType, DefaultIndexType);

        template<typename MatrixType>
        constexpr decltype(auto) determinant(MatrixType const &);
        
        template<typename MatrixType>
        constexpr decltype(auto) es_eigenvectors_matrix(MatrixType const &arg);
        
        template<typename MatrixType>
        constexpr decltype(auto) evaluate(MatrixType const &arg);

        template<typename... ScalarTypes>
        constexpr decltype(auto) fill_column_matrix(ScalarTypes &&...);

        template<DefaultIndexType RowsAtCompileTime, DefaultIndexType MaxRowsAtCompileTime, typename IteratorType, typename... ExtraScalarTypes>
        constexpr decltype(auto) fill_column_matrix_using_iterator(IteratorType, IteratorType, ExtraScalarTypes &&...);

        template<typename MatrixType>
        constexpr decltype(auto) inverse(MatrixType const &);

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod(FirstMatrixType const &, SecondMatrixType const &);

        template<DefaultIndexType FirstBlockRowsAtCompileTime, DefaultIndexType SecondBlockRowsAtCompileTime, DefaultIndexType SecondBlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, SecondMatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType);

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType, SecondMatrixType const &);

        template<DefaultIndexType BlockRowsAtCompileTime, DefaultIndexType BlockColsAtCompileTime, typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) prod_block(FirstMatrixType const &, SecondMatrixType const &, DefaultIndexType, DefaultIndexType, DefaultIndexType, DefaultIndexType);

        template<typename MatrixType>
        constexpr decltype(auto) qr_orthogonal_matrix(MatrixType const &);

        template<typename MatrixType>
        constexpr decltype(auto) reverse_columns(MatrixType const &);

        template<typename MatrixType>
        constexpr decltype(auto) split_columns_and_swap(MatrixType const &, DefaultIndexType);

        template<typename FirstMatrixType, typename SecondMatrixType>
        constexpr decltype(auto) subtract(FirstMatrixType const &arg1, SecondMatrixType const &arg2);

        template<typename MatrixType>
        constexpr decltype(auto) transpose(MatrixType const &);

    }

}

#endif // __TBGAL_MATRIX_DECLARATIONS_HPP__
