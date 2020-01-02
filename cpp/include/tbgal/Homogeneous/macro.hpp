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

#ifndef __TBGAL_HOMOGENEOUS_MACRO_HPP__
#define __TBGAL_HOMOGENEOUS_MACRO_HPP__

#define _TBGAL_OVERLOAD_HOMOGENEOUS_UTILS(METRIC_SPACE) \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) direction(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 0); \
    } \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 0); \
    } \
    \
    template<typename... ScalarTypes> \
    constexpr decltype(auto) point(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)..., 1); \
    }

//TODO [NEXT] Inclue other helpers.
//TODO [NEXT] Inclue euclidean_vector with iterators.

#endif // __TBGAL_HOMOGENEOUS_MACRO_HPP__
