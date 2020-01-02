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

#ifndef __TBGAL_MACRO_HPP__
#define __TBGAL_MACRO_HPP__

#define _TBGAL_OVERLOAD_UTILS(METRIC_SPACE) \
    \
    decltype(auto) e(tbgal::DefaultIndexType index) noexcept { \
        return tbgal::e<tbgal::DefaultScalarType>(METRIC_SPACE, index); \
    } \
    \
    template<typename ScalarType> \
    decltype(auto) scalar(ScalarType &&scalar) noexcept { \
        return tbgal::scalar(METRIC_SPACE, std::move(scalar)); \
    } \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > > \
    decltype(auto) vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::vector(METRIC_SPACE, std::move(coords)...); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType> > > \
    decltype(auto) vector(IteratorType begin, IteratorType end) noexcept { \
        return tbgal::vector(METRIC_SPACE, begin, end); \
    }

#endif // __TBGAL_MACRO_HPP__
