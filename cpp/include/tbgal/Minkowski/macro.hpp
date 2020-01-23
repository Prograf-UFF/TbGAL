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

#ifndef __TBGAL_MINKOWSKI_MACRO_HPP__
#define __TBGAL_MINKOWSKI_MACRO_HPP__

#define _TBGAL_OVERLOAD_MINKOWSKI_UTILS(METRIC_SPACE) \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > > \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) noexcept { \
        return tbgal::detail::make_vector(&METRIC_SPACE, std::move(coords)..., 0, 0); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType> > > \
    constexpr decltype(auto) euclidean_vector(IteratorType begin, IteratorType end) noexcept { \
        return tbgal::detail::make_vector_using_iterator(&METRIC_SPACE, begin, end, 0, 0); \
    } \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > > \
    constexpr decltype(auto) point(ScalarTypes &&... coords) noexcept { \
        auto aux = ((std::move(coords) * std::move(coords)) + ... + 0); \
        return tbgal::detail::make_vector(&METRIC_SPACE, std::move(coords)..., (aux - 1) / 2, (aux + 1) / 2); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType> > > \
    constexpr decltype(auto) point(IteratorType begin, IteratorType end) noexcept { \
        std::remove_cv_t<std::remove_reference_t<typename std::iterator_traits<IteratorType>::value_type> > aux = 0; \
        for (IteratorType itr = begin; itr != end; ++itr) { \
            aux += (*itr) * (*itr); \
        } \
        return tbgal::detail::make_vector_using_iterator(&METRIC_SPACE, begin, end, (aux - 1) / 2, (aux + 1) / 2); \
    }

//TODO [NEXT] Inclue other helpers.

#endif // __TBGAL_MINKOWSKI_MACRO_HPP__
