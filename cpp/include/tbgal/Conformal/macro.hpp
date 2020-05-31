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

#ifndef __TBGAL_CONFORMAL_MACRO_HPP__
#define __TBGAL_CONFORMAL_MACRO_HPP__

#define TBGAL_OVERLOAD_CONFORMAL_UTILS(METRIC_SPACE) \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...>, int> > \
    constexpr decltype(auto) euclidean_vector(ScalarTypes &&... coords) { \
        return tbgal::detail::make_vector(&METRIC_SPACE, std::move(coords)..., 0, 0); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType>, int> > \
    constexpr decltype(auto) euclidean_vector(IteratorType begin, IteratorType end) { \
        return tbgal::detail::make_vector_using_iterator(&METRIC_SPACE, begin, end, 0, 0); \
    } \
    \
    template<typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...>, int> > \
    constexpr decltype(auto) point(ScalarTypes &&... coords) { \
        return tbgal::detail::make_vector(&METRIC_SPACE, std::move(coords)..., 1, ((std::move(coords) * std::move(coords)) + ... + 0) / 2); \
    } \
    \
    template<typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType>, int> > \
    constexpr decltype(auto) point(IteratorType begin, IteratorType end) { \
        std::remove_cv_t<std::remove_reference_t<typename std::iterator_traits<IteratorType>::value_type> > aux = 0; \
        for (IteratorType itr = begin; itr != end; ++itr) { \
            aux += (*itr) * (*itr); \
        } \
        return tbgal::detail::make_vector_using_iterator(&METRIC_SPACE, begin, end, 1, aux / 2); \
    }

//TODO [NEXT] Inclue other helpers.

#endif // __TBGAL_CONFORMAL_MACRO_HPP__
