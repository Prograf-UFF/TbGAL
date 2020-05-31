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

#ifndef __TBGAL_REVERSE_NORM_HPP__
#define __TBGAL_REVERSE_NORM_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) rnorm_sqr(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) {
        if (arg.factors_count() > 0) {
            return arg.scalar() * arg.scalar() * detail::metric_factor(arg.space_ptr(), arg.factors_in_signed_metric());
        }
        else {
            return arg.scalar() * arg.scalar();
        }
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) rnorm_sqr(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) {
        if (arg.factors_count() > 0) {
            auto factors_tuple = detail::from_outer_to_geometric_factors(arg.space_ptr(), arg.factors_in_signed_metric());
            return arg.scalar() * arg.scalar() * std::get<0>(factors_tuple) * std::get<0>(factors_tuple) * detail::metric_factor(arg.space_ptr(), std::get<1>(factors_tuple));
        }
        else {
            return arg.scalar() * arg.scalar();
        }
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type>, int> >
    constexpr decltype(auto) rnorm_sqr(Type const &arg) {
        return arg * arg;
    }

    template<typename Type>
    constexpr decltype(auto) rnorm(Type const &arg) {
        return sqrt(rnorm_sqr(arg));
    }

}

#endif // __TBGAL_REVERSE_NORM_HPP__
