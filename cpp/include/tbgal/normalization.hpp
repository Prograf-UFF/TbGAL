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

#ifndef __TBGAL_NORMALIZATION_HPP__
#define __TBGAL_NORMALIZATION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) unit(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = std::common_type_t<ScalarType, typename MetricSpaceType::ScalarType>;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        if (arg.factors_count() > 0) {
            return ResultingFactoredMultivectorType(arg.space_ptr(), arg.scalar() / sqrt(abs(arg.scalar() * arg.scalar() * detail::metric_factor(arg.space_ptr(), arg.factors_in_signed_metric()))), arg.factors_in_signed_metric());
        }
        else {
            return ResultingFactoredMultivectorType(arg.space_ptr(), arg.scalar() / abs(arg.scalar()), arg.factors_in_signed_metric());
        }
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) unit(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = std::common_type_t<ScalarType, typename MetricSpaceType::ScalarType>;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        if (arg.factors_count() > 0) {
            auto factors_tuple = detail::from_outer_to_geometric_factors(arg.space_ptr(), arg.factors_in_signed_metric());
            return ResultingFactoredMultivectorType(arg.space_ptr(), arg.scalar() / sqrt(abs(arg.scalar() * arg.scalar() * std::get<0>(factors_tuple) * std::get<0>(factors_tuple) * detail::metric_factor(arg.space_ptr(), std::get<1>(factors_tuple)))), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
        }
        else {
            return ResultingFactoredMultivectorType(arg.space_ptr(), arg.scalar() / abs(arg.scalar()), arg.factors_and_complement_in_signed_metric(), 0);
        }
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type unit(Type const &arg) noexcept {
        return arg / abs(arg);
    }

}

#endif // __TBGAL_NORMALIZATION_HPP__
