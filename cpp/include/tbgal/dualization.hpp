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

#ifndef __TBGAL_DUALIZATION_HPP__
#define __TBGAL_DUALIZATION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) dual(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg);
    //TODO [FUTURE] Implement dualization for the general case (when the input multivector is not a blade).

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) dual(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(
            arg.space_ptr(),
            (((arg.factors_count() * (arg.factors_count() - 1) + arg.space().dimensions() * (arg.space().dimensions() - 1)) & 2) ? -arg.scalar() : arg.scalar()) * detail::determinant(arg.factors_and_complement_in_signed_metric()),
            detail::apply_signed_metric(arg.space_ptr(), detail::split_columns_and_swap(arg.factors_and_complement_in_signed_metric(), arg.factors_count())),
            arg.space().dimensions() - arg.factors_count()
        );
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) undual(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg);
    //TODO [FUTURE] Implement undualization for the general case (when the input multivector is not a blade).

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) undual(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto factors_and_complement_in_signed_metric = detail::evaluate(detail::apply_signed_metric(arg.space_ptr(), arg.factors_and_complement_in_signed_metric()));
        return ResultingFactoredMultivectorType(
            arg.space_ptr(),
            (((arg.factors_count() * (arg.factors_count() - 1)) & 2) ? -arg.scalar() : arg.scalar()) * detail::determinant(factors_and_complement_in_signed_metric),
            detail::split_columns_and_swap(factors_and_complement_in_signed_metric, arg.factors_count()),
            arg.space().dimensions() - arg.factors_count()
        );
    }

}

#endif // __TBGAL_DUALIZATION_HPP__
