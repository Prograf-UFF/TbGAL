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

#ifndef __TBGAL_UNARY_MINUS_HPP__
#define __TBGAL_UNARY_MINUS_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) unary_minus(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space_ptr(), -arg.scalar(), arg.factors_in_signed_metric());
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) unary_minus(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space_ptr(), -arg.scalar(), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type>, int> >
    constexpr decltype(auto) unary_minus(Type const &arg) {
        return -arg;
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) operator-(FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return unary_minus(arg);
    }

}

#endif // __TBGAL_UNARY_MINUS_HPP__
