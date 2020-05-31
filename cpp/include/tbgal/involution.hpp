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

#ifndef __TBGAL_INVOLUTION_HPP__
#define __TBGAL_INVOLUTION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) involute(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space_ptr(), (arg.factors_count() & 1) ? -arg.scalar() : arg.scalar(), arg.factors_in_signed_metric());
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) involute(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space_ptr(), (arg.factors_count() & 1) ? -arg.scalar() : arg.scalar(), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type>, int> >
    constexpr Type involute(Type const &arg) {
        return arg;
    }

}

#endif // __TBGAL_INVOLUTION_HPP__
