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

#ifndef __TBGAL_RIGHT_CONTRACTION_HPP__
#define __TBGAL_RIGHT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoredMultivectorType = decltype(undual(op(arg2, dual(arg1))));
        if (!is_blade(arg1) || !is_blade(arg2) || (arg1.factors_count() >= arg2.factors_count())) {
            if ((arg2.factors_count() * (arg1.factors_count() + 1)) & 1) {
                return -undual(op(arg2, dual(arg1)));
            }
            return undual(op(arg2, dual(arg1)));
        }
        return ResultingFactoredMultivectorType(detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstMetricSpaceType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType>, int> >
    constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, GeometricProduct<FirstMetricSpaceType> > const &arg1, SecondScalarType const &arg2) {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = GeometricProduct<FirstMetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1.scalar() * arg2;
        if (!is_zero(resulting_scalar)) {
            return ResultingFactoredMultivectorType(arg1.space_ptr(), resulting_scalar, arg1.factors_in_signed_metric());
        }
        return ResultingFactoredMultivectorType(arg1.space_ptr(), 0);
    }

    template<typename FirstScalarType, typename FirstMetricSpaceType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType>, int> >
    constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, OuterProduct<FirstMetricSpaceType> > const &arg1, SecondScalarType const &arg2) {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = OuterProduct<FirstMetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1.scalar() * arg2;
        if (!is_zero(resulting_scalar)) {
            return ResultingFactoredMultivectorType(arg1.space_ptr(), resulting_scalar, arg1.factors_and_complement_in_signed_metric(), arg1.factors_count());
        }
        return ResultingFactoredMultivectorType(arg1.space_ptr(), 0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType>, int> >
    constexpr decltype(auto) rcont(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        using ResultingType = std::common_type_t<FirstScalarType, SecondScalarType>;
        return arg2.factors_count() == 0 ? arg1 * arg2.scalar() : ResultingType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>), int> >
    constexpr decltype(auto) rcont(FirstScalarType const &arg1, SecondScalarType const &arg2) {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_RIGHT_CONTRACTION_HPP__
