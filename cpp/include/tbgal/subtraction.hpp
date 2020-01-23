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

#ifndef __TBGAL_SUBTRACTION_HPP__
#define __TBGAL_SUBTRACTION_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement associativity.

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) subtraction(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = OuterProduct<typename FirstFactoringProductType::MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        static_assert(std::is_same_v<typename FirstFactoringProductType::MetricSpaceType, typename SecondFactoringProductType::MetricSpaceType>, "The multivectors must have the same metric space.");
        assert(arg1.space_ptr() == arg2.space_ptr());
        if (arg1.factors_count() == arg2.factors_count()) {
            if (arg1.factors_count() == 0) {
                return ResultingFactoredMultivectorType(arg1.space_ptr(), arg1.scalar() - arg2.scalar());
            }
            //TODO [FUTURE] Handle vectors.
            //TODO [FUTURE] Handle pseudovetors.
            //TODO [FUTURE] Handle pseudoscalars.
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) subtraction(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        if (arg1.factors_count() == 0) {
            return scalar(arg1.space_ptr(), arg1.scalar() - arg2);
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) subtraction(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        if (arg2.factors_count() == 0) {
            return scalar(arg2.space_ptr(), arg1 - arg2.scalar());
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) subtraction(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 - arg2;
    }

    template<typename FirstType, typename SecondType>
    constexpr decltype(auto) sub(FirstType const &arg1, SecondType const &arg2) noexcept {
        return subtraction(arg1, arg2);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator-(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return subtraction(arg1, arg2);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator-(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return subtraction(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator-(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return subtraction(arg1, arg2);
    }

}

#endif // __TBGAL_SUBTRACTION_HPP__
