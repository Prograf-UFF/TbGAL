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

#ifndef __TBGAL_SCALAR_PRODUCT_HPP__
#define __TBGAL_SCALAR_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement scalar product for the general case (when the input multivector is not a blade).

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) sp(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        static_assert(std::is_same_v<typename FirstFactoringProductType::MetricSpaceType, typename SecondFactoringProductType::MetricSpaceType>, "The multivectors must have the same metric space.");
        using MetricSpaceType = typename FirstFactoringProductType::MetricSpaceType;
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType, typename MetricSpaceType::ScalarType>;
        if (!is_blade(arg1) || !is_blade(arg2)) {
            throw NotSupportedError("The scalar product is only implemented for blades.");
        }
        if (arg1.factors_count() == arg2.factors_count()) {
            return (((arg1.factors_count() * (arg1.factors_count() - 1)) & 2) ? -arg1.scalar() : arg1.scalar()) * arg2.scalar() * detail::determinant(detail::prod(detail::transpose(arg1.factors_in_signed_metric()), detail::apply_signed_metric(arg2.space_ptr(), arg2.factors_in_signed_metric())));
        }
        return ResultingScalarType(0);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType>, int> >
    constexpr decltype(auto) sp(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        if (!is_blade(arg1)) {
            throw NotSupportedError("The scalar product is only implemented for blades.");
        }
        if (arg1.factors_count() == 0) {
            return arg1.scalar() * arg2;
        }
        return ResultingScalarType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType>, int> >
    constexpr decltype(auto) sp(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        if (!is_blade(arg2)) {
            throw NotSupportedError("The scalar product is only implemented for blades.");
        }
        if (arg2.factors_count() == 0) {
            return arg1 * arg2.scalar();
        }
        return ResultingScalarType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>), int> >
    constexpr decltype(auto) sp(FirstScalarType const &arg1, SecondScalarType const &arg2) {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_SCALAR_PRODUCT_HPP__
