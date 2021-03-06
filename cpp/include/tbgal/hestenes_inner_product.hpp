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

#ifndef __TBGAL_HESTENES_INNER_PRODUCT_HPP__
#define __TBGAL_HESTENES_INNER_PRODUCT_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) hip(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        using ResultingFactoredMultivectorType = decltype(lcont(arg1, arg2));
        if (arg1.factors_count() != 0 && arg2.factors_count() != 0) {
            if (arg1.factors_count() <= arg2.factors_count()) {
                return lcont(arg1, arg2);
            }
            return rcont(arg1, arg2);
        }
        return ResultingFactoredMultivectorType(detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType>, int> >
    constexpr decltype(auto) hip(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, SecondScalarType const &) {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType>, int> >
    constexpr decltype(auto) hip(FirstScalarType const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>), int> >
    constexpr decltype(auto) hip(FirstScalarType const &, SecondScalarType const &) {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

}

#endif // __TBGAL_HESTENES_INNER_PRODUCT_HPP__
