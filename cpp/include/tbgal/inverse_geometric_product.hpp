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

#ifndef __TBGAL_INVERSE_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_INVERSE_GEOMETRIC_PRODUCT_HPP__

namespace tbgal {

    template<typename FirstType, typename SecondType>
    constexpr decltype(auto) igp(FirstType const &arg1, SecondType const &arg2) noexcept {
        return gp(arg1, inverse(arg2));
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator/(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return igp(arg1, arg2);
    }

    template<typename FirstScsalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator/(FactoredMultivector<FirstScsalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return igp(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator/(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return igp(arg1, arg2);
    }

}

#endif // __TBGAL_INVERSE_GEOMETRIC_PRODUCT_HPP__
