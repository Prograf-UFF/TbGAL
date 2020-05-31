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

#ifndef __TBGAL_APPLY_VERSOR_HPP__
#define __TBGAL_APPLY_VERSOR_HPP__

namespace tbgal {

    template<typename VersorScalarType, typename VersorFactoringProductType, typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) apply_even_versor(FactoredMultivector<VersorScalarType, VersorFactoringProductType> const &versor, FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return gp(versor, arg, inverse(versor));
    }

    template<typename VersorScalarType, typename VersorFactoringProductType, typename Type, typename = std::enable_if_t<!is_multivector_v<Type>, int> >
    constexpr Type apply_even_versor(FactoredMultivector<VersorScalarType, VersorFactoringProductType> const &, Type const &arg) {
        return arg;
    }

    template<typename VersorType, typename ScalarType, typename FactoringProductType, typename = std::enable_if_t<!is_multivector_v<VersorType>, int> >
    constexpr FactoredMultivector<ScalarType, FactoringProductType> apply_even_versor(VersorType const &, FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return arg;
    }

    template<typename VersorType, typename Type, typename = std::enable_if_t<!(is_multivector_v<VersorType> || is_multivector_v<Type>), int> >
    constexpr Type apply_even_versor(VersorType const &, Type const &arg) {
        return arg;
    }

    template<typename VersorScalarType, typename VersorFactoringProductType, typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) apply_odd_versor(FactoredMultivector<VersorScalarType, VersorFactoringProductType> const &versor, FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return gp(versor, involute(arg), inverse(versor));
    }

    template<typename VersorScalarType, typename VersorFactoringProductType, typename Type, typename = std::enable_if_t<!is_multivector_v<Type>, int> >
    constexpr Type apply_odd_versor(FactoredMultivector<VersorScalarType, VersorFactoringProductType> const &, Type const &arg) {
        return arg;
    }

    template<typename VersorType, typename ScalarType, typename FactoringProductType, typename = std::enable_if_t<!is_multivector_v<VersorType>, int> >
    constexpr FactoredMultivector<ScalarType, FactoringProductType> apply_odd_versor(VersorType const &, FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return arg;
    }

    template<typename VersorType, typename Type, typename = std::enable_if_t<!(is_multivector_v<VersorType> || is_multivector_v<Type>), int> >
    constexpr Type apply_odd_versor(VersorType const &, Type const &arg) {
        return arg;
    }

    template<typename RotorScalarType, typename RotorFactoringProductType, typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) apply_rotor(FactoredMultivector<RotorScalarType, RotorFactoringProductType> const &rotor, FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return gp(rotor, arg, reverse(rotor));
    }

    template<typename RotorScalarType, typename RotorFactoringProductType, typename Type, typename = std::enable_if_t<!is_multivector_v<Type>, int> >
    constexpr Type apply_rotor(FactoredMultivector<RotorScalarType, RotorFactoringProductType> const &, Type const &arg) {
        return arg;
    }

    template<typename RotorType, typename ScalarType, typename FactoringProductType, typename = std::enable_if_t<!is_multivector_v<RotorType>, int> >
    constexpr FactoredMultivector<ScalarType, FactoringProductType> apply_rotor(RotorType const &, FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return arg;
    }

    template<typename RotorType, typename Type, typename = std::enable_if_t<!(is_multivector_v<RotorType> || is_multivector_v<Type>), int> >
    constexpr Type apply_rotor(RotorType const &, Type const &arg) {
        return arg;
    }

}

#endif // __TBGAL_APPLY_VERSOR_HPP__
