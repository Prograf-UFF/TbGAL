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

    //TODO [Parei aqui]
    /*
    template<typename VersorCoefficientType, typename VersorExpression, typename CoefficientType, typename Expression>
    constexpr decltype(auto) apply_even_versor(clifford_expression<VersorCoefficientType, VersorExpression> const &versor, clifford_expression<CoefficientType, Expression> const &arg) {
        auto const lazy = make_lazy_context(versor, arg);
        return lazy.eval(gp(gp(lazy.template argument<0>(), lazy.template argument<1>(), mtr), inv(lazy.template argument<0>(), mtr), mtr));
    }

    template<typename VersorCoefficientType, typename VersorExpression, typename Type>
    constexpr decltype(auto) apply_even_versor(clifford_expression<VersorCoefficientType, VersorExpression> const &versor, Type const &arg) {
        return apply_even_versor(versor, scalar(arg), mtr);
    }

    template<typename VersorType, typename CoefficientType, typename Expression>
    constexpr decltype(auto) apply_even_versor(VersorType const &versor, clifford_expression<CoefficientType, Expression> const &arg) {
        return apply_even_versor(scalar(versor), arg, mtr);
    }

    template<typename VersorType, typename Type>
    constexpr decltype(auto) apply_even_versor(VersorType const &versor, Type const &arg) {
        return apply_even_versor(scalar(versor), scalar(arg), mtr);
    }

    template<typename VersorCoefficientType, typename VersorExpression, typename CoefficientType, typename Expression>
    constexpr decltype(auto) apply_odd_versor(clifford_expression<VersorCoefficientType, VersorExpression> const &versor, clifford_expression<CoefficientType, Expression> const &arg) {
        auto const lazy = make_lazy_context(versor, arg);
        return lazy.eval(gp(gp(lazy.template argument<0>(), involution(lazy.template argument<1>()), mtr), inv(lazy.template argument<0>(), mtr), mtr));
    }

    template<typename VersorCoefficientType, typename VersorExpression, typename Type>
    constexpr decltype(auto) apply_odd_versor(clifford_expression<VersorCoefficientType, VersorExpression> const &versor, Type const &arg) {
        return apply_odd_versor(versor, scalar(arg), mtr);
    }

    template<typename VersorType, typename CoefficientType, typename Expression>
    constexpr decltype(auto) apply_odd_versor(VersorType const &versor, clifford_expression<CoefficientType, Expression> const &arg) {
        return apply_odd_versor(scalar(versor), arg, mtr);
    }

    template<typename VersorType, typename Type>
    constexpr decltype(auto) apply_odd_versor(VersorType const &versor, Type const &arg) {
        return apply_odd_versor(scalar(versor), scalar(arg), mtr);
    }

    template<typename RotorCoefficientType, typename RotorExpression, typename CoefficientType, typename Expression>
    constexpr decltype(auto) apply_rotor(clifford_expression<RotorCoefficientType, RotorExpression> const &rotor, clifford_expression<CoefficientType, Expression> const &arg) {
        auto const lazy = make_lazy_context(rotor, arg);
        return lazy.eval(gp(gp(lazy.template argument<0>(), lazy.template argument<1>(), mtr), reversion(lazy.template argument<0>()), mtr));
    }

    template<typename RotorCoefficientType, typename RotorExpression, typename Type>
    constexpr decltype(auto) apply_rotor(clifford_expression<RotorCoefficientType, RotorExpression> const &rotor, Type const &arg) {
        return apply_rotor(rotor, scalar(arg), mtr);
    }

    template<typename RotorType, typename CoefficientType, typename Expression>
    constexpr decltype(auto) apply_rotor(RotorType const &rotor, clifford_expression<CoefficientType, Expression> const &arg) {
        return apply_rotor(scalar(rotor), arg, mtr);
    }

    template<typename RotorType, typename Type>
    constexpr decltype(auto) apply_rotor(RotorType const &rotor, Type const &arg) {
        return apply_rotor(scalar(rotor), scalar(arg), mtr);
    }
    */

}

#endif // __TBGAL_APPLY_VERSOR_HPP__
