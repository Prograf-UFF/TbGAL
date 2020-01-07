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

#ifndef __TBGAL_CORE_HPP__
#define __TBGAL_CORE_HPP__

#ifndef TBGAL_USING_MATRIX_DEFINITIONS
    #error "Before assuming some model of geometry you have to use some matrix definition by calling the command '#include <using_{SomeMatrixAlgebraLibrary}.hpp>.'"
#endif // TBGAL_USING_MATRIX_DEFINITIONS

#include <cassert>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace tbgal {

    template<typename ScalarType, typename FactoringProductType> class FactoredMultivector;

    template<typename MetricSpaceType> struct GeometricProduct;
    template<typename MetricSpaceType> struct OuterProduct;

    namespace detail {

        template<typename T, typename... Rest>
        constexpr bool is_any_v = std::disjunction_v<std::bool_constant<std::is_same_v<T, Rest> >...>;

        template <typename T, typename = void>
        struct is_iterator : std::false_type {
        };

        template <class T>
        struct is_iterator<T, std::void_t<typename std::iterator_traits<T>::iterator_category> > : std::true_type {
        };

        template<typename T>
        constexpr bool is_iterator_v = is_iterator<T>::value;

        struct gp_impl;
        struct op_impl;

        template<typename MetricSpaceType> struct from_actual_to_signed_metric_impl;
        template<typename MetricSpaceType> struct from_signed_to_actual_metric_impl;

        template<typename FirstType, typename... NextTypes>
        struct common_scalar_type;

        template<typename... Types>
        using common_scalar_type_t = typename common_scalar_type<Types...>::type;

        template<typename FirstScalarType, typename... NextTypes>
        struct common_scalar_type {
            using type = std::common_type_t<FirstScalarType, common_scalar_type_t<NextTypes...> >;
        };

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        struct common_scalar_type<FactoredMultivector<FirstScalarType, FirstFactoringProductType>, NextTypes...> {
            using type = std::common_type_t<FirstScalarType, common_scalar_type_t<NextTypes...> >;
        };

        template<typename ScalarType>
        struct common_scalar_type<ScalarType> {
            using type = ScalarType;
        };

        template<typename ScalarType, typename FactoringProductType>
        struct common_scalar_type<FactoredMultivector<ScalarType, FactoringProductType> > {
            using type = ScalarType;
        };

        template<typename FirstType, typename... NextTypes>
        struct metric_space_type;

        template<typename... Types>
        using metric_space_type_t = typename metric_space_type<Types...>::type;

        template<typename FirstType, typename... NextTypes>
        struct metric_space_type :
            metric_space_type<NextTypes...> {
        };

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        struct metric_space_type<FactoredMultivector<FirstScalarType, FirstFactoringProductType>, NextTypes...> {
        private:

            using NextSpaceType = metric_space_type_t<NextTypes...>;
            static_assert(std::is_same_v<NextSpaceType, std::nullptr_t> || std::is_same_v<NextSpaceType, typename FirstFactoringProductType::MetricSpaceType>, "The multivectors must have the same metric space.");

        public:

            using type = typename FirstFactoringProductType::MetricSpaceType;
        };

        template<typename Type>
        struct metric_space_type<Type> {
            using type = std::nullptr_t;
        };

        template<typename ScalarType, typename FactoringProductType>
        struct metric_space_type<FactoredMultivector<ScalarType, FactoringProductType> > {
            using type = typename FactoringProductType::MetricSpaceType;
        };

        constexpr decltype(auto) space_ptr() noexcept {
            return nullptr;
        }

        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) space_ptr(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return &arg.space();
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr decltype(auto) space_ptr(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            assert(space_ptr(args...) == &arg1.space() || space_ptr(args...) == nullptr);
            return &arg1.space();
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr decltype(auto) space_ptr(FirstScalarType const &, NextTypes const &... args) noexcept {
            return space_ptr(args...);
        }

        template<typename MetricSpaceType, typename... ScalarTypes>
        constexpr decltype(auto) make_vector(MetricSpaceType const &, ScalarTypes &&...) noexcept;

        template<typename MetricSpaceType, typename IteratorType, typename... ExtraScalarTypes>
        constexpr decltype(auto) make_vector_using_iterator(MetricSpaceType const &, IteratorType, IteratorType, ExtraScalarTypes &&...) noexcept;

    }

    template<typename Type>
    struct is_multivector :
        std::false_type {
    };

    template<typename Type>
    using is_multivector_t = typename is_multivector<Type>::type;

    template<typename Type>
    constexpr bool is_multivector_v = is_multivector<Type>::value;

}

#include "exception.hpp"

#include "metric_space.hpp"

#include "factoring_product.hpp"
#include "factored_multivector.hpp"

#include "utils.hpp"

#include "unary_plus.hpp"
#include "unary_minus.hpp"

#include "addition.hpp"
#include "subtraction.hpp"

#include "outer_product.hpp"
#include "geometric_product.hpp"
#include "left_contraction.hpp"
#include "right_contraction.hpp"
#include "dot_product.hpp"
#include "hestenes_inner_product.hpp"
#include "scalar_product.hpp"

#include "conjugation.hpp"
#include "reversion.hpp"
#include "involution.hpp"

#include "reverse_norm.hpp"
#include "inversion.hpp"
#include "dualization.hpp"

#include "apply_versor.hpp"
#include "inverse_geometric_product.hpp"
#include "normalization.hpp"

#include "write.hpp"
#include "macro.hpp"

#include "Conformal/metric_space.hpp"
#include "Conformal/macro.hpp"

#include "Euclidean/metric_space.hpp"
#include "Euclidean/macro.hpp"

#include "Homogeneous/metric_space.hpp"
#include "Homogeneous/macro.hpp"

#include "Minkowski/metric_space.hpp"
#include "Minkowski/macro.hpp"

#include "Signed/metric_space.hpp"
#include "Signed/macro.hpp"

#endif // __TBGAL_CORE_HPP__
