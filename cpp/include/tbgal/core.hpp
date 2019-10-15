#ifndef __TBGAL_CORE_HPP__
#define __TBGAL_CORE_HPP__

#ifndef TBGAL_USING_MATRIX_DEFINITIONS
    #error "Before assuming some model of geometry you have to use some matrix definition by calling the command '#include <using_{SomeMatrixAlgebraLibrary}.hpp>.'"
#endif // TBGAL_USING_MATRIX_DEFINITIONS

#include <array>
#include <cassert>
#include <iostream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "matrix_declarations.hpp"

namespace tbgal {

    template<typename MetricSpaceType> class MetricSpace;
    template<typename ScalarType, typename FactoringProductType> class FactoredMultivector;

    template<typename MetricSpaceType> struct GeometricProduct;
    template<typename MetricSpaceType> struct OuterProduct;

    namespace detail {

        template<typename T, typename... Rest>
        constexpr bool is_any_v = std::disjunction_v<std::bool_constant<std::is_same_v<T, Rest> >...>;

        template<bool AnyMultivectorType> struct GP_impl;
        template<bool AnyMultivectorType> struct OP_impl;

        template<typename MetricSpaceType> struct apply_signed_metric_impl;
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

        template<typename ResultingMatrixType>
        constexpr static decltype(auto) fill_matrix_with_factors(ResultingMatrixType &result) noexcept {
            return std::make_tuple(result, cols(result));
        }

        template<typename ResultingMatrixType, typename FirstScalarType, typename... NextTypes>
        constexpr static decltype(auto) fill_matrix_with_factors(ResultingMatrixType &result, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
            return fill_matrix_with_factors(result, args...);
        }

        template<typename ResultingMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr static decltype(auto) fill_matrix_with_factors(ResultingMatrixType &result, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            assert(is_blade(arg1));
            auto end_column_index = std::get<1>(fill_matrix_with_factors(result, args...));
            return std::make_tuple(copy_columns<Dynamic>(arg1.factors_in_signed_metric(), 0, result, end_column_index - arg1.factors_count(), arg1.factors_count()), end_column_index - arg1.factors_count());
        }

        constexpr static decltype(auto) multiply_scalars() noexcept {
            return 1;
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr static decltype(auto) multiply_scalars(FirstScalarType const &arg1, NextTypes const &... args) noexcept {
            return arg1 * multiply_scalars(args...);
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr static decltype(auto) multiply_scalars(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            return arg1.scalar() * multiply_scalars(args...);
        }

        constexpr static decltype(auto) space_ptr() noexcept {
            return nullptr;
        }

        template<typename ScalarType, typename FactoringProductType>
        constexpr static decltype(auto) space_ptr(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return &arg.space();
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr static decltype(auto) space_ptr(FirstScalarType const &, NextTypes const &... args) noexcept {
            return space_ptr(args...);
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr static decltype(auto) space_ptr(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            assert(space_ptr(args...) == &arg1.space() || space_ptr(args...) == nullptr);
            return &arg1.space();
        }

        constexpr static decltype(auto) sum_factors_count() noexcept {
            return 0;
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr static decltype(auto) sum_factors_count(FirstScalarType const &, NextTypes const &... args) noexcept {
            return sum_factors_count(args...);
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr static decltype(auto) sum_factors_count(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            assert(is_blade(arg1));
            return arg1.factors_count() + sum_factors_count(args...);
        }

    }

    template<typename Type>
    struct is_factoring_product :
        std::false_type {
    };

    template<typename Type>
    using is_factoring_product_t = typename is_factoring_product<Type>::type;

    template<typename Type>
    constexpr bool is_factoring_product_v = is_factoring_product<Type>::value;

    template<typename Type>
    struct is_metric_space :
        std::false_type {
    };

    template<typename Type>
    using is_metric_space_t = typename is_metric_space<Type>::type;

    template<typename Type>
    constexpr bool is_metric_space_v = is_metric_space<Type>::value;

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

#include "geometric_product.hpp"
#include "outer_product.hpp"
#include "left_contraction.hpp"

#include "reversion.hpp"

#include "reverse_norm.hpp"
#include "inversion.hpp"
#include "dualization.hpp"

#include "write.hpp"
#include "macro.hpp"

#include "Conformal/metric_space.hpp"
#include "Conformal/macro.hpp"

#include "Euclidean/metric_space.hpp"
#include "Euclidean/macro.hpp"

#include "Homogeneous/metric_space.hpp"
#include "Homogeneous/macro.hpp"

#endif // __TBGAL_CORE_HPP__
