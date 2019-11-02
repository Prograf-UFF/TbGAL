#ifndef __TBGAL_CORE_HPP__
#define __TBGAL_CORE_HPP__

#ifndef TBGAL_USING_MATRIX_DEFINITIONS
    #error "Before assuming some model of geometry you have to use some matrix definition by calling the command '#include <using_{SomeMatrixAlgebraLibrary}.hpp>.'"
#endif // TBGAL_USING_MATRIX_DEFINITIONS

#include <cassert>
#include <iostream>
#include <iterator>
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

        template <typename T, typename = void>
        struct is_iterator : std::false_type {
        };

        template <class T>
        struct is_iterator<T, std::void_t<typename std::iterator_traits<T>::iterator_category> > : std::true_type {
        };

        template<typename T>
        constexpr bool is_iterator_v = is_iterator<T>::value;

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
        constexpr void fill_matrix_with_first_factors(ResultingMatrixType &) noexcept {
            // Nothing to be done.
        }

        template<typename ResultingMatrixType, typename ScalarType>
        constexpr void fill_matrix_with_first_factors(ResultingMatrixType &, ScalarType const &) noexcept {
            // Nothing to be done.
        }

        template<typename ResultingMatrixType, typename ScalarType, typename FactoringProductType>
        constexpr void fill_matrix_with_first_factors(ResultingMatrixType &result, FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            if (arg.factors_count() > 0) {
                assign_block<rows_at_compile_time_v<ResultingMatrixType>, Dynamic>(arg.factors_in_signed_metric(), 0, 0, result, 0, 0, rows(result), arg.factors_count());
            }
        }

        template<typename ResultingMatrixType, typename FirstScalarType, typename... NextTypes>
        constexpr void fill_matrix_with_first_factors(ResultingMatrixType &result, FirstScalarType const &, NextTypes const &... args) noexcept {
            fill_matrix_with_first_factors(result, args...);
        }

        template<typename ResultingMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr void fill_matrix_with_first_factors(ResultingMatrixType &result, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            if (arg1.factors_count() > 0) {
                assign_block<rows_at_compile_time_v<ResultingMatrixType>, Dynamic>(arg1.factors_in_signed_metric(), 0, 0, result, 0, 0, rows(result), arg1.factors_count());
            }
            else {
                fill_matrix_with_first_factors(result, args...);
            }
        }

        template<typename ResultingMatrixType>
        constexpr std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix_with_tail_factors(ResultingMatrixType &result) noexcept {
            return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, cols(result));
        }

        template<typename ResultingMatrixType, typename ScalarType>
        constexpr std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix_with_tail_factors(ResultingMatrixType &result, ScalarType const &) noexcept {
            return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, cols(result));
        }

        template<typename ResultingMatrixType, typename ScalarType, typename FactoringProductType>
        constexpr std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix_with_tail_factors(ResultingMatrixType &result, FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            if (arg.factors_count() > 0) {
                auto end_column_index = cols(result);
                if (end_column_index >= arg.factors_count()) {
                    assign_block<rows_at_compile_time_v<ResultingMatrixType>, Dynamic>(arg.factors_in_signed_metric(), 0, 0, result, 0, end_column_index - arg.factors_count(), rows(result), arg.factors_count());
                    return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, end_column_index - arg.factors_count());
                }
                assert(end_column_index == 0);
                return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, 0);
            }
            return fill_matrix_with_tail_factors(result);
        }

        template<typename ResultingMatrixType, typename FirstScalarType, typename... NextTypes>
        constexpr std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix_with_tail_factors(ResultingMatrixType &result, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
            return fill_matrix_with_tail_factors(result, args...);
        }

        template<typename ResultingMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix_with_tail_factors(ResultingMatrixType &result, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            if (arg1.factors_count() > 0) {
                auto end_column_index = std::get<1>(fill_matrix_with_tail_factors(result, args...));
                if (end_column_index >= arg1.factors_count()) {
                    assign_block<rows_at_compile_time_v<ResultingMatrixType>, Dynamic>(arg1.factors_in_signed_metric(), 0, 0, result, 0, end_column_index - arg1.factors_count(), rows(result), arg1.factors_count());
                    return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, end_column_index - arg1.factors_count());
                }
                assert(end_column_index == 0);
                return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, 0);
            }
            return fill_matrix_with_tail_factors(result, args...);
        }

        constexpr decltype(auto) first_factors_count() noexcept {
            return 0;
        }

        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) first_factors_count(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return arg.factors_count();
        }

        template<typename ScalarType>
        constexpr decltype(auto) first_factors_count(ScalarType const &) noexcept {
            return 0;
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr decltype(auto) first_factors_count(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            if (arg1.factors_count() > 0) {
                return arg1.factors_count();
            }
            else {
                return first_factors_count(args...);
            }
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr decltype(auto) first_factors_count(FirstScalarType const &, NextTypes const &... args) noexcept {
            return first_factors_count(args...);
        }

        constexpr decltype(auto) multiply_scalars() noexcept {
            return 1;
        }

        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) multiply_scalars(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return arg.scalar();
        }

        template<typename ScalarType>
        constexpr decltype(auto) multiply_scalars(ScalarType const &arg) noexcept {
            return arg;
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr decltype(auto) multiply_scalars(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            return arg1.scalar() * multiply_scalars(args...);
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr decltype(auto) multiply_scalars(FirstScalarType const &arg1, NextTypes const &... args) noexcept {
            return arg1 * multiply_scalars(args...);
        }

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

        constexpr decltype(auto) sum_factors_count() noexcept {
            return 0;
        }

        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) sum_factors_count(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            assert(is_blade(arg));
            return arg.factors_count();
        }

        template<typename ScalarType>
        constexpr decltype(auto) sum_factors_count(ScalarType const &) noexcept {
            return 0;
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr decltype(auto) sum_factors_count(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            assert(is_blade(arg1));
            return arg1.factors_count() + sum_factors_count(args...);
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr decltype(auto) sum_factors_count(FirstScalarType const &, NextTypes const &... args) noexcept {
            return sum_factors_count(args...);
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
