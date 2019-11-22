#ifndef __TBGAL_CORE_HPP__
#define __TBGAL_CORE_HPP__

#ifndef TBGAL_USING_MATRIX_DEFINITIONS
    #error "Before assuming some model of geometry you have to use some matrix definition by calling the command '#include <using_{SomeMatrixAlgebraLibrary}.hpp>.'"
#endif // TBGAL_USING_MATRIX_DEFINITIONS

#include <cassert>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>

#include "matrix_declarations.hpp"

namespace tbgal {

    //TODO {DEBUG}
    template<typename MatrixType>
    void print_matrix(MatrixType const &arg) {
        std::cout << "{";
        for (int row = 0; row != arg.rows(); ++row) {
            if (row != 0) std::cout << ", ";
            std::cout << "{";
            for (int col = 0; col != arg.cols(); ++col) {
                if (col != 0) std::cout << ", ";
                std::cout << std::setprecision(20) << arg(row, col);
            }
            std::cout << "}";
        }
        std::cout << "}";
    }

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

#include "Signed/metric_space.hpp"
#include "Signed/macro.hpp"

#endif // __TBGAL_CORE_HPP__
