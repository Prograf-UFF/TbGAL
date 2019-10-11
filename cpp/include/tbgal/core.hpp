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
    template<typename FactoringProductType, typename SquareMatrixType> class FactoredMultivector;

    template<typename MetricSpaceType> struct GeometricProduct;
    template<typename MetricSpaceType> struct OuterProduct;

    namespace detail {

        template<typename T, typename... Rest>
        constexpr bool is_any_v = std::disjunction_v<std::bool_constant<std::is_same_v<T, Rest> >...>;

        template<bool AnyMultivectorType> struct OP_impl;

        template<typename MetricSpaceType> struct apply_signed_metric_impl;
        template<typename MetricSpaceType> struct from_actual_to_signed_metric_impl;
        template<typename MetricSpaceType> struct from_signed_to_actual_metric_impl;

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
