#ifndef __TBGAL_CORE_HPP__
#define __TBGAL_CORE_HPP__

#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <tuple>
#include <vector>

#ifndef TBGAL_DEFAULT_SCALAR_TYPE
    #define TBGAL_DEFAULT_SCALAR_TYPE std::double_t
#endif // TBGAL_DEFAULT_SCALAR_TYPE

#ifndef TBGAL_DEFAULT_INDEX_TYPE
    #define TBGAL_DEFAULT_INDEX_TYPE std::int64_t
#endif // TBGAL_DEFAULT_INDEX_TYPE

namespace tbgal {

    using DefaultScalarType = TBGAL_DEFAULT_SCALAR_TYPE;
    using DefaultIndexType = TBGAL_DEFAULT_INDEX_TYPE;

    constexpr static DefaultIndexType Dynamic = -1;

    namespace detail {

        template<bool AnyMultivectorType> struct OP_impl;

    }

}

#include "exception.hpp"

#include "traits.hpp"
#include "matrix_traits.hpp"

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

#include "Euclidean/metric_space.hpp"
#include "Euclidean/macro.hpp"

#include "Homogeneous/metric_space.hpp"
#include "Homogeneous/macro.hpp"

#endif // __TBGAL_CORE_HPP__
