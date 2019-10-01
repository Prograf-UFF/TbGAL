#ifndef __TBGAL_TOOLS_TEST_COMMON_HPP__
#define __TBGAL_TOOLS_TEST_COMMON_HPP__

#include <random>
#include <gtest/gtest.h>

using scalar_factor_t = std::double_t;
using vector_factor_t = std::array<scalar_factor_t, TESTING_VECTOR_SPACE_DIMENSIONS>;

std::default_random_engine random_engine{ static_cast<long unsigned int>(32) };
std::uniform_real_distribution<scalar_factor_t> uniform_distribution(0, 1);

auto tbgal_Dualization = [](auto const &arg) noexcept {
    return tbgal::DUAL(arg);
};

auto tbgal_LeftContraction = [](auto const &lhs, auto const &rhs) noexcept {
    return tbgal::LCONT(lhs, rhs);
};

auto tbgal_OuterProduct = [](auto const &lhs, auto const &rhs) noexcept {
    return tbgal::OP(lhs, rhs);
};

auto tbgal_ReverseNorm = [](auto const &arg) noexcept {
    return tbgal::RNORM(arg);
};

auto tbgal_Reversion = [](auto const &arg) noexcept {
    return tbgal::REVERSE(arg);
};

auto tbgal_SquaredReverseNorm = [](auto const &arg) noexcept {
    return tbgal::RNORM_SQR(arg);
};

auto tbgal_UnaryMinus = [](auto const &arg) noexcept {
    return -arg;
};

auto tbgal_UnaryPlus = [](auto const &arg) noexcept {
    return +arg;
};

auto tbgal_Undualization = [](auto const &arg) noexcept {
    return tbgal::UNDUAL(arg);
};

auto gatl_Dualization = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return dual(arg);
};

auto gatl_LeftContraction = [&](auto const &lhs, auto const &rhs) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return lcont(lhs, rhs);
};

auto gatl_OuterProduct = [&](auto const &lhs, auto const &rhs) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return op(lhs, rhs);
};

auto gatl_ReverseNorm = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return rnorm(arg);
};

auto gatl_Reversion = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return reversion(arg);
};

auto gatl_SquaredReverseNorm = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return rnorm_sqr(arg);
};

auto gatl_UnaryMinus = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return -arg;
};

auto gatl_UnaryPlus = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return +arg;
};

auto gatl_Undualization = [](auto const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return undual(arg);
};

template<std::size_t EndIndex>
struct _make_tbgal_vector_impl {
    template<typename... Types>
    constexpr static decltype(auto) eval(vector_factor_t const &vector_factor, Types const &... coords) noexcept {
        return _make_tbgal_vector_impl<EndIndex - 1>::eval(vector_factor, vector_factor[EndIndex - 1], coords...);
    }
};

template<>
struct _make_tbgal_vector_impl<0> {
    template<typename... Types>
    constexpr static decltype(auto) eval(vector_factor_t const &, Types const &... coords) noexcept {
        using namespace TESTING_TBGAL_MODEL_NAMESPACE;
        return vector(coords...);
    }
};

decltype(auto) make_tbgal_vector(vector_factor_t const &arg) noexcept {
    return _make_tbgal_vector_impl<std::tuple_size_v<vector_factor_t> >::eval(arg);
}

template<std::size_t EndIndex>
struct _make_gatl_vector_impl {
    template<typename... Types>
    constexpr static decltype(auto) eval(vector_factor_t const &vector_factor, Types const &... coords) noexcept {
        return _make_gatl_vector_impl<EndIndex - 1>::eval(vector_factor, vector_factor[EndIndex - 1], coords...);
    }
};

template<>
struct _make_gatl_vector_impl<0> {
    template<typename... Types>
    constexpr static decltype(auto) eval(vector_factor_t const &, Types const &... coords) noexcept {
        using namespace TESTING_GATL_MODEL_NAMESPACE;
        return vector(coords...);
    }
};

constexpr decltype(auto) make_gatl_vector(vector_factor_t const &arg) noexcept {
    return _make_gatl_vector_impl<std::tuple_size_v<vector_factor_t> >::eval(arg);
}

template<std::size_t EndIndex>
struct _make_tbgal_multivector_using_outer_product_impl {
    template<typename VectorFactorType, std::size_t K, typename... Types>
    constexpr static decltype(auto) eval(std::array<VectorFactorType, K> const &vector_factors, Types const &... factors) noexcept {
        return _make_tbgal_multivector_using_outer_product_impl<EndIndex - 1>::eval(vector_factors, make_tbgal_vector(vector_factors[EndIndex - 1]), factors...);
    }
};

template<>
struct _make_tbgal_multivector_using_outer_product_impl<0> {
    template<typename VectorFactorType, std::size_t K, typename... Types>
    constexpr static decltype(auto) eval(std::array<VectorFactorType, K> const &, Types const &... factors) noexcept {
        return tbgal::OP(factors...);
    }
};

template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
constexpr decltype(auto) make_tbgal_multivector_using_outer_product(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
    using namespace TESTING_TBGAL_MODEL_NAMESPACE;
    return _make_tbgal_multivector_using_outer_product_impl<K>::eval(vector_factors, scalar(scalar_factor));
}

template<std::size_t EndIndex>
struct _make_gatl_multivector_using_outer_product_impl {
    template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
    constexpr static decltype(auto) eval(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
        return gatl_OuterProduct(_make_gatl_multivector_using_outer_product_impl<EndIndex - 1>::eval(scalar_factor, vector_factors), make_gatl_vector(vector_factors[EndIndex - 1]));
    }
};

template<>
struct _make_gatl_multivector_using_outer_product_impl<0> {
    template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
    constexpr static decltype(auto) eval(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &) noexcept {
        return scalar_factor;
    }
};

template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
constexpr decltype(auto) make_gatl_multivector_using_outer_product(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
    return _make_gatl_multivector_using_outer_product_impl<K>::eval(scalar_factor, vector_factors);
}

template<typename MetricSpaceType, typename SquareMatrixType>
decltype(auto) from_tbgal_to_gatl(tbgal::FactoredMultivector<tbgal::OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
    using IndexType = typename tbgal::FactoredMultivector<tbgal::OuterProduct<MetricSpaceType>, SquareMatrixType>::IndexType;
    ga::full_multivector_t<scalar_factor_t, std::tuple_size_v<vector_factor_t> > result;
    ga::trivial_copy(arg.scalar(), result);
    for (IndexType col = 0; col != arg.factors_count(); ++col) {
        vector_factor_t factor;
        for (std::size_t row = 0; row != std::tuple_size_v<vector_factor_t>; ++row) {
            factor[row] = tbgal::detail::coeff(arg.factors(), row, col);
        }
        ga::trivial_copy(gatl_OuterProduct(result, make_gatl_vector(factor)), result);
    }
    return result;
}

template<typename ScalarType, typename = std::enable_if_t<!tbgal::is_multivector_v<ScalarType> > >
decltype(auto) from_tbgal_to_gatl(ScalarType const &arg) noexcept {
    return arg;
}

template<std::size_t K>
std::tuple<scalar_factor_t, std::array<vector_factor_t, K> > make_random_factors() noexcept {
    scalar_factor_t scalar_factor = uniform_distribution(random_engine);
    std::array<vector_factor_t, K> vector_factors;
    for (std::size_t k = 0; k != K; ++k) {
        for (std::size_t i = 0; i != k; ++i) {
            vector_factors[k][i] = 0;
        }
        for (std::size_t i = k; i != std::tuple_size_v<vector_factor_t>; ++i) {
            vector_factors[k][i] = uniform_distribution(random_engine);
        }
    }
    return std::make_tuple(scalar_factor, vector_factors);
}

template<typename LeftCoefficientType, typename LeftExpression, typename RightCoefficientType, typename RightExpression>
constexpr bool same_multivector(ga::clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, ga::clifford_expression<RightCoefficientType, RightExpression> const &rhs) noexcept {
    return ga::is_zero(lhs - rhs);
}

template<typename LeftType, typename RightCoefficientType, typename RightExpression, typename = std::enable_if_t<!ga::is_clifford_expression_v<LeftType> > >
constexpr bool same_multivector(LeftType const &lhs, ga::clifford_expression<RightCoefficientType, RightExpression> const &rhs) noexcept {
    return same_multivector(ga::scalar(lhs), rhs);
}

template<typename LeftCoefficientType, typename LeftExpression, typename RightType, typename = std::enable_if_t<!ga::is_clifford_expression_v<RightType> > >
constexpr bool same_multivector(ga::clifford_expression<LeftCoefficientType, LeftExpression> const &lhs, RightType const &rhs) noexcept {
    return same_multivector(lhs, ga::scalar(rhs));
}

template<typename LeftType, typename RightType, typename = std::enable_if_t<!(ga::is_clifford_expression_v<LeftType> || ga::is_clifford_expression_v<RightType>)> >
constexpr bool same_multivector(LeftType const &lhs, RightType const &rhs) noexcept {
    return same_multivector(ga::scalar(lhs), ga::scalar(rhs));
}

#define BINARY_OPERATION_TESTS_FOR(OPERATION, LHS_K, RHS_K) \
    TEST(OPERATION##Test, Outer##LHS_K##_Outer##RHS_K) { \
        scalar_factor_t lhs_scalar_factor, rhs_scalar_factor; \
        std::array<vector_factor_t, LHS_K> lhs_vector_factors; \
        std::array<vector_factor_t, RHS_K> rhs_vector_factors; \
        std::tie(lhs_scalar_factor, lhs_vector_factors) = make_random_factors<LHS_K>(); \
        std::tie(rhs_scalar_factor, rhs_vector_factors) = make_random_factors<RHS_K>(); \
        EXPECT_TRUE(same_multivector(from_tbgal_to_gatl(tbgal_##OPERATION(make_tbgal_multivector_using_outer_product(lhs_scalar_factor, lhs_vector_factors), make_tbgal_multivector_using_outer_product(rhs_scalar_factor, rhs_vector_factors))), gatl_##OPERATION(make_gatl_multivector_using_outer_product(lhs_scalar_factor, lhs_vector_factors), make_gatl_multivector_using_outer_product(rhs_scalar_factor, rhs_vector_factors)))); \
    }

#define BINARY_OPERATION_TESTS(LHS_K, RHS_K) \
    BINARY_OPERATION_TESTS_FOR(LeftContraction, LHS_K, RHS_K) \
    BINARY_OPERATION_TESTS_FOR(OuterProduct, LHS_K, RHS_K)

    //TODO [TEST] Addition
    //TODO [TEST] GeometricProduct
    //TODO [TEST] Subtraction

#define UNARY_OPERATION_TESTS_FOR(OPERATION, ARG_K) \
    TEST(OPERATION##Test, Outer##ARG_K) { \
        scalar_factor_t scalar_factor; \
        std::array<vector_factor_t, ARG_K> vector_factors; \
        std::tie(scalar_factor, vector_factors) = make_random_factors<ARG_K>(); \
        EXPECT_TRUE(same_multivector(from_tbgal_to_gatl(tbgal_##OPERATION(make_tbgal_multivector_using_outer_product(scalar_factor, vector_factors))), gatl_##OPERATION(make_gatl_multivector_using_outer_product(scalar_factor, vector_factors)))); \
    }

#define UNARY_OPERATION_TESTS(ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Dualization, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(ReverseNorm, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Reversion, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(SquaredReverseNorm, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(UnaryMinus, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(UnaryPlus, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Undualization, ARG_K)

    //TODO [TEST] Inversion

#if TESTING_VECTOR_SPACE_DIMENSIONS >= 0
    BINARY_OPERATION_TESTS(0, 0)
    UNARY_OPERATION_TESTS(0)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 1
    BINARY_OPERATION_TESTS(0, 1)
    BINARY_OPERATION_TESTS(1, 0)
    BINARY_OPERATION_TESTS(1, 1)
    UNARY_OPERATION_TESTS(1)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 2
    BINARY_OPERATION_TESTS(0, 2)
    BINARY_OPERATION_TESTS(1, 2)
    BINARY_OPERATION_TESTS(2, 0)
    BINARY_OPERATION_TESTS(2, 1)
    BINARY_OPERATION_TESTS(2, 2)
    UNARY_OPERATION_TESTS(2)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 3
    BINARY_OPERATION_TESTS(0, 3)
    BINARY_OPERATION_TESTS(1, 3)
    BINARY_OPERATION_TESTS(2, 3)
    BINARY_OPERATION_TESTS(3, 0)
    BINARY_OPERATION_TESTS(3, 1)
    BINARY_OPERATION_TESTS(3, 2)
    BINARY_OPERATION_TESTS(3, 3)
    UNARY_OPERATION_TESTS(3)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 4
    BINARY_OPERATION_TESTS(0, 4)
    BINARY_OPERATION_TESTS(1, 4)
    BINARY_OPERATION_TESTS(2, 4)
    BINARY_OPERATION_TESTS(3, 4)
    BINARY_OPERATION_TESTS(4, 0)
    BINARY_OPERATION_TESTS(4, 1)
    BINARY_OPERATION_TESTS(4, 2)
    BINARY_OPERATION_TESTS(4, 3)
    BINARY_OPERATION_TESTS(4, 4)
    UNARY_OPERATION_TESTS(4)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 5
    BINARY_OPERATION_TESTS(0, 5)
    BINARY_OPERATION_TESTS(1, 5)
    BINARY_OPERATION_TESTS(2, 5)
    BINARY_OPERATION_TESTS(3, 5)
    BINARY_OPERATION_TESTS(4, 5)
    BINARY_OPERATION_TESTS(5, 0)
    BINARY_OPERATION_TESTS(5, 1)
    BINARY_OPERATION_TESTS(5, 2)
    BINARY_OPERATION_TESTS(5, 3)
    BINARY_OPERATION_TESTS(5, 4)
    BINARY_OPERATION_TESTS(5, 5)
    UNARY_OPERATION_TESTS(5)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 6
    #error "The testing tool is not prepared for too many dimensions."
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 6
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 5
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 4
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 3
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 2
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 1
#endif // TESTING_VECTOR_SPACE_DIMENSIONS >= 0

#endif // __TBGAL_TOOLS_TEST_COMMON_HPP__
