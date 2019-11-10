#ifndef __TBGAL_TOOLS_TEST_COMMON_HPP__
#define __TBGAL_TOOLS_TEST_COMMON_HPP__

#include <random>
#include <gtest/gtest.h>

using scalar_factor_t = std::double_t;
using vector_factor_t = std::array<scalar_factor_t, TESTING_VECTOR_SPACE_DIMENSIONS>;

std::default_random_engine random_engine{ static_cast<long unsigned int>(32) };
std::uniform_real_distribution<scalar_factor_t> uniform_distribution(0, 1);

template<typename FirstType, typename SecondType>
constexpr decltype(auto) tbgal_DotProduct(FirstType const &arg1, SecondType const &arg2) noexcept {
    return tbgal::dot(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) tbgal_Dualization(Type const &arg) noexcept {
    return tbgal::dual(arg);
}

template<typename FirstType, typename... NextTypes>
constexpr decltype(auto) tbgal_GeometricProduct(FirstType const &arg1, NextTypes const &... args) noexcept {
    return tbgal::gp(arg1, args...);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) tbgal_HestenesInnerProduct(FirstType const &arg1, SecondType const &arg2) noexcept {
    return tbgal::hip(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) tbgal_Inversion(Type const &arg) noexcept {
    return tbgal::inverse(arg);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) tbgal_LeftContraction(FirstType const &arg1, SecondType const &arg2) noexcept {
    return tbgal::lcont(arg1, arg2);
}

template<typename FirstType, typename... NextTypes>
constexpr decltype(auto) tbgal_OuterProduct(FirstType const &arg1, NextTypes const &... args) noexcept {
    return tbgal::op(arg1, args...);
}

template<typename Type>
constexpr decltype(auto) tbgal_ReverseNorm(Type const &arg) noexcept {
    return tbgal::rnorm(arg);
}

template<typename Type>
constexpr decltype(auto) tbgal_Reversion(Type const &arg) noexcept {
    return tbgal::reverse(arg);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) tbgal_RightContraction(FirstType const &arg1, SecondType const &arg2) noexcept {
    return tbgal::rcont(arg1, arg2);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) tbgal_ScalarProduct(FirstType const &arg1, SecondType const &arg2) noexcept {
    return tbgal::sp(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) tbgal_SquaredReverseNorm(Type const &arg) noexcept {
    return tbgal::rnorm_sqr(arg);
}

template<typename Type>
constexpr decltype(auto) tbgal_UnaryMinus(Type const &arg) noexcept {
    return -arg;
}

template<typename Type>
constexpr decltype(auto) tbgal_UnaryPlus(Type const &arg) noexcept {
    return +arg;
}

template<typename Type>
constexpr decltype(auto) tbgal_Undualization(Type const &arg) noexcept {
    return tbgal::undual(arg);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) gatl_DotProduct(FirstType const &arg1, SecondType const &arg2) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return dot(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) gatl_Dualization(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return dual(arg);
}

template<typename Type>
constexpr decltype(auto) gatl_GeometricProduct(Type const &arg) noexcept {
    return arg;
}

template<typename FirstType, typename... NextTypes>
constexpr decltype(auto) gatl_GeometricProduct(FirstType const &arg1, NextTypes const &... args) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return gp(arg1, gatl_GeometricProduct(args...));
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) gatl_HestenesInnerProduct(FirstType const &arg1, SecondType const &arg2) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return hip(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) gatl_Inversion(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return inv(arg);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) gatl_LeftContraction(FirstType const &arg1, SecondType const &arg2) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return lcont(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) gatl_OuterProduct(Type const &arg) noexcept {
    return arg;
}

template<typename FirstType, typename... NextTypes>
constexpr decltype(auto) gatl_OuterProduct(FirstType const &arg1, NextTypes const &... args) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return op(arg1, gatl_OuterProduct(args...));
}

template<typename Type>
constexpr decltype(auto) gatl_ReverseNorm(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return rnorm(arg);
}

template<typename Type>
constexpr decltype(auto) gatl_Reversion(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return reversion(arg);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) gatl_RightContraction(FirstType const &arg1, SecondType const &arg2) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return rcont(arg1, arg2);
}

template<typename FirstType, typename SecondType>
constexpr decltype(auto) gatl_ScalarProduct(FirstType const &arg1, SecondType const &arg2) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return sp(arg1, arg2);
}

template<typename Type>
constexpr decltype(auto) gatl_SquaredReverseNorm(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return rnorm_sqr(arg);
}

template<typename Type>
constexpr decltype(auto) gatl_UnaryMinus(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return -arg;
}

template<typename Type>
constexpr decltype(auto) gatl_UnaryPlus(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return +arg;
}

template<typename Type>
constexpr decltype(auto) gatl_Undualization(Type const &arg) noexcept {
    using namespace TESTING_GATL_MODEL_NAMESPACE;
    return undual(arg);
}

decltype(auto) make_tbgal_vector(vector_factor_t const &arg) noexcept {
    using namespace TESTING_TBGAL_MODEL_NAMESPACE;
    return vector(arg.begin(), arg.end());
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
struct _make_tbgal_multivector_using_geometric_product_impl {
    template<typename VectorFactorType, std::size_t K, typename... Types>
    constexpr static decltype(auto) eval(std::array<VectorFactorType, K> const &vector_factors, Types const &... factors) noexcept {
        return _make_tbgal_multivector_using_geometric_product_impl<EndIndex - 1>::eval(vector_factors, make_tbgal_vector(vector_factors[EndIndex - 1]), factors...);
    }
};

template<>
struct _make_tbgal_multivector_using_geometric_product_impl<0> {
    template<typename VectorFactorType, std::size_t K, typename... Types>
    constexpr static decltype(auto) eval(std::array<VectorFactorType, K> const &, Types const &... factors) noexcept {
        return tbgal::gp(factors...);
    }
};

template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
constexpr decltype(auto) make_tbgal_multivector_using_geometric_product(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
    using namespace TESTING_TBGAL_MODEL_NAMESPACE;
    return _make_tbgal_multivector_using_geometric_product_impl<K>::eval(vector_factors, scalar(scalar_factor));
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
        return tbgal::op(factors...);
    }
};

template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
constexpr decltype(auto) make_tbgal_multivector_using_outer_product(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
    using namespace TESTING_TBGAL_MODEL_NAMESPACE;
    return _make_tbgal_multivector_using_outer_product_impl<K>::eval(vector_factors, scalar(scalar_factor));
}

template<std::size_t EndIndex>
struct _make_gatl_multivector_using_geometric_product_impl {
    template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
    constexpr static decltype(auto) eval(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
        return gatl_GeometricProduct(_make_gatl_multivector_using_geometric_product_impl<EndIndex - 1>::eval(scalar_factor, vector_factors), make_gatl_vector(vector_factors[EndIndex - 1]));
    }
};

template<>
struct _make_gatl_multivector_using_geometric_product_impl<0> {
    template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
    constexpr static decltype(auto) eval(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &) noexcept {
        return scalar_factor;
    }
};

template<typename ScalarFactorType, typename VectorFactorType, std::size_t K>
constexpr decltype(auto) make_gatl_multivector_using_geometric_product(ScalarFactorType const &scalar_factor, std::array<VectorFactorType, K> const &vector_factors) noexcept {
    return _make_gatl_multivector_using_geometric_product_impl<K>::eval(scalar_factor, vector_factors);
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

template<typename ScalarType, typename MetricSpaceType>
decltype(auto) from_tbgal_to_gatl(tbgal::FactoredMultivector<ScalarType, tbgal::GeometricProduct<MetricSpaceType> > const &arg) noexcept {
    using IndexType = typename tbgal::FactoredMultivector<ScalarType, tbgal::GeometricProduct<MetricSpaceType> >::IndexType;
    ga::full_multivector_t<scalar_factor_t, std::tuple_size_v<vector_factor_t> > result;
    ga::trivial_copy(arg.scalar(), result);
    auto factors = arg.factors_in_actual_metric();
    for (IndexType col = 0; col != arg.factors_count(); ++col) {
        vector_factor_t factor;
        for (std::size_t row = 0; row != std::tuple_size_v<vector_factor_t>; ++row) {
            factor[row] = tbgal::detail::coeff(factors, row, col);
        }
        ga::trivial_copy(gatl_GeometricProduct(result, make_gatl_vector(factor)), result);
    }
    return result;
}

template<typename ScalarType, typename MetricSpaceType>
decltype(auto) from_tbgal_to_gatl(tbgal::FactoredMultivector<ScalarType, tbgal::OuterProduct<MetricSpaceType> > const &arg) noexcept {
    using IndexType = typename tbgal::FactoredMultivector<ScalarType, tbgal::OuterProduct<MetricSpaceType> >::IndexType;
    ga::full_multivector_t<scalar_factor_t, std::tuple_size_v<vector_factor_t> > result;
    ga::trivial_copy(arg.scalar(), result);
    auto factors = arg.factors_in_actual_metric();
    for (IndexType col = 0; col != arg.factors_count(); ++col) {
        vector_factor_t factor;
        for (std::size_t row = 0; row != std::tuple_size_v<vector_factor_t>; ++row) {
            factor[row] = tbgal::detail::coeff(factors, row, col);
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

    return std::tuple<scalar_factor_t, std::array<vector_factor_t, K> >(scalar_factor, vector_factors);
}

template<typename LeftCoefficientType, typename LeftExpression, typename RightCoefficientType, typename RightExpression>
constexpr bool same_multivector(ga::clifford_expression<LeftCoefficientType, LeftExpression> const &arg1, ga::clifford_expression<RightCoefficientType, RightExpression> const &arg2) noexcept {
    return ga::is_zero(arg1 - arg2);
}

template<typename LeftType, typename RightCoefficientType, typename RightExpression, typename = std::enable_if_t<!ga::is_clifford_expression_v<LeftType> > >
constexpr bool same_multivector(LeftType const &arg1, ga::clifford_expression<RightCoefficientType, RightExpression> const &arg2) noexcept {
    return same_multivector(ga::scalar(arg1), arg2);
}

template<typename LeftCoefficientType, typename LeftExpression, typename RightType, typename = std::enable_if_t<!ga::is_clifford_expression_v<RightType> > >
constexpr bool same_multivector(ga::clifford_expression<LeftCoefficientType, LeftExpression> const &arg1, RightType const &arg2) noexcept {
    return same_multivector(arg1, ga::scalar(arg2));
}

template<typename LeftType, typename RightType, typename = std::enable_if_t<!(ga::is_clifford_expression_v<LeftType> || ga::is_clifford_expression_v<RightType>)> >
constexpr bool same_multivector(LeftType const &arg1, RightType const &arg2) noexcept {
    return same_multivector(ga::scalar(arg1), ga::scalar(arg2));
}

#define BINARY_OPERATION_TESTS_FOR(OPERATION, ARG1_K, ARG2_K) \
    TEST(OPERATION##Test, Outer##ARG1_K##_Outer##ARG2_K) { \
        scalar_factor_t arg1_scalar_factor, arg2_scalar_factor; \
        std::array<vector_factor_t, ARG1_K> arg1_vector_factors; \
        std::array<vector_factor_t, ARG2_K> arg2_vector_factors; \
        std::tie(arg1_scalar_factor, arg1_vector_factors) = make_random_factors<ARG1_K>(); \
        std::tie(arg2_scalar_factor, arg2_vector_factors) = make_random_factors<ARG2_K>(); \
        auto tbgal_result = from_tbgal_to_gatl(tbgal_##OPERATION(make_tbgal_multivector_using_outer_product(arg1_scalar_factor, arg1_vector_factors), make_tbgal_multivector_using_outer_product(arg2_scalar_factor, arg2_vector_factors))); \
        auto gatl_result = gatl_##OPERATION(make_gatl_multivector_using_outer_product(arg1_scalar_factor, arg1_vector_factors), make_gatl_multivector_using_outer_product(arg2_scalar_factor, arg2_vector_factors)); \
        bool is_same = same_multivector(tbgal_result, gatl_result); \
        if (!is_same) { \
            using namespace TESTING_GATL_MODEL_NAMESPACE; \
            std::cout << "  TbGAL = " << tbgal_result << std::endl; \
            std::cout << "   GATL = " << gatl_result << std::endl; \
        } \
        EXPECT_TRUE(is_same); \
    }

#define BINARY_OPERATION_TESTS(ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(DotProduct, ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(GeometricProduct, ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(HestenesInnerProduct, ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(LeftContraction, ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(OuterProduct, ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(ScalarProduct, ARG1_K, ARG2_K) \
    BINARY_OPERATION_TESTS_FOR(RightContraction, ARG1_K, ARG2_K)

    //TODO [TEST] Addition
    //TODO [TEST] Subtraction

#define TERNARY_OPERATION_TESTS_FOR(OPERATION, ARG1_K, ARG2_K, ARG3_K) \
    TEST(OPERATION##Test, Outer##ARG1_K##_Outer##ARG2_K##_Outer##ARG3_K) { \
        scalar_factor_t arg1_scalar_factor, arg2_scalar_factor, arg3_scalar_factor; \
        std::array<vector_factor_t, ARG1_K> arg1_vector_factors; \
        std::array<vector_factor_t, ARG2_K> arg2_vector_factors; \
        std::array<vector_factor_t, ARG3_K> arg3_vector_factors; \
        std::tie(arg1_scalar_factor, arg1_vector_factors) = make_random_factors<ARG1_K>(); \
        std::tie(arg2_scalar_factor, arg2_vector_factors) = make_random_factors<ARG2_K>(); \
        std::tie(arg3_scalar_factor, arg3_vector_factors) = make_random_factors<ARG3_K>(); \
        auto tbgal_result = from_tbgal_to_gatl(tbgal_##OPERATION(make_tbgal_multivector_using_outer_product(arg1_scalar_factor, arg1_vector_factors), make_tbgal_multivector_using_outer_product(arg2_scalar_factor, arg2_vector_factors), make_tbgal_multivector_using_outer_product(arg3_scalar_factor, arg3_vector_factors))); \
        auto gatl_result = gatl_##OPERATION(make_gatl_multivector_using_outer_product(arg1_scalar_factor, arg1_vector_factors), make_gatl_multivector_using_outer_product(arg2_scalar_factor, arg2_vector_factors), make_gatl_multivector_using_outer_product(arg3_scalar_factor, arg3_vector_factors)); \
        bool is_same = same_multivector(tbgal_result, gatl_result); \
        if (!is_same) { \
            using namespace TESTING_GATL_MODEL_NAMESPACE; \
            std::cout << "  TbGAL = " << tbgal_result << std::endl; \
            std::cout << "   GATL = " << gatl_result << std::endl; \
        } \
        EXPECT_TRUE(is_same); \
    }

#define TERNARY_OPERATION_TESTS(ARG1_K, ARG2_K, ARG3_K) \
    /*TODO TERNARY_OPERATION_TESTS_FOR(GeometricProduct, ARG1_K, ARG2_K, ARG3_K)*/ \
    TERNARY_OPERATION_TESTS_FOR(OuterProduct, ARG1_K, ARG2_K, ARG3_K)

    //TODO [TEST] Addition
    //TODO [TEST] Subtraction

#define UNARY_OPERATION_TESTS_FOR(OPERATION, ARG_K) \
    /*TODO [TEST] Someday*/ \
    /*TEST(OPERATION##Test, Geometric##ARG_K) {*/ \
    /*    scalar_factor_t scalar_factor;*/ \
    /*    std::array<vector_factor_t, ARG_K> vector_factors;*/ \
    /*    std::tie(scalar_factor, vector_factors) = make_random_factors<ARG_K>();*/ \
    /*   auto tbgal_result = from_tbgal_to_gatl(tbgal_##OPERATION(make_tbgal_multivector_using_geometric_product(scalar_factor, vector_factors)));*/ \
    /*    auto gatl_result = gatl_##OPERATION(make_gatl_multivector_using_geometric_product(scalar_factor, vector_factors));*/ \
    /*    bool is_same = same_multivector(tbgal_result, gatl_result);*/ \
    /*    if (!is_same) {*/ \
    /*        using namespace TESTING_GATL_MODEL_NAMESPACE;*/ \
    /*        std::cout << "  TbGAL = " << tbgal_result << std::endl;*/ \
    /*        std::cout << "   GATL = " << gatl_result << std::endl;*/ \
    /*    }*/ \
    /*    EXPECT_TRUE(is_same);*/ \
    /*}*/ \
    \
    TEST(OPERATION##Test, Outer##ARG_K) { \
        scalar_factor_t scalar_factor; \
        std::array<vector_factor_t, ARG_K> vector_factors; \
        std::tie(scalar_factor, vector_factors) = make_random_factors<ARG_K>(); \
        auto tbgal_result = from_tbgal_to_gatl(tbgal_##OPERATION(make_tbgal_multivector_using_outer_product(scalar_factor, vector_factors))); \
        auto gatl_result = gatl_##OPERATION(make_gatl_multivector_using_outer_product(scalar_factor, vector_factors)); \
        bool is_same = same_multivector(tbgal_result, gatl_result); \
        if (!is_same) { \
            using namespace TESTING_GATL_MODEL_NAMESPACE; \
            std::cout << "  TbGAL = " << tbgal_result << std::endl; \
            std::cout << "   GATL = " << gatl_result << std::endl; \
        } \
        EXPECT_TRUE(is_same); \
    }

#define UNARY_OPERATION_TESTS(ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Dualization, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Inversion, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Reversion, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(SquaredReverseNorm, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(UnaryMinus, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(UnaryPlus, ARG_K) \
    UNARY_OPERATION_TESTS_FOR(Undualization, ARG_K)

#if TESTING_VECTOR_SPACE_DIMENSIONS >= 0
    BINARY_OPERATION_TESTS(0, 0)
    TERNARY_OPERATION_TESTS(0, 0, 0)
    UNARY_OPERATION_TESTS(0)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 1
    BINARY_OPERATION_TESTS(0, 1)
    BINARY_OPERATION_TESTS(1, 0)
    BINARY_OPERATION_TESTS(1, 1)
    TERNARY_OPERATION_TESTS(0, 0, 1)
    TERNARY_OPERATION_TESTS(0, 1, 0)
    TERNARY_OPERATION_TESTS(0, 1, 1)
    TERNARY_OPERATION_TESTS(1, 0, 0)
    TERNARY_OPERATION_TESTS(1, 0, 1)
    TERNARY_OPERATION_TESTS(1, 1, 0)
    TERNARY_OPERATION_TESTS(1, 1, 1)
    UNARY_OPERATION_TESTS(1)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 2
    BINARY_OPERATION_TESTS(0, 2)
    BINARY_OPERATION_TESTS(1, 2)
    BINARY_OPERATION_TESTS(2, 0)
    BINARY_OPERATION_TESTS(2, 1)
    BINARY_OPERATION_TESTS(2, 2)
    TERNARY_OPERATION_TESTS(0, 0, 2)
    TERNARY_OPERATION_TESTS(0, 1, 2)
    TERNARY_OPERATION_TESTS(0, 2, 0)
    TERNARY_OPERATION_TESTS(0, 2, 1)
    TERNARY_OPERATION_TESTS(0, 2, 2)
    TERNARY_OPERATION_TESTS(1, 0, 2)
    TERNARY_OPERATION_TESTS(1, 1, 2)
    TERNARY_OPERATION_TESTS(1, 2, 0)
    TERNARY_OPERATION_TESTS(1, 2, 1)
    TERNARY_OPERATION_TESTS(1, 2, 2)
    TERNARY_OPERATION_TESTS(2, 0, 0)
    TERNARY_OPERATION_TESTS(2, 0, 1)
    TERNARY_OPERATION_TESTS(2, 0, 2)
    TERNARY_OPERATION_TESTS(2, 1, 0)
    TERNARY_OPERATION_TESTS(2, 1, 1)
    TERNARY_OPERATION_TESTS(2, 1, 2)
    TERNARY_OPERATION_TESTS(2, 2, 0)
    TERNARY_OPERATION_TESTS(2, 2, 1)
    TERNARY_OPERATION_TESTS(2, 2, 2)
    UNARY_OPERATION_TESTS(2)
#if TESTING_VECTOR_SPACE_DIMENSIONS >= 3
    BINARY_OPERATION_TESTS(0, 3)
    BINARY_OPERATION_TESTS(1, 3)
    BINARY_OPERATION_TESTS(2, 3)
    BINARY_OPERATION_TESTS(3, 0)
    BINARY_OPERATION_TESTS(3, 1)
    BINARY_OPERATION_TESTS(3, 2)
    BINARY_OPERATION_TESTS(3, 3)
    TERNARY_OPERATION_TESTS(0, 0, 3)
    TERNARY_OPERATION_TESTS(0, 1, 3)
    TERNARY_OPERATION_TESTS(0, 2, 3)
    TERNARY_OPERATION_TESTS(0, 3, 0)
    TERNARY_OPERATION_TESTS(0, 3, 1)
    TERNARY_OPERATION_TESTS(0, 3, 2)
    TERNARY_OPERATION_TESTS(0, 3, 3)
    TERNARY_OPERATION_TESTS(1, 0, 3)
    TERNARY_OPERATION_TESTS(1, 1, 3)
    TERNARY_OPERATION_TESTS(1, 2, 3)
    TERNARY_OPERATION_TESTS(1, 3, 0)
    TERNARY_OPERATION_TESTS(1, 3, 1)
    TERNARY_OPERATION_TESTS(1, 3, 2)
    TERNARY_OPERATION_TESTS(1, 3, 3)
    TERNARY_OPERATION_TESTS(2, 0, 3)
    TERNARY_OPERATION_TESTS(2, 1, 3)
    TERNARY_OPERATION_TESTS(2, 2, 3)
    TERNARY_OPERATION_TESTS(2, 3, 0)
    TERNARY_OPERATION_TESTS(2, 3, 1)
    TERNARY_OPERATION_TESTS(2, 3, 2)
    TERNARY_OPERATION_TESTS(2, 3, 3)
    TERNARY_OPERATION_TESTS(3, 0, 0)
    TERNARY_OPERATION_TESTS(3, 0, 1)
    TERNARY_OPERATION_TESTS(3, 0, 2)
    TERNARY_OPERATION_TESTS(3, 0, 3)
    TERNARY_OPERATION_TESTS(3, 1, 0)
    TERNARY_OPERATION_TESTS(3, 1, 1)
    TERNARY_OPERATION_TESTS(3, 1, 2)
    TERNARY_OPERATION_TESTS(3, 1, 3)
    TERNARY_OPERATION_TESTS(3, 2, 0)
    TERNARY_OPERATION_TESTS(3, 2, 1)
    TERNARY_OPERATION_TESTS(3, 2, 2)
    TERNARY_OPERATION_TESTS(3, 2, 3)
    TERNARY_OPERATION_TESTS(3, 3, 0)
    TERNARY_OPERATION_TESTS(3, 3, 1)
    TERNARY_OPERATION_TESTS(3, 3, 2)
    TERNARY_OPERATION_TESTS(3, 3, 3)
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
    TERNARY_OPERATION_TESTS(0, 0, 4)
    TERNARY_OPERATION_TESTS(0, 1, 4)
    TERNARY_OPERATION_TESTS(0, 2, 4)
    TERNARY_OPERATION_TESTS(0, 3, 4)
    TERNARY_OPERATION_TESTS(0, 4, 0)
    TERNARY_OPERATION_TESTS(0, 4, 1)
    TERNARY_OPERATION_TESTS(0, 4, 2)
    TERNARY_OPERATION_TESTS(0, 4, 3)
    TERNARY_OPERATION_TESTS(0, 4, 4)
    TERNARY_OPERATION_TESTS(1, 0, 4)
    TERNARY_OPERATION_TESTS(1, 1, 4)
    TERNARY_OPERATION_TESTS(1, 2, 4)
    TERNARY_OPERATION_TESTS(1, 3, 4)
    TERNARY_OPERATION_TESTS(1, 4, 0)
    TERNARY_OPERATION_TESTS(1, 4, 1)
    TERNARY_OPERATION_TESTS(1, 4, 2)
    TERNARY_OPERATION_TESTS(1, 4, 3)
    TERNARY_OPERATION_TESTS(1, 4, 4)
    TERNARY_OPERATION_TESTS(2, 0, 4)
    TERNARY_OPERATION_TESTS(2, 1, 4)
    TERNARY_OPERATION_TESTS(2, 2, 4)
    TERNARY_OPERATION_TESTS(2, 3, 4)
    TERNARY_OPERATION_TESTS(2, 4, 0)
    TERNARY_OPERATION_TESTS(2, 4, 1)
    TERNARY_OPERATION_TESTS(2, 4, 2)
    TERNARY_OPERATION_TESTS(2, 4, 3)
    TERNARY_OPERATION_TESTS(2, 4, 4)
    TERNARY_OPERATION_TESTS(3, 0, 4)
    TERNARY_OPERATION_TESTS(3, 1, 4)
    TERNARY_OPERATION_TESTS(3, 2, 4)
    TERNARY_OPERATION_TESTS(3, 3, 4)
    TERNARY_OPERATION_TESTS(3, 4, 0)
    TERNARY_OPERATION_TESTS(3, 4, 1)
    TERNARY_OPERATION_TESTS(3, 4, 2)
    TERNARY_OPERATION_TESTS(3, 4, 3)
    TERNARY_OPERATION_TESTS(3, 4, 4)
    TERNARY_OPERATION_TESTS(4, 0, 0)
    TERNARY_OPERATION_TESTS(4, 0, 1)
    TERNARY_OPERATION_TESTS(4, 0, 2)
    TERNARY_OPERATION_TESTS(4, 0, 3)
    TERNARY_OPERATION_TESTS(4, 0, 4)
    TERNARY_OPERATION_TESTS(4, 1, 0)
    TERNARY_OPERATION_TESTS(4, 1, 1)
    TERNARY_OPERATION_TESTS(4, 1, 2)
    TERNARY_OPERATION_TESTS(4, 1, 3)
    TERNARY_OPERATION_TESTS(4, 1, 4)
    TERNARY_OPERATION_TESTS(4, 2, 0)
    TERNARY_OPERATION_TESTS(4, 2, 1)
    TERNARY_OPERATION_TESTS(4, 2, 2)
    TERNARY_OPERATION_TESTS(4, 2, 3)
    TERNARY_OPERATION_TESTS(4, 2, 4)
    TERNARY_OPERATION_TESTS(4, 3, 0)
    TERNARY_OPERATION_TESTS(4, 3, 1)
    TERNARY_OPERATION_TESTS(4, 3, 2)
    TERNARY_OPERATION_TESTS(4, 3, 3)
    TERNARY_OPERATION_TESTS(4, 3, 4)
    TERNARY_OPERATION_TESTS(4, 4, 0)
    TERNARY_OPERATION_TESTS(4, 4, 1)
    TERNARY_OPERATION_TESTS(4, 4, 2)
    TERNARY_OPERATION_TESTS(4, 4, 3)
    TERNARY_OPERATION_TESTS(4, 4, 4)
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
    TERNARY_OPERATION_TESTS(0, 0, 5)
    TERNARY_OPERATION_TESTS(0, 1, 5)
    TERNARY_OPERATION_TESTS(0, 2, 5)
    TERNARY_OPERATION_TESTS(0, 3, 5)
    TERNARY_OPERATION_TESTS(0, 4, 5)
    TERNARY_OPERATION_TESTS(0, 5, 0)
    TERNARY_OPERATION_TESTS(0, 5, 1)
    TERNARY_OPERATION_TESTS(0, 5, 2)
    TERNARY_OPERATION_TESTS(0, 5, 3)
    TERNARY_OPERATION_TESTS(0, 5, 4)
    TERNARY_OPERATION_TESTS(0, 5, 5)
    TERNARY_OPERATION_TESTS(1, 0, 5)
    TERNARY_OPERATION_TESTS(1, 1, 5)
    TERNARY_OPERATION_TESTS(1, 2, 5)
    TERNARY_OPERATION_TESTS(1, 3, 5)
    TERNARY_OPERATION_TESTS(1, 4, 5)
    TERNARY_OPERATION_TESTS(1, 5, 0)
    TERNARY_OPERATION_TESTS(1, 5, 1)
    TERNARY_OPERATION_TESTS(1, 5, 2)
    TERNARY_OPERATION_TESTS(1, 5, 3)
    TERNARY_OPERATION_TESTS(1, 5, 4)
    TERNARY_OPERATION_TESTS(1, 5, 5)
    TERNARY_OPERATION_TESTS(2, 0, 5)
    TERNARY_OPERATION_TESTS(2, 1, 5)
    TERNARY_OPERATION_TESTS(2, 2, 5)
    TERNARY_OPERATION_TESTS(2, 3, 5)
    TERNARY_OPERATION_TESTS(2, 4, 5)
    TERNARY_OPERATION_TESTS(2, 5, 0)
    TERNARY_OPERATION_TESTS(2, 5, 1)
    TERNARY_OPERATION_TESTS(2, 5, 2)
    TERNARY_OPERATION_TESTS(2, 5, 3)
    TERNARY_OPERATION_TESTS(2, 5, 4)
    TERNARY_OPERATION_TESTS(2, 5, 5)
    TERNARY_OPERATION_TESTS(3, 0, 5)
    TERNARY_OPERATION_TESTS(3, 1, 5)
    TERNARY_OPERATION_TESTS(3, 2, 5)
    TERNARY_OPERATION_TESTS(3, 3, 5)
    TERNARY_OPERATION_TESTS(3, 4, 5)
    TERNARY_OPERATION_TESTS(3, 5, 0)
    TERNARY_OPERATION_TESTS(3, 5, 1)
    TERNARY_OPERATION_TESTS(3, 5, 2)
    TERNARY_OPERATION_TESTS(3, 5, 3)
    TERNARY_OPERATION_TESTS(3, 5, 4)
    TERNARY_OPERATION_TESTS(3, 5, 5)
    TERNARY_OPERATION_TESTS(4, 0, 5)
    TERNARY_OPERATION_TESTS(4, 1, 5)
    TERNARY_OPERATION_TESTS(4, 2, 5)
    TERNARY_OPERATION_TESTS(4, 3, 5)
    TERNARY_OPERATION_TESTS(4, 4, 5)
    TERNARY_OPERATION_TESTS(4, 5, 0)
    TERNARY_OPERATION_TESTS(4, 5, 1)
    TERNARY_OPERATION_TESTS(4, 5, 2)
    TERNARY_OPERATION_TESTS(4, 5, 3)
    TERNARY_OPERATION_TESTS(4, 5, 4)
    TERNARY_OPERATION_TESTS(4, 5, 5)
    TERNARY_OPERATION_TESTS(5, 0, 0)
    TERNARY_OPERATION_TESTS(5, 0, 1)
    TERNARY_OPERATION_TESTS(5, 0, 2)
    TERNARY_OPERATION_TESTS(5, 0, 3)
    TERNARY_OPERATION_TESTS(5, 0, 4)
    TERNARY_OPERATION_TESTS(5, 0, 5)
    TERNARY_OPERATION_TESTS(5, 1, 0)
    TERNARY_OPERATION_TESTS(5, 1, 1)
    TERNARY_OPERATION_TESTS(5, 1, 2)
    TERNARY_OPERATION_TESTS(5, 1, 3)
    TERNARY_OPERATION_TESTS(5, 1, 4)
    TERNARY_OPERATION_TESTS(5, 1, 5)
    TERNARY_OPERATION_TESTS(5, 2, 0)
    TERNARY_OPERATION_TESTS(5, 2, 1)
    TERNARY_OPERATION_TESTS(5, 2, 2)
    TERNARY_OPERATION_TESTS(5, 2, 3)
    TERNARY_OPERATION_TESTS(5, 2, 4)
    TERNARY_OPERATION_TESTS(5, 2, 5)
    TERNARY_OPERATION_TESTS(5, 3, 0)
    TERNARY_OPERATION_TESTS(5, 3, 1)
    TERNARY_OPERATION_TESTS(5, 3, 2)
    TERNARY_OPERATION_TESTS(5, 3, 3)
    TERNARY_OPERATION_TESTS(5, 3, 4)
    TERNARY_OPERATION_TESTS(5, 3, 5)
    TERNARY_OPERATION_TESTS(5, 4, 0)
    TERNARY_OPERATION_TESTS(5, 4, 1)
    TERNARY_OPERATION_TESTS(5, 4, 2)
    TERNARY_OPERATION_TESTS(5, 4, 3)
    TERNARY_OPERATION_TESTS(5, 4, 4)
    TERNARY_OPERATION_TESTS(5, 4, 5)
    TERNARY_OPERATION_TESTS(5, 5, 0)
    TERNARY_OPERATION_TESTS(5, 5, 1)
    TERNARY_OPERATION_TESTS(5, 5, 2)
    TERNARY_OPERATION_TESTS(5, 5, 3)
    TERNARY_OPERATION_TESTS(5, 5, 4)
    TERNARY_OPERATION_TESTS(5, 5, 5)
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
