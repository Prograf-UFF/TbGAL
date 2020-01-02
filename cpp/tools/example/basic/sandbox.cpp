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

#include "../../../include/tbgal/using_Eigen.hpp" //TODO [DEBUG]
#include "../../../include/tbgal/assuming_Conformal1.hpp" //TODO [DEBUG]

namespace tbgal_mdl = tbgal::Conformal1;

#include <gatl/ga1m.hpp>

namespace gatl_mdl = ga1m;

template<typename ScalarType, std::size_t K>
using factors_row_t = std::array<ScalarType, K>;

template<typename ScalarType, std::size_t N, std::size_t K>
using factors_matrix_t = std::array<factors_row_t<ScalarType, K>, N>;

template<typename ScalarType, std::size_t N, std::size_t K, std::size_t... Indices>
constexpr decltype(auto) _make_tbgal_vector_impl(factors_matrix_t<ScalarType, N, K> const &factors, std::size_t k, std::index_sequence<Indices...>) noexcept {
    return tbgal_mdl::vector(factors[Indices][k]...);
}

template<typename ScalarType, std::size_t N, std::size_t K>
constexpr decltype(auto) make_tbgal_vector(factors_matrix_t<ScalarType, N, K> const &factors, std::size_t k) noexcept {
    return _make_tbgal_vector_impl(factors, k, std::make_index_sequence<N>{});
}

template<typename ScalarType, std::size_t N, std::size_t K, std::size_t... Indices>
constexpr decltype(auto) _make_gatl_vector_impl(factors_matrix_t<ScalarType, N, K> const &factors, std::size_t k, std::index_sequence<Indices...>) noexcept {
    return gatl_mdl::vector(factors[Indices][k]...);
}

template<typename ScalarType, std::size_t N, std::size_t K>
constexpr decltype(auto) make_gatl_vector(factors_matrix_t<ScalarType, N, K> const &factors, std::size_t k) noexcept {
    return _make_gatl_vector_impl(factors, k, std::make_index_sequence<N>{});
}

template<typename ScalarType, std::size_t N, std::size_t K, std::size_t... Indices>
constexpr decltype(auto) _make_tbgal_multivector_using_GeometricProduct_impl(factors_matrix_t<ScalarType, N, K> const &factors, std::index_sequence<Indices...>) noexcept {
    return tbgal::gp(make_tbgal_vector(factors, Indices)...);
}
 
template<typename ScalarType, std::size_t N, std::size_t K, typename Indices = std::make_index_sequence<K> >
constexpr decltype(auto) make_tbgal_multivector_using_GeometricProduct(factors_matrix_t<ScalarType, N, K> const &factors) noexcept {
    return _make_tbgal_multivector_using_GeometricProduct_impl(factors, Indices{});
}

template<typename Type>
constexpr decltype(auto) gatl_GeometricProduct(Type const &arg) noexcept {
    return arg;
}

template<typename FirstType, typename... NextTypes>
constexpr decltype(auto) gatl_GeometricProduct(FirstType const &arg1, NextTypes const &... args) noexcept {
    return gatl_mdl::gp(arg1, gatl_GeometricProduct(args...));
}

template<typename ScalarType, std::size_t N, std::size_t K, std::size_t... Indices>
constexpr decltype(auto) _make_gatl_multivector_using_GeometricProduct_impl(factors_matrix_t<ScalarType, N, K> const &factors, std::index_sequence<Indices...>) noexcept {
    return gatl_GeometricProduct(make_gatl_vector(factors, Indices)...);
}
 
template<typename ScalarType, std::size_t N, std::size_t K, typename Indices = std::make_index_sequence<K> >
constexpr decltype(auto) make_gatl_multivector_using_GeometricProduct(factors_matrix_t<ScalarType, N, K> const &factors) noexcept {
    return _make_gatl_multivector_using_GeometricProduct_impl(factors, Indices{});
}

template<typename Type>
constexpr decltype(auto) gatl_OuterProduct(Type const &arg) noexcept {
    return arg;
}

template<typename FirstType, typename... NextTypes>
constexpr decltype(auto) gatl_OuterProduct(FirstType const &arg1, NextTypes const &... args) noexcept {
    return gatl_mdl::op(arg1, gatl_OuterProduct(args...));
}

template<typename ScalarType, typename MetricSpaceType>
decltype(auto) from_tbgal_to_gatl(tbgal::FactoredMultivector<ScalarType, tbgal::GeometricProduct<MetricSpaceType> > const &arg) noexcept {
    using IndexType = typename tbgal::FactoredMultivector<ScalarType, tbgal::GeometricProduct<MetricSpaceType> >::IndexType;

    constexpr std::size_t N = MetricSpaceType::DimensionsAtCompileTime;

    ga::full_multivector_t<ScalarType, N> result;
    ga::trivial_copy(arg.scalar(), result);
    auto factors = arg.factors_in_actual_metric();
    for (IndexType col = 0; col != arg.factors_count(); ++col) {
        factors_matrix_t<ScalarType, N, 1> factor;
        for (std::size_t row = 0; row != N; ++row) {
            factor[row][0] = tbgal::detail::coeff(factors, row, col);
        }
        ga::trivial_copy(gatl_GeometricProduct(result, make_gatl_vector(factor, 0)), result);
    }
    return result;
}

template<typename ScalarType, typename MetricSpaceType>
decltype(auto) from_tbgal_to_gatl(tbgal::FactoredMultivector<ScalarType, tbgal::OuterProduct<MetricSpaceType> > const &arg) noexcept {
    using IndexType = typename tbgal::FactoredMultivector<ScalarType, tbgal::OuterProduct<MetricSpaceType> >::IndexType;

    constexpr std::size_t N = MetricSpaceType::DimensionsAtCompileTime;

    ga::full_multivector_t<ScalarType, N> result;
    ga::trivial_copy(arg.scalar(), result);
    auto factors = arg.factors_in_actual_metric();
    for (IndexType col = 0; col != arg.factors_count(); ++col) {
        factors_matrix_t<ScalarType, N, 1> factor;
        for (std::size_t row = 0; row != N; ++row) {
            factor[row][0] = tbgal::detail::coeff(factors, row, col);
        }
        ga::trivial_copy(gatl_OuterProduct(result, make_gatl_vector(factor, 0)), result);
    }
    return result;
}

template<typename ScalarType, typename = std::enable_if_t<!tbgal::is_multivector_v<ScalarType> > >
decltype(auto) from_tbgal_to_gatl(ScalarType const &arg) noexcept {
    return arg;
}

int main(int argc, char *argv[]) {
    using namespace gatl_mdl;

    constexpr std::size_t N = tbgal_mdl::MetricSpaceType::DimensionsAtCompileTime;

    // b0
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;

    auto b0 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{-0.68847728521057272033}, factors_row_t<double, 1>{-0.27595604617189428698}, factors_row_t<double, 1>{-0.67070655903327403013}});
    auto b0_ = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{-0.68847728521057272033}, factors_row_t<double, 1>{-0.27595604617189428698}, factors_row_t<double, 1>{-0.67070655903327403013}});

    auto A0_before_0 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{-0.37622305113100784624, 0, 0}, factors_row_t<double, 3>{-0.62867855215218748643, 0.12707459166396728456, 0.70710678118654746172}, factors_row_t<double, 3>{-0.680602302274613713, -0.99189316367915147943, -0.70710678118654746172}});
    auto A0_after_0 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{0.19888768175555027184, -0.68847728521057183215, 0}, factors_row_t<double, 3>{0.32442418766603348113, -1.3634181472902262566, 0.70710678118654746172}, factors_row_t<double, 3>{0.42163414109225028081, 0.41675554208505749543, -0.70710678118654746172}});

    auto A0_before_1 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{0.19888768175555027184, -0.68847728521057183215, 0}, factors_row_t<double, 3>{0.32442418766603348113, -1.3634181472902262566, 0.70710678118654746172}, factors_row_t<double, 3>{0.42163414109225028081, 0.41675554208505749543, -0.70710678118654746172}});
    auto A0_after_1 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{0.19888768175555027184, 9.1891516952123897255, -0.6884772852105730534}, factors_row_t<double, 3>{0.32442418766603348113, 2.9760967521137651204, -0.2759560461718963964}, factors_row_t<double, 3>{0.42163414109225028081, 9.6590714230492604742, -0.67070655903327414116}});

    auto c0 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{-0.6884772852105730534}, factors_row_t<double, 1>{-0.2759560461718963964}, factors_row_t<double, 1>{-0.67070655903327414116}});
    
    auto A0_resized = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 2>{factors_row_t<double, 2>{0.19888768175555027184, 9.1891516952123897255}, factors_row_t<double, 2>{0.32442418766603348113, 2.9760967521137651204}, factors_row_t<double, 2>{0.42163414109225028081, 9.6590714230492604742}});

    auto A0_result = gp(sp(c0, b0), A0_resized);

    auto A0_target = gp(A0_before_0, b0);

    std::cout << "b0 = " << b0 << std::endl;
    std::cout << "b0_ = " << b0_ << std::endl;
    std::cout << std::endl;

    std::cout << "A0_before_0 = " << A0_before_0 << std::endl;
    std::cout << "A0_after_0  = " << A0_after_0 << std::endl;
    std::cout << std::endl;

    std::cout << "A0_before_1 = " << A0_before_1 << std::endl;
    std::cout << "A0_after_1  = " << A0_after_1 << std::endl;
    std::cout << std::endl;

    std::cout << "c0 = " << c0 << std::endl;
    std::cout << std::endl;

    std::cout << "A0_resized = " << A0_resized << std::endl;
    std::cout << std::endl;

    std::cout << "A0_result = " << A0_result << std::endl;
    std::cout << "A0_target = " << A0_target << std::endl;
    std::cout << std::endl;

    // b1
    /**/
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    auto b1 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{0}, factors_row_t<double, 1>{-0.32485272479837540294}, factors_row_t<double, 1>{-0.94576461510836351554}});

    auto A1_resized = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{0.19888768175555027184, 9.1891516952123897255, 0}, factors_row_t<double, 3>{0.32442418766603348113, 2.9760967521137651204, -0.32485272479837540294}, factors_row_t<double, 3>{0.42163414109225028081, 9.6590714230492604742, -0.94576461510836351554}});
    
    auto A1_result = gp(sp(c0, b0), A1_resized);

    auto A1_target = gp(A0_target, b1);

    std::cout << "A1_resized = " << A1_resized << std::endl;
    std::cout << std::endl;

    std::cout << "A1_result = " << A1_result << std::endl;
    std::cout << "A1_target = " << A1_target << std::endl;
    std::cout << std::endl;
    /**/

    // b2
    /**/
    std::cout << "-----------------------------------------------------------" << std::endl;
    std::cout << std::endl;
    
    auto b2 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{0}, factors_row_t<double, 1>{0.70710678118654746172}, factors_row_t<double, 1>{-0.70710678118654746172}});
    auto b2_ = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{0}, factors_row_t<double, 1>{0.70710678118654746172}, factors_row_t<double, 1>{0.70710678118654746172}});

    auto A2_before_0 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{0.19888768175555027184, 9.1891516952123897255, 0}, factors_row_t<double, 3>{0.32442418766603348113, 2.9760967521137651204, -0.32485272479837540294}, factors_row_t<double, 3>{0.42163414109225028081, 9.6590714230492604742, -0.94576461510836351554}});
    auto A2_after_0 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{-0.51139345980684169923, -3.3306690738754696213e-16, 0}, factors_row_t<double, 3>{-1.8174080096434812592, 0.64549617335332853951, -0.32485272479837540294}, factors_row_t<double, 3>{-1.8879870614624163494, 0.52773585905516440242, -0.94576461510836351554}});

    auto A2_before_1 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{-0.51139345980684169923, -3.3306690738754696213e-16, 0}, factors_row_t<double, 3>{-1.8174080096434812592, 0.64549617335332853951, -0.32485272479837540294}, factors_row_t<double, 3>{-1.8879870614624163494, 0.52773585905516440242, -0.94576461510836351554}});
    auto A2_after_1 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 3>{factors_row_t<double, 3>{-0.51139345980684169923, -2.4263048870542632274e-16, -3.3306690738754696213e-16}, factors_row_t<double, 3>{-1.8174080096434812592, 0.51510870663168084604, 0.7071067811865479058}, factors_row_t<double, 3>{-1.8879870614624163494, 0.51510870663168084604, 0.70710678118654768376}});

    auto c2 = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 1>{factors_row_t<double, 1>{-3.3306690738754696213e-16}, factors_row_t<double, 1>{0.7071067811865479058}, factors_row_t<double, 1>{0.70710678118654768376}});
    
    auto A2_resized = make_gatl_multivector_using_GeometricProduct(factors_matrix_t<double, 3, 2>{factors_row_t<double, 2>{-0.51139345980684169923, -2.4263048870542632274e-16}, factors_row_t<double, 2>{-1.8174080096434812592, 0.51510870663168084604}, factors_row_t<double, 2>{-1.8879870614624163494, 0.51510870663168084604}});

    auto A2_result = gp(sp(c2, b2), A2_resized);

    auto A2_target = gp(A2_before_0, b2);

    std::cout << "b2 = " << b2 << std::endl;
    std::cout << "b2_ = " << b2_ << std::endl;
    std::cout << std::endl;

    std::cout << "A2_before_0 = " << A2_before_0 << std::endl;
    std::cout << "A2_after_0  = " << A2_after_0 << std::endl;
    std::cout << std::endl;

    std::cout << "A2_before_1 = " << A2_before_1 << std::endl;
    std::cout << "A2_after_1  = " << A2_after_1 << std::endl;
    std::cout << std::endl;

    std::cout << "c2 = " << c2 << std::endl;
    std::cout << std::endl;

    std::cout << "A2_resized = " << A2_resized << std::endl;
    std::cout << std::endl;

    std::cout << "A2_result = " << A2_result << std::endl;
    std::cout << "A2_target = " << A2_target << std::endl;
    std::cout << std::endl;
    /**/

    return EXIT_SUCCESS;
}
