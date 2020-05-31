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

#ifndef __TBGAL_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_GEOMETRIC_PRODUCT_HPP__

namespace tbgal {

    namespace detail {

        struct gp_impl {
        private:

            template<typename ScalarType, typename MetricSpaceType>
            constexpr static decltype(auto) to_geometric_factors(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) {
                return std::make_tuple(arg.scalar(), arg.factors_in_signed_metric());
            }
            
            template<typename ScalarType, typename MetricSpaceType>
            constexpr static decltype(auto) to_geometric_factors(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) {
                auto factors_tuple = from_outer_to_geometric_factors(arg.space_ptr(), arg.factors_in_signed_metric());
                return std::make_tuple(arg.scalar() * std::get<0>(factors_tuple), std::get<1>(factors_tuple));
            }
            
            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename ColumnMatrixType>
            constexpr static void update_factors(MetricSpaceType const *, ScalarType &, FactorsMatrixType const &, DynamicSquareMatrixType const &, FactorsMatrixType const &, TwoByTwoMatrixType const &, TwoByTwoMatrixType const &, DynamicColumnMatrixType const &, ColumnMatrixType const &) {
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename ColumnMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const *space_ptr, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, FactorsMatrixType &QA, TwoByTwoMatrixType &T, TwoByTwoMatrixType &inv_T, DynamicColumnMatrixType &x, ColumnMatrixType &z, FirstScalarType const &arg1, NextTypes const &... args) {
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;

                alpha *= arg1;
                if (is_zero(alpha)) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space_ptr->dimensions(), 0);
                    return;
                }
                update_factors(space_ptr, alpha, A, MA, QA, T, inv_T, x, z, args...);
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename ColumnMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const *space_ptr, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, FactorsMatrixType &QA, TwoByTwoMatrixType &T, TwoByTwoMatrixType &inv_T, DynamicColumnMatrixType &x, ColumnMatrixType &z, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) {
                using IndexType = typename MetricSpaceType::IndexType;

                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;

                IndexType const n = space_ptr->dimensions();

                auto column = [&](auto const &F, IndexType const col) -> auto {
                    return block<MetricSpaceType::DimensionsAtCompileTime, 1>(F, 0, col, n, 1);
                };

                auto inner_product_em = [&](auto const &F1, IndexType const col1, auto const &F2, IndexType const col2) -> auto {
                    return coeff(prod(transpose(column(F1, col1)), column(F2, col2)), 0, 0);
                };

                IndexType s = arg1.factors_count();

                if (s != 0) {
                    auto B_tuple = to_geometric_factors(arg1);

                    alpha *= std::get<0>(B_tuple);
                    if (is_zero(alpha)) {
                        A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                        return;
                    }

                    auto const &B = std::get<1>(B_tuple);
                    IndexType r = cols(A);

                    if (r != 0) {
                        for (IndexType k = 0; k != s; ++k) {
                            auto b = column(B, k);
                            auto b_metric = apply_signed_metric(space_ptr, b);

                            x = prod(transpose(QA), b_metric); // Projection of b onto A using reciprocal frame vectors (thus, b = A . x)
                            z = subtract(b, prod(A, x));

                            if (is_zero(coeff(prod(transpose(z), z), 0, 0))) {
                                if (r > 1) {
                                    for (IndexType i = 0; i != (r - 1); ++i) {
                                        auto xi1 = coeff(x, i, 0);
                                        if (!is_zero(xi1)) {
                                            auto xi2 = coeff(x, i + 1, 0);

                                            auto mu11 = coeff(MA, i, i);
                                            auto mu12 = coeff(MA, i, i + 1);
                                            auto mu22 = coeff(MA, i + 1, i + 1);

                                            auto mu11_em = inner_product_em(A, i, A, i);
                                            auto mu12_em = inner_product_em(A, i, A, i + 1);
                                            auto mu22_em = inner_product_em(A, i + 1, A, i + 1);

                                            auto delta1 = 2 * mu12 * xi1 + mu22 * xi2;
                                            auto delta2 = -mu11 * xi1;
                                            
                                            auto rnorm_r1_em = sqrt(delta1 * delta1 * mu11_em + 2 * delta1 * delta2 * mu12_em + delta2 * delta2 * mu22_em); // Reverse norm under Euclidean metric
                                            auto rnorm_r2_sqr = xi1 * xi1 * mu11 + 2 * xi1 * xi2 * mu12 + xi2 * xi2 * mu22;

                                            bool r2_is_null = is_zero(rnorm_r2_sqr);

                                            auto gamma1 = 1 / rnorm_r1_em;
                                            auto gamma2 = r2_is_null ? rnorm_r1_em : (rnorm_r1_em / rnorm_r2_sqr);

                                            coeff(T, 0, 0) = gamma1 * delta1;
                                            coeff(T, 1, 0) = gamma1 * delta2;
                                            coeff(T, 0, 1) = gamma2 * xi1;
                                            coeff(T, 1, 1) = gamma2 * xi2;

                                            assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, T), A, 0, i, n, 2);

                                            assign_block<2, Dynamic>(prod_block<2, Dynamic>(transpose(T), MA, i, 0, 2, r), MA, i, 0, 2, r);
                                            assign_block<Dynamic, 2>(prod_block<Dynamic, 2>(MA, 0, i, r, 2, T), MA, 0, i, r, 2);

                                            auto inv_det_T = 1 / (coeff(T, 0, 0) * coeff(T, 1, 1) - coeff(T, 0, 1) * coeff(T, 1, 0));
                                            
                                            coeff(inv_T, 0, 0) = coeff(T, 1, 1) * inv_det_T;
                                            coeff(inv_T, 1, 0) = -coeff(T, 1, 0) * inv_det_T;
                                            coeff(inv_T, 0, 1) = -coeff(T, 0, 1) * inv_det_T;
                                            coeff(inv_T, 1, 1) = coeff(T, 0, 0) * inv_det_T;

                                            assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(QA, 0, i, n, 2, transpose(inv_T)), QA, 0, i, n, 2);
                                            
                                            assign_block<2, 1>(prod_block<2, 1>(inv_T, x, i, 0, 2, 1), x, i, 0, 2, 1);

                                            if (r2_is_null) {
                                                assign_block<Dynamic, 1>(prod_block<DimensionsAtCompileTime, 1>(transpose(A), A, 0, i + 1, n, 1), MA, 0, i + 1, r, 1);
                                                assign_block<1, Dynamic>(transpose(MA), i + 1, 0, MA, i + 1, 0, 1, r);

                                                assign_block<DimensionsAtCompileTime, 1>(A, 0, i + 1, QA, 0, i + 1, n, 1);
                                                assign_block<DimensionsAtCompileTime, 1>(apply_signed_metric(space_ptr, column(A, i + 1)), A, 0, i + 1, n, 1);
                                            }
                                        }
                                    }
                                }

                                alpha *= coeff(prod(transpose(column(A, r - 1)), b_metric), 0, 0); // Inner product
                                if (is_zero(alpha)) {
                                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                                    return;
                                }

                                conservative_resize(A, n, r - 1);
                                conservative_resize(MA, r - 1, r - 1);
                                conservative_resize(QA, n, r - 1);
                                --r;
                            }
                            else {
                                conservative_resize(A, n, r + 1);
                                assign_block<DimensionsAtCompileTime, 1>(B, 0, k, A, 0, r, n, 1);

                                conservative_resize(MA, r + 1, r + 1);
                                assign_block<Dynamic, 1>(prod(transpose(A), b_metric), MA, 0, r, r + 1, 1);
                                assign_block<1, Dynamic>(transpose(MA), r, 0, MA, r, 0, 1, r);

                                QA = prod(A, inverse(MA)); //TODO Use blockwise inversion

                                ++r;
                            }
                        }
                    }
                    else {
                        A = B;
                        MA = prod(transpose(A), apply_signed_metric(space_ptr, A)); // Metric matrix of A
                        QA = prod(A, inverse(MA)); // Reciprocal frame vectors of A
                    }
                }
                else {
                    alpha *= arg1.scalar();
                    if (is_zero(alpha)) {
                        A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                        return;
                    }
                }
                update_factors(space_ptr, alpha, A, MA, QA, T, inv_T, x, z, args...);
            }

        public:

            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) {
                using ResultingScalarType = common_scalar_type_t<Types...>;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = GeometricProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;

                constexpr DefaultIndexType DimensionsAtCompileTime = ResultingMetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = ResultingMetricSpaceType::MaxDimensionsAtCompileTime;

                using DynamicColumnMatrixType = matrix_type_t<ResultingScalarType, Dynamic, 1, MaxDimensionsAtCompileTime, 1>;
                using DynamicSquareMatrixType = matrix_type_t<ResultingScalarType, Dynamic, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using ColumnMatrixType = matrix_type_t<ResultingScalarType, DimensionsAtCompileTime, 1, MaxDimensionsAtCompileTime, 1>;
                using FactorsMatrixType = matrix_type_t<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using TwoByTwoMatrixType = matrix_type_t<ResultingScalarType, 2, 2, 2, 2>;

                auto const *space_ptr = detail::space_ptr(args...);

                ResultingScalarType alpha = 1;
                FactorsMatrixType A = make_matrix<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space_ptr->dimensions(), 0);
                
                FactorsMatrixType QA;
                TwoByTwoMatrixType T, inv_T;
                DynamicSquareMatrixType MA;
                DynamicColumnMatrixType x;
                ColumnMatrixType z;
                update_factors(space_ptr, alpha, A, MA, QA, T, inv_T, x, z, args...);

                return ResultingFactoredMultivectorType(space_ptr, alpha, A);
            }
        };

        template<typename... Types>
        struct gp_impl_reduces_to_op_impl;

        template<typename... Types>
        constexpr bool gp_impl_reduces_to_op_impl_v = gp_impl_reduces_to_op_impl<Types...>::value;

        template<>
        struct gp_impl_reduces_to_op_impl<> {
            constexpr static bool value = true;
        };

        template<typename FirstScalarType, typename... NextTypes>
        struct gp_impl_reduces_to_op_impl<FirstScalarType, NextTypes...> :
            gp_impl_reduces_to_op_impl<NextTypes...> {
        };

        template<typename FirstScalarType, typename FirstMetricSpaceType, typename... NextTypes>
        struct gp_impl_reduces_to_op_impl<FactoredMultivector<FirstScalarType, GeometricProduct<FirstMetricSpaceType> >, NextTypes...> {
            constexpr static bool value = false;
        };

        template<typename FirstScalarType, typename FirstMetricSpaceType, typename... NextTypes>
        struct gp_impl_reduces_to_op_impl<FactoredMultivector<FirstScalarType, OuterProduct<FirstMetricSpaceType> >, NextTypes...> {
            constexpr static bool value = !is_any_v<std::true_type, is_multivector_t<NextTypes>...>;
        };

    }

    template<typename FirstType, typename... NextTypes>
    constexpr decltype(auto) gp(FirstType const &arg1, NextTypes const &... args) {
        return std::conditional_t<detail::gp_impl_reduces_to_op_impl_v<FirstType, NextTypes...>, detail::op_impl, detail::gp_impl>::eval(arg1, args...);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator*(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        return gp(arg1, arg2);
    }

    template<typename FirstScsalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType>, int> >
    constexpr decltype(auto) operator*(FactoredMultivector<FirstScsalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) {
        return gp(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType>, int> >
    constexpr decltype(auto) operator*(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) {
        return gp(arg1, arg2);
    }

}

#endif // __TBGAL_GEOMETRIC_PRODUCT_HPP__
