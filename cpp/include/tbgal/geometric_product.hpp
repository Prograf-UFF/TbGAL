#ifndef __TBGAL_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_GEOMETRIC_PRODUCT_HPP__

namespace tbgal {

    namespace detail {

        struct gp_impl {
        private:

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType>
            constexpr static void update_factors(MetricSpaceType const &, ScalarType &, FactorsMatrixType const &, DynamicSquareMatrixType const &, TwoByTwoMatrixType const &, DynamicColumnMatrixType const &, DynamicColumnMatrixType const &) noexcept {
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, TwoByTwoMatrixType &R, DynamicColumnMatrixType &u, DynamicColumnMatrixType &x, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;

                alpha *= arg1;
                if (is_zero(alpha)) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);
                    return;
                }

                update_factors(space, alpha, A, MA, R, u, x, args...);
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, TwoByTwoMatrixType &R, DynamicColumnMatrixType &u, DynamicColumnMatrixType &x, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                using IndexType = typename MetricSpaceType::IndexType;

                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;

                IndexType const n = space.dimensions();

                auto column = [&](auto const &F, IndexType const col) -> auto {
                    return block_view<DimensionsAtCompileTime, 1>(F, 0, col, n, 1);
                };

                auto inner_product = [&](auto const &F1, IndexType const col1, auto const &F2, IndexType const col2) -> auto {
                    return coeff(prod(transpose(column(F1, col1)), apply_signed_metric(space, column(F2, col2))), 0, 0);
                };

                auto inner_product_em = [&](auto const &F1, IndexType const col1, auto const &F2, IndexType const col2) -> auto {
                    return coeff(prod(transpose(column(F1, col1)), column(F2, col2)), 0, 0);
                };

                alpha *= arg1.scalar();
                if (is_zero(alpha)) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                    return;
                }

                auto const &B = arg1.factors_in_signed_metric();
                IndexType r = cols(A);
                IndexType s = cols(B);

                if (r != 0) {
                    if (s != 0) {
                        for (IndexType k = 0; k != s; ++k) {
                            auto qr_tuple_A = qr_orthogonal_matrix(A);
                            assert(r == std::get<1>(qr_tuple_A));

                            u = prod_block<Dynamic, DimensionsAtCompileTime, 1>(transpose(std::get<0>(qr_tuple_A)), 0, 0, r, B, 0, k, n, 1);

                            if (is_zero(coeff(prod(transpose(u), u), 0, 0) - 1)) {
                                if (r > 1) {
                                    for (IndexType i = 0; i != (r - 1); ++i) {
                                        MA = prod(transpose(A), apply_signed_metric(space, A));
                                        x = prod(prod(inverse(MA), transpose(A)), apply_signed_metric(space, column(B, k)));

                                        auto xi1 = coeff(x, i, 0);
                                        if (!is_zero(xi1)) {
                                            auto xi2 = coeff(x, i + 1, 0);

                                            auto mu11 = coeff(MA, i, i);
                                            auto mu12 = coeff(MA, i, i + 1);
                                            auto mu22 = coeff(MA, i + 1, i + 1);

                                            auto delta1 = 2 * mu12 * xi1 + mu22 * xi2;
                                            auto delta2 = -mu11 * xi1;
                                            auto gamma1 = 1 / sqrt(delta1 * delta1 * inner_product_em(A, i, A, i) + 2 * delta1 * delta2 * inner_product_em(A, i, A, i + 1) + delta2 * delta2); // Normalization under Euclidean metric

                                            coeff(R, 0, 0) = gamma1 * delta1;
                                            coeff(R, 1, 0) = gamma1 * delta2;
                                            coeff(R, 0, 1) = xi1;
                                            coeff(R, 1, 1) = xi2;

                                            assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, R), A, 0, i, n, 2);
                                        }
                                    }
                                }

                                alpha *= inner_product(A, r - 1, B, k);
                                conservative_resize(A, n, r - 1);
                                --r;
                            }
                            else {
                                conservative_resize(A, n, r + 1);
                                assign_block<DimensionsAtCompileTime, 1>(B, 0, k, A, 0, r, n, 1);
                                ++r;
                            }
                        }
                    }
                }
                else {
                    A = B;
                }
                update_factors(space, alpha, A, MA, R, u, x, args...);
            }

        public:

            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) noexcept {
                using ResultingScalarType = common_scalar_type_t<Types...>;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = GeometricProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;

                constexpr DefaultIndexType DimensionsAtCompileTime = ResultingMetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = ResultingMetricSpaceType::MaxDimensionsAtCompileTime;

                using DynamicColumnMatrixType = matrix_type_t<ResultingScalarType, Dynamic, 1, MaxDimensionsAtCompileTime, 1>;
                using DynamicSquareMatrixType = matrix_type_t<ResultingScalarType, Dynamic, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using TwoByTwoMatrixType = matrix_type_t<ResultingScalarType, 2, 2, 2, 2>;

                auto const &space = *space_ptr(args...);

                ResultingScalarType alpha = 1;
                auto A = make_matrix<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);
                
                TwoByTwoMatrixType R;
                DynamicSquareMatrixType MA;
                DynamicColumnMatrixType u, x;
                update_factors(space, alpha, A, MA, R, u, x, args...);

                return ResultingFactoredMultivectorType(space, alpha, A);
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
    constexpr decltype(auto) gp(FirstType const &arg1, NextTypes const &... args) noexcept {
        return std::conditional_t<detail::gp_impl_reduces_to_op_impl_v<FirstType, NextTypes...>, detail::op_impl<detail::is_any_v<std::true_type, is_multivector_t<FirstType>, is_multivector_t<NextTypes>...> >, detail::gp_impl>::eval(arg1, args...);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator*(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return gp(arg1, arg2);
    }

    template<typename FirstScsalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator*(FactoredMultivector<FirstScsalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return gp(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator*(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return gp(arg1, arg2);
    }

}

#endif // __TBGAL_GEOMETRIC_PRODUCT_HPP__
