#ifndef __TBGAL_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_GEOMETRIC_PRODUCT_HPP__

namespace tbgal {

    namespace detail {

        //TODO {DEBUG}
        template<typename MatrixType>
        void print_entry(int level, std::string const &name, MatrixType const &arg) {
            for (int i = 0; i != level; ++i) std::cout << "    ";
            std::cout << name << std::endl;
            for (int i = 0; i != level; ++i) std::cout << "    ";
            std::cout << "----------------------------------------------" << std::endl;
            for (int i = 0; i != level; ++i) std::cout << "    ";
            print_matrix(arg);
            std::cout << std::endl << std::endl;
        }

        //TODO {DEBUG}
        void print_entry(int level, std::string const &name, double const &arg) {
            for (int i = 0; i != level; ++i) std::cout << "    ";
            std::cout << name << std::endl;
            for (int i = 0; i != level; ++i) std::cout << "    ";
            std::cout << "----------------------------------------------" << std::endl;
            for (int i = 0; i != level; ++i) std::cout << "    ";
            std::cout << arg;
            std::cout << std::endl << std::endl;
        }

        struct gp_impl {
        private:

            template<typename ScalarType, typename MetricSpaceType>
            constexpr static decltype(auto) to_geometric_factors(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
                return std::make_tuple(arg.scalar(), arg.factors_in_signed_metric());
            }
            
            template<typename ScalarType, typename MetricSpaceType>
            constexpr static decltype(auto) to_geometric_factors(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
                auto factors_tuple = from_outer_to_geometric_factors(arg.space(), arg.factors_in_signed_metric());
                return std::make_tuple(arg.scalar() * std::get<0>(factors_tuple), std::get<1>(factors_tuple));
            }
            
            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType>
            constexpr static void update_factors(MetricSpaceType const &, ScalarType &, FactorsMatrixType const &, DynamicSquareMatrixType const &, TwoByTwoMatrixType const &, DynamicColumnMatrixType const &, DynamicColumnMatrixType const &) noexcept {
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, TwoByTwoMatrixType &T, DynamicColumnMatrixType &u, DynamicColumnMatrixType &x, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;

                alpha *= arg1;
                if (is_zero(alpha)) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);
                    return;
                }
                update_factors(space, alpha, A, MA, T, u, x, args...);
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename TwoByTwoMatrixType, typename DynamicColumnMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, TwoByTwoMatrixType &T, DynamicColumnMatrixType &u, DynamicColumnMatrixType &x, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                using IndexType = typename MetricSpaceType::IndexType;

                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;

                IndexType const n = space.dimensions();

                auto column = [&](auto const &F, IndexType const col) -> auto {
                    return block<DimensionsAtCompileTime, 1>(F, 0, col, n, 1);
                };

                auto inner_product = [&](auto const &F1, IndexType const col1, auto const &F2, IndexType const col2) -> auto {
                    return coeff(prod(transpose(column(F1, col1)), apply_signed_metric(space, column(F2, col2))), 0, 0);
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
                        //print_entry(0, "A", A);
                        //print_entry(0, "B", B);
                        for (IndexType k = 0; k != s; ++k) {
                            auto qr_tuple_A = qr_orthogonal_matrix(A);
                            assert(r == std::get<1>(qr_tuple_A));

                            auto const &OA = std::get<0>(qr_tuple_A);

                            u = prod_block<Dynamic, DimensionsAtCompileTime, 1>(transpose(OA), 0, 0, r, B, 0, k, n, 1);

                            if (is_zero(coeff(prod(transpose(u), u), 0, 0) - 1)) {
                                if (r > 1) {
                                    auto b = column(B, k); //TODO Embutir no código

                                    //print_entry(0, "b" + std::to_string(k), column(B, k));
                                    //print_entry(0, "b'" + std::to_string(k), b);

                                    for (IndexType i = 0; i != (r - 1); ++i) {
                                        //print_entry(1, "i = " + std::to_string(i) + " | A - before", A);
                                        MA = prod(transpose(A), apply_signed_metric(space, A)); // Metric matrix of A
                                        auto QA = prod(A, inverse(MA)); // Reciprocal frame vectors of A
                                        x = prod(transpose(QA), apply_signed_metric(space, b)); // Projection of b onto A using reciprocal frame vectors (thus, b = A . x)
                                        //print_entry(1, "i = " + std::to_string(i) + " | x - before", x);

                                        //TODO {DEBUG}
                                        //auto x_rec = prod(transpose(A), apply_signed_metric(space, b));
                                        //print_entry(2, "x", x);
                                        //print_entry(2, "x_rec", x_rec);
                                        //print_entry(2, "b", prod(A, x));
                                        //print_entry(2, "b_rec", prod(QA, x_rec));

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
                                            auto gamma1 = 1 / rnorm_r1_em;

                                            coeff(T, 0, 0) = gamma1 * delta1;
                                            coeff(T, 1, 0) = gamma1 * delta2;

                                            auto rnorm_r2_sqr = xi1 * xi1 * mu11 + 2 * xi1 * xi2 * mu12 + xi2 * xi2 * mu22;
                                            if (is_zero(rnorm_r2_sqr)) {
                                                auto gamma2 = rnorm_r1_em;

                                                coeff(T, 0, 1) = gamma2 * xi1;
                                                coeff(T, 1, 1) = gamma2 * xi2;

                                                assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, T), A, 0, i, n, 2);
                                                assign_block<DimensionsAtCompileTime, 1>(apply_signed_metric(space, column(A, i + 1)), A, 0, i + 1, n, 1);
                                            }
                                            else {
                                                auto gamma2 = rnorm_r1_em / rnorm_r2_sqr;

                                                coeff(T, 0, 1) = gamma2 * xi1;
                                                coeff(T, 1, 1) = gamma2 * xi2;

                                                assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, T), A, 0, i, n, 2);
                                            }
                                        }
                                        //print_entry(1, "i = " + std::to_string(i) + " | A - after", A);
                                        //print_entry(1, "i = " + std::to_string(i) + " | x - after", evaluate(prod(prod(inverse(prod(transpose(A), apply_signed_metric(space, A))), transpose(A)), apply_signed_metric(space, b))));
                                    }
                                }

                                alpha *= inner_product(A, r - 1, B, k);
                                if (is_zero(alpha)) {
                                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                                    return;
                                }

                                conservative_resize(A, n, r - 1);
                                --r;
                                //print_entry(0, "A - resized", A);
                            }
                            else {
                                conservative_resize(A, n, r + 1);
                                assign_block<DimensionsAtCompileTime, 1>(B, 0, k, A, 0, r, n, 1);
                                ++r;
                            }
                        }
                    }
                    else {
                        A = B;
                    }
                }
                else {
                    alpha *= arg1.scalar();
                    if (is_zero(alpha)) {
                        A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                        return;
                    }
                }
                update_factors(space, alpha, A, MA, T, u, x, args...);
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
                
                TwoByTwoMatrixType T;
                DynamicSquareMatrixType MA;
                DynamicColumnMatrixType u, x;
                update_factors(space, alpha, A, MA, T, u, x, args...);

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
        return std::conditional_t<detail::gp_impl_reduces_to_op_impl_v<FirstType, NextTypes...>, detail::op_impl, detail::gp_impl>::eval(arg1, args...);
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
