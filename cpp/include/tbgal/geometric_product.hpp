#ifndef __TBGAL_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_GEOMETRIC_PRODUCT_HPP__

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

namespace tbgal {

    namespace detail {

        struct gp_impl {
        private:

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename DynamicColumnMatrixType, typename TwoByTwoMatrixType>
            constexpr static void update_factors(MetricSpaceType const &, ScalarType const &, FactorsMatrixType const &, DynamicSquareMatrixType const &, DynamicColumnMatrixType const &, TwoByTwoMatrixType const &) noexcept {
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename DynamicColumnMatrixType, typename TwoByTwoMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, DynamicColumnMatrixType &x, TwoByTwoMatrixType &R, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;
                constexpr ScalarType ZeroTolerance = 100 * std::numeric_limits<ScalarType>::epsilon(); //TODO Adequado?

                alpha *= arg1;
                if (abs(alpha) <= ZeroTolerance) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);
                    return;
                }

                update_factors(space, alpha, A, MA, x, R, args...);
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename DynamicSquareMatrixType, typename DynamicColumnMatrixType, typename TwoByTwoMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, DynamicSquareMatrixType &MA, DynamicColumnMatrixType &x, TwoByTwoMatrixType &R, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                using IndexType = typename MetricSpaceType::IndexType;

                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;
                constexpr ScalarType ZeroTolerance = 100 * std::numeric_limits<ScalarType>::epsilon(); //TODO Adequado?

                IndexType const n = space.dimensions();

                alpha *= arg1.scalar();
                if (abs(alpha) <= ZeroTolerance) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                    return;
                }

                auto const &B = arg1.factors_in_signed_metric();
                IndexType r = cols(A);
                IndexType s = cols(B);

                if (r != 0) {
                    if (s != 0) {
                        std::cout << std::endl;
                        std::cout << "**************************************************" << std::endl;

                        for (IndexType k = 0; k != s; ++k) {
                            //TODO Como atualizar?
                            auto qr_tuple_A = qr_orthogonal_matrix(A);
                            assert(r == std::get<1>(qr_tuple_A));

                            auto const &QA = std::get<0>(qr_tuple_A);
                            auto OA = block_view<DimensionsAtCompileTime, Dynamic>(QA, 0, 0, n, r);

                            //TODO Forcei aqui
                            //alpha *= determinant(prod_block<Dynamic, DimensionsAtCompileTime>(transpose(QA), 0, 0, r, n, A));
                            
                            std::cout << "--- A --------" << std::endl; print_matrix(A); std::cout << std::endl << std::endl;
                            std::cout << "--- QA -------" << std::endl; print_matrix(QA); std::cout << std::endl << std::endl;
                            std::cout << "--- OA -------" << std::endl; print_matrix(OA); std::cout << std::endl << std::endl;

                            std::cout << "    --- b" << k << " -------" << std::endl << "    "; print_matrix(evaluate(block_view<DimensionsAtCompileTime, 1>(B, 0, k, n, 1))); std::cout << std::endl << std::endl;

                            auto u = prod_block<DimensionsAtCompileTime, 1>(transpose(OA), B, 0, k, n, 1); //TODO Explorar o fato de que OA é block de QA

                            std::cout << "    --- u" << k << " -------" << std::endl << "    "; print_matrix(u); std::cout << std::endl << std::endl;

                            if (abs(coeff(prod(transpose(u), u), 0, 0) - 1) <= ZeroTolerance) {
                                /**
                                if (r > 2) {
                                    x = prod_block<DimensionsAtCompileTime, 1>(transpose(A), B, 0, k, n, 1); //TODO A métrica entra aqui
                                    std::cout << "    --- x -before-" << std::endl << "    "; print_matrix(evaluate(x)); std::cout << std::endl << std::endl;

                                    for (IndexType i = 0; i != (r - 2); ++i) {
                                        auto xi1 = coeff(x, i, 0);
                                        if (abs(xi1) > ZeroTolerance) {
                                            auto xi2 = coeff(x, i + 1, 0);

                                            auto dot11 = dot_product_column(A, i, A, i); //TODO A métrica vai aqui
                                            auto dot12 = dot_product_column(A, i, A, i + 1); //TODO A métrica vai aqui
                                            auto dot22 = dot_product_column(A, i + 1, A, i + 1); //TODO A métrica vai aqui

                                            auto inv_norm = 1 / sqrt(xi1 * xi1 - 2 * xi1 * xi2 * dot_product_column(A, i, A, i + 1) + xi2 * xi2); // Normalization under Eiclidean metric
                                            auto rho1 = -inv_norm * xi2;
                                            auto rho2 = inv_norm * xi1;

                                            auto inv_norm_sqr_r1 = 1 / (rho1 * rho1 * dot11 + 2 * rho1 * rho2 * dot12 + rho2 * rho2 * dot22);

                                            coeff(R, 0, 0) = rho1;
                                            coeff(R, 1, 0) = rho2;
                                            coeff(R, 0, 1) = -inv_norm_sqr_r1 * rho2 * dot22;
                                            coeff(R, 1, 1) = inv_norm_sqr_r1 * (rho1 + 2 * rho2 * dot12);

                                            assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, R), A, 0, i, n, 2);
                                            std::cout << "    --- A -updated" << std::endl << "    "; print_matrix(A); std::cout << std::endl << std::endl;

                                            //TODO Como Atualizar x de forma mais rápida?
                                            x = prod_block<DimensionsAtCompileTime, 1>(transpose(A), B, 0, k, n, 1); //TODO A métrica entra aqui, em dois lugares
                                            std::cout << "    --- x -updated" << std::endl << "    "; print_matrix(x); std::cout << std::endl << std::endl;
                                        }
                                    }
                                }

                                if (r > 1) {
                                    x = prod(inverse(prod(transpose(A), A)), prod_block<DimensionsAtCompileTime, 1>(transpose(A), B, 0, k, n, 1)); //TODO A métrica entra aqui, em dois lugares
                                    std::cout << "    --- x -before-" << std::endl << "    "; print_matrix(evaluate(x)); std::cout << std::endl << std::endl;

                                    IndexType i = r - 2;
                                    auto xi1 = coeff(x, i, 0);
                                    if (abs(xi1) > ZeroTolerance) {
                                        auto xi2 = coeff(x, i + 1, 0);

                                        auto dot11 = dot_product_column(A, i, A, i); //TODO A métrica vai aqui
                                        auto dot12 = dot_product_column(A, i, A, i + 1); //TODO A métrica vai aqui
                                        auto dot22 = dot_product_column(A, i + 1, A, i + 1); //TODO A métrica vai aqui

                                        auto delta1 = 2 * xi1 + dot12 + xi2 * dot22;
                                        auto delta2 = -xi1 * dot11;
                                        auto inv_norm = 1 / sqrt(delta1 * delta1 + 2 * delta1 * delta2 * dot_product_column(A, i, A, i + 1) + delta2 * delta2); // Normalization under Eiclidean metric

                                        coeff(R, 0, 0) = inv_norm * delta1;
                                        coeff(R, 1, 0) = inv_norm * delta2;
                                        coeff(R, 0, 1) = xi1;
                                        coeff(R, 1, 1) = xi2;

                                        assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, R), A, 0, i, n, 2);
                                        std::cout << "    --- A -updated" << std::endl << "    "; print_matrix(A); std::cout << std::endl << std::endl;

                                        //TODO Como Atualizar x de forma mais rápida?
                                        x = prod(inverse(prod(transpose(A), A)), prod_block<DimensionsAtCompileTime, 1>(transpose(A), B, 0, k, n, 1)); //TODO A métrica entra aqui, em dois lugares
                                        std::cout << "    --- x -updated" << std::endl << "    "; print_matrix(x); std::cout << std::endl << std::endl;
                                    }
                                }
                                /*/
                                if (r > 1) {
                                    std::cout << "    --- A -before-" << std::endl << "    "; print_matrix(A); std::cout << std::endl << std::endl;

                                    for (IndexType i = 0; i != (r - 1); ++i) {
                                        //TODO Como atualizar MA de forma mais rápida?
                                        MA = prod(transpose(A), A); //TODO A métrica entra aqui
                                        
                                        //TODO Como Atualizar x de forma mais rápida?
                                        x = prod(inverse(MA), prod_block<DimensionsAtCompileTime, 1>(transpose(A), B, 0, k, n, 1)); //TODO A métrica entra aqui
                                        std::cout << "    --- x -current" << std::endl << "    "; print_matrix(x); std::cout << std::endl << std::endl;

                                        auto xi1 = coeff(x, i, 0);
                                        if (abs(xi1) > ZeroTolerance) {
                                            auto xi2 = coeff(x, i + 1, 0);

                                            auto dot11 = coeff(MA, i, i);
                                            auto dot12 = coeff(MA, i, i + 1);
                                            auto dot22 = coeff(MA, i + 1, i + 1);

                                            auto delta1 = 2 * xi1 + dot12 + xi2 * dot22;
                                            auto delta2 = -xi1 * dot11;
                                            auto inv_norm = 1 / sqrt(delta1 * delta1 * dot_product_column(A, i, A, i) + 2 * delta1 * delta2 * dot_product_column(A, i, A, i + 1) + delta2 * delta2); // Normalization under Eiclidean metric

                                            coeff(R, 0, 0) = inv_norm * delta1;
                                            coeff(R, 1, 0) = inv_norm * delta2;
                                            coeff(R, 0, 1) = xi1;
                                            coeff(R, 1, 1) = xi2;

                                            assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, R), A, 0, i, n, 2);
                                            std::cout << "    --- A -updated" << std::endl << "    "; print_matrix(A); std::cout << std::endl << std::endl;
                                        }
                                    }
                                    
                                    //TODO Mudei aqui!
                                    //alpha /= sqrt(dot_product_column(A, r - 1, A, r - 1));
                                }
                                /**/

                                //TODO Mudei aqui
                                alpha *= dot_product_column(A, r - 1, B, k); //TODO A métrica vai aqui
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

                std::cout << "--- C --------" << std::endl; print_matrix(A); std::cout << std::endl << std::endl;

                update_factors(space, alpha, A, MA, x, R, args...);
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

                using TwoByTwoMatrixType = matrix_type_t<ResultingScalarType, 2, 2, 2, 2>;
                using DynamicSquareMatrixType = matrix_type_t<ResultingScalarType, Dynamic, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using DynamicColumnMatrixType = matrix_type_t<ResultingScalarType, Dynamic, 1, MaxDimensionsAtCompileTime, 1>;

                auto const &space = *space_ptr(args...);

                ResultingScalarType alpha = 1;
                auto A = make_matrix<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);
                
                TwoByTwoMatrixType R;
                DynamicSquareMatrixType MA;
                DynamicColumnMatrixType x;

                update_factors(space, alpha, A, MA, x, R, args...);

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

