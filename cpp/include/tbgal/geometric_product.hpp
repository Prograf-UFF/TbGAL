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

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType>
            constexpr static void update_factors(MetricSpaceType const &, ScalarType const &, FactorsMatrixType const &) noexcept {
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;
                constexpr ScalarType ZeroTolerance = 100 * std::numeric_limits<ScalarType>::epsilon(); //TODO Adequado?

                alpha *= arg1;
                if (abs(alpha) <= ZeroTolerance) {
                    A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);
                    return;
                }

                update_factors(space, alpha, A, args...);
            }

            template<typename MetricSpaceType, typename ScalarType, typename FactorsMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ScalarType &alpha, FactorsMatrixType &A, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                using IndexType = typename MetricSpaceType::IndexType;

                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = MetricSpaceType::MaxDimensionsAtCompileTime;
                constexpr ScalarType ZeroTolerance = 100 * std::numeric_limits<ScalarType>::epsilon(); //TODO Adequado?

                using DynamicColumnMatrixType = matrix_type_t<ScalarType, Dynamic, 1, MaxDimensionsAtCompileTime, 1>;
                using TwoByTwoMatrixType = matrix_type_t<ScalarType, 2, 2, 2, 2>;
                
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

                        auto qr_tuple_A = qr_orthogonal_matrix(A);
                        assert(r == std::get<1>(qr_tuple_A));

                        auto const &QA = std::get<0>(qr_tuple_A);
                        auto OA = block_view<DimensionsAtCompileTime, Dynamic>(QA, 0, 0, n, r);
                        
                        std::cout << "--- A --------" << std::endl; print_matrix(A); std::cout << std::endl << std::endl;
                        std::cout << "--- QA -------" << std::endl; print_matrix(QA); std::cout << std::endl << std::endl;
                        std::cout << "--- OA -------" << std::endl; print_matrix(OA); std::cout << std::endl << std::endl;

                        auto qr_tuple_B = qr_orthogonal_matrix(evaluate(B));
                        assert(s == std::get<1>(qr_tuple_B));

                        auto const &QB = std::get<0>(qr_tuple_B);
                        auto OB = block_view<DimensionsAtCompileTime, Dynamic>(QB, 0, 0, n, s);

                        std::cout << "--- B --------" << std::endl; print_matrix(B); std::cout << std::endl << std::endl;
                        std::cout << "--- QB -------" << std::endl; print_matrix(QB); std::cout << std::endl << std::endl;
                        std::cout << "--- OB -------" << std::endl; print_matrix(OB); std::cout << std::endl << std::endl;

                        auto svd_tuple = singular_value_decomposition(prod(transpose(OA), OB));
                        auto const &lambda = std::get<0>(svd_tuple);
                        auto const &U = std::get<1>(svd_tuple);
                        auto const &V = std::get<2>(svd_tuple);
                        auto const rank = std::get<3>(svd_tuple);

                        std::cout << "--- lambda ---" << std::endl << lambda << std::endl << std::endl;
                        std::cout << "--- U --------" << std::endl; print_matrix(U); std::cout << std::endl << std::endl;
                        std::cout << "--- V --------" << std::endl; print_matrix(V); std::cout << std::endl << std::endl;
                        std::cout << "--- rank -----" << std::endl << rank << std::endl << std::endl;

                        IndexType t = 0;
                        while (t != rank && abs(abs(lambda[t]) - 1) <= ZeroTolerance) {
                            ++t;
                        }

                        if (t > 0) {
                            auto K = evaluate(prod_block<Dynamic, Dynamic>(OA, U, 0, 0, rows(U), t));

                            std::cout << "--- K --------" << std::endl; print_matrix(K); std::cout << std::endl << std::endl;

                            alpha *= metric_factor(space, K); //TODO Essa é o fator de métrica que devo utilizar aqui?
                            if (abs(alpha) <= ZeroTolerance) {
                                A = make_zero_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(n, 0);
                                return;
                            }

                            FactorsMatrixType B_ = B; //TODO Mudarei aqui!
                            TwoByTwoMatrixType R;
                            DynamicColumnMatrixType x;
                            
                            for (IndexType k = 0; k != t; ++k) {
                                std::cout << "    --- k" << k << " -------" << std::endl << "    "; print_matrix(evaluate(block_view<DimensionsAtCompileTime, 1>(K, 0, k, n, 1))); std::cout << std::endl << std::endl;

                                std::cout << "    --- A -before-" << std::endl << "    "; print_matrix(A); std::cout << std::endl << std::endl;

                                x = prod_block<DimensionsAtCompileTime, 1>(transpose(A), K, 0, k, n, 1); //TODO A métrica entra aqui
                                std::cout << "    --- x -before-" << std::endl << "    "; print_matrix(evaluate(x)); std::cout << std::endl << std::endl;

                                for (IndexType i = 0; i != (r - 1); ++i) {
                                    auto xi1 = coeff(x, i, 0);
                                    if (abs(xi1) > ZeroTolerance) {
                                        auto xi2 = coeff(x, i + 1, 0);

                                        auto inv_norm = 1 / sqrt(xi1 * xi1 - 2 * xi1 * xi2 * dot_product_column(A, i, A, i + 1) + xi2 * xi2); // Normalization under Euclidean metric
                                        
                                        auto dot11 = dot_product_column(A, i, A, i); //TODO A métrica vai aqui
                                        auto dot12 = dot_product_column(A, i, A, i + 1); //TODO A métrica vai aqui
                                        auto dot22 = dot_product_column(A, i + 1, A, i + 1); //TODO A métrica vai aqui

                                        auto rho11 = -inv_norm * xi2;
                                        auto rho21 = inv_norm * xi1;
                                        auto inv_norm_sqr_r1 = 1 / (rho11 * rho11 * dot11 + 2 * rho11 * rho21 * dot12 + rho21 * rho21 * dot22);
                                        
                                        coeff(R, 0, 0) = rho11;
                                        coeff(R, 1, 0) = rho21;
                                        coeff(R, 0, 1) = -inv_norm_sqr_r1 * (rho21 * dot22);
                                        coeff(R, 1, 1) = inv_norm_sqr_r1 * (rho11 * dot11 + 2 * rho21 * dot12);

                                        assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(A, 0, i, n, 2, R), A, 0, i, n, 2);
                                        assign_block<2, 1>(prod_block<2, DimensionsAtCompileTime, 1>(transpose(A), i, 0, 2, K, 0, k, n, 1), x, i, 0, 2, 1); //TODO A métrica entra aqui
                                        std::cout << "    --- A -updated" << std::endl << "    "; print_matrix(A); std::cout << std::endl << std::endl;
                                        std::cout << "    --- x -updated" << std::endl << "    "; print_matrix(x); std::cout << std::endl << std::endl;
                                    }
                                }
                                
                                alpha /= sqrt(dot_product_column(A, r - 1, A, r - 1));
                                conservative_resize(A, n, --r);
                                
                                std::cout << "    --- B -before-" << std::endl << "    "; print_matrix(B_); std::cout << std::endl << std::endl;

                                x = prod_block<DimensionsAtCompileTime, 1>(transpose(B_), K, 0, k, n, 1); //TODO A métrica entra aqui
                                std::cout << "    --- y -before-" << std::endl << "    "; print_matrix(evaluate(x)); std::cout << std::endl << std::endl;

                                for (IndexType i = s - 1; i != 0; --i) {
                                    auto xi2 = coeff(x, i, 0);
                                    if (abs(xi2) > ZeroTolerance) {
                                        auto xi1 = coeff(x, i - 1, 0);

                                        auto inv_norm = 1 / sqrt(xi1 * xi1 - 2 * xi1 * xi2 * dot_product_column(B_, i - 1, B_, i) + xi2 * xi2); // Normalization under Euclidean metric

                                        auto dot11 = dot_product_column(B_, i - 1, B_, i - 1); //TODO A métrica vai aqui
                                        auto dot12 = dot_product_column(B_, i - 1, B_, i); //TODO A métrica vai aqui
                                        auto dot22 = dot_product_column(B_, i, B_, i); //TODO A métrica vai aqui

                                        auto rho12 = -inv_norm * xi2;
                                        auto rho22 = inv_norm * xi1;
                                        auto inv_norm_sqr_r2 = 1 / (rho12 * rho12 * dot11 + 2 * rho12 * rho22 * dot12 + rho22 * rho22 * dot22);
                                        
                                        coeff(R, 0, 0) = -inv_norm_sqr_r2 * (rho22 * dot22);
                                        coeff(R, 1, 0) = inv_norm_sqr_r2 * (rho12 * dot11 + 2 * rho22 * dot12);
                                        coeff(R, 0, 1) = rho12;
                                        coeff(R, 1, 1) = rho22;

                                        assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(B_, 0, i - 1, n, 2, R), B_, 0, i - 1, n, 2);
                                        assign_block<2, 1>(prod_block<2, DimensionsAtCompileTime, 1>(transpose(B_), i - 1, 0, 2, K, 0, k, n, 1), x, i - 1, 0, 2, 1); //TODO A métrica entra aqui
                                        std::cout << "    --- B -updated" << std::endl << "    "; print_matrix(B_); std::cout << std::endl << std::endl;
                                        std::cout << "    --- y -updated" << std::endl << "    "; print_matrix(x); std::cout << std::endl << std::endl;
                                    }
                                }

                                alpha /= sqrt(dot_product_column(B_, 0, B_, 0));
                                assign_block<DimensionsAtCompileTime, Dynamic>(B_, 0, 1, B_, n, s - 1);
                                conservative_resize(B_, n, --s);
                            }

                            conservative_resize(A, n, r + s);
                            assign_block<DimensionsAtCompileTime, Dynamic>(B_, A, 0, r, n, s);
                        }
                        else {
                            conservative_resize(A, n, r + s);
                            assign_block<DimensionsAtCompileTime, Dynamic>(B, A, 0, r, n, s);
                        }
                    }
                }
                else {
                    A = B;
                }

                std::cout << "--- C --------" << std::endl; print_matrix(A); std::cout << std::endl << std::endl;

                update_factors(space, alpha, A, args...);
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

                auto const &space = *space_ptr(args...);
                ResultingScalarType alpha = 1;
                auto A = make_matrix<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);

                update_factors(space, alpha, A, args...);

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

