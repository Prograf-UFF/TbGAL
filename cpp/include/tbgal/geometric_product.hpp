#ifndef __TBGAL_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_GEOMETRIC_PRODUCT_HPP__

namespace tbgal {

    namespace detail {

        struct gp_impl {
        private:

            template<typename MetricSpaceType, typename ResultingScalarType, typename NxKMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &, ResultingScalarType const &, NxKMatrixType const &) noexcept {
            }

            template<typename MetricSpaceType, typename ResultingScalarType, typename NxKMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ResultingScalarType &alpha, NxKMatrixType &A, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                alpha *= arg1.scalar();
                
                if (cols(A) != 0) {
                    if (arg1.factors_count() != 0) {
                        auto B = arg1.factors_in_signed_metric();
                        //TODO [Parei aqui!] prod(transpose(B), B);
                    }
                }
                else {
                    A = arg1.factors_in_signed_metric();
                }

                if (alpha != 0) {
                    update_factors(space, alpha, A, args...);
                }
            }

            template<typename MetricSpaceType, typename ResultingScalarType, typename NxKMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static void update_factors(MetricSpaceType const &space, ResultingScalarType &alpha, NxKMatrixType &A, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                alpha *= arg1;
                if (alpha != 0) {
                    update_factors(space, alpha, A, args...);
                }
            }

        public:

            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) noexcept {
                /**
                using ResultingScalarType = common_scalar_type_t<Types...>;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = GeometricProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;

                constexpr DefaultIndexType DimensionsAtCompileTime = ResultingMetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = ResultingMetricSpaceType::MaxDimensionsAtCompileTime;

                using MatrixNxKType = matrix_type_t<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;

                auto const &space = *space_ptr(args...);
                ScalarType alpha = 1;
                MatrixNxKType A = make_matrix<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), 0);

                update_factors(space, alpha, A, args...);

                return ResultingFactoredMultivectorType(space, alpha, A);
                //TODO Remover funções criadas para serem utilizadas na implementação antiga.
                /*/
                using ScalarType = common_scalar_type_t<Types...>;

                using ResultingScalarType = ScalarType;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = GeometricProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;

                constexpr ScalarType ZeroTolerance = std::numeric_limits<ScalarType>::epsilon();
                constexpr DefaultIndexType DimensionsAtCompileTime = ResultingMetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = ResultingMetricSpaceType::MaxDimensionsAtCompileTime;

                using MatrixNxKType = matrix_type_t<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using MatrixNxDType = matrix_type_t<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, Dynamic>;
                using MatrixKxKType = matrix_type_t<ScalarType, Dynamic, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using MatrixNx1Type = matrix_type_t<ScalarType, DimensionsAtCompileTime, 1, MaxDimensionsAtCompileTime, 1>;
                using Matrix1xNType = matrix_type_t<ScalarType, 1, DimensionsAtCompileTime, 1, MaxDimensionsAtCompileTime>;
                using MatrixKx1Type = matrix_type_t<ScalarType, Dynamic, 1, MaxDimensionsAtCompileTime, 1>;
                using Matrix2x2Type = matrix_type_t<ScalarType, 2, 2, 2, 2>;
                
                using IndexType = index_type_t<MatrixNxKType>;

                auto& space = *space_ptr(args...);
                IndexType factors_count = sum_factors_count(args...);

                IndexType F_cols = first_factors_count(args...);
                MatrixNxKType F = make_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), F_cols);
                fill_matrix_with_first_factors(F, args...);

                ScalarType scalar = multiply_scalars(args...);

                if (factors_count != F_cols) {
                    MatrixKxKType M = prod(transpose(F), F);
                    MatrixKxKType W = eigen_eigenvectors(M);
                    MatrixNxKType O = prod(F, W);

                    IndexType A_cols = factors_count - F_cols;
                    MatrixNxDType A = make_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, Dynamic>(space.dimensions(), A_cols);
                    fill_matrix_with_tail_factors(A, args...);

                    MatrixKx1Type c, t;
                    Matrix2x2Type R = make_matrix<ScalarType, 2, 2, 2, 2>(2, 2);
                    MatrixNx1Type r, a = make_matrix<ScalarType, DimensionsAtCompileTime, 1, MaxDimensionsAtCompileTime, 1>(space.dimensions(), 1);

                    for (IndexType i = 0; i != A_cols; ++i) {
                        // Take the current vector factor a and compute c as a in the coordinate system defines by the column of F, and r as the residual vector.
                        assign_block<DimensionsAtCompileTime, 1>(A, 0, i, a, space.dimensions(), 1);
                        c = prod(prod(W, transpose(O)), a);
                        r = sub(a, prod(F, c));

                        ScalarType r_sqr_norm = dot_product_column(r, 0, r, 0);

                        // Check whether a is included in the subspace spanned by the column of F.
                        if (r_sqr_norm <= ZeroTolerance) {
                            for (IndexType j = 0; j != (F_cols - 1); ++j) {
                                if (abs(coeff(c, j, 0)) > ZeroTolerance) {
                                    // R^{curr} = [phi * alpha, gamma * c_{j}; -phi * beta, gamma * c_{j+1}]  ---  R^{curr} is in the block diagonal matrix Q^{curr} at row j and column j. The other portions of Q^{curr} define an identity matrix.
                                    ScalarType delta = dot_product_column(F, j, F, j + 1);
                                    ScalarType gamma = 1 / sqrt(coeff(c, j, 0) * coeff(c, j, 0) + 2 * coeff(c, j, 0) * coeff(c, j + 1, 0) * delta + coeff(c, j + 1, 0) * coeff(c, j + 1, 0));
                                    ScalarType rho = sqrt(1 / (1 - delta * delta));
                                    ScalarType alpha = coeff(c, j, 0) * delta * (1 + rho) + coeff(c, j + 1, 0);
                                    ScalarType beta = coeff(c, j, 0) * rho;
                                    ScalarType phi = 1 / sqrt(alpha * alpha - 2 * alpha * beta * delta + beta * beta);

                                    coeff(R, 0, 0) = phi * alpha;
                                    coeff(R, 0, 1) = gamma * coeff(c, j, 0);
                                    coeff(R, 1, 0) = -phi * beta;
                                    coeff(R, 1, 1) = gamma * coeff(c, j + 1, 0);

                                    // F^{next} = F^{curr} . Q^{curr}
                                    assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(F, 0, j, space.dimensions(), 2, R), F, 0, j, space.dimensions(), 2);

                                    // W^{next} = inv(Q^{curr}) . W^{curr}
                                    assign_block<2, Dynamic>(prod_block<2, Dynamic>(inverse(R), W, j, 0, 2, F_cols), W, j, 0, 2, F_cols);

                                    // M^{next} = transpose(Q^{curr}) . M^{curr} . Q^{curr}
                                    assign_block<2, Dynamic>(prod_block<2, Dynamic>(transpose(R), M, j, 0, 2, F_cols), M, j, 0, 2, F_cols);
                                    assign_block<Dynamic, 2>(prod_block<Dynamic, 2>(M, 0, j, F_cols, 2, R), M, 0, j, F_cols, 2);

                                    // c^{next} = (W^{next} . transpose(O)) . a
                                    assign_block<2, 1>(prod(prod_block<2, Dynamic>(W, j, 0, 2, F_cols, transpose(O)), a), c, j, 0, 2, 1);
                                }
                            }

                            // Update the resulting scalar factor.
                            scalar *= dot_product_column(F, F_cols - 1, a, 0);
                            
                            // Remove the rightmost last column of F and update matrices M, O, and W accordingly.
                            conservative_resize(F, space.dimensions(), F_cols - 1);
                            conservative_resize(M, F_cols - 1, F_cols - 1);
                            if (F_cols > 1) {
                                //TODO Posso simplificar as contas para obter a última linha ou coluna desejada.
                                //auto O_ = evaluate(prod(compute_reflection_matrix_from_rows(O, F_cols - 1, prod(W, transpose(O)), F_cols - 1), O));
                                //conservative_resize(O_, space.dimensions(), F_cols - 1);
                                //
                                //auto W_ = evaluate(prod(compute_reflection_matrix_from_cols(transpose(W), F_cols - 1, prod(W, transpose(W)), F_cols - 1), W));
                                //conservative_resize(W_, F_cols - 1, F_cols - 1);
                                
                                W = eigen_eigenvectors(M); //TODO [FUTURE] It could be faster (see Fontijine's Thesis, Section 5.4.4).
                                O = prod(F, W);            //TODO [FUTURE] It could be faster (see Fontijine's Thesis, Section 5.4.4).
                            }
                            else {
                                W = make_zero_matrix<ScalarType, 0, 0, 0, 0>(0, 0);
                                O = make_zero_matrix<ScalarType, DimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime, 0>(space.dimensions(), 0);
                            }
                            --F_cols;
                        }
                        else {
                            ScalarType inv_r_norm = 1 / sqrt(r_sqr_norm);

                            // M = [M (transpose(F) . a); (transpose(a) . F) 1]
                            t = prod(transpose(F), a);

                            conservative_resize(M, F_cols + 1, F_cols + 1);
                            assign_block<Dynamic, 1>(t, M, 0, F_cols, F_cols, 1);
                            assign_block<1, Dynamic>(transpose(t), M, F_cols, 0, 1, F_cols);
                            coeff(M, F_cols, F_cols) = 1;

                            // F = [F a]
                            conservative_resize(F, space.dimensions(), F_cols + 1);
                            assign_block<DimensionsAtCompileTime, 1>(a, F, 0, F_cols, space.dimensions(), 1);

                            // O = [O r/r_norm]
                            for (IndexType k = 0; k != space.dimensions(); ++k) {
                                coeff(r, k, 0) *= inv_r_norm;
                            }

                            conservative_resize(O, space.dimensions(), F_cols + 1);
                            assign_block<DimensionsAtCompileTime, 1>(r, O, 0, F_cols, space.dimensions(), 1);
                            
                            // W = [W -c/r_norm; 0 1/r_norm]
                            for (IndexType k = 0; k != F_cols; ++k) {
                                coeff(c, k, 0) *= -inv_r_norm;
                            }

                            conservative_resize(W, F_cols + 1, F_cols + 1);
                            assign_block<Dynamic, 1>(c, W, 0, F_cols, F_cols, 1);
                            assign_block<1, Dynamic>(make_zero_matrix<ScalarType, 1, Dynamic, 1, MaxDimensionsAtCompileTime>(1, F_cols), W, F_cols, 0, 1, F_cols);
                            coeff(W, F_cols, F_cols) = inv_r_norm;

                            ++F_cols;
                        }
                    }
                }

                return ResultingFactoredMultivectorType(space, scalar, F);
                /**/
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

