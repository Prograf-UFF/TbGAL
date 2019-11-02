#ifndef __TBGAL_GEOMETRIC_PRODUCT_HPP__
#define __TBGAL_GEOMETRIC_PRODUCT_HPP__

namespace tbgal {

    namespace detail {

        //TODO {DEBUG}
        template<typename MatrixType>
        void print_matrix(MatrixType const &arg) noexcept {
            std::cout << "[";
            for (DefaultIndexType row = 0; row != rows(arg); ++row) {
                for (DefaultIndexType col = 0; col != cols(arg); ++col) {
                    std::cout << coeff(arg, row, col);
                    if (col + 1 != cols(arg)) {
                        std::cout << ", ";
                    }
                }
                if (row + 1 != rows(arg)) {
                    std::cout << "; ";
                }
            }
            std::cout << "]";
        }

        template<bool AnyMultivectorType>
        struct GP_impl {
            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) noexcept {
                using ScalarType = common_scalar_type_t<Types...>;
                
                using ResultingScalarType = ScalarType;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = GeometricProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
                
                constexpr ScalarType ZeroTolerance = std::numeric_limits<ScalarType>::epsilon();
                constexpr DefaultIndexType DimensionsAtCompileTime = ResultingMetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = ResultingMetricSpaceType::MaxDimensionsAtCompileTime;

                using TxT_MatrixType = matrix_type_t<ScalarType, 2, 2, 2, 2>;
                using NxK_MatrixType = matrix_type_t<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using NxD_MatrixType = matrix_type_t<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, Dynamic>;
                using KxK_MatrixType = matrix_type_t<ScalarType, Dynamic, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
                using Nx1_MatrixType = matrix_type_t<ScalarType, DimensionsAtCompileTime, 1, MaxDimensionsAtCompileTime, 1>;
                using Kx1_MatrixType = matrix_type_t<ScalarType, Dynamic, 1, MaxDimensionsAtCompileTime, 1>;
                
                using IndexType = index_type_t<NxK_MatrixType>;

                auto& space = *space_ptr(args...);
                IndexType factors_count = sum_factors_count(args...);

                IndexType F_cols = first_factors_count(args...);
                NxK_MatrixType F = make_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime>(space.dimensions(), F_cols);
                fill_matrix_with_first_factors(F, args...);

                if (factors_count != F_cols) {
                    KxK_MatrixType Mf = prod(transpose(F), F);
                    KxK_MatrixType W = svd_left_eigenvectors(Mf);
                    NxK_MatrixType O = prod(F, W);

                    IndexType A_cols = factors_count - F_cols;
                    NxD_MatrixType A = make_matrix<ScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, Dynamic>(space.dimensions(), A_cols);
                    fill_matrix_with_tail_factors(A, args...);

                    TxT_MatrixType R;
                    Kx1_MatrixType c, t;
                    Nx1_MatrixType r, a = make_matrix<ScalarType, DimensionsAtCompileTime, 1, MaxDimensionsAtCompileTime, 1>(space.dimensions(), 1);
                    
                    for (IndexType ind = 0; ind != A_cols; ++ind) {
                        assign_block<DimensionsAtCompileTime, 1>(A, 0, ind, a, 0, 0, space.dimensions(), 1);
                        c = prod(prod(W, transpose(O)), a);
                        r = sub(a, prod(F, c));

                        ScalarType r_sqr_norm = 0;
                        for (IndexType k = 0; k != space.dimensions(); ++k) {
                            r_sqr_norm += coeff(r, k, 0) * coeff(r, k, 0);
                        }

                        if (r_sqr_norm <= ZeroTolerance) {
                            for (IndexType j = 0; j != (F_cols - 1); ++j) {
                                if (abs(coeff(c, j, 0)) > ZeroTolerance) {
                                    // R^{curr} = [alpha c_{[j]}^{curr}; beta c_{[j+1]}^{curr}]
                                    // Q^{curr} = [I_{j-1} 0 0; 0 R^{curr} 0; 0 0 I_{k - j - 1}] --- Q does not matter, since only two columns or rows of each affected matrix are modied
                                    // inv(Q^{curr}) = [I_{j-1} 0 0; 0 inv(R^{curr}) 0; 0 0 I_{k - j - 1}]
                                    ScalarType dot = 0;
                                    for (IndexType k = 0; k != F_cols; ++k) {
                                        dot += coeff(F, k, j) * coeff(F, k, j + 1);
                                    }
                                    ScalarType inv_weight = 1 / (coeff(c, j, 0) * coeff(c, j, 0) + 2 * coeff(c, j, 0) * coeff(c, j + 1, 0) * dot + coeff(c, j + 1, 0) * coeff(c, j + 1, 0));

                                    coeff(R, 0, 0) = inv_weight * coeff(c, j, 0);
                                    coeff(R, 0, 1) = coeff(c, j, 0);
                                    coeff(R, 1, 0) = inv_weight * coeff(c, j + 1, 0);
                                    coeff(R, 1, 1) = coeff(c, j + 1, 0);

                                    // F^{next} = F^{curr} . Q^{curr}

                                    //TODO {DEBUG}
                                    std::cout << "---- j = " << j << " ------------------------------------------" << std::endl;
                                    std::cout << "    a = "; print_matrix(a); std::cout << std::endl;
                                    std::cout << "    c = "; print_matrix(c); std::cout << std::endl;
                                    std::cout << "    r = "; print_matrix(r); std::cout << std::endl;
                                    std::cout << "    F = "; print_matrix(F); std::cout << std::endl;
                                    std::cout << "    R = "; print_matrix(R); std::cout << std::endl;

                                    assign_block<DimensionsAtCompileTime, 2>(prod_block<DimensionsAtCompileTime, 2>(F, 0, j, space.dimensions(), 2, R), F, 0, j, space.dimensions(), 2);

                                    //TODO {DEBUG}
                                    std::cout << "F . Q = "; print_matrix(F); std::cout << std::endl;
                                    std::cout << "-----------------------------------------------------" << std::endl << std::endl;

                                    // W^{next} = inv_Q^{curr} . W^{curr}
                                    assign_block<2, Dynamic>(prod_block<2, Dynamic>(inverse(R), W, j, 0, 2, F_cols), W, j, 0, 2, F_cols);

                                    // Mf^{next} = transpose(Q^{curr}) . Mf^{curr} . Q^{curr}
                                    assign_block<2, Dynamic>(prod_block<2, Dynamic>(transpose(R), Mf, j, 0, 2, F_cols), Mf, j, 0, 2, F_cols);
                                    assign_block<Dynamic, 2>(prod_block<Dynamic, 2>(Mf, 0, j, F_cols, 2, R), Mf, 0, j, F_cols, 2);
                                }
                            }

                            conservative_resize(Mf, F_cols - 1, F_cols - 1);
                            conservative_resize(F, space.dimensions(), F_cols - 1);
                            if (F_cols > 1) {
                                W = svd_left_eigenvectors(Mf); //TODO [FUTURE] It could be faster (see Fontijine's Thesis, Section 5.4.4).
                                O = prod(F, W);                //TODO [FUTURE] It could be faster (see Fontijine's Thesis, Section 5.4.4).
                            }
                            else {
                                W = make_zero_matrix<ScalarType, 0, 0, 0, 0>(0, 0);
                                O = make_zero_matrix<ScalarType, DimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime, 0>(space.dimensions(), 0);
                            }
                            --F_cols;
                        }
                        else {
                            ScalarType inv_r_norm = 1 / sqrt(r_sqr_norm);

                            // Mf = [Mf (transpose(F) . a); (transpose(a) . F) 1]
                            t = prod(transpose(F), a);

                            conservative_resize(Mf, F_cols + 1, F_cols + 1);
                            assign_block<Dynamic, 1>(t, Mf, 0, F_cols, F_cols, 1);
                            assign_block<1, Dynamic>(transpose(t), Mf, F_cols, 0, 1, F_cols);
                            coeff(Mf, F_cols, F_cols) = 1;

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
                            for (IndexType k = 0; k != F_cols; ++k) {
                                coeff(W, F_cols, k) = 0;
                            }
                            coeff(W, F_cols, F_cols) = inv_r_norm;

                            ++F_cols;
                        }
                    }
                }

                return ResultingFactoredMultivectorType(space, multiply_scalars(args...), F);
            }
        };

        template<>
        struct GP_impl<false> {

            template<typename FirstScalarType, typename... NextScalarTypes>
            constexpr static decltype(auto) eval(FirstScalarType const &arg1, NextScalarTypes const &... args) noexcept {
                return arg1 * eval(args...);
            }

            template<typename ScalarType>
            constexpr static decltype(auto) eval(ScalarType const &arg) noexcept {
                return arg;
            }
        };

    }

    template<typename FirstType, typename... NextTypes>
    constexpr decltype(auto) GP(FirstType const &arg1, NextTypes const &... args) noexcept {
        return detail::GP_impl<detail::is_any_v<std::true_type, is_multivector_t<FirstType>, is_multivector_t<NextTypes>...> >::eval(arg1, args...);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator*(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return GP(arg1, arg2);
    }

    template<typename FirstScsalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator*(FactoredMultivector<FirstScsalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return GP(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator*(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return GP(arg1, arg2);
    }

}

#endif // __TBGAL_GEOMETRIC_PRODUCT_HPP__

