#ifndef __TBGAL_OUTER_PRODUCT_HPP__
#define __TBGAL_OUTER_PRODUCT_HPP__

namespace tbgal {

    namespace detail {

        template<bool AllNativeScalars>
        struct OP_impl {

            template<typename FirstType, typename... NextTypes>
            constexpr decltype(auto) eval(FirstType const &arg1, NextTypes const &... args) noexcept {
                return arg1 * eval(args...);
            }

            template<typename FirstType, typename SecondType>
            constexpr decltype(auto) eval(FirstType const &arg1, SecondType const &arg2) noexcept {
                return arg1 * arg2;
            }
        };

        template<>
        struct OP_impl<false> {
        private:

            template<typename FirstType, typename... NextTypes>
            struct factoring_product_type;

            template<typename FirstType, typename... NextTypes>
            using factoring_product_type_t = typename factoring_product_type<FirstType, NextTypes...>::type;

            template<typename FirstType, typename... NextTypes>
            struct factoring_product_type :
                factoring_product_type<NextTypes...> {
            };

            template<typename FactoringProductType, typename SquareMatrixType, typename... NextTypes>
            struct factoring_product_type<FactoredMultivector<FactoringProductType, SquareMatrixType>, NextTypes...> {
                using type = std::conditional_t<
                    detail::is_any_v<std::true_type, detail::is_multivector_t<NextTypes>...>,
                    FactoredMultivector<OuterProduct<typename FactoringProductType::SpaceType>, SquareMatrixType>,
                    FactoredMultivector<FactoringProductType, SquareMatrixType>
                >;
            };

            template<typename Type>
            struct factoring_product_type {
                using type = std::nullptr_t;
            };

            template<typename FactoringProductType, typename SquareMatrixType>
            struct factoring_product_type<FactoredMultivector<FactoringProductType, SquareMatrixType> > {
                using type = FactoringProductType;
            };

            template<typename ResultingMatrixType, typename FirstType, typename... NextTypes>
            constexpr std::tuple<ResultingMatrixType&, std::size_t> build_input_matrix(ResultingMatrixType &result, FirstType const &arg1, NextTypes const &... args) noexcept {
                return build_input_matrix(result, args...);
            }

            template<typename ResultingMatrixType, typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr std::tuple<ResultingMatrixType&, std::size_t> build_input_matrix(ResultingMatrixType &result, FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                constexpr std::size_t end_column_index = std::get<1>(build_input_matrix(result, args...));
                return std::make_tuple(copy_first_columns(result, end_column_index - arg1.factors_count(), arg1, arg1.factors_count()), end_column_index - arg1.factors_count());
            }

            template<typename ResultingMatrixType, typename Type>
            constexpr std::tuple<ResultingMatrixType&, std::size_t> build_input_matrix(ResultingMatrixType &result, Type const &) noexcept {
                return std::make_tuple(result, cols(result));
            }

            template<typename ResultingMatrixType, typename FactoringProductType, typename SquareMatrixType>
            constexpr std::tuple<ResultingMatrixType&, std::size_t> build_input_matrix(ResultingMatrixType &result, FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                constexpr std::size_t end_column_index = cols(result);
                return std::make_tuple(copy_first_columns(result, end_column_index - arg.factors_count(), arg, arg.factors_count()), end_column_index - arg.factors_count());
            }

            template<typename FirstType, typename... NextTypes>
            constexpr decltype(auto) multiply_scalars(FirstType const &arg1, NextTypes const &... args) noexcept {
                return arg1 * multiply_scalars(args...);
            }

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr decltype(auto) multiply_scalars(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                return arg1.scalar() + multiply_scalars(args...);
            }

            template<typename Type>
            constexpr decltype(auto) multiply_scalars(Type const &arg) noexcept {
                return arg;
            }

            template<typename FactoringProductType, typename SquareMatrixType>
            constexpr decltype(auto) multiply_scalars(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                return arg.scalar();
            }

            template<typename FirstType, typename... NextTypes>
            constexpr decltype(auto) space_ptr(FirstType const &, NextTypes const &... args) noexcept {
                return space_ptr(args...);
            }

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr decltype(auto) space_ptr(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                return &arg1.space();
            }

            template<typename Type>
            constexpr decltype(auto) space_ptr(Type const &) noexcept {
                return nullptr;
            }

            template<typename FactoringProductType, typename SquareMatrixType>
            constexpr decltype(auto) space_ptr(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                return &arg.space();
            }

            template<typename FirstType, typename... NextTypes>
            constexpr decltype(auto) sum_factors_count(FirstType const &, NextTypes const &... args) noexcept {
                return sum_factors_count(args...);
            }

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr decltype(auto) sum_factors_count(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                return arg1.factors_count() + sum_factors_count(args...);
            }

            template<typename Type>
            constexpr decltype(auto) sum_factors_count(Type const &) noexcept {
                return 0;
            }

            template<typename FactoringProductType, typename SquareMatrixType>
            constexpr decltype(auto) sum_factors_count(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                return arg.factors_count();
            }

        public:

            template<typename FirstType, typename... NextTypes>
            constexpr decltype(auto) eval(FirstType const &arg1, NextTypes const &... args) noexcept {
                using InputFactorsMatrixType = input_factors_matrix_type_t<FirstType, NextTypes...>;
                using ResultingFactoringProductType = factoring_product_type_t<FirstType, NextTypes...>;
                using ResultingSquareMatrixType = squared_matrix_type_t<FirstType, NextTypes...>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
                auto& space = *space_ptr(arg1, args...);
                auto factors_count = sum_factors_count(arg1, args...);
                if (factors_count <= space.dimensions()) {
                    auto prod_scalar = multiply_scalars(arg1, args...);
                    if (prod_scalar != 0) {
                        auto qr = qr_decomposition(std::get<0>(build_input_matrix(make_matrix<InputFactorsMatrixType>(space.dimensions(), factors_count), arg1, args...)));
                        if (factors_count == qr.rank()) {
                            return ResultingFactoredMultivectorType(space, prod_scalar * determinant_triangular(qr.matrix_r()), qr.matrix_q(), factors_count);
                        }
                    }
                }
                return ResultingFactoredMultivectorType(space, 0);
            }
        };

    }

    template<typename FirstType, typename... NextTypes>
    constexpr decltype(auto) OP(FirstType const &arg1, NextTypes const &... args) noexcept {
        return detail::OP_impl<!detail::is_any_v<std::true_type, detail::is_multivector_t<FirstType>, detail::is_multivector_t<NextTypes>...> >::eval(arg1, args...);
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) operator^(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return OP(arg1, arg2);
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) operator^(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        return OP(arg1, arg2);
    }

    template<typename FirstType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) operator^(FirstType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return OP(arg1, arg2);
    }

}

#endif // __TBGAL_OUTER_PRODUCT_HPP__
