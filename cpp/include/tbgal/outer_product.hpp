#ifndef __TBGAL_OUTER_PRODUCT_HPP__
#define __TBGAL_OUTER_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement the outer product for the general case (when input multivectors may not be blades).

    namespace detail {

        template<bool AnyMultivectorType>
        struct OP_impl {
            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) noexcept {
                using ResultingScalarType = common_scalar_type_t<Types...>;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = OuterProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
                auto& space = *space_ptr(args...);
                auto factors_count = sum_factors_count(args...);
                if (factors_count <= space.dimensions()) {
                    auto prod_scalar = multiply_scalars(args...);
                    if (factors_count > 0 && prod_scalar != 0) {
                        auto input = make_matrix<ResultingScalarType, ResultingMetricSpaceType::DimensionsAtCompileTime, Dynamic>(space.dimensions(), factors_count);
                        auto qr_tuple = qr_decomposition(std::get<0>(fill_matrix_with_factors(input, args...)));
                        if (std::get<2>(qr_tuple) == factors_count) {
                            auto const &matrix_q = std::get<0>(qr_tuple);
                            return ResultingFactoredMultivectorType(
                                space,
                                prod_scalar * determinant(prod(transpose(left_columns(matrix_q, factors_count)), input)),
                                matrix_q,
                                factors_count
                            );
                        }
                    }
                    else {
                        return ResultingFactoredMultivectorType(space, prod_scalar);
                    }
                }
                return ResultingFactoredMultivectorType(space, 0);
            }
        };

        template<>
        struct OP_impl<false> {

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
    constexpr decltype(auto) OP(FirstType const &arg1, NextTypes const &... args) noexcept {
        return detail::OP_impl<detail::is_any_v<std::true_type, is_multivector_t<FirstType>, is_multivector_t<NextTypes>...> >::eval(arg1, args...);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator^(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return OP(arg1, arg2);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator^(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return OP(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator^(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return OP(arg1, arg2);
    }

}

#endif // __TBGAL_OUTER_PRODUCT_HPP__
