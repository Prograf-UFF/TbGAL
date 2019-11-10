#ifndef __TBGAL_OUTER_PRODUCT_HPP__
#define __TBGAL_OUTER_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement the outer product for the general case (when input multivectors may not be blades).

    namespace detail {

        template<bool AnyFactoredMultivector>
        struct op_impl {
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
                        auto input = make_matrix<ResultingScalarType, ResultingMetricSpaceType::DimensionsAtCompileTime, Dynamic, ResultingMetricSpaceType::MaxDimensionsAtCompileTime, Dynamic>(space.dimensions(), factors_count);
                        auto qr_tuple = qr_orthogonal_matrix(std::get<0>(fill_matrix_with_tail_factors(input, args...)));
                        if (std::get<1>(qr_tuple) == factors_count) {
                            auto const &matrix_q = std::get<0>(qr_tuple);
                            return ResultingFactoredMultivectorType(
                                space,
                                prod_scalar * determinant(prod_block<Dynamic, ResultingMetricSpaceType::DimensionsAtCompileTime>(transpose(matrix_q), 0, 0, factors_count, space.dimensions(), input)),
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
        struct op_impl<false> {
            template<typename... ScalarTypes>
            constexpr static decltype(auto) eval(ScalarTypes const &... args) noexcept {
                return multiply_scalars(args...);
            }
        };

    }

    template<typename FirstType, typename... NextTypes>
    constexpr decltype(auto) op(FirstType const &arg1, NextTypes const &... args) noexcept {
        return detail::op_impl<detail::is_any_v<std::true_type, is_multivector_t<FirstType>, is_multivector_t<NextTypes>...> >::eval(arg1, args...);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator^(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return op(arg1, arg2);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator^(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return op(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator^(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return op(arg1, arg2);
    }

}

#endif // __TBGAL_OUTER_PRODUCT_HPP__
