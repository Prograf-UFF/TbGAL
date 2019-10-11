#ifndef __TBGAL_OUTER_PRODUCT_HPP__
#define __TBGAL_OUTER_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement the outer product for the general case (when input multivectors may not be blades).

    namespace detail {

        template<bool AnyMultivectorType>
        struct OP_impl {
        private:

            template<typename FirstType, typename... NextTypes>
            struct common_scalar_type;

            template<typename... Types>
            using common_scalar_type_t = typename common_scalar_type<Types...>::type;

            template<typename FirstScalarType, typename... NextTypes>
            struct common_scalar_type {
                using type = std::common_type_t<FirstScalarType, common_scalar_type_t<NextTypes...> >;
            };

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            struct common_scalar_type<FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType>, NextTypes...> {
                using type = std::common_type_t<typename FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType>::ScalarType, common_scalar_type_t<NextTypes...> >;
            };

            template<typename ScalarType>
            struct common_scalar_type<ScalarType> {
                using type = ScalarType;
            };

            template<typename FactoringProductType, typename SquareMatrixType>
            struct common_scalar_type<FactoredMultivector<FactoringProductType, SquareMatrixType> > {
                using type = typename FactoredMultivector<FactoringProductType, SquareMatrixType>::ScalarType;
            };

            template<typename FirstType, typename... NextTypes>
            struct metric_space_type;

            template<typename... Types>
            using metric_space_type_t = typename metric_space_type<Types...>::type;

            template<typename FirstType, typename... NextTypes>
            struct metric_space_type :
                metric_space_type<NextTypes...> {
            };

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            struct metric_space_type<FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType>, NextTypes...> {
            private:

                using NextSpaceType = metric_space_type_t<NextTypes...>;
                static_assert(std::is_same_v<NextSpaceType, std::nullptr_t> || std::is_same_v<NextSpaceType, typename FirstFactoringProductType::MetricSpaceType>, "The multivectors must have the same metric space.");

            public:

                using type = typename FirstFactoringProductType::MetricSpaceType;
            };

            template<typename Type>
            struct metric_space_type<Type> {
                using type = std::nullptr_t;
            };

            template<typename FactoringProductType, typename SquareMatrixType>
            struct metric_space_type<FactoredMultivector<FactoringProductType, SquareMatrixType> > {
                using type = typename FactoringProductType::MetricSpaceType;
            };

            template<typename ResultingMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static decltype(auto) fill_input_matrix(ResultingMatrixType &result, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                return fill_input_matrix(result, args...);
            }

            template<typename ResultingMatrixType, typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr static decltype(auto) fill_input_matrix(ResultingMatrixType &result, FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                assert(is_blade(arg1));
                auto end_column_index = std::get<1>(fill_input_matrix(result, args...));
                return std::make_tuple(copy_columns<Dynamic>(arg1.factors_in_signed_metric(), 0, result, end_column_index - arg1.factors_count(), arg1.factors_count()), end_column_index - arg1.factors_count());
            }

            template<typename ResultingMatrixType, typename Type>
            constexpr static decltype(auto) fill_input_matrix(ResultingMatrixType &result, Type const &) noexcept {
                return std::make_tuple(result, cols(result));
            }

            template<typename ResultingMatrixType, typename FactoringProductType, typename SquareMatrixType>
            constexpr static decltype(auto) fill_input_matrix(ResultingMatrixType &result, FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                assert(is_blade(arg));
                return std::make_tuple(copy_columns<Dynamic>(arg.factors_in_signed_metric(), 0, result, cols(result) - arg.factors_count(), arg.factors_count()), cols(result) - arg.factors_count());
            }

            template<typename FirstScalarType, typename... NextTypes>
            constexpr static decltype(auto) multiply_scalars(FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                return arg1 * multiply_scalars(args...);
            }

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr static decltype(auto) multiply_scalars(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                return arg1.scalar() * multiply_scalars(args...);
            }

            template<typename Type>
            constexpr static decltype(auto) multiply_scalars(Type const &arg) noexcept {
                return arg;
            }

            template<typename FactoringProductType, typename SquareMatrixType>
            constexpr static decltype(auto) multiply_scalars(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                return arg.scalar();
            }

            template<typename FirstScalarType, typename... NextTypes>
            constexpr static decltype(auto) space_ptr(FirstScalarType const &, NextTypes const &... args) noexcept {
                return space_ptr(args...);
            }

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr static decltype(auto) space_ptr(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                assert(space_ptr(args...) == &arg1.space() || space_ptr(args...) == nullptr);
                return &arg1.space();
            }

            template<typename Type>
            constexpr static decltype(auto) space_ptr(Type const &) noexcept {
                return nullptr;
            }

            template<typename FactoringProductType, typename SquareMatrixType>
            constexpr static decltype(auto) space_ptr(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                return &arg.space();
            }

            template<typename FirstScalarType, typename... NextTypes>
            constexpr static decltype(auto) sum_factors_count(FirstScalarType const &, NextTypes const &... args) noexcept {
                return sum_factors_count(args...);
            }

            template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename... NextTypes>
            constexpr static decltype(auto) sum_factors_count(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, NextTypes const &... args) noexcept {
                assert(is_blade(arg1));
                return arg1.factors_count() + sum_factors_count(args...);
            }

            template<typename Type>
            constexpr static decltype(auto) sum_factors_count(Type const &) noexcept {
                return 0;
            }

            template<typename FactoringProductType, typename SquareMatrixType>
            constexpr static decltype(auto) sum_factors_count(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
                assert(is_blade(arg));
                return arg.factors_count();
            }

        public:

            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) noexcept {
                using ResultingScalarType = common_scalar_type_t<Types...>;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = OuterProduct<ResultingMetricSpaceType>;
                using ResultingSquareMatrixType = matrix_type_t<ResultingScalarType, ResultingMetricSpaceType::DimensionsAtCompileTime, ResultingMetricSpaceType::DimensionsAtCompileTime>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
                auto& space = *space_ptr(args...);
                auto factors_count = sum_factors_count(args...);
                if (factors_count <= space.dimensions()) {
                    auto prod_scalar = multiply_scalars(args...);
                    if (factors_count > 0 && prod_scalar != 0) {
                        auto input = make_matrix<ResultingScalarType, ResultingMetricSpaceType::DimensionsAtCompileTime, Dynamic>(space.dimensions(), factors_count);
                        auto qr_tuple = qr_decomposition(std::get<0>(fill_input_matrix(input, args...)));
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

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) operator^(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return OP(arg1, arg2);
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator^(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondScalarType const &arg2) noexcept {
        return OP(arg1, arg2);
    }

    template<typename FirstType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!is_multivector_v<FirstType> > >
    constexpr decltype(auto) operator^(FirstType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return OP(arg1, arg2);
    }

}

#endif // __TBGAL_OUTER_PRODUCT_HPP__
