#ifndef __TBGAL_OUTER_PRODUCT_HPP__
#define __TBGAL_OUTER_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement the outer product for the general case (when input multivectors may not be blades).

    namespace detail {

        //TODO Pode ser embutido desde que encontre uma solução packed para *
        constexpr decltype(auto) multiply_scalars() noexcept {
            return 1;
        }

        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) multiply_scalars(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return arg.scalar();
        }

        template<typename ScalarType>
        constexpr decltype(auto) multiply_scalars(ScalarType const &arg) noexcept {
            return arg;
        }

        template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
        constexpr decltype(auto) multiply_scalars(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
            return arg1.scalar() * multiply_scalars(args...);
        }

        template<typename FirstScalarType, typename... NextTypes>
        constexpr decltype(auto) multiply_scalars(FirstScalarType const &arg1, NextTypes const &... args) noexcept {
            return arg1 * multiply_scalars(args...);
        }
        
        template<bool AnyFactoredMultivector>
        struct op_impl {
        private:

            template<typename ResultingMatrixType>
            constexpr static std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix(ResultingMatrixType &result) noexcept {
                return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, cols(result));
            }

            template<typename ResultingMatrixType, typename ScalarType>
            constexpr static std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix(ResultingMatrixType &result, ScalarType const &) noexcept {
                return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, cols(result));
            }

            template<typename ResultingMatrixType, typename ScalarType, typename FactoringProductType>
            constexpr static std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix(ResultingMatrixType &result, FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
                using MetricSpaceType = typename FactoringProductType::MetricSpaceType;
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                if (arg.factors_count() > 0) {
                    auto end_column_index = cols(result);
                    if (end_column_index >= arg.factors_count()) {
                        assign_block<DimensionsAtCompileTime, Dynamic>(arg.factors_in_signed_metric(), result, 0, end_column_index - arg.factors_count(), rows(result), arg.factors_count());
                        return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, end_column_index - arg.factors_count());
                    }
                    assert(end_column_index == 0);
                    return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, 0);
                }
                return fill_matrix(result);
            }

            template<typename ResultingMatrixType, typename FirstScalarType, typename... NextTypes>
            constexpr static std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix(ResultingMatrixType &result, FirstScalarType const &arg1, NextTypes const &... args) noexcept {
                return fill_matrix(result, args...);
            }

            template<typename ResultingMatrixType, typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static std::tuple<ResultingMatrixType&, DefaultIndexType> fill_matrix(ResultingMatrixType &result, FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                using MetricSpaceType = typename FirstFactoringProductType::MetricSpaceType;
                constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;
                if (arg1.factors_count() > 0) {
                    auto end_column_index = std::get<1>(fill_matrix(result, args...));
                    if (end_column_index >= arg1.factors_count()) {
                        assign_block<DimensionsAtCompileTime, Dynamic>(arg1.factors_in_signed_metric(), result, 0, end_column_index - arg1.factors_count(), rows(result), arg1.factors_count());
                        return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, end_column_index - arg1.factors_count());
                    }
                    assert(end_column_index == 0);
                    return std::tuple<ResultingMatrixType&, DefaultIndexType>(result, 0);
                }
                return fill_matrix(result, args...);
            }

            //TODO Pode ser embutido, desde de que encontre uma versão packed para +
            constexpr static decltype(auto) sum_factors_count() noexcept {
                return 0;
            }

            template<typename ScalarType, typename FactoringProductType>
            constexpr static decltype(auto) sum_factors_count(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
                assert(is_blade(arg));
                return arg.factors_count();
            }

            template<typename ScalarType>
            constexpr static decltype(auto) sum_factors_count(ScalarType const &) noexcept {
                return 0;
            }

            template<typename FirstScalarType, typename FirstFactoringProductType, typename... NextTypes>
            constexpr static decltype(auto) sum_factors_count(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, NextTypes const &... args) noexcept {
                assert(is_blade(arg1));
                return arg1.factors_count() + sum_factors_count(args...);
            }

            template<typename FirstScalarType, typename... NextTypes>
            constexpr static decltype(auto) sum_factors_count(FirstScalarType const &, NextTypes const &... args) noexcept {
                return sum_factors_count(args...);
            }

        public:

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
                    if (factors_count > 0 && !is_zero(prod_scalar)) {
                        auto input = make_matrix<ResultingScalarType, ResultingMetricSpaceType::DimensionsAtCompileTime, Dynamic, ResultingMetricSpaceType::MaxDimensionsAtCompileTime, Dynamic>(space.dimensions(), factors_count);
                        auto qr_tuple = qr_orthogonal_matrix(std::get<0>(fill_matrix(input, args...)));
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
