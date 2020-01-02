/* Copyright (C) Eduardo Vera Sousa and Leandro Augusto Frata Fernandes
 * 
 * authors    : Sousa, Eduardo V.
 *              Fernandes, Leandro A. F.
 * repository : https://github.com/Prograf-UFF/TbGAL
 * 
 * This file is part of the Tensor-based Geometric Algebra Library (TbGAL).
 * 
 * TbGAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * TbGAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with TbGAL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef __TBGAL_OUTER_PRODUCT_HPP__
#define __TBGAL_OUTER_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement the outer product for the general case (when input multivectors may not be blades).

    namespace detail {

        template<typename ScalarType>
        constexpr decltype(auto) _op_impl_get_factors_count(ScalarType const &arg) noexcept {
            return 0;
        }
        
        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) _op_impl_get_factors_count(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return arg.factors_count();
        }
        
        template<typename ScalarType>
        constexpr decltype(auto) _op_impl_get_scalar(ScalarType const &arg) noexcept {
            return arg;
        }
        
        template<typename ScalarType, typename FactoringProductType>
        constexpr decltype(auto) _op_impl_get_scalar(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            return arg.scalar();
        }
        
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

        public:

            template<typename... Types>
            constexpr static decltype(auto) eval(Types const &... args) noexcept {
                using ResultingScalarType = common_scalar_type_t<Types...>;
                using ResultingMetricSpaceType = metric_space_type_t<Types...>;
                using ResultingFactoringProductType = OuterProduct<ResultingMetricSpaceType>;
                using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;

                constexpr DefaultIndexType DimensionsAtCompileTime = ResultingMetricSpaceType::DimensionsAtCompileTime;
                constexpr DefaultIndexType MaxDimensionsAtCompileTime = ResultingMetricSpaceType::MaxDimensionsAtCompileTime;

                auto& space = *space_ptr(args...);
                auto factors_count = (_op_impl_get_factors_count(args) + ... + 0);
                if (factors_count <= space.dimensions()) {
                    auto prod_scalar = (_op_impl_get_scalar(args) * ... * 1);
                    if (factors_count > 0 && !is_zero(prod_scalar)) {
                        auto input = make_matrix<ResultingScalarType, DimensionsAtCompileTime, Dynamic, MaxDimensionsAtCompileTime, Dynamic>(space.dimensions(), factors_count);
                        auto qr_tuple = qr_orthogonal_matrix(std::get<0>(fill_matrix(input, args...)));
                        if (std::get<1>(qr_tuple) == factors_count) {
                            auto const &matrix_q = std::get<0>(qr_tuple);
                            return ResultingFactoredMultivectorType(
                                space,
                                prod_scalar * determinant(prod_block<Dynamic, DimensionsAtCompileTime>(transpose(matrix_q), 0, 0, factors_count, space.dimensions(), input)),
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

    }

    template<typename FirstType, typename... NextTypes>
    constexpr decltype(auto) op(FirstType const &arg1, NextTypes const &... args) noexcept {
        return detail::op_impl::eval(arg1, args...);
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
