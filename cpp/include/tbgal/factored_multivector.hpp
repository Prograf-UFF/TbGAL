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

#ifndef __TBGAL_MULTIVECTOR_HPP__
#define __TBGAL_MULTIVECTOR_HPP__

namespace tbgal {

    template<typename ScalarType_, typename MetricSpaceType_>
    class FactoredMultivector<ScalarType_, GeometricProduct<MetricSpaceType_> > final {
    private:

        using FactorsMatrixType = detail::matrix_type_t<ScalarType_, MetricSpaceType_::DimensionsAtCompileTime, Dynamic, MetricSpaceType_::MaxDimensionsAtCompileTime, MetricSpaceType_::MaxDimensionsAtCompileTime>;

    public:

        using MetricSpaceType = MetricSpaceType_;

        using IndexType = detail::index_type_t<FactorsMatrixType>;
        using ScalarType = ScalarType_;

        constexpr FactoredMultivector(FactoredMultivector const &) = default;
        constexpr FactoredMultivector(FactoredMultivector &&) = default;

        constexpr FactoredMultivector(MetricSpaceType const &space) noexcept :
            FactoredMultivector(space, 0, detail::make_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, Dynamic, MetricSpaceType::MaxDimensionsAtCompileTime, MetricSpaceType::MaxDimensionsAtCompileTime>(space.dimensions(), 0)) {
        }

        template<typename OtherScalarType>
        constexpr FactoredMultivector(MetricSpaceType const &space, OtherScalarType &&scalar) noexcept :
            FactoredMultivector(space, std::move(scalar), detail::make_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, Dynamic, MetricSpaceType::MaxDimensionsAtCompileTime, MetricSpaceType::MaxDimensionsAtCompileTime>(space.dimensions(), 0)) {
        }

        //TODO [FUTURE] Specialization.
        //template<typename OtherScalarType, typename OtherFactoringProductType>
        //constexpr FactoredMultivector(FactoredMultivector<OtherScalarType, OtherFactoringProductType> const &other) noexcept;

        constexpr FactoredMultivector & operator=(FactoredMultivector const &) = default;
        constexpr FactoredMultivector & operator=(FactoredMultivector &&) = default;

        //TODO [FUTURE] Specialization.
        //template<typename OtherScalarType, typename OtherFactoringProductType>
        //constexpr FactoredMultivector & operator=(FactoredMultivector<OtherScalarType, OtherFactoringProductType> const &) = default;

        constexpr MetricSpaceType const & space() const noexcept {
            return space_;
        }

        constexpr ScalarType const & scalar() const noexcept {
            return scalar_;
        }

        constexpr decltype(auto) factors_in_actual_metric() const noexcept {
            return detail::evaluate(detail::from_signed_to_actual_metric(space_, factors_in_signed_metric_));
        }

        constexpr auto const & factors_in_signed_metric() const noexcept {
            return factors_in_signed_metric_;
        }

        constexpr IndexType factors_count() const noexcept {
            return detail::cols(factors_in_signed_metric_);
        }

    private:

        template<typename OtherScalarType, typename OtherFactorsMatrixType>
        constexpr FactoredMultivector(MetricSpaceType const &space, OtherScalarType &&scalar, OtherFactorsMatrixType &&factors_in_signed_metric) noexcept :
            space_(space),
            scalar_(std::move(scalar)),
            factors_in_signed_metric_(std::move(factors_in_signed_metric)) {
        }

    private:

        MetricSpaceType const &space_;

        ScalarType scalar_;

        FactorsMatrixType factors_in_signed_metric_;

        friend struct detail::gp_impl;
        friend struct detail::op_impl;

        template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType> friend constexpr decltype(auto) addition(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) inverse(FactoredMultivector<SomeScalarType, GeometricProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename SecondScalarType, typename SecondMetricSpaceType, typename> friend constexpr decltype(auto) lcont(FirstScalarType const &, FactoredMultivector<SecondScalarType, GeometricProduct<SecondMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename FirstMetricSpaceType, typename SecondScalarType, typename> friend constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, GeometricProduct<FirstMetricSpaceType> > const &, SecondScalarType const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) reverse(FactoredMultivector<SomeScalarType, GeometricProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType> friend constexpr decltype(auto) subtraction(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) unary_minus(FactoredMultivector<SomeScalarType, GeometricProduct<SomeMetricSpaceType> > const &) noexcept;
    };

    template<typename ScalarType_, typename MetricSpaceType_>
    class FactoredMultivector<ScalarType_, OuterProduct<MetricSpaceType_> > final {
    private:

        using FactorsAndComplementMatrixType = detail::matrix_type_t<ScalarType_, MetricSpaceType_::DimensionsAtCompileTime, MetricSpaceType_::DimensionsAtCompileTime, MetricSpaceType_::MaxDimensionsAtCompileTime, MetricSpaceType_::MaxDimensionsAtCompileTime>;

    public:

        using MetricSpaceType = MetricSpaceType_;

        using IndexType = detail::index_type_t<FactorsAndComplementMatrixType>;
        using ScalarType = ScalarType_;

        constexpr FactoredMultivector(FactoredMultivector const &) = default;
        constexpr FactoredMultivector(FactoredMultivector &&) = default;

        constexpr FactoredMultivector(MetricSpaceType const &space) noexcept :
            FactoredMultivector(space, 0, detail::make_identity_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::MaxDimensionsAtCompileTime>(space.dimensions()), 0) {
        }

        template<typename OtherScalarType>
        constexpr FactoredMultivector(MetricSpaceType const &space, OtherScalarType &&scalar) noexcept :
            FactoredMultivector(space, std::move(scalar), detail::make_identity_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::MaxDimensionsAtCompileTime>(space.dimensions()), 0) {
        }

        //TODO [FUTURE] Specialization.
        //template<typename OtherScalarType, typename OtherFactoringProductType>
        //constexpr FactoredMultivector(FactoredMultivector<OtherScalarType, OtherFactoringProductType> const &other) noexcept;

        constexpr FactoredMultivector & operator=(FactoredMultivector const &) = default;
        constexpr FactoredMultivector & operator=(FactoredMultivector &&) = default;

        //TODO [FUTURE] Specialization.
        //template<typename OtherScalarType, typename OtherFactoringProductType>
        //constexpr FactoredMultivector & operator=(FactoredMultivector<OtherScalarType, OtherFactoringProductType> const &) = default;

        constexpr MetricSpaceType const & space() const noexcept {
            return space_;
        }

        constexpr ScalarType const & scalar() const noexcept {
            return scalar_;
        }

        constexpr decltype(auto) factors_in_actual_metric() const noexcept {
            return detail::evaluate(detail::from_signed_to_actual_metric(space_, factors_in_signed_metric()));
        }

        constexpr decltype(auto) factors_in_signed_metric() const noexcept {
            return detail::block<MetricSpaceType::DimensionsAtCompileTime, Dynamic>(factors_and_complement_in_signed_metric_, 0, 0, space_.dimensions(), factors_count_);
        }

        constexpr auto const & factors_and_complement_in_signed_metric() const noexcept {
            return factors_and_complement_in_signed_metric_;
        }

        constexpr IndexType factors_count() const noexcept {
            return factors_count_;
        }

    private:

        template<typename OtherScalarType, typename OtherFactorsAndComplementMatrixType, typename OtherIndexType>
        constexpr FactoredMultivector(MetricSpaceType const &space, OtherScalarType &&scalar, OtherFactorsAndComplementMatrixType &&factors_and_complement_in_signed_metric, OtherIndexType &&factors_count) noexcept :
            space_(space),
            scalar_(std::move(scalar)),
            factors_and_complement_in_signed_metric_(std::move(factors_and_complement_in_signed_metric)),
            factors_count_(std::move(factors_count)) {
        }

    private:

        MetricSpaceType const &space_;

        ScalarType scalar_;

        FactorsAndComplementMatrixType factors_and_complement_in_signed_metric_;
        IndexType factors_count_;

        friend struct detail::gp_impl;
        friend struct detail::op_impl;

        template<typename SomeScalarType, typename SomeMetricSpaceType> friend decltype(auto) e(SomeMetricSpaceType const &, DefaultIndexType) noexcept;
        template<typename SomeMetricSpaceType, typename SomeScalarType, typename> friend decltype(auto) scalar(SomeMetricSpaceType const &, SomeScalarType const &) noexcept;
        template<typename SomeMetricSpaceType, typename... ScalarTypes, typename> friend decltype(auto) vector(SomeMetricSpaceType const &, ScalarTypes &&...) noexcept;
        template<typename SomeMetricSpaceType, typename IteratorType, typename> friend decltype(auto) vector(SomeMetricSpaceType const &, IteratorType, IteratorType) noexcept;

        template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType> friend constexpr decltype(auto) addition(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) dual(FactoredMultivector<SomeScalarType, OuterProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) inverse(FactoredMultivector<SomeScalarType, OuterProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename SecondScalarType, typename SecondMetricSpaceType, typename> friend constexpr decltype(auto) lcont(FirstScalarType const &, FactoredMultivector<SecondScalarType, OuterProduct<SecondMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename FirstMetricSpaceType, typename SecondScalarType, typename> friend constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, OuterProduct<FirstMetricSpaceType> > const &, SecondScalarType const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) reverse(FactoredMultivector<SomeScalarType, OuterProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType> friend constexpr decltype(auto) subtraction(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) unary_minus(FactoredMultivector<SomeScalarType, OuterProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) undual(FactoredMultivector<SomeScalarType, OuterProduct<SomeMetricSpaceType> > const &) noexcept;
    };

    template<typename ScalarType, typename FactoringProductType>
    struct is_multivector<FactoredMultivector<ScalarType, FactoringProductType> > :
        std::true_type {
    };

}

#endif // __TBGAL_MULTIVECTOR_HPP__
