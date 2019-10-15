#ifndef __TBGAL_MULTIVECTOR_HPP__
#define __TBGAL_MULTIVECTOR_HPP__

namespace tbgal {

    template<typename ScalarType_, typename MetricSpaceType_>
    class FactoredMultivector<ScalarType_, OuterProduct<MetricSpaceType_> > final {
    public:

        using MetricSpaceType = MetricSpaceType_;

        using SquareMatrixType = detail::matrix_type_t<ScalarType_, MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::DimensionsAtCompileTime>;

        using IndexType = detail::index_type_t<SquareMatrixType>;
        using ScalarType = ScalarType_;

        constexpr FactoredMultivector(FactoredMultivector const &) = default;
        constexpr FactoredMultivector(FactoredMultivector &&) = default;

        constexpr FactoredMultivector(MetricSpaceType const &space) noexcept :
            space_(space),
            scalar_(0),
            factors_in_signed_metric_(detail::make_identity_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime>(space.dimensions())),
            factors_count_(0) {
        }

        template<typename OtherScalarType>
        constexpr FactoredMultivector(MetricSpaceType const &space, OtherScalarType &&scalar) noexcept :
            space_(space),
            scalar_(std::move(scalar)),
            factors_in_signed_metric_(detail::make_identity_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime>(space.dimensions())),
            factors_count_(0) {
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
            return detail::from_signed_to_actual_metric(space_, factors_in_signed_metric_);
        }

        constexpr SquareMatrixType const & factors_in_signed_metric() const noexcept {
            return factors_in_signed_metric_;
        }

        constexpr IndexType const & factors_count() const noexcept {
            return factors_count_;
        }

    private:

        template<typename OtherScalarType, typename OtherSquareMatrixType, typename OtherIndexType>
        constexpr FactoredMultivector(MetricSpaceType const &space, OtherScalarType &&scalar, OtherSquareMatrixType &&factors_in_signed_metric, OtherIndexType &&factors_count) noexcept :
            space_(space),
            scalar_(std::move(scalar)),
            factors_in_signed_metric_(std::move(factors_in_signed_metric)),
            factors_count_(std::move(factors_count)) {
        }

    private:

        MetricSpaceType const &space_;

        ScalarType scalar_;

        SquareMatrixType factors_in_signed_metric_;
        IndexType factors_count_;

        template<bool AnyMultivectorType> friend struct detail::GP_impl;
        template<bool AnyMultivectorType> friend struct detail::OP_impl;

        template<typename SomeScalarType, typename SomeMetricSpaceType> friend decltype(auto) e(SomeMetricSpaceType const &, DefaultIndexType) noexcept;
        template<typename SomeMetricSpaceType, typename SomeScalarType, typename> friend decltype(auto) scalar(SomeMetricSpaceType const &, SomeScalarType const &) noexcept;
        template<typename SomeMetricSpaceType, typename... ScalarTypes> friend decltype(auto) vector(SomeMetricSpaceType const &, ScalarTypes &&...) noexcept;

        template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType> friend constexpr decltype(auto) ADD(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) DUAL(FactoredMultivector<SomeScalarType,OuterProduct<SomeMetricSpaceType> > const &) noexcept;
        template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProduct, typename> friend constexpr decltype(auto) LCONT(FirstScalarType const &, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &) noexcept;
        template<typename SomeScalarType, typename SomeFactoringProductType> friend constexpr decltype(auto) REVERSE(FactoredMultivector<SomeScalarType, SomeFactoringProductType> const &) noexcept;
        template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType> friend constexpr decltype(auto) SUB(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeFactoringProductType> friend constexpr FactoredMultivector<SomeScalarType, SomeFactoringProductType> UNARY_MINUS(FactoredMultivector<SomeScalarType, SomeFactoringProductType> const &) noexcept;
        template<typename SomeScalarType, typename SomeMetricSpaceType> friend constexpr decltype(auto) UNDUAL(FactoredMultivector<SomeScalarType, OuterProduct<SomeMetricSpaceType> > const &) noexcept;
    };

    template<typename ScalarType, typename FactoringProductType>
    struct is_multivector<FactoredMultivector<ScalarType, FactoringProductType> > :
        std::true_type {
    };

}

#endif // __TBGAL_MULTIVECTOR_HPP__
