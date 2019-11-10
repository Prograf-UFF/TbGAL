#ifndef __TBGAL_DUALIZATION_HPP__
#define __TBGAL_DUALIZATION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) dual(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg);
    //TODO [FUTURE] Implement dualization for the general case (when the input multivector is not a blade).

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) dual(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(
            arg.space(),
            (((arg.factors_count() * (arg.factors_count() - 1) + arg.space().dimensions() * (arg.space().dimensions() - 1)) & 2) ? -arg.scalar() : arg.scalar()) * detail::determinant(arg.factors_and_complement_in_signed_metric()),
            detail::apply_signed_metric(arg.space(), detail::split_columns_and_swap(arg.factors_and_complement_in_signed_metric(), arg.factors_count())),
            arg.space().dimensions() - arg.factors_count()
        );
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) undual(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg);
    //TODO [FUTURE] Implement undualization for the general case (when the input multivector is not a blade).

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) undual(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto factors_and_complement_in_signed_metric = detail::evaluate(detail::apply_signed_metric(arg.space(), arg.factors_and_complement_in_signed_metric()));
        return ResultingFactoredMultivectorType(
            arg.space(),
            (((arg.factors_count() * (arg.factors_count() - 1)) & 2) ? -arg.scalar() : arg.scalar()) * detail::determinant(factors_and_complement_in_signed_metric),
            detail::split_columns_and_swap(factors_and_complement_in_signed_metric, arg.factors_count()),
            arg.space().dimensions() - arg.factors_count()
        );
    }

}

#endif // __TBGAL_DUALIZATION_HPP__
