#ifndef __TBGAL_DUALIZATION_HPP__
#define __TBGAL_DUALIZATION_HPP__

namespace tbgal {

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) DUAL(FactoredMultivector<GeometricProduct<MetricSpaceType>, SquareMatrixType> const &arg);
    //TODO [FUTURE] Implement dualization for the general case (when the input multivector is not a blade).

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) DUAL(FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = decltype(detail::split_columns_and_swap(arg.factors(), arg.factors_count()));
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        return ResultingFactoredMultivectorType(
            arg.space(),
            (((arg.factors_count() * (arg.factors_count() - 1) + arg.space().dimensions() * (arg.space().dimensions() - 1)) & 2) ? -arg.scalar() : arg.scalar()) * detail::determinant(arg.factors()) * detail::orthogonal_metric_factor(arg.space(), arg.factors(), arg.factors_count()),
            detail::split_columns_and_swap(arg.factors(), arg.factors_count()),
            arg.space().dimensions() - arg.factors_count()
        );
    }

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) UNDUAL(FactoredMultivector<GeometricProduct<MetricSpaceType>, SquareMatrixType> const &arg);
    //TODO [FUTURE] Implement undualization for the general case (when the input multivector is not a blade).

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) UNDUAL(FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = decltype(detail::split_columns_and_swap(arg.factors(), arg.factors_count()));
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        return ResultingFactoredMultivectorType(
            arg.space(),
            (((arg.factors_count() * (arg.factors_count() - 1)) & 2) ? -arg.scalar() : arg.scalar()) * detail::determinant(arg.factors()) * detail::orthogonal_metric_factor(arg.space(), arg.factors(), arg.factors_count()),
            detail::split_columns_and_swap(arg.factors(), arg.factors_count()),
            arg.space().dimensions() - arg.factors_count()
        );
    }

}

#endif // __TBGAL_DUALIZATION_HPP__
