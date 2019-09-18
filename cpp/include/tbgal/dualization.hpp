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
            (((arg.factors_count() * (arg.space().dimensions() - 1)) & 1) ? -arg.scalar() : arg.scalar()) * detail::orthogonal_metric_factor(arg.space(), arg.factors(), arg.factors_count()),
            detail::split_columns_and_swap(arg.factors(), arg.factors_count()),
            arg.space().dimensions() - arg.factors_count()
        );
        //TODO [CHECK] Sign change.
    }

    template<typename ScalarType>
    constexpr decltype(auto) DUAL(ScalarType const &arg);
    //TODO [FUTURE] Implement dualization for the general case (when the input multivector is a native scalar value).

    template<typename Type>
    constexpr decltype(auto) UNDUAL(Type const &arg) noexcept {
        return DUAL(arg);
        //TODO [CHECK] Is the UNDUAL necessary?
    }

}

#endif // __TBGAL_DUALIZATION_HPP__
