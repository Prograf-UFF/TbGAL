#ifndef __TBGAL_REVERSION_HPP__
#define __TBGAL_REVERSION_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg);

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<GeometricProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = SquareMatrixType;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        if (is_blade(arg)) {
            return ResultingFactoredMultivectorType(arg.space(), (((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -arg.scalar() ? arg.scalar(), arg.factors(), arg.factors_count());
        }
        //TODO [FUTURE] Implement reversion for the general case (when the input multivector is not a blade).
        throw NotSupportedError("The general case of reversion of factored multivectors is not supported yet.");
    }

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = SquareMatrixType;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        return ResultingFactoredMultivectorType(arg.space(), (((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -arg.scalar() ? arg.scalar(), arg.factors(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!detail::is_multivector_v<Type> > >
    constexpr Type REVERSE(Type const &arg) noexcept {
        return arg;
    }

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) operator~(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return REVERSE(arg);
    }

}

#endif // __TBGAL_REVERSION_HPP__
