#ifndef __TBGAL_UNARY_MINUS_HPP__
#define __TBGAL_UNARY_MINUS_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) unary_minus(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), -arg.scalar(), arg.factors_in_signed_metric());
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) unary_minus(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), -arg.scalar(), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) operator-(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return unary_minus(arg);
    }

}

#endif // __TBGAL_UNARY_MINUS_HPP__
