#ifndef __TBGAL_REVERSION_HPP__
#define __TBGAL_REVERSION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) reverse(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), arg.scalar(), detail::reverse_columns(arg.factors_in_signed_metric()));
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) reverse(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), ((arg.factors_count() * (arg.factors_count() - 1)) & 2) ? -arg.scalar() : arg.scalar(), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type reverse(Type const &arg) noexcept {
        return arg;
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) operator~(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return reverse(arg);
    }

}

#endif // __TBGAL_REVERSION_HPP__
