#ifndef __TBGAL_INVERSION_HPP__
#define __TBGAL_INVERSION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) inverse(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = std::common_type_t<ScalarType, typename MetricSpaceType::ScalarType>;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), ((((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -1 : 1) / (arg.scalar() * detail::metric_factor(arg.space(), arg.factors_in_signed_metric())), arg.factors_in_signed_metric());
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) inverse(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = std::common_type_t<ScalarType, typename MetricSpaceType::ScalarType>;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), ((((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -1 : 1) / (arg.scalar() * detail::metric_factor(arg.space(), arg.factors_in_signed_metric())), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type inverse(Type const &arg) noexcept {
        return Type(1) / arg;
    }

}

#endif // __TBGAL_INVERSION_HPP__
