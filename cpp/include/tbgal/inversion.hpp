#ifndef __TBGAL_INVERSION_HPP__
#define __TBGAL_INVERSION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) inverse(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = std::common_type_t<ScalarType, typename MetricSpaceType::ScalarType>;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        if (arg.factors_count() > 0) {
            return ResultingFactoredMultivectorType(arg.space(), 1 / (arg.scalar() * detail::metric_factor(arg.space(), arg.factors_in_signed_metric())), detail::reverse_columns(arg.factors_in_signed_metric()));
        }
        else {
            return ResultingFactoredMultivectorType(arg.space(), 1 / arg.scalar(), arg.factors_in_signed_metric());
        }
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) inverse(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = std::common_type_t<ScalarType, typename MetricSpaceType::ScalarType>;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        if (arg.factors_count() > 0) {
            auto factors_tuple = detail::from_outer_to_geometric_factors(arg.space(), arg.factors_in_signed_metric());
            return ResultingFactoredMultivectorType(arg.space(), (((arg.factors_count() * (arg.factors_count() - 1)) & 2) ? -1 : 1) / (arg.scalar() * std::get<0>(factors_tuple) * std::get<0>(factors_tuple) * detail::metric_factor(arg.space(), std::get<1>(factors_tuple))), arg.factors_and_complement_in_signed_metric(), arg.factors_count());
        }
        else {
            return ResultingFactoredMultivectorType(arg.space(), 1 / arg.scalar(), arg.factors_and_complement_in_signed_metric(), 0);
        }
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type inverse(Type const &arg) noexcept {
        return Type(1) / arg;
    }

}

#endif // __TBGAL_INVERSION_HPP__
