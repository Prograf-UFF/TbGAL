#ifndef __TBGAL_REVERSE_NORM_HPP__
#define __TBGAL_REVERSE_NORM_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) rnorm_sqr(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) {
        if (arg.factors_count() > 0) {
            return arg.scalar() * arg.scalar() * detail::metric_factor(arg.space(), arg.factors_in_signed_metric());
        }
        else {
            return arg.scalar() * arg.scalar();
        }
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) rnorm_sqr(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) {
        if (arg.factors_count() > 0) {
            auto factors_tuple = detail::from_outer_to_geometric_factors(arg.space(), arg.factors_in_signed_metric());
            return arg.scalar() * arg.scalar() * std::get<0>(factors_tuple) * std::get<0>(factors_tuple) * detail::metric_factor(arg.space(), std::get<1>(factors_tuple));
        }
        else {
            return arg.scalar() * arg.scalar();
        }
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr decltype(auto) rnorm_sqr(Type const &arg) {
        return arg * arg;
    }

    template<typename Type>
    constexpr decltype(auto) rnorm(Type const &arg) {
        return sqrt(rnorm_sqr(arg));
    }

}

#endif // __TBGAL_REVERSE_NORM_HPP__
