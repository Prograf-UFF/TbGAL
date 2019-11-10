#ifndef __TBGAL_REVERSE_NORM_HPP__
#define __TBGAL_REVERSE_NORM_HPP__

namespace tbgal {

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) rnorm_sqr(FactoredMultivector<ScalarType, FactoringProductType> const &arg) {
        return arg.scalar() * arg.scalar() * detail::metric_factor(arg.space(), arg.factors_in_signed_metric());
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
