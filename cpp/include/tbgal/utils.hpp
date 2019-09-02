#ifndef __TBGAL_UTILS_HPP__
#define __TBGAL_UTILS_HPP__

namespace tbgal {

    template<typename MetricSpaceType, typename ScalarType, typename = std::enable_if_t<!detail::is_multivector_v<ScalarType> > >
    constexpr decltype(auto) scalar(MetricSpaceType const &space, ScalarType const &scalar) noexcept {
        return FactoredMultivector<OuterProduct<MetricSpaceType>, detail::identity_matrix_type_t<ScalarType, MetricSpaceType> >(space, scalar);
    }

    template<typename MetricSpaceType, typename FirstScalarType, typename... NextScalarTypes>
    constexpr decltype(auto) vector(MetricSpaceType const &space, FirstScalarType const &arg1, NextScalarTypes const &... args) noexcept {
        //TODO Implementar
    }

}

#endif // __TBGAL_UTILS_HPP__
