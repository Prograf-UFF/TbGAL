#ifndef __TBGAL_UTILS_HPP__
#define __TBGAL_UTILS_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType, typename = std::enable_if_t<!detail::is_multivector_v<ScalarType> > >
    constexpr decltype(auto) scalar(ScalarType const &scalar, MetricSpaceType const &space) noexcept {
        return FactoredMultivector<OuterProduct<MetricSpaceType>, detail::identity_matrix_t<ScalarType, MetricSpaceType> >(scalar, space);
    }

}

#endif // __TBGAL_UTILS_HPP__
