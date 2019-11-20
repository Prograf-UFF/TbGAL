#ifndef __TBGAL_SCALAR_PRODUCT_HPP__
#define __TBGAL_SCALAR_PRODUCT_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement scalar product for the general case (when the input multivector is not a blade).

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) sp(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        static_assert(std::is_same_v<typename FirstFactoringProductType::MetricSpaceType, typename SecondFactoringProductType::MetricSpaceType>, "The multivectors must have the same metric space.");
        using MetricSpaceType = typename FirstFactoringProductType::MetricSpaceType;
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType, typename MetricSpaceType::ScalarType>;
        assert(is_blade(arg1) && is_blade(arg2));
        if (arg1.factors_count() == arg2.factors_count()) {
            return (((arg1.factors_count() * (arg1.factors_count() - 1)) & 2) ? -arg1.scalar() : arg1.scalar()) * arg2.scalar() * detail::determinant(detail::prod(detail::transpose(arg1.factors_in_signed_metric()), detail::apply_signed_metric(arg2.space(), arg2.factors_in_signed_metric())));
        }
        return ResultingScalarType(0);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) sp(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        assert(is_blade(arg1));
        if (arg1.factors_count() == 0) {
            return arg1.scalar() * arg2;
        }
        return ResultingScalarType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) sp(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        assert(is_blade(arg2));
        if (arg2.factors_count() == 0) {
            return arg1 * arg2.scalar();
        }
        return ResultingScalarType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) sp(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_SCALAR_PRODUCT_HPP__
