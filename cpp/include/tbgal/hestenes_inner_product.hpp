#ifndef __TBGAL_HESTENES_INNER_PRODUCT_HPP__
#define __TBGAL_HESTENES_INNER_PRODUCT_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProduct, typename SecondScalarType, typename SecondFactoringProduct>
    constexpr decltype(auto) HIP(FactoredMultivector<FirstScalarType, FirstFactoringProduct> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &arg2) noexcept {
        using ResultingFactoredMultivectorType = decltype(LCONT(arg1, arg2));
        if (arg1.factors_count() != 0 && arg2.factors_count() != 0) {
            if (arg1.factors_count() <= arg2.factors_count()) {
                return LCONT(arg1, arg2);
            }
            return RCONT(arg1, arg2);
        }
        return ResultingFactoredMultivectorType(*detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstFactoringProduct, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) HIP(FactoredMultivector<FirstScalarType, FirstFactoringProduct> const &, SecondScalarType const &) noexcept {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProduct, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) HIP(FirstScalarType const &, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &) noexcept {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) HIP(FirstScalarType const &, SecondScalarType const &) noexcept {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

}

#endif // __TBGAL_HESTENES_INNER_PRODUCT_HPP__
