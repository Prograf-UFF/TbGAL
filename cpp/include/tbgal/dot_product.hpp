#ifndef __TBGAL_DOT_PRODUCT_HPP__
#define __TBGAL_DOT_PRODUCT_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProduct, typename SecondScalarType, typename SecondFactoringProduct>
    constexpr decltype(auto) DOT(FactoredMultivector<FirstScalarType, FirstFactoringProduct> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &arg2) noexcept {
        if (arg1.factors_count() <= arg2.factors_count()) {
            return LCONT(arg1, arg2);
        }
        return RCONT(arg1, arg2);
    }

    template<typename FirstScalarType, typename FirstFactoringProduct, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) DOT(FactoredMultivector<FirstScalarType, FirstFactoringProduct> const &arg1, SecondScalarType const &arg2) noexcept {
        return RCONT(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProduct, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) DOT(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &arg2) noexcept {
        return LCONT(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) DOT(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_DOT_PRODUCT_HPP__
