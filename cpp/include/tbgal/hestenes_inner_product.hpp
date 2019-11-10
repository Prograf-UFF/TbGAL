#ifndef __TBGAL_HESTENES_INNER_PRODUCT_HPP__
#define __TBGAL_HESTENES_INNER_PRODUCT_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) hip(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingFactoredMultivectorType = decltype(lcont(arg1, arg2));
        if (arg1.factors_count() != 0 && arg2.factors_count() != 0) {
            if (arg1.factors_count() <= arg2.factors_count()) {
                return lcont(arg1, arg2);
            }
            return rcont(arg1, arg2);
        }
        return ResultingFactoredMultivectorType(*detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) hip(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &, SecondScalarType const &) noexcept {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) hip(FirstScalarType const &, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &) noexcept {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) hip(FirstScalarType const &, SecondScalarType const &) noexcept {
        return std::common_type_t<FirstScalarType, SecondScalarType>(0);
    }

}

#endif // __TBGAL_HESTENES_INNER_PRODUCT_HPP__
