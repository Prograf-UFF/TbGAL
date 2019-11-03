#ifndef __TBGAL_LEFT_CONTRACTION_HPP__
#define __TBGAL_LEFT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProduct, typename SecondScalarType, typename SecondFactoringProduct>
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstScalarType, FirstFactoringProduct> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &arg2) noexcept {
        using ResultingFactoredMultivectorType = decltype(UNDUAL(OP(arg1, DUAL(arg2))));
        if (!is_blade(arg1) || !is_blade(arg2) || (arg1.factors_count() <= arg2.factors_count())) {
            return UNDUAL(OP(arg1, DUAL(arg2)));
        }
        return ResultingFactoredMultivectorType(*detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstFactoringProduct, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstScalarType, FirstFactoringProduct> const &arg1, SecondScalarType const &arg2) noexcept {
        using ResultingType = std::common_type_t<FirstScalarType, SecondScalarType>;
        return arg1.factors_count() == 0 ? arg1.scalar() * arg2 : ResultingType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProduct, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) LCONT(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProduct> const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = SecondFactoringProduct;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1 * arg2.scalar();
        if (resulting_scalar != 0) {
            return ResultingFactoredMultivectorType(arg2.space(), resulting_scalar, arg2.factors_in_signed_metric(), arg2.factors_count());
        }
        return ResultingFactoredMultivectorType(arg2.space(), 0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) LCONT(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_LEFT_CONTRACTION_HPP__
