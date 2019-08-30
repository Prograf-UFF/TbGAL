#ifndef __TBGAL_LEFT_CONTRACTION_HPP__
#define __TBGAL_LEFT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstFactoringProduct, typename FirstSquareMatrixType, typename SecondFactoringProduct, typename SecondSquareMatrixType>
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &arg2) noexcept {
        return UNDUAL(OP(arg1, DUAL(arg2)));
    }

    template<typename FirstFactoringProduct, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        return scalar(arg1.factors_count() == 0 ? arg1.scalar() * arg2 : 0, arg2.space());
    }

    template<typename FirstType, typename SecondFactoringProduct, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) LCONT(FirstType const &arg1, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &arg2) noexcept {
        using ResultingSquareMatrixType = detail::promote_to_common_scalar_type_t<SecondSquareMatrixType, FirstType>;
        auto resulting_scalar = arg1 * arg2.scalar();
        if (resulting_scalar != 0) {
            return FactoredMultivector<SecondFactoringProduct, ResultingSquareMatrixType>(resulting_scalar, arg2.factors(), arg2.factors_count(), arg2.space());
        }
        else {
            return FactoredMultivector<SecondFactoringProduct, ResultingSquareMatrixType>(0, arg2.space());
        }
    }

    template<typename FirstType, typename SecondType, typename = std::enable_if_t<!(detail::is_multivector_v<FirstType> || detail::is_multivector_v<SecondType>)> >
    constexpr decltype(auto) LCONT(FirstType const &arg1, SecondType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_LEFT_CONTRACTION_HPP__
