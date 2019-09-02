#ifndef __TBGAL_LEFT_CONTRACTION_HPP__
#define __TBGAL_LEFT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstFactoringProduct, typename FirstSquareMatrixType, typename SecondFactoringProduct, typename SecondSquareMatrixType>
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &arg2) noexcept {
        return UNDUAL(OP(arg1, DUAL(arg2)));
    }

    template<typename FirstFactoringProduct, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        return scalar(arg2.space(), arg1.factors_count() == 0 ? arg1.scalar() * arg2 : 0);
    }

    template<typename FirstType, typename SecondFactoringProduct, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) LCONT(FirstType const &arg1, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &arg2) noexcept {
        using ResultingFactoringProductType = SecondFactoringProduct;
        using ResultingSquareMatrixType = detail::common_type_t<FirstType, SecondSquareMatrixType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        auto resulting_scalar = arg1 * arg2.scalar();
        if (resulting_scalar != 0) {
            return ResultingFactoredMultivectorType(arg2.space(), resulting_scalar, arg2.factors(), arg2.factors_count());
        }
        else {
            return ResultingFactoredMultivectorType(arg2.space(), 0);
        }
    }

    template<typename FirstType, typename SecondType, typename = std::enable_if_t<!(detail::is_multivector_v<FirstType> || detail::is_multivector_v<SecondType>)> >
    constexpr decltype(auto) LCONT(FirstType const &arg1, SecondType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_LEFT_CONTRACTION_HPP__
