#ifndef __TBGAL_LEFT_CONTRACTION_HPP__
#define __TBGAL_LEFT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstFactoringProduct, typename FirstSquareMatrixType, typename SecondFactoringProduct, typename SecondSquareMatrixType>
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &arg2) noexcept {
        return UNDUAL(OP(arg1, DUAL(arg2)));
    }

    template<typename FirstFactoringProduct, typename FirstSquareMatrixType, typename SecondScalarType, typename = std::enable_if_t<!detail::is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) LCONT(FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType> const &arg1, SecondScalarType const &arg2) noexcept {
        using ResultingType = std::common_type_t<typename FactoredMultivector<FirstFactoringProduct, FirstSquareMatrixType>::ScalarType, SecondScalarType>;
        return arg1.factors_count() == 0 ? arg1.scalar() * arg2 : ResultingType(0);
    }

    template<typename FirstScalarType, typename SecondFactoringProduct, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) LCONT(FirstScalarType const &arg1, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &arg2) noexcept {
        using ResultingFactoringProductType = SecondFactoringProduct;
        using ResultingSquareMatrixType = detail::common_type_t<FirstScalarType, SecondSquareMatrixType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        auto resulting_scalar = arg1 * arg2.scalar();
        if (resulting_scalar != 0) {
            return ResultingFactoredMultivectorType(arg2.space(), resulting_scalar, arg2.factors(), arg2.factors_count());
        }
        else {
            return ResultingFactoredMultivectorType(arg2.space(), 0);
        }
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(detail::is_multivector_v<FirstScalarType> || detail::is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) LCONT(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_LEFT_CONTRACTION_HPP__
