#ifndef __TBGAL_ADDITION_HPP__
#define __TBGAL_ADDITION_HPP__

namespace tbgal {

    //TODO Implementar associatividade.

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) ADD(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        using ResultingFactoringProductType = OuterProduct<typename FirstFactoringProductType::SpaceType>;
        using ResultingSquareMatrixType = common_type_t<FirstSquareMatrixType, SecondSquareMatrixType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        if (arg1.factors_count() == arg2.factors_count()) {
            if (arg1.factors_count() == 0) {
                return ResultingFactoredMultivectorType(arg1.space(), arg1.scalar() + arg2.scalar());
            }
            //TODO Tratar vetor, pseudovetor e pseudoscalar. Lembrar de normalizar os fatores.
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) ADD(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        if (arg1.factors_count() == 0) {
            return scalar(arg1.space(), arg1.scalar() + arg2);
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) ADD(FirstType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        if (arg2.factors_count() == 0) {
            return scalar(arg2.space(), arg1 + arg2.scalar());
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstType, typename SecondType, typename = std::enable_if_t<!(detail::is_multivector_v<FirstType> || detail::is_multivector_v<SecondType>)> >
    constexpr decltype(auto) ADD(FirstType const &arg1, SecondType const &arg2) noexcept {
        return arg1 + arg2;
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) operator+(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return ADD(arg1, arg2);
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) operator+(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        return ADD(arg1, arg2);
    }

    template<typename FirstType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) operator+(FirstType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return ADD(arg1, arg2);
    }

}

#endif // __TBGAL_ADDITION_HPP__
