#ifndef __TBGAL_SUBTRACTION_HPP__
#define __TBGAL_SUBTRACTION_HPP__

namespace tbgal {

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) SUB(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        using ResultingFactoringProductType = OuterProduct<typename FirstFactoringProductType::SpaceType>;
        using ResultingSquareMatrixType = common_matrix_type_t<FirstSquareMatrixType, SecondSquareMatrixType>;
        std::assert(&arg1.space() == &arg2.space());
        if (arg1.factors_count() == arg2.factors_count()) {
            if (arg1.factors_count() == 0) {
                return FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>(arg1.scalar() - arg2.scalar(), arg1.space());
            }
            //TODO Tratar vetor, pseudovetor e pseudoscalar. Lembrar de normalizar os fatores.
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) SUB(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        if (arg1.factors_count() != 0) {
            throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
        }
        return scalar(arg1.scalar() - arg2, arg1.space());
    }

    template<typename FirstType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) SUB(FirstType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        if (arg2.factors_count() != 0) {
            throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
        }
        return scalar(arg1 - arg2.scalar(), arg2.space());
    }

    template<typename FirstType, typename SecondType, typename = std::enable_if_t<!(detail::is_multivector_v<FirstType> || detail::is_multivector_v<SecondType>)> >
    constexpr decltype(auto) SUB(FirstType const &arg1, SecondType const &arg2) noexcept {
        return arg1 - arg2;
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) operator-(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return SUB(arg1, arg2);
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondType, typename = std::enable_if_t<!detail::is_multivector_v<SecondType> > >
    constexpr decltype(auto) operator-(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondType const &arg2) noexcept {
        return SUB(arg1, arg2);
    }

    template<typename FirstType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstType> > >
    constexpr decltype(auto) operator-(FirstType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return SUB(arg1, arg2);
    }

}

#endif // __TBGAL_SUBTRACTION_HPP__
