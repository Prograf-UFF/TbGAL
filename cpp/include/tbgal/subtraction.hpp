#ifndef __TBGAL_SUBTRACTION_HPP__
#define __TBGAL_SUBTRACTION_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement associativity.

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) SUB(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        using ResultingFactoringProductType = OuterProduct<typename FirstFactoringProductType::SpaceType>;
        using ResultingSquareMatrixType = detail::common_type_t<FirstSquareMatrixType, SecondSquareMatrixType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        static_assert(std::is_same_v<typename FirstFactoringProductType::SpaceType, typename SecondFactoringProductType::SpaceType>, "The multivectors must have the same metric space.");
        assert(&arg1.space() == &arg2.space());
        if (arg1.factors_count() == arg2.factors_count()) {
            if (arg1.factors_count() == 0) {
                return ResultingFactoredMultivectorType(arg1.space(), arg1.scalar() - arg2.scalar());
            }
            //TODO [FUTURE] Handle vectors.
            //TODO [FUTURE] Handle pseudovetors.
            //TODO [FUTURE] Handle pseudoscalars.
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondScalarType, typename = std::enable_if_t<!detail::is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) SUB(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondScalarType const &arg2) noexcept {
        if (arg1.factors_count() == 0) {
            return scalar(arg1.space(), arg1.scalar() - arg2);
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) SUB(FirstScalarType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        if (arg2.factors_count() == 0) {
            return scalar(arg2.space(), arg1 - arg2.scalar());
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(detail::is_multivector_v<FirstScalarType> || detail::is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) SUB(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 - arg2;
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType>
    constexpr decltype(auto) operator-(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return SUB(arg1, arg2);
    }

    template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondScalarType, typename = std::enable_if_t<!detail::is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator-(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &arg1, SecondScalarType const &arg2) noexcept {
        return SUB(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondFactoringProductType, typename SecondSquareMatrixType, typename = std::enable_if_t<!detail::is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator-(FirstScalarType const &arg1, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &arg2) noexcept {
        return SUB(arg1, arg2);
    }

}

#endif // __TBGAL_SUBTRACTION_HPP__
