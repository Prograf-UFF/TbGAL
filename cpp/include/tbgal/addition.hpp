#ifndef __TBGAL_ADDITION_HPP__
#define __TBGAL_ADDITION_HPP__

namespace tbgal {

    //TODO [FUTURE] Implement associativity.

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) addition(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = OuterProduct<typename FirstFactoringProductType::MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        static_assert(std::is_same_v<typename FirstFactoringProductType::MetricSpaceType, typename SecondFactoringProductType::MetricSpaceType>, "The multivectors must have the same metric space.");
        assert(&arg1.space() == &arg2.space());
        if (arg1.factors_count() == arg2.factors_count()) {
            if (arg1.factors_count() == 0) {
                return ResultingFactoredMultivectorType(arg1.space(), arg1.scalar() + arg2.scalar());
            }
            //TODO [FUTURE] Handle vectors.
            //TODO [FUTURE] Handle pseudovetors.
            //TODO [FUTURE] Handle pseudoscalars.
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) addition(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        if (arg1.factors_count() == 0) {
            return scalar(arg1.space(), arg1.scalar() + arg2);
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) addition(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        if (arg2.factors_count() == 0) {
            return scalar(arg2.space(), arg1 + arg2.scalar());
        }
        throw NotSupportedError("The general case of summation of factored multivectors is not supported.");
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) addition(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 + arg2;
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) operator+(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return addition(arg1, arg2);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) operator+(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        return addition(arg1, arg2);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) operator+(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        return addition(arg1, arg2);
    }

}

#endif // __TBGAL_ADDITION_HPP__
