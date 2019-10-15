#ifndef __TBGAL_REVERSION_HPP__
#define __TBGAL_REVERSION_HPP__

namespace tbgal {

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = FactoringProductType;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(
            arg.space(),
            (((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -arg.scalar() : arg.scalar(),
            arg.factors_in_signed_metric(),
            arg.factors_count()
        );
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type REVERSE(Type const &arg) noexcept {
        return arg;
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) operator~(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return REVERSE(arg);
    }

}

#endif // __TBGAL_REVERSION_HPP__
