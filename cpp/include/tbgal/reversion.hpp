#ifndef __TBGAL_REVERSION_HPP__
#define __TBGAL_REVERSION_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        using ResultingFactoringProductType = FactoringProductType;
        using ResultingSquareMatrixType = SquareMatrixType;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        return ResultingFactoredMultivectorType(arg.space(), (((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -arg.scalar() : arg.scalar(), arg.factors(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type REVERSE(Type const &arg) noexcept {
        return arg;
    }

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) operator~(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return REVERSE(arg);
    }

}

#endif // __TBGAL_REVERSION_HPP__
