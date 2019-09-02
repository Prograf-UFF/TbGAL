#ifndef __TBGAL_REVERSION_HPP__
#define __TBGAL_REVERSION_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg);

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        return FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType>(arg.space(), (((arg.factors() * (arg.factors() - 1)) >> 1) & 1) ? -arg.scalar() ? arg.scalar(), arg.factors(), arg.factors_count());
    }

    template<typename Type, typename = std::enable_if_t<!detail::is_multivector_v<Type> > >
    constexpr Type REVERSE(Type const &arg) noexcept {
        return arg;
    }

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) operator~(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return REVERSE(arg);
    }

}

#endif // __TBGAL_REVERSION_HPP__
