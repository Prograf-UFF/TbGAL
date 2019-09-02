#ifndef __TBGAL_DUALIZATION_HPP__
#define __TBGAL_DUALIZATION_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) DUAL(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg);

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr decltype(auto) DUAL(FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        //TODO Verificar troca seletiva de sinal
        //TODO Verificar tipo de retorno
        return FactoredMultivector<OuterProduct<MetricSpaceType>, decltype(detail::split_columns_and_swap(arg.factors(), arg.factors_count()))>(
            arg.space(),
            ((arg.factors_count() * (arg.space().dimensions() - 1)) & 1) ? -arg.scalar() ? arg.scalar(),
            detail::split_columns_and_swap(arg.factors(), arg.factors_count()),
            arg.space().dimensions() - arg.factors_count()
        );
    }

    //TODO Especializar para FactoredMultivector<GeometricProduct<MetricSpaceType>, SquareMatrixType>

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) UNDUAL(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return DUAL(arg); //TODO Acho que não precisa do UNDUAL
    }

}

#endif // __TBGAL_DUALIZATION_HPP__
