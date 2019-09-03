#ifndef __TBGAL_UNARY_MINUS_HPP__
#define __TBGAL_UNARY_MINUS_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr FactoredMultivector<FactoringProductType, SquareMatrixType> UNARY_MINUS(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        using ResultingFactoringProductType = FactoringProductType;
        using ResultingSquareMatrixType = SquareMatrixType;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        return ResultingFactoredMultivectorType(arg.space(), -arg.scalar(), arg.factors(), arg.factors_count());
    }

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) operator-(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return UNARY_MINUS(arg);
    }

}

#endif // __TBGAL_UNARY_MINUS_HPP__
