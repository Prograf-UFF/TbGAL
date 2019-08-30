#ifndef __TBGAL_UNARY_PLUS_HPP__
#define __TBGAL_UNARY_PLUS_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr FactoredMultivector<FactoringProductType, SquareMatrixType> UNARY_PLUS(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return arg;
    }

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr FactoredMultivector<FactoringProductType, SquareMatrixType> operator+(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
        return arg;
    }

}

#endif // __TBGAL_UNARY_PLUS_HPP__
