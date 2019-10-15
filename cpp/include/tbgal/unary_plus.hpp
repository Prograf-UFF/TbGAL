#ifndef __TBGAL_UNARY_PLUS_HPP__
#define __TBGAL_UNARY_PLUS_HPP__

namespace tbgal {

    template<typename ScalarType, typename FactoringProductType>
    constexpr FactoredMultivector<ScalarType, FactoringProductType> UNARY_PLUS(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return arg;
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr FactoredMultivector<ScalarType, FactoringProductType> operator+(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return arg;
    }

}

#endif // __TBGAL_UNARY_PLUS_HPP__
