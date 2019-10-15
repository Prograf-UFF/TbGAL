#ifndef __TBGAL_UNARY_MINUS_HPP__
#define __TBGAL_UNARY_MINUS_HPP__

namespace tbgal {

    template<typename ScalarType, typename FactoringProductType>
    constexpr FactoredMultivector<ScalarType, FactoringProductType> UNARY_MINUS(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = FactoringProductType;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(arg.space(), -arg.scalar(), arg.factors_in_signed_metric(), arg.factors_count());
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) operator-(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return UNARY_MINUS(arg);
    }

}

#endif // __TBGAL_UNARY_MINUS_HPP__
