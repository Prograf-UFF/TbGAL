#ifndef __TBGAL_LEFT_CONTRACTION_HPP__
#define __TBGAL_LEFT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) lcont(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingFactoredMultivectorType = decltype(undual(op(arg1, dual(arg2))));
        if (!is_blade(arg1) || !is_blade(arg2) || (arg1.factors_count() <= arg2.factors_count())) {
            return undual(op(arg1, dual(arg2)));
        }
        return ResultingFactoredMultivectorType(*detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) lcont(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, SecondScalarType const &arg2) noexcept {
        using ResultingType = std::common_type_t<FirstScalarType, SecondScalarType>;
        return arg1.factors_count() == 0 ? arg1.scalar() * arg2 : ResultingType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondMetricSpaceType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) lcont(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, GeometricProduct<SecondMetricSpaceType> > const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = GeometricProduct<SecondMetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1 * arg2.scalar();
        if (resulting_scalar != 0) {
            return ResultingFactoredMultivectorType(arg2.space(), resulting_scalar, arg2.factors_in_signed_metric());
        }
        return ResultingFactoredMultivectorType(arg2.space(), 0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondMetricSpaceType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) lcont(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, OuterProduct<SecondMetricSpaceType> > const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = OuterProduct<SecondMetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1 * arg2.scalar();
        if (resulting_scalar != 0) {
            return ResultingFactoredMultivectorType(arg2.space(), resulting_scalar, arg2.factors_and_complement_in_signed_metric(), arg2.factors_count());
        }
        return ResultingFactoredMultivectorType(arg2.space(), 0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) lcont(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_LEFT_CONTRACTION_HPP__
