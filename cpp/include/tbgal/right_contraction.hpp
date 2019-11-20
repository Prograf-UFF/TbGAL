#ifndef __TBGAL_RIGHT_CONTRACTION_HPP__
#define __TBGAL_RIGHT_CONTRACTION_HPP__

namespace tbgal {

    template<typename FirstScalarType, typename FirstFactoringProductType, typename SecondScalarType, typename SecondFactoringProductType>
    constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, FirstFactoringProductType> const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoredMultivectorType = decltype(undual(op(arg2, dual(arg1))));
        if (!is_blade(arg1) || !is_blade(arg2) || (arg1.factors_count() >= arg2.factors_count())) {
            if ((arg2.factors_count() * (arg1.factors_count() + 1)) & 1) {
                return -undual(op(arg2, dual(arg1)));
            }
            return undual(op(arg2, dual(arg1)));
        }
        return ResultingFactoredMultivectorType(*detail::space_ptr(arg1, arg2), 0);
    }

    template<typename FirstScalarType, typename FirstMetricSpaceType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, GeometricProduct<FirstMetricSpaceType> > const &arg1, SecondScalarType const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = GeometricProduct<FirstMetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1.scalar() * arg2;
        if (!is_zero(resulting_scalar)) {
            return ResultingFactoredMultivectorType(arg1.space(), resulting_scalar, arg1.factors_in_signed_metric(), arg1.factors_count());
        }
        return ResultingFactoredMultivectorType(arg1.space(), 0);
    }

    template<typename FirstScalarType, typename FirstMetricSpaceType, typename SecondScalarType, typename = std::enable_if_t<!is_multivector_v<SecondScalarType> > >
    constexpr decltype(auto) rcont(FactoredMultivector<FirstScalarType, OuterProduct<FirstMetricSpaceType> > const &arg1, SecondScalarType const &arg2) noexcept {
        using ResultingScalarType = std::common_type_t<FirstScalarType, SecondScalarType>;
        using ResultingFactoringProductType = OuterProduct<FirstMetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        auto resulting_scalar = arg1.scalar() * arg2;
        if (!is_zero(resulting_scalar)) {
            return ResultingFactoredMultivectorType(arg1.space(), resulting_scalar, arg1.factors_and_complement_in_signed_metric(), arg1.factors_count());
        }
        return ResultingFactoredMultivectorType(arg1.space(), 0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename SecondFactoringProductType, typename = std::enable_if_t<!is_multivector_v<FirstScalarType> > >
    constexpr decltype(auto) rcont(FirstScalarType const &arg1, FactoredMultivector<SecondScalarType, SecondFactoringProductType> const &arg2) noexcept {
        using ResultingType = std::common_type_t<FirstScalarType, SecondScalarType>;
        return arg2.factors_count() == 0 ? arg1 * arg2.scalar() : ResultingType(0);
    }

    template<typename FirstScalarType, typename SecondScalarType, typename = std::enable_if_t<!(is_multivector_v<FirstScalarType> || is_multivector_v<SecondScalarType>)> >
    constexpr decltype(auto) rcont(FirstScalarType const &arg1, SecondScalarType const &arg2) noexcept {
        return arg1 * arg2;
    }

}

#endif // __TBGAL_RIGHT_CONTRACTION_HPP__
