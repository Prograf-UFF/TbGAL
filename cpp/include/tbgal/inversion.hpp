#ifndef __TBGAL_INVERSION_HPP__
#define __TBGAL_INVERSION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) INVERSE(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        //TODO [FUTURE] It could be faster!
        return GP(ScalarType(1) / GP(REVERSE(arg), arg).scalar(), REVERSE(arg));
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) INVERSE(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        //TODO [FUTURE] It could be faster!
        return OP(ScalarType(1) / LCONT(REVERSE(arg), arg).scalar(), REVERSE(arg));
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type INVERSE(Type const &arg) noexcept {
        return Type(1) / arg;
    }

}

#endif // __TBGAL_INVERSION_HPP__
