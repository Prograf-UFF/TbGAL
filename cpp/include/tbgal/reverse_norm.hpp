#ifndef __TBGAL_REVERSE_NORM_HPP__
#define __TBGAL_REVERSE_NORM_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    constexpr decltype(auto) RNORM_SQR(FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) {
        return LCONT(arg, REVERSE(arg)).scalar();
    }

    template<typename Type, typename = std::enable_if_t<!detail::is_multivector_v<Type> > >
    constexpr decltype(auto) RNORM_SQR(Type const &arg) {
        return arg * arg;
    }

    template<typename Type>
    constexpr decltype(auto) RNORM(Type const &arg) {
        return sqrt(RNORM_SQR(arg));
    }

}

#endif // __TBGAL_REVERSE_NORM_HPP__
