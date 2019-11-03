#ifndef __TBGAL_REVERSION_HPP__
#define __TBGAL_REVERSION_HPP__

namespace tbgal {

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = GeometricProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        using ResultingFactorsMatrixType = typename ResultingFactoringProductType::FactorsMatrixType;
        using IndexType = detail::index_type_t<ResultingFactorsMatrixType>;
        auto const &input_factors = arg.factors_in_signed_metric();
        IndexType const rows = arg.space().dimensions();
        IndexType const cols = arg.factors_count();
        ResultingFactorsMatrixType resulting_factors = detail::make_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, Dynamic, MetricSpaceType::MaxDimensionsAtCompileTime, MetricSpaceType::MaxDimensionsAtCompileTime>(rows, cols);
        for (IndexType col = 0; col != cols; ++col) {
            detail::assign_block<MetricSpaceType::DimensionsAtCompileTime, 1>(input_factors, 0, col, resulting_factors, 0, cols - (col + 1));
        }
        return ResultingFactoredMultivectorType(arg.space(), arg.scalar(), resulting_factors, cols);
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr decltype(auto) REVERSE(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(
            arg.space(),
            (((arg.factors_count() * (arg.factors_count() - 1)) >> 1) & 1) ? -arg.scalar() : arg.scalar(),
            arg.factors_in_signed_metric(),
            arg.factors_count()
        );
    }

    template<typename Type, typename = std::enable_if_t<!is_multivector_v<Type> > >
    constexpr Type REVERSE(Type const &arg) noexcept {
        return arg;
    }

    template<typename ScalarType, typename FactoringProductType>
    constexpr decltype(auto) operator~(FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
        return REVERSE(arg);
    }

}

#endif // __TBGAL_REVERSION_HPP__
