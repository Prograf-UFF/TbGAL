#ifndef __TBGAL_UTILS_HPP__
#define __TBGAL_UTILS_HPP__

namespace tbgal {

    template<typename Type>
    constexpr bool is_blade(Type const &arg) noexcept {
        return true;
    }

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr bool is_blade(FactoredMultivector<GeometricProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        return arg.factors_count() == 0 || arg.factors_count() == 1 || arg.scalar() == 0;
        //TODO [FUTURE] Implement the general case.
    }

    template<typename MetricSpaceType, typename SquareMatrixType>
    constexpr bool is_blade(FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        return true;
    }

    template<typename MetricSpaceType, typename ScalarType, typename = std::enable_if_t<!is_multivector_v<std::remove_cv_t<std::remove_reference_t<ScalarType> > > > >
    decltype(auto) scalar(MetricSpaceType const &space, ScalarType &&value) noexcept {
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = detail::matrix_type_t<std::remove_cv_t<std::remove_reference_t<ScalarType> >, MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::DimensionsAtCompileTime>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        return ResultingFactoredMultivectorType(space, std::move(value));
    }

    template<typename MetricSpaceType, typename... ScalarTypes>
    decltype(auto) vector(MetricSpaceType const &space, ScalarTypes &&... coords) noexcept {
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = detail::matrix_type_t<std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>, MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::DimensionsAtCompileTime>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        static_assert(MetricSpaceType::DimensionsAtCompileTime == Dynamic || MetricSpaceType::DimensionsAtCompileTime == sizeof...(ScalarTypes), "Invalid number of coordinates.");
        assert(space.dimensions() == sizeof...(ScalarTypes));
        auto input = detail::from_actual_to_signed_metric(space, detail::fill_column_matrix(std::move(coords)...));
        auto qr_tuple = detail::qr_decomposition(input);
        if (std::get<2>(qr_tuple) == 1) {
            auto const &matrix_q = std::get<0>(qr_tuple);
            return ResultingFactoredMultivectorType(
                space,
                detail::determinant(detail::prod(detail::transpose(detail::left_columns(matrix_q, 1)), input)),
                matrix_q,
                1
            );
        }
        return ResultingFactoredMultivectorType(space, 0);
    }

    template<typename ScalarType, typename MetricSpaceType>
    decltype(auto) e(MetricSpaceType const &space, DefaultIndexType index) noexcept {
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingSquareMatrixType = detail::matrix_type_t<ScalarType, MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::DimensionsAtCompileTime>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingFactoringProductType, ResultingSquareMatrixType>;
        assert(1 <= index && index <= space.dimensions());
        auto aux = detail::make_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, 1>(space.dimensions(), 1);
        auto input = detail::from_actual_to_signed_metric(space, detail::copy_columns<1>(detail::make_identity_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime>(space.dimensions()), index - 1, aux, 0, 1));
        auto qr_tuple = detail::qr_decomposition(input);
        auto const &matrix_q = std::get<0>(qr_tuple);
        return ResultingFactoredMultivectorType(
            space,
            detail::determinant(detail::prod(detail::transpose(detail::left_columns(matrix_q, 1)), input)),
            matrix_q,
            1
        );
    }

}

#endif // __TBGAL_UTILS_HPP__
