#ifndef __TBGAL_UTILS_HPP__
#define __TBGAL_UTILS_HPP__

namespace tbgal {

    template<typename Type>
    constexpr bool is_blade(Type const &arg) noexcept {
        return true;
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr bool is_blade(FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType> > const &arg) noexcept {
        return arg.factors_count() == 0 || arg.factors_count() == 1 || arg.scalar() == 0;
        //TODO [FUTURE] Implement the general case.
    }

    template<typename ScalarType, typename MetricSpaceType>
    constexpr bool is_blade(FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType> > const &arg) noexcept {
        return true;
    }

    template<typename MetricSpaceType, typename ScalarType, typename = std::enable_if_t<!is_multivector_v<std::remove_cv_t<std::remove_reference_t<ScalarType> > > > >
    decltype(auto) scalar(MetricSpaceType const &space, ScalarType &&value) noexcept {
        using ResultingScalarType = std::remove_cv_t<std::remove_reference_t<ScalarType> >;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        return ResultingFactoredMultivectorType(space, std::move(value));
    }

    template<typename MetricSpaceType, typename... ScalarTypes, typename = std::enable_if_t<std::disjunction_v<std::bool_constant<!detail::is_iterator_v<ScalarTypes> >...> > >
    decltype(auto) vector(MetricSpaceType const &space, ScalarTypes &&... coords) noexcept {
        using ResultingScalarType = std::common_type_t<std::remove_cv_t<std::remove_reference_t<ScalarTypes> >...>;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        static_assert(MetricSpaceType::DimensionsAtCompileTime == Dynamic || MetricSpaceType::DimensionsAtCompileTime == sizeof...(ScalarTypes), "Invalid number of coordinates.");
        assert(space.dimensions() == sizeof...(ScalarTypes));
        auto input = detail::evaluate(detail::from_actual_to_signed_metric(space, detail::fill_column_matrix(std::move(coords)...)));
        auto qr_tuple = detail::qr_orthogonal_matrix(input);
        if (std::get<1>(qr_tuple) == 1) {
            auto const &matrix_q = std::get<0>(qr_tuple);
            return ResultingFactoredMultivectorType(
                space,
                detail::coeff(detail::prod_block<1, MetricSpaceType::DimensionsAtCompileTime>(detail::transpose(matrix_q), 0, 0, 1, space.dimensions(), input), 0, 0),
                matrix_q,
                1
            );
        }
        return ResultingFactoredMultivectorType(space, 0);
    }

    template<typename MetricSpaceType, typename IteratorType, typename = std::enable_if_t<detail::is_iterator_v<IteratorType> > >
    decltype(auto) vector(MetricSpaceType const &space, IteratorType begin, IteratorType end) noexcept {
        using ResultingScalarType = std::remove_cv_t<std::remove_reference_t<typename std::iterator_traits<IteratorType>::value_type> >;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        assert(space.dimensions() == std::distance(begin, end));
        auto input = detail::evaluate(detail::from_actual_to_signed_metric(space, detail::fill_column_matrix_using_iterators<MetricSpaceType::DimensionsAtCompileTime, MetricSpaceType::MaxDimensionsAtCompileTime>(begin, end)));
        auto qr_tuple = detail::qr_orthogonal_matrix(input);

        if (std::get<1>(qr_tuple) == 1) {
            auto const &matrix_q = std::get<0>(qr_tuple);
            return ResultingFactoredMultivectorType(
                space,
                detail::coeff(detail::prod_block<1, MetricSpaceType::DimensionsAtCompileTime>(detail::transpose(matrix_q), 0, 0, 1, space.dimensions(), input), 0, 0),
                matrix_q,
                1
            );
        }
        return ResultingFactoredMultivectorType(space, 0);
    }

    template<typename ScalarType, typename MetricSpaceType>
    decltype(auto) e(MetricSpaceType const &space, DefaultIndexType index) noexcept {
        using ResultingScalarType = ScalarType;
        using ResultingFactoringProductType = OuterProduct<MetricSpaceType>;
        using ResultingFactoredMultivectorType = FactoredMultivector<ResultingScalarType, ResultingFactoringProductType>;
        assert(1 <= index && index <= space.dimensions());
        auto aux = detail::make_matrix<ScalarType, MetricSpaceType::DimensionsAtCompileTime, 1, MetricSpaceType::MaxDimensionsAtCompileTime, 1>(space.dimensions(), 1);
        detail::assign_block<Dynamic, 1>(detail::make_zero_matrix<ScalarType, Dynamic, 1, MetricSpaceType::DimensionsAtCompileTime, 1>(index - 1, 1), aux, 0, 0, index - 1, 1);
        detail::coeff(aux, index - 1, 0) = 1;
        detail::assign_block<Dynamic, 1>(detail::make_zero_matrix<ScalarType, Dynamic, 1, MetricSpaceType::DimensionsAtCompileTime, 1>(space.dimensions() - index, 1), aux, index, 0, space.dimensions() - index, 1);
        auto input = detail::evaluate(detail::from_actual_to_signed_metric(space, aux));
        auto qr_tuple = detail::qr_orthogonal_matrix(input);
        auto const &matrix_q = std::get<0>(qr_tuple);
        return ResultingFactoredMultivectorType(
            space,
            detail::coeff(detail::prod_block<1, MetricSpaceType::DimensionsAtCompileTime>(detail::transpose(matrix_q), 0, 0, 1, space.dimensions(), input), 0, 0),
            matrix_q,
            1
        );
    }

}

#endif // __TBGAL_UTILS_HPP__
