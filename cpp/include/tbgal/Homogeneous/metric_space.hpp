#ifndef __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
#define __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename ScalarType, std::int32_t BaseDimensionsAtCompileTime>
    class HomogeneousMetricSpace : public MetricSpace<detail::identity_matrix_type_t<ScalarType, BaseDimensionsAtCompileTime + 1> > {
    private:

        static_assert(BaseDimensionsAtCompileTime > 0, "Invalid number of base dimensions.");

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, BaseDimensionsAtCompileTime + 1> >;
    
    public:

        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        constexpr HomogeneousMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, BaseDimensionsAtCompileTime + 1>(BaseDimensionsAtCompileTime + 1)) {
        }
    };

    template<typename ScalarType>
    class HomogeneousMetricSpace<ScalarType, Dynamic> : public MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> > {
    private:

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> >;
    
    public:

        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        constexpr HomogeneousMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(1)) {
        }

        constexpr HomogeneousMetricSpace(std::int32_t base_dimensions) noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(base_dimensions + 1)) {
            std::assert(base_dimensions > 0);
        }

        constexpr void set_base_dimensions(std::int32_t base_dimensions) noexcept {
            std::assert(base_dimensions > 0);
            Super::metric_matrix_ = detail::make_identity_matrix<ScalarType, Dynamic>(base_dimensions + 1);
        }
    };

    namespace detail {

        template<typename ScalarType, std::int32_t BaseDimensionsAtCompileTime>
        struct is_metric_space<HomogeneousMetricSpace<ScalarType, BaseDimensionsAtCompileTime> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
