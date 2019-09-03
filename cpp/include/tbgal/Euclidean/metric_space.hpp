#ifndef __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
#define __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename ScalarType, std::int32_t DimensionsAtCompileTime>
    class EuclideanMetricSpace : public MetricSpace<detail::identity_matrix_type_t<ScalarType, DimensionsAtCompileTime> > {
    private:

        static_assert(DimensionsAtCompileTime > 0, "Invalid number of dimensions.");

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, DimensionsAtCompileTime> >;
    
    public:

        constexpr EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        constexpr EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        constexpr EuclideanMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, DimensionsAtCompileTime>(DimensionsAtCompileTime)) {
        }
    };

    template<typename ScalarType>
    class EuclideanMetricSpace<ScalarType, Dynamic> : public MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> > {
    private:

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> >;
    
    public:

        constexpr EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        constexpr EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        constexpr EuclideanMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(0)) {
        }

        constexpr EuclideanMetricSpace(std::int32_t dimensions) noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(dimensions)) {
            std::assert(dimensions > 0);
        }

        constexpr void set_dimensions(std::int32_t dimensions) noexcept {
            std::assert(dimensions > 0);
            Super::metric_matrix_ = detail::make_identity_matrix<ScalarType, Dynamic>(dimensions);
        }
    };

    namespace detail {

        template<typename ScalarType, std::int32_t DimensionsAtCompileTime>
        struct is_metric_space<EuclideanMetricSpace<ScalarType, DimensionsAtCompileTime> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
