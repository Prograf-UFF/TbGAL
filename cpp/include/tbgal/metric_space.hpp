#ifndef __TBGAL_METRIC_SPACE_HPP__
#define __TBGAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename MetricSpaceType>
    class MetricSpace {
    public:

        using IndexType = DefaultIndexType;

        virtual std::string const & basis_vector_str(IndexType) const noexcept = 0;

        virtual IndexType dimensions() const noexcept = 0;

    protected:

        constexpr MetricSpace() = default;
        constexpr MetricSpace(MetricSpace const &) = default;
        constexpr MetricSpace(MetricSpace &&) = default;

        constexpr MetricSpace & operator=(MetricSpace const &) = default;
        constexpr MetricSpace & operator=(MetricSpace &&) = default;
    };

    namespace detail {

        template<typename MetricSpaceType>
        struct is_metric_space<MetricSpace<MetricSpaceType> > :
            std::true_type {
        };

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) from_actual_to_orthogonal_metric(MetricSpaceType const &space, MatrixType const &factors) noexcept {
            return from_actual_to_orthogonal_metric_impl<MetricSpaceType>::eval(space, factors);
        }

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) from_orthogonal_to_actual_metric(MetricSpaceType const &space, MatrixType const &factors) noexcept {
            return from_orthogonal_to_actual_metric_impl<MetricSpaceType>::eval(space, factors);
        }

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) orthogonal_metric_factor(MetricSpaceType const &space, MatrixType const &factors, DefaultIndexType factors_count) noexcept {
            return orthogonal_metric_factor_impl<MetricSpaceType>::eval(space, factors, factors_count);
        }

    }

}

#endif // __TBGAL_METRIC_SPACE_HPP__
