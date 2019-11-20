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

    template<typename MetricSpaceType>
    struct is_metric_space<MetricSpace<MetricSpaceType> > :
        std::true_type {
    };

    namespace detail {

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) apply_signed_metric(MetricSpaceType const &space, MatrixType const &factors_in_signed_metric) noexcept {
            return apply_signed_metric_impl<MetricSpaceType>::eval(space, factors_in_signed_metric);
        }

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) from_actual_to_signed_metric(MetricSpaceType const &space, MatrixType &&factors_in_actual_metric) noexcept {
            return from_actual_to_signed_metric_impl<MetricSpaceType>::eval(space, std::move(factors_in_actual_metric));
        }

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) from_signed_to_actual_metric(MetricSpaceType const &space, MatrixType &&factors_in_signed_metric) noexcept {
            return from_signed_to_actual_metric_impl<MetricSpaceType>::eval(space, std::move(factors_in_signed_metric));
        }

        template<typename MetricSpaceType, typename MatrixType>
        constexpr decltype(auto) metric_factor(MetricSpaceType const &space, MatrixType const &factors_in_signed_metric) noexcept {
            return metric_factor_impl<MetricSpaceType>::eval(space, factors_in_signed_metric);
        }

    }

}

#endif // __TBGAL_METRIC_SPACE_HPP__
