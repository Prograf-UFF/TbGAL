#ifndef __TBGAL_METRIC_SPACE_HPP__
#define __TBGAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename SymmetricMatrixType>
    class MetricSpace {
    public:

        using IndexType = detail::index_type_t<SymmetricMatrixType>;

        constexpr MetricSpace(MetricSpace const &) = default;
        constexpr MetricSpace(MetricSpace &&) = default;

        constexpr MetricSpace(SymmetricMatrixType &&metric_matrix) noexcept :
            metric_matrix_{ std::move(metric_matrix) } {
            std::assert(detail::is_symmetric(metric_matrix_));
        }

        constexpr IndexType dimensions() const noexcept {
            return detail::rows_count(metric_matrix_);
        }

    private:

        constexpr SymmetricMatrixType metric_matrix_;
    };

    namespace detail {

        template<typename SymmetricMatrixType>
        struct is_metric_space<MetricSpace<SymmetricMatrixType> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_METRIC_SPACE_HPP__
