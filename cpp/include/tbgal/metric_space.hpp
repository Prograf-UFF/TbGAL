#ifndef __TBGAL_METRIC_SPACE_HPP__
#define __TBGAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename SymmetricMatrixType>
    class MetricSpace final {
    public:

        using IndexType = detail::index_type_t<SymmetricMatrixType>;
        
        constexpr IndexType DimensionsAtCompileTime = /*TODO Implementar*/;

        constexpr MetricSpace(MetricSpace const &) = default;
        constexpr MetricSpace(MetricSpace &&) = default;

        constexpr MetricSpace(SymmetricMatrixType &&metric_matrix) noexcept :
            metric_matrix_{ std::move(metric_matrix) } {
        }

        constexpr IndexType dimensions() const noexcept {
            return detail::rows(metric_matrix_);
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
