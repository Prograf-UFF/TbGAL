#ifndef __TBGAL_METRIC_SPACE_HPP__
#define __TBGAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename MetricSpaceType>
    class MetricSpace {
    public:

        using ActualType = MetricSpaceType;
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

    }

}

#endif // __TBGAL_METRIC_SPACE_HPP__
