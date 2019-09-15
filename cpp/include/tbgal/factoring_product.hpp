#ifndef __TBGAL_FACTORING_PRODUCT_HPP__
#define __TBGAL_FACTORING_PRODUCT_HPP__

namespace tbgal {

    template<typename MetricSpaceType_>
    struct GeometricProduct final {
    public:

        static_assert(detail::is_metric_space_v<MetricSpaceType_>, "Invalid MetricSpaceType.");

        using MetricSpaceType = MetricSpaceType_;

    private:

        constexpr GeometricProduct() {
            // Nothing to do. It just avoid the instantiation of GeometricProduct<MetricSpaceType>.
        }
    };

    template<typename MetricSpaceType_>
    struct OuterProduct final {
    public:

        static_assert(detail::is_metric_space_v<MetricSpaceType_>, "Invalid MetricSpaceType.");

        using MetricSpaceType = MetricSpaceType_;

    private:

        constexpr OuterProduct() {
            // Nothing to do. It just avoid the instantiation of OuterProduct<MetricSpaceType>.
        }
    };

    namespace detail {

        template<typename MetricSpaceType>
        struct is_factoring_product<GeometricProduct<MetricSpaceType> > :
            std::true_type {
        };

        template<typename MetricSpaceType>
        struct is_factoring_product<OuterProduct<MetricSpaceType> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_FACTORING_PRODUCT_HPP__
