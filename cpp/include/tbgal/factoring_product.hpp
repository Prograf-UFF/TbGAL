#ifndef __TBGAL_FACTORING_PRODUCT_HPP__
#define __TBGAL_FACTORING_PRODUCT_HPP__

namespace tbgal {

    template<typename MetricSpaceType_>
    struct GeometricProduct final {
    public:

        using MetricSpaceType = MetricSpaceType_;

    private:

        constexpr GeometricProduct() {
            // Nothing to do. It just avoid the instantiation of GeometricProduct<MetricSpaceType>.
        }
    };

    template<typename MetricSpaceType_>
    struct OuterProduct final {
    public:

        using MetricSpaceType = MetricSpaceType_;

    private:

        constexpr OuterProduct() {
            // Nothing to do. It just avoid the instantiation of OuterProduct<MetricSpaceType>.
        }
    };

}

#endif // __TBGAL_FACTORING_PRODUCT_HPP__
