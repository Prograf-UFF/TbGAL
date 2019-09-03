#ifndef __TBGAL_ASSUMING_EUCLIDEAN2_HPP__
#define __TBGAL_ASSUMING_EUCLIDEAN2_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Euclidean2 {
        
        using MetricSpaceType = EuclideanMetricSpace<DefaultScalarType, 2>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(SPACE)

        static auto const e1 = e(1);
        static auto const e2 = e(2);

    }

}

#endif // __TBGAL_ASSUMING_EUCLIDEAN2_HPP__
