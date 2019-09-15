#ifndef __TBGAL_ASSUMING_EUCLIDEAN3_HPP__
#define __TBGAL_ASSUMING_EUCLIDEAN3_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Euclidean3 {
        
        using MetricSpaceType = EuclideanMetricSpace<3>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(SPACE)

        static auto const e1 = e(1);
        static auto const e2 = e(2);
        static auto const e3 = e(3);

    }

}

#endif // __TBGAL_ASSUMING_EUCLIDEAN3_HPP__
