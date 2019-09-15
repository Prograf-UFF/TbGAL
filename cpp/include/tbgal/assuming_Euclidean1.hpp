#ifndef __TBGAL_ASSUMING_EUCLIDEAN1_HPP__
#define __TBGAL_ASSUMING_EUCLIDEAN1_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Euclidean1 {
        
        using MetricSpaceType = EuclideanMetricSpace<1>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(SPACE)

        static auto const e1 = e(1);

    }

}

#endif // __TBGAL_ASSUMING_EUCLIDEAN1_HPP__
