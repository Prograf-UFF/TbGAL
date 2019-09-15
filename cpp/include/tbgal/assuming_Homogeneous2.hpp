#ifndef __TBGAL_ASSUMING_HOMOGENEOUS2_HPP__
#define __TBGAL_ASSUMING_HOMOGENEOUS2_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Homogeneous2 {
        
        using MetricSpaceType = HomogeneousMetricSpace<2>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_HOMOGENEOUS_UTILS(SPACE)

        static auto const e1 = e(1);
        static auto const e2 = e(2);
        static auto const ep = e(3);

    }

}

#endif // __TBGAL_ASSUMING_HOMOGENEOUS2_HPP__
