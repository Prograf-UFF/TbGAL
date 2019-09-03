#ifndef __TBGAL_ASSUMING_HOMOGENEOUS3_HPP__
#define __TBGAL_ASSUMING_HOMOGENEOUS3_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Homogeneous3 {
        
        using MetricSpaceType = HomogeneousMetricSpace<DefaultScalarType, 3>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_HOMOGENEOUS_UTILS(SPACE)

        static auto const e1 = e(1);
        static auto const e2 = e(2);
        static auto const e3 = e(3);
        static auto const ep = e(4);

    }

}

#endif // __TBGAL_ASSUMING_HOMOGENEOUS3_HPP__
