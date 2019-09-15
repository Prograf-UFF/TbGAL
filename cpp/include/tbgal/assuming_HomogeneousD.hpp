#ifndef __TBGAL_ASSUMING_HOMOGENEOUSD_HPP__
#define __TBGAL_ASSUMING_HOMOGENEOUSD_HPP__

#include "core.hpp"

namespace tbgal {

    namespace HomogeneousD {
        
        using MetricSpaceType = HomogeneousMetricSpace<Dynamic>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_HOMOGENEOUS_UTILS(SPACE)

        constexpr decltype(auto) ep() noexcept {
            return e(SPACE.dimensions());
        }

    }

}

#endif // __TBGAL_ASSUMING_HOMOGENEOUSD_HPP__
