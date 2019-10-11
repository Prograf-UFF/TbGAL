#ifndef __TBGAL_ASSUMING_CONFORMAL2_HPP__
#define __TBGAL_ASSUMING_CONFORMAL2_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Conformal2 {
        
        using MetricSpaceType = ConformalMetricSpace<2>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_CONFORMAL_UTILS(SPACE)

        static auto const e1 = e(1);
        static auto const e2 = e(2);
        static auto const no = e(3);
        static auto const ni = e(4);

    }

}

#endif // __TBGAL_ASSUMING_CONFORMAL2_HPP__
