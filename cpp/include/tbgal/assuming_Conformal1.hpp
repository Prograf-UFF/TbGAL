#ifndef __TBGAL_ASSUMING_CONFORMAL1_HPP__
#define __TBGAL_ASSUMING_CONFORMAL1_HPP__

#include "core.hpp"

namespace tbgal {

    namespace Conformal1 {
        
        using MetricSpaceType = ConformalMetricSpace<1>;

        static MetricSpaceType const SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_CONFORMAL_UTILS(SPACE)

        static auto const e1 = e(1);
        static auto const no = e(2);
        static auto const ni = e(3);

    }

}

#endif // __TBGAL_ASSUMING_CONFORMAL1_HPP__
