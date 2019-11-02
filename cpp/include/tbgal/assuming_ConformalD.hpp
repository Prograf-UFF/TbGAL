#ifndef __TBGAL_ASSUMING_CONFORMALD_HPP__
#define __TBGAL_ASSUMING_CONFORMALD_HPP__

#include "core.hpp"

#ifndef TBGAL_ConformalD_MaxBaseSpaceDimensions
    #define TBGAL_ConformalD_MaxBaseSpaceDimensions Dynamic
#endif // TBGAL_ConformalD_MaxBaseSpaceDimensions

namespace tbgal {

    namespace ConformalD {
        
        using MetricSpaceType = ConformalMetricSpace<Dynamic, (TBGAL_ConformalD_MaxBaseSpaceDimensions)>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_CONFORMAL_UTILS(SPACE)

        constexpr decltype(auto) no() noexcept {
            return e(SPACE.dimensions() - 1);
        }

        constexpr decltype(auto) ni() noexcept {
            return e(SPACE.dimensions());
        }

    }

}

#endif // __TBGAL_ASSUMING_CONFORMALD_HPP__
