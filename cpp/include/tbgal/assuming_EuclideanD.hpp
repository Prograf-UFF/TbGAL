#ifndef __TBGAL_ASSUMING_EUCLIDEAND_HPP__
#define __TBGAL_ASSUMING_EUCLIDEAND_HPP__

#include "core.hpp"

namespace tbgal {

    namespace EuclideanD {
        
        using MetricSpaceType = EuclideanMetricSpace<DefaultScalarType, Dynamic>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(SPACE)

    }

}

#endif // __TBGAL_ASSUMING_EUCLIDEAND_HPP__
