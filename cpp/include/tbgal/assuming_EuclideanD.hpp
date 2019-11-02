#ifndef __TBGAL_ASSUMING_EUCLIDEAND_HPP__
#define __TBGAL_ASSUMING_EUCLIDEAND_HPP__

#include "core.hpp"

#ifndef TBGAL_EuclideanD_MaxBaseSpaceDimensions
    #define TBGAL_EuclideanD_MaxBaseSpaceDimensions Dynamic
#endif // TBGAL_EuclideanD_MaxBaseSpaceDimensions

namespace tbgal {

    namespace EuclideanD {
        
        using MetricSpaceType = EuclideanMetricSpace<Dynamic, (TBGAL_EuclideanD_MaxBaseSpaceDimensions)>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(SPACE)

    }

}

#endif // __TBGAL_ASSUMING_EUCLIDEAND_HPP__
