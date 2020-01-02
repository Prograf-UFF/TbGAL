/* Copyright (C) Eduardo Vera Sousa and Leandro Augusto Frata Fernandes
 * 
 * authors    : Sousa, Eduardo V.
 *              Fernandes, Leandro A. F.
 * repository : https://github.com/Prograf-UFF/TbGAL
 * 
 * This file is part of the Tensor-based Geometric Algebra Library (TbGAL).
 * 
 * TbGAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * TbGAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * 
 * You should have received a copy of the GNU General Public License
 * along with TbGAL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef __TBGAL_ASSUMING_EUCLIDEAND_HPP__
#define __TBGAL_ASSUMING_EUCLIDEAND_HPP__

#include "core.hpp"

#ifndef TBGAL_EuclideanD_BaseSpaceDimensions
    #define TBGAL_EuclideanD_BaseSpaceDimensions Dynamic
#endif // TBGAL_EuclideanD_BaseSpaceDimensions

#ifndef TBGAL_EuclideanD_MaxBaseSpaceDimensions
    #define TBGAL_EuclideanD_MaxBaseSpaceDimensions (TBGAL_EuclideanD_BaseSpaceDimensions)
#endif // TBGAL_EuclideanD_MaxBaseSpaceDimensions

namespace tbgal {

    namespace EuclideanD {
        
        using MetricSpaceType = EuclideanMetricSpace<(TBGAL_EuclideanD_BaseSpaceDimensions), (TBGAL_EuclideanD_MaxBaseSpaceDimensions)>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_EUCLIDEAN_UTILS(SPACE)

    }

}

#endif // __TBGAL_ASSUMING_EUCLIDEAND_HPP__
