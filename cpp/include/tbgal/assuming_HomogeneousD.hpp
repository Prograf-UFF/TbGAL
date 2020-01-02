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

#ifndef __TBGAL_ASSUMING_HOMOGENEOUSD_HPP__
#define __TBGAL_ASSUMING_HOMOGENEOUSD_HPP__

#include "core.hpp"

#ifndef TBGAL_HomogeneousD_BaseSpaceDimensions
    #define TBGAL_HomogeneousD_BaseSpaceDimensions Dynamic
#endif // TBGAL_HomogeneousD_BaseSpaceDimensions

#ifndef TBGAL_HomogeneousD_MaxBaseSpaceDimensions
    #define TBGAL_HomogeneousD_MaxBaseSpaceDimensions (TBGAL_HomogeneousD_BaseSpaceDimensions)
#endif // TBGAL_HomogeneousD_MaxBaseSpaceDimensions

namespace tbgal {

    namespace HomogeneousD {
        
        using MetricSpaceType = HomogeneousMetricSpace<(TBGAL_HomogeneousD_BaseSpaceDimensions), (TBGAL_HomogeneousD_MaxBaseSpaceDimensions)>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_HOMOGENEOUS_UTILS(SPACE)

        constexpr decltype(auto) ep() noexcept {
            return e(SPACE.dimensions());
        }

    }

}

#endif // __TBGAL_ASSUMING_HOMOGENEOUSD_HPP__
