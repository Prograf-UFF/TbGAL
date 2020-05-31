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

#ifndef __TBGAL_ASSUMING_CONFORMALD_HPP__
#define __TBGAL_ASSUMING_CONFORMALD_HPP__

#include "core.hpp"

#ifndef TBGAL_ConformalD_BaseSpaceDimensions
    #define TBGAL_ConformalD_BaseSpaceDimensions Dynamic
#endif // TBGAL_ConformalD_BaseSpaceDimensions

#ifndef TBGAL_ConformalD_MaxBaseSpaceDimensions
    #define TBGAL_ConformalD_MaxBaseSpaceDimensions (TBGAL_ConformalD_BaseSpaceDimensions)
#endif // TBGAL_ConformalD_MaxBaseSpaceDimensions

namespace tbgal {

    namespace ConformalD {
        
        using MetricSpaceType = ConformalMetricSpace<(TBGAL_ConformalD_BaseSpaceDimensions), (TBGAL_ConformalD_MaxBaseSpaceDimensions)>;

        static MetricSpaceType SPACE;
        
        TBGAL_OVERLOAD_UTILS(SPACE)
        TBGAL_OVERLOAD_CONFORMAL_UTILS(SPACE)

        inline decltype(auto) no() {
            return e(SPACE.dimensions() - 1);
        }

        inline decltype(auto) ni() {
            return e(SPACE.dimensions());
        }

    }

}

#endif // __TBGAL_ASSUMING_CONFORMALD_HPP__
