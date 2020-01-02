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

#ifndef __TBGAL_ASSUMING_MINKOWSID_HPP__
#define __TBGAL_ASSUMING_MINKOWSID_HPP__

#include "core.hpp"

#ifndef TBGAL_MinkowsiD_BaseSpaceDimensions
    #define TBGAL_MinkowsiD_BaseSpaceDimensions Dynamic
#endif // TBGAL_MinkowsiD_BaseSpaceDimensions

#ifndef TBGAL_MinkowsiD_MaxBaseSpaceDimensions
    #define TBGAL_MinkowsiD_MaxBaseSpaceDimensions (TBGAL_MinkowsiD_BaseSpaceDimensions)
#endif // TBGAL_MinkowsiD_MaxBaseSpaceDimensions

namespace tbgal {

    namespace MinkowsiD {
        
        using MetricSpaceType = MinkowsiMetricSpace<(TBGAL_MinkowsiD_BaseSpaceDimensions), (TBGAL_MinkowsiD_MaxBaseSpaceDimensions)>;

        static MetricSpaceType SPACE;
        
        _TBGAL_OVERLOAD_UTILS(SPACE)
        _TBGAL_OVERLOAD_MINKOWSI_UTILS(SPACE)

        constexpr decltype(auto) ep() noexcept {
            return e(SPACE.dimensions() - 1);
        }

        constexpr decltype(auto) em() noexcept {
            return e(SPACE.dimensions());
        }

    }

}

#endif // __TBGAL_ASSUMING_MINKOWSID_HPP__
