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

#ifndef __TBGAL_ASSUMING_SIGNEDPQ_HPP__
#define __TBGAL_ASSUMING_SIGNEDPQ_HPP__

#include "core.hpp"

#ifndef TBGAL_SignedPQ_PDimensions
    #define TBGAL_SignedPQ_PDimensions Dynamic
#endif // TBGAL_SignedPQ_PDimensions

#ifndef TBGAL_SignedPQ_QDimensions
    #define TBGAL_SignedPQ_QDimensions Dynamic
#endif // TBGAL_SignedPQ_QDimensions

#ifndef TBGAL_SignedPQ_MaxDimensions
    #define TBGAL_SignedPQ_MaxDimensions ((TBGAL_SignedPQ_PDimensions) != Dynamic && (TBGAL_SignedPQ_QDimensions) != Dynamic ? (TBGAL_SignedPQ_PDimensions) + (TBGAL_SignedPQ_QDimensions) : Dynamic)
#endif // TBGAL_SignedPQ_MaxDimensions

namespace tbgal {

    namespace SignedPQ {
        
        using MetricSpaceType = SignedMetricSpace<(TBGAL_SignedPQ_PDimensions), (TBGAL_SignedPQ_QDimensions), (TBGAL_SignedPQ_MaxDimensions)>;

        static MetricSpaceType SPACE;
        
        TBGAL_OVERLOAD_UTILS(SPACE)
        TBGAL_OVERLOAD_SIGNED_UTILS(SPACE)

    }

}

#endif // __TBGAL_ASSUMING_SIGNEDPQ_HPP__
