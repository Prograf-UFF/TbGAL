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

#ifndef __TBGAL_FACTORING_PRODUCT_HPP__
#define __TBGAL_FACTORING_PRODUCT_HPP__

namespace tbgal {

    template<typename MetricSpaceType_>
    struct GeometricProduct final {
    public:

        using MetricSpaceType = MetricSpaceType_;

    private:

        constexpr GeometricProduct() {
            // Nothing to do. It just avoid the instantiation of GeometricProduct<MetricSpaceType>.
        }
    };

    template<typename MetricSpaceType_>
    struct OuterProduct final {
    public:

        using MetricSpaceType = MetricSpaceType_;

    private:

        constexpr OuterProduct() {
            // Nothing to do. It just avoid the instantiation of OuterProduct<MetricSpaceType>.
        }
    };

}

#endif // __TBGAL_FACTORING_PRODUCT_HPP__
