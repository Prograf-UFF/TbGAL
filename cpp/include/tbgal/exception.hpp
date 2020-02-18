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

#ifndef __TBGAL_EXCEPTION_HPP__
#define __TBGAL_EXCEPTION_HPP__

namespace tbgal {

    class NotSupportedError : public std::logic_error {
    public:

        explicit NotSupportedError(const std::string &what_arg) :
            std::logic_error(what_arg) {
        }

        explicit NotSupportedError(const char *what_arg) :
            std::logic_error(what_arg) {
        }
    };

}

#endif // __TBGAL_EXCEPTION_HPP__
