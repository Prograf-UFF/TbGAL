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

#include <tbgal/using_Eigen.hpp>
#include <tbgal/assuming_Euclidean3.hpp>

using namespace tbgal;
using namespace tbgal::Euclidean3;

int main(int argc, char *argv[]) {
    auto a = vector(0.5, 0.0, 0.5);    // a = 0.5 * (e1 + e3)

    auto M = 3.0 * (e1 ^ e2);          // M = 3.0 * e1^e2
    auto v = dual(M);                  // v = 3.0 * e3

    auto a_ = -gp(v, a, inverse(v));   // a_ = 0.5 * (e1 - e3)

    std::cout << "a = " << a << std::endl;
    std::cout << "M = " << M << std::endl;
    std::cout << "a_ = " << a_ << std::endl;

    return EXIT_SUCCESS;
}
