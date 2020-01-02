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

#include "../../../include/tbgal/using_Eigen.hpp" //TODO [DEBUG]
#include "../../../include/tbgal/assuming_Homogeneous3.hpp" //TODO [DEBUG]

using namespace tbgal;

int main(int argc, char *argv[]) {
    double x, y, z;

    auto func = [](auto a, auto b) {
        return a + b;
    };

    std::cout << "-- Input" << std::endl;
    std::cout << std::endl;
    std::cout << "x = "; std::cin >> x;
    std::cout << "y = "; std::cin >> y;
    std::cout << "z = "; std::cin >> z;
    std::cout << std::endl;

    auto p = Homogeneous3::point(x, y, z);
    auto d = Homogeneous3::direction(x, y, z);
    auto l = p ^ d;

    std::cout << "-- Result" << std::endl;
    std::cout << std::endl;
    std::cout << "p = " << p << std::endl;
    std::cout << "d = " << d << std::endl;
    std::cout << "l = p ^ d = " << l << std::endl;
    std::cout << std::endl;

    return EXIT_SUCCESS;
}
