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
#include <tbgal/assuming_ConformalD.hpp>

using namespace tbgal;
using namespace tbgal::ConformalD;

int main(int argc, char *argv[]) {
    SPACE.set_base_space_dimensions(10);
    std::cout << "Conformal model of " << SPACE.base_space_dimensions() << "-D Euclidean space (n =" << SPACE.dimensions() << ")" << std::endl;
    std::cout << std::endl;

    auto x = euclidean_vector(10.0, 20.0, 30.0, 40.0, 50.0, 60.0, 70.0, 80.0, 90.0, 100.0);
    auto x1 = sp(x, e(1));
    auto x3 = sp(x, e(3));                     // We have to use the e(i), no(), and ni()
                                               // functions instead of the ei, no, and ni
    auto Y = x^ni();                           // constants because the number of dimensions
                                               // of the base space is defined at runtime
    std::cout << "x = " << x << std::endl;
    std::cout << "x1 = " << x1 << std::endl;
    std::cout << "x3 = " << x3 << std::endl;
    std::cout << std::endl;
    std::cout << "Y = " << Y << std::endl;

    return EXIT_SUCCESS;
}
