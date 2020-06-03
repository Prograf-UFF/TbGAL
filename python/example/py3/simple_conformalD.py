# Copyright (C) Eduardo Vera Sousa and Leandro Augusto Frata Fernandes
# 
# authors    : Sousa, Eduardo V.
#              Fernandes, Leandro A. F.
# repository : https://github.com/Prograf-UFF/TbGAL
# 
# This file is part of the Tensor-based Geometric Algebra Library (TbGAL).
# 
# TbGAL is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# TbGAL is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with TbGAL. If not, see <https://www.gnu.org/licenses/>.

from tbgal.conformalD import *

space.set_base_space_dimensions(10)
print('Conformal model of {}-D Euclidean space (n = {})'.format(space.base_space_dimensions(), space.dimensions()))
print()

x = euclidean_vector(10, 20, 30, 40 ,50, 60, 70, 80, 90, 100)
x1 = sp(x, e(1))
x3 = sp(x, e(3))    # We have to use the e(i), no(), and ni()
                    # functions instead of the ei, no, and ni
Y = x^ni()          # constants because the number of dimensions
                    # of the base space is defined at runtime
print('x =', x)
print('x1 =', x1)
print('x3 =', x3)
print()
print('Y =', Y)
