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

from tbgal.homogeneous3 import *

print '-- Input'
x = float(raw_input('x = '))
y = float(raw_input('y = '))
z = float(raw_input('z = '))
print

p = point(x, y, z);
d = direction(x, y, z);
l = p ^ d;

print '-- Result'
print 'p =', p
print 'd =', d
print 'l = p ^ d', l
