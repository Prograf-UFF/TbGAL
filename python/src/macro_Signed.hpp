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

#ifndef __TBGAL_PYTHON_MACRO_SIGNED_HPP__
#define __TBGAL_PYTHON_MACRO_SIGNED_HPP__

#define PY_TBGAL_EXPOSE_SIGNED_METRIC_SPACE(METRIC_SPACE_TYPE) \
    py::class_<METRIC_SPACE_TYPE>("MetricSpace") \
        .def("dimensions", &METRIC_SPACE_TYPE::dimensions) \
        .def("p_dimensions", &METRIC_SPACE_TYPE::p_dimensions) \
        .def("q_dimensions", &METRIC_SPACE_TYPE::q_dimensions) \
        .def("set_dimensions", &METRIC_SPACE_TYPE::set_dimensions, py::args("p_dimensions", "q_dimensions"))

#define PY_TBGAL_EXPOSE_SIGNED_UTILS()

#endif // __TBGAL_PYTHON_MACRO_SIGNED_HPP__
