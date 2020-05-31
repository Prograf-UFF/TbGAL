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

#ifndef __TBGAL_PYTHON_MACRO_EUCLIDEAN_HPP__
#define __TBGAL_PYTHON_MACRO_EUCLIDEAN_HPP__

#define PY_TBGAL_EXPOSE_EUCLIDEAN_METRIC_SPACE(METRIC_SPACE_TYPE) \
    py::class_<METRIC_SPACE_TYPE>("MetricSpace") \
        .def("dimensions", &METRIC_SPACE_TYPE::dimensions) \
        .def("set_dimensions", &METRIC_SPACE_TYPE::set_dimensions)

#define PY_TBGAL_EXPOSE_EUCLIDEAN_UTILS() \
    py::def("euclidean_vector", py::raw_function(+[](py::tuple const &args, py::dict const &) { return euclidean_vector(py::stl_input_iterator<tbgal::DefaultScalarType>(args), py::stl_input_iterator<tbgal::DefaultScalarType>()); }))

#endif // __TBGAL_PYTHON_MACRO_EUCLIDEAN_HPP__
