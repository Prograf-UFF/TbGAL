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

#ifndef __TBGAL_PYTHON_MACRO_MINKOWSKI_HPP__
#define __TBGAL_PYTHON_MACRO_MINKOWSKI_HPP__

#define PY_TBGAL_EXPOSE_MINKOWSKI_METRIC_SPACE(METRIC_SPACE_TYPE) \
    py::class_<METRIC_SPACE_TYPE>("MetricSpace") \
        .def("dimensions", &METRIC_SPACE_TYPE::dimensions) \
        .def("base_space_dimensions", &METRIC_SPACE_TYPE::base_space_dimensions) \
        .def("set_base_space_dimensions", &METRIC_SPACE_TYPE::set_base_space_dimensions)

#define PY_TBGAL_EXPOSE_MINKOWSKI_UTILS() \
    py::def("euclidean_vector", py::raw_function(+[](py::tuple const &args, py::dict const &) { \
        std::vector<DefaultScalarType> dummy; \
        dummy.assign(py::stl_input_iterator<DefaultScalarType>(args), py::stl_input_iterator<DefaultScalarType>()); \
        return vector(dummy.begin(), dummy.end()); \
    })); \
    py::def("point", py::raw_function(+[](py::tuple const &args, py::dict const &) { \
        std::vector<DefaultScalarType> dummy; \
        dummy.assign(py::stl_input_iterator<DefaultScalarType>(args), py::stl_input_iterator<DefaultScalarType>()); \
        return vector(dummy.begin(), dummy.end()); \
    }));

#endif // __TBGAL_PYTHON_MACRO_MINKOWSKI_HPP__
