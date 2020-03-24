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

#include "../cpp/include/tbgal/using_Eigen.hpp"
#include "../cpp/include/tbgal/assuming_HomogeneousD.hpp"
#include "wrapper_functions.hpp"
#include "macros.hpp"

namespace python = boost::python;

auto py_vector(python::tuple args, python::dict kwargs) {
    tbgal::HomogeneousD::SPACE = HomogeneousMetricSpace<Dynamic, Dynamic>(len(args));
    std::vector<double> args_as_container = py_list_to_std_vector<double>(args);
    auto a = tbgal::HomogeneousD::vector(args_as_container.begin(), args_as_container.end());
    return a;
}

BOOST_PYTHON_MODULE(homogeneous) {

    _BUILD_FACTORED_MULTIVECTOR_PYTHON(tbgal::OuterProduct, tbgal::GeometricProduct, tbgal::HomogeneousMetricSpace, double, int);
    _BUILD_FACTORED_MULTIVECTOR_PYTHON(tbgal::GeometricProduct, tbgal::OuterProduct, tbgal::HomogeneousMetricSpace, double, int);

    python::def("vector", python::raw_function(py_vector) );

    python::def("hip", &py_hip< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("dot", &py_dot< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("lcont", &py_lcont< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("rcont", &py_rcont< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("reverse", &py_reverse< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);
    python::def("reverse", &py_reverse< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("inverse", &py_inverse< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);
    python::def("inverse", &py_inverse< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("rnorm", &py_rnorm< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);
    python::def("rnorm", &py_rnorm< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("rnorm_sqr", &py_rnorm_sqr< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);
    python::def("rnorm_sqr", &py_rnorm_sqr< tbgal::FactoredMultivector<double, tbgal::GeometricProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("dual", &py_dual< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("undual", &py_undual< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("sp", &py_sp< tbgal::FactoredMultivector<double, tbgal::OuterProduct<tbgal::HomogeneousMetricSpace<Dynamic, Dynamic>>> >);

    python::def("ep", &tbgal::HomogeneousD::ep);

}