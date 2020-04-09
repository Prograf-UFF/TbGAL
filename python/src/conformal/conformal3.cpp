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

#include "../../../cpp/include/tbgal/using_Eigen.hpp"
#include "../../../cpp/include/tbgal/assuming_Conformal3.hpp"
#include "../wrapper_functions.hpp"
#include "../macros.hpp"

namespace py_tbgal {

    namespace python = boost::python;

    auto py_vector(python::tuple args, python::dict kwargs) {
        std::vector<std::double_t> args_as_container = py_list_to_std_vector<std::double_t>(args);
        auto a = tbgal::Conformal3::vector(args_as_container.begin(), args_as_container.end());
        return a;
    }

    auto py_euclidean_vector(python::tuple args, python::dict kwargs) {
        std::vector<std::double_t> args_as_container = py_list_to_std_vector<std::double_t>(args);
        auto a = tbgal::Conformal3::euclidean_vector(args_as_container.begin(), args_as_container.end());
        return a;
    }

    auto py_point(python::tuple args, python::dict kwargs) {
        std::vector<std::double_t> args_as_container = py_list_to_std_vector<std::double_t>(args);
        auto a = tbgal::Conformal3::point(args_as_container.begin(), args_as_container.end());
        return a;
    }

    auto no() {
        return tbgal::Conformal3::no;
    }

    auto ni() {
        return tbgal::Conformal3::ni;
    }


    BOOST_PYTHON_MODULE(conformal3) {
        using namespace tbgal::Conformal3;

        using OP = tbgal::OuterProduct<tbgal::ConformalMetricSpace<3>>;
        using GP = tbgal::GeometricProduct<tbgal::ConformalMetricSpace<3>>;

        _DECLARE_FACTORED_MULTIVECTOR_PYTHON(std::double_t, OP, std::double_t, GP);
        _DECLARE_ALL_OPERATIONS(std::double_t, OP, std::double_t, GP);

        python::def("vector", python::raw_function(py_vector) );
        python::def("euclidean_vector", python::raw_function(py_euclidean_vector));
        python::def("point", python::raw_function(py_point));

        python::def("e", &tbgal::Conformal3::e);
        python::def("no", &no);
        python::def("ni", &ni);


    }
}