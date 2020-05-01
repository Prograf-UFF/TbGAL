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
#include "../../../cpp/include/tbgal/assuming_Conformal2.hpp"
#include "../wrapper_functions.hpp"
#include "../macros.hpp"
#include "../../../cpp/include/tbgal/exception.hpp"
#include <boost/python/exception_translator.hpp>

namespace py_tbgal {

    namespace python = boost::python;

    auto py_vector(python::tuple args, python::dict kwargs) {
        std::vector<std::double_t> args_as_container = py_list_to_std_vector<std::double_t>(args);
        auto a = tbgal::Conformal2::vector(args_as_container.begin(), args_as_container.end());
        return a;
    }

    auto py_euclidean_vector(python::tuple args, python::dict kwargs) {
        std::vector<std::double_t> args_as_container = py_list_to_std_vector<std::double_t>(args);
        auto a = tbgal::Conformal2::euclidean_vector(args_as_container.begin(), args_as_container.end());
        return a;
    }

    auto py_point(python::tuple args, python::dict kwargs) {
        std::vector<std::double_t> args_as_container = py_list_to_std_vector<std::double_t>(args);
        auto a = tbgal::Conformal2::point(args_as_container.begin(), args_as_container.end());
        return a;
    }

    auto no() {
        return tbgal::Conformal2::no;
    }

    auto ni() {
        return tbgal::Conformal2::ni;
    }


    BOOST_PYTHON_MODULE(conformal2) {
        using namespace tbgal::Conformal2;

        using OP = tbgal::OuterProduct<tbgal::ConformalMetricSpace<2>>;
        using GP = tbgal::GeometricProduct<tbgal::ConformalMetricSpace<2>>;

        _DECLARE_FACTORED_MULTIVECTOR_PYTHON(std::double_t, OP, std::double_t, GP);
        _DECLARE_ALL_OPERATIONS(std::double_t, OP, std::double_t, GP);

        python::def("vector", python::raw_function(py_vector) );
        python::def("euclidean_vector", python::raw_function(py_euclidean_vector));
        python::def("point", python::raw_function(py_point));

        python::def("e", &tbgal::Conformal2::e);
        python::def("no", &no);
        python::def("ni", &ni);


    }
}