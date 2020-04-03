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

#include "../../cpp/include/tbgal/using_Eigen.hpp"
#include "../../cpp/include/tbgal/assuming_EuclideanD.hpp"
#include "wrapper_functions.hpp"
#include "macros.hpp"

namespace py_tbgal {

    namespace python = boost::python;

    auto py_vector(python::tuple args, python::dict kwargs) {
        tbgal::EuclideanD::SPACE = EuclideanMetricSpace<Dynamic, Dynamic>(len(args));
        std::vector<double> args_as_container = py_list_to_std_vector<double>(args);
        auto a = tbgal::EuclideanD::vector(args_as_container.begin(), args_as_container.end());
        return a;
    }

    BOOST_PYTHON_MODULE(euclidean) {

        using OP = tbgal::OuterProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>;
        using GP = tbgal::GeometricProduct<tbgal::EuclideanMetricSpace<Dynamic, Dynamic>>;

        _DECLARE_FACTORED_MULTIVECTOR_PYTHON(std::double_t, OP, std::double_t, GP);

        python::def("vector", python::raw_function(py_vector) );

        python::def("hip", +[](
            const tbgal::FactoredMultivector<double, OP> &lhs, 
            const tbgal::FactoredMultivector<double, OP> &rhs) {
                return hip(lhs, rhs); 
                }
            );

        python::def("dot", +[](
            const tbgal::FactoredMultivector<double, OP> &lhs, 
            const tbgal::FactoredMultivector<double, OP> &rhs) {
                return dot(lhs, rhs); 
                }
            );
        python::def("lcont", +[](
            const tbgal::FactoredMultivector<double, OP> &lhs, 
            const tbgal::FactoredMultivector<double, OP> &rhs) {
                return lcont(lhs, rhs); 
                }
            );

        python::def("rcont", +[](
            const tbgal::FactoredMultivector<double, OP> &lhs, 
            const tbgal::FactoredMultivector<double, OP> &rhs) {
                return rcont(lhs, rhs); 
                }
            );

        python::def("reverse", +[](
            const tbgal::FactoredMultivector<double, OP> &mv) {
                return reverse(mv);
                }
            );

        python::def("reverse", +[](
            const tbgal::FactoredMultivector<double, GP> &mv) {
                return reverse(mv);
                }
            );

        python::def("inverse", +[](
            const tbgal::FactoredMultivector<double, OP> &mv) {
                return inverse(mv);
                }
            );

        python::def("reverse", +[](
            const tbgal::FactoredMultivector<double, GP> &mv) {
                return inverse(mv);
                }
            );

        python::def("rnorm", +[](
            const tbgal::FactoredMultivector<double, OP> &mv) {
                return rnorm(mv);
                }
            );

        python::def("rnorm", +[](
            const tbgal::FactoredMultivector<double, GP> &mv) {
                return rnorm(mv);
                }
            );

        python::def("rnorm_sqr", +[](
            const tbgal::FactoredMultivector<double, OP> &mv) {
                return rnorm_sqr(mv);
                }
            );

        python::def("rnorm_sqr", +[](
            const tbgal::FactoredMultivector<double, GP> &mv) {
                return rnorm_sqr(mv);
                }
            );

        python::def("dual", +[](
            const tbgal::FactoredMultivector<double, OP> &mv) {
                return dual(mv);
                }
            );

        python::def("undual", +[](
            const tbgal::FactoredMultivector<double, OP> &mv) {
                return undual(mv);
                }
            );

        python::def("sp", +[](
            const tbgal::FactoredMultivector<double, OP> &lhs, 
            const tbgal::FactoredMultivector<double, OP> &rhs) {
                return sp(lhs, rhs); 
                }
            );

    }

}