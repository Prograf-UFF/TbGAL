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

#ifndef __TBGAL_PYTHON_MACROS_HPP__
#define __TBGAL_PYTHON_MACROS_HPP__


#include <boost/python/module.hpp>
#include <boost/python/def.hpp>  
#include <boost/python.hpp>  
#include <boost/python/object.hpp>  
#include <boost/python/args.hpp>  
#include <boost/python/list.hpp>  
#include <boost/python/tuple.hpp> 
#include <boost/python/dict.hpp>  
#include "boost/python/stl_iterator.hpp" 
#include <boost/python/raw_function.hpp> 

namespace py_tbgal {

    #define _DECLARE_FACTORED_MULTIVECTOR_SCALAR_OPERATIONS(SCALAR)\
            .def(python::self ^ python::other<SCALAR>())\
            .def(python::other<SCALAR>() ^ python::self)\
            .def(python::self * python::other<SCALAR>())\
            .def(python::other<SCALAR>() * python::self)\

    #define _DECLARE_FACTORED_MULTIVECTOR_SELF_OPERATIONS(SCALAR_TYPE, FACTORING_PRODUCT_TYPE)\
            .def(python::self ^ python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >())\
            .def(python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>() ^ python::self)\
            .def(python::self * python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>())\
            .def(python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>() * python::self)\
            .def(python::self + python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>())\
            .def(python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>() + python::self)\
            .def(python::self - python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>())\
            .def(python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>() - python::self)\
            .def(python::self / python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>())\
            .def(python::other<tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>>() / python::self)\

    #define _DECLARE_FACTORED_MULTIVECTOR_CLASS(SCALAR_TYPE, FACTORING_PRODUCT_TYPE)\
        python::class_< tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >("FactoredMultivector")\
            .def("__repr__", &tbgal::FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>::repr)\
            _DECLARE_FACTORED_MULTIVECTOR_SCALAR_OPERATIONS(std::double_t)\
            _DECLARE_FACTORED_MULTIVECTOR_SCALAR_OPERATIONS(std::int16_t)\
            _DECLARE_FACTORED_MULTIVECTOR_SCALAR_OPERATIONS(std::int32_t)\
            _DECLARE_FACTORED_MULTIVECTOR_SCALAR_OPERATIONS(std::int64_t)\
            _DECLARE_FACTORED_MULTIVECTOR_SCALAR_OPERATIONS(std::float_t)\
            .def(python::self ^ python::self)\
            .def(python::self * python::self)\
            .def(python::self / python::self)\
            .def(+python::self)\
            .def(-python::self)\
            .def(python::self + python::self)\
            .def(python::self - python::self)\
            .def(~python::self)\



    #define _DECLARE_FACTORED_MULTIVECTOR_PYTHON(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE)\
        namespace python = boost::python;\
        _DECLARE_FACTORED_MULTIVECTOR_CLASS(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE)\
        _DECLARE_FACTORED_MULTIVECTOR_SELF_OPERATIONS(SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE)\
        _DECLARE_FACTORED_MULTIVECTOR_SELF_OPERATIONS(std::double_t, SECOND_FACTORING_PRODUCT_TYPE)\
        ;\
        _DECLARE_FACTORED_MULTIVECTOR_CLASS(SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE)\
        _DECLARE_FACTORED_MULTIVECTOR_SELF_OPERATIONS(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE)\
        _DECLARE_FACTORED_MULTIVECTOR_SELF_OPERATIONS(std::double_t, FIRST_FACTORING_PRODUCT_TYPE)\
        ;\

    #define _DECLARE_ALL_OPERATIONS(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE)\
            namespace python = boost::python;\
            python::def("gp", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return gp(lhs, rhs);  \
                    } \
                ); \
            python::def("gp", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return gp(lhs, rhs);  \
                    } \
                ); \
            python::def("gp", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return gp(lhs, rhs);  \
                    } \
                ); \
            python::def("gp", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return gp(lhs, rhs);  \
                    } \
                ); \
            python::def("igp", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return igp(lhs, rhs);  \
                    } \
                ); \
            python::def("igp", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return igp(lhs, rhs);  \
                    } \
                ); \
            python::def("igp", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return igp(lhs, rhs);  \
                    } \
                ); \
            python::def("igp", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return igp(lhs, rhs);  \
                    } \
                ); \
            python::def("op", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return op(lhs, rhs);  \
                    } \
                ); \
            python::def("op", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return op(lhs, rhs);  \
                    } \
                ); \
            python::def("op", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return op(lhs, rhs);  \
                    } \
                ); \
            python::def("op", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return op(lhs, rhs);  \
                    } \
                ); \
            python::def("hip", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return hip(lhs, rhs);  \
                    } \
                ); \
            python::def("dot", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return dot(lhs, rhs);  \
                    } \
                ); \
            python::def("lcont", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return lcont(lhs, rhs);  \
                    } \
                ); \
            python::def("rcont", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return rcont(lhs, rhs);  \
                    } \
                ); \
            python::def("reverse", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &mv) { \
                    return reverse(mv); \
                    } \
                ); \
            python::def("reverse", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &mv) { \
                    return reverse(mv); \
                    } \
                ); \
            python::def("inverse", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &mv) { \
                    return inverse(mv); \
                    } \
                ); \
            python::def("reverse", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &mv) { \
                    return inverse(mv); \
                    } \
                ); \
            python::def("rnorm", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &mv) { \
                    return rnorm(mv); \
                    } \
                ); \
            python::def("rnorm", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &mv) { \
                    return rnorm(mv); \
                    } \
                ); \
            python::def("rnorm_sqr", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &mv) { \
                    return rnorm_sqr(mv); \
                    } \
                ); \
            python::def("rnorm_sqr", +[]( \
                const tbgal::FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> &mv) { \
                    return rnorm_sqr(mv); \
                    } \
                ); \
            python::def("dual", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &mv) { \
                    return dual(mv); \
                    } \
                ); \
            python::def("undual", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &mv) { \
                    return undual(mv); \
                    } \
                ); \
            python::def("sp", +[]( \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &lhs,  \
                const tbgal::FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> &rhs) { \
                    return sp(lhs, rhs);  \
                    } \
                );



}

    #endif