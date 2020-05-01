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

#ifndef __TBGAL_PYTHON_WRAPPER_FUNCTIONS_HPP__
#define __TBGAL_PYTHON_WRAPPER_FUNCTIONS_HPP__

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

    using namespace tbgal;
    namespace python = boost::python;

    template<typename T>
    void translate_exception(T const& e)
    {
        PyErr_SetString(PyExc_RuntimeError, e.what());
    }

    template<typename T>
    inline
    std::vector< T > py_list_to_std_vector( const boost::python::object& iterable )
    {
        return std::vector< T >( boost::python::stl_input_iterator< T >( iterable ),
                                boost::python::stl_input_iterator< T >( ) );
    }

}

#endif