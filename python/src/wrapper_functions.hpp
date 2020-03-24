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

using namespace tbgal;
namespace python = boost::python;

template<typename T>
inline
std::vector< T > py_list_to_std_vector( const boost::python::object& iterable )
{
    return std::vector< T >( boost::python::stl_input_iterator< T >( iterable ),
                             boost::python::stl_input_iterator< T >( ) );
}

template<typename T>
T py_hip(const T &lhs, const T &rhs) {
    return tbgal::hip(lhs, rhs);
}

template<typename T>
T py_dot(const T &lhs, const T &rhs) {
    return tbgal::dot(lhs, rhs);
}

template<typename T>
T py_lcont(const T &lhs, const T &rhs) {
    return tbgal::lcont(lhs, rhs);
}

template<typename T>
T py_rcont(const T &lhs, const T &rhs) {
    return tbgal::rcont(lhs, rhs);
}

template<typename T>
T py_reverse(const T &mv) {
    return tbgal::reverse(mv);
}

template<typename T>
T py_inverse(const T &mv) {
    return tbgal::inv(mv);
}

template<typename T>
double py_rnorm(const T &mv) {
    return tbgal::rnorm(mv);
}

template<typename T>
double py_rnorm_sqr(const T &mv) {
    return tbgal::rnorm_sqr(mv);
}

template<typename T>
T py_dual(const T &mv) {
    return tbgal::dual(mv);
}

template<typename T>
T py_undual(const T &mv) {
    return tbgal::undual(mv);
}

template<typename T>
double py_sp(const T &lhs, const T &rhs) {
    return tbgal::sp(lhs, rhs);
}

#endif