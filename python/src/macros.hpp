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

#define _BUILD_FACTORED_MULTIVECTOR_PYTHON(PRODUCT, OTHER_PRODUCT, METRIC, SCALAR_MAIN, SCALAR_EXTRA) \
    namespace python = boost::python; \
    python::class_< tbgal::FactoredMultivector<SCALAR_MAIN, PRODUCT<METRIC<Dynamic, Dynamic>>> >("FactoredMultivector")  \ 
        .def("__repr__", &tbgal::FactoredMultivector<SCALAR_MAIN, PRODUCT<METRIC<Dynamic, Dynamic>> >::repr)  \ 
        .def(python::self ^ python::other<tbgal::FactoredMultivector<SCALAR_MAIN, OTHER_PRODUCT<METRIC<Dynamic, Dynamic>> >>())  \ 
        .def(python::self ^ python::self)  \ 
        .def(python::self ^ python::other<SCALAR_EXTRA>())  \ 
        .def(python::other<SCALAR_EXTRA>() ^ python::self)  \ 
        .def(python::self ^ python::other<SCALAR_MAIN>())  \ 
        .def(python::other<SCALAR_MAIN>() ^ python::self)  \ 
        .def(python::self * python::other<tbgal::FactoredMultivector<SCALAR_MAIN, OTHER_PRODUCT<METRIC<Dynamic, Dynamic>> >>())  \ 
        .def(python::self * python::self)  \ 
        .def(python::self * python::other<SCALAR_EXTRA>())  \ 
        .def(python::other<SCALAR_EXTRA>() * python::self)  \ 
        .def(python::self * python::other<SCALAR_MAIN>())  \ 
        .def(python::other<SCALAR_MAIN>() * python::self)  \ 
		.def(+python::self)  \ 
		.def(-python::self)  \ 
		.def(python::self + python::self)  \ 
        .def(python::self + python::other<tbgal::FactoredMultivector<SCALAR_MAIN, OTHER_PRODUCT<METRIC<Dynamic, Dynamic>> > >())  \ 
        .def(~python::self);  \ 

#endif