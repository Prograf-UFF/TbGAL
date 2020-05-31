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

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

#include "../../cpp/include/tbgal/using_Eigen.hpp"
#include "../../cpp/include/tbgal/assuming_ConformalD.hpp"

#include "macro.hpp"
#include "macro_Conformal.hpp"

namespace py = boost::python;
namespace np = boost::python::numpy;

BOOST_PYTHON_MODULE(conformalD) {
    using namespace tbgal;
    using namespace tbgal::ConformalD;

    PY_TBGAL_INITIALIZE();
    
    PY_TBGAL_EXPOSE_CORE(MetricSpaceType);
    PY_TBGAL_EXPOSE_UTILS();

    PY_TBGAL_EXPOSE_CONFORMAL_METRIC_SPACE(MetricSpaceType);
    PY_TBGAL_EXPOSE_CONFORMAL_UTILS();
    
    PY_TBGAL_EXPOSE_VARIABLE("space", SPACE);
    PY_TBGAL_EXPOSE_FUNCTION("no", no);
    PY_TBGAL_EXPOSE_FUNCTION("ni", ni);
}
