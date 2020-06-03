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
#include "../../cpp/include/tbgal/assuming_Conformal3.hpp"

#include "common.hpp"
#include "macro_Conformal.hpp"

BOOST_PYTHON_MODULE(conformal3) {
    using namespace tbgal;
    using namespace tbgal::Conformal3;

    PY_TBGAL_INITIALIZE();
    
    PY_TBGAL_EXPOSE_CORE(MetricSpaceType);
    PY_TBGAL_EXPOSE_UTILS();

    PY_TBGAL_EXPOSE_CONFORMAL_METRIC_SPACE(MetricSpaceType);
    PY_TBGAL_EXPOSE_CONFORMAL_UTILS();
    
    PY_TBGAL_EXPOSE_GLOBAL_VARIABLE("space", SPACE);
    PY_TBGAL_EXPOSE_GLOBAL_VARIABLE("e1", e1);
    PY_TBGAL_EXPOSE_GLOBAL_VARIABLE("e2", e2);
    PY_TBGAL_EXPOSE_GLOBAL_VARIABLE("e3", e3);
    PY_TBGAL_EXPOSE_GLOBAL_VARIABLE("no", no);
    PY_TBGAL_EXPOSE_GLOBAL_VARIABLE("ni", ni);
}
