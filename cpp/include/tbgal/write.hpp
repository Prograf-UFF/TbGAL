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

#ifndef __TBGAL_WRITE_HPP__
#define __TBGAL_WRITE_HPP__

namespace tbgal {

    namespace detail {

        template<typename ScalarType, typename FactoringProductType>
        std::ostream & write(std::ostream &os, FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            using IndexType = typename FactoredMultivector<ScalarType, FactoringProductType>::IndexType;
            os << arg.scalar();
            if (arg.factors_count() > 0) {
                auto factors = arg.factors_in_actual_metric();
                for (IndexType ind = 0; ind != arg.factors_count(); ++ind) {
                    os << ", (";
                    for (IndexType dim = 0; dim != arg.space().dimensions(); ++dim) {
                        auto value = coeff(factors, dim, ind);
                        if (dim > 0) {
                            os << (value >= 0 ? " + " : " - ") << std::abs(value);
                        }
                        else {
                            os << value;
                        }
                        os << " * " << arg.space().basis_vector_str(dim);
                    }
                    os << ")";
                }
            }
            return os;
        }

    }

    template<typename ScalarType, typename MetricSpaceType>
    std::ostream & operator <<(std::ostream &os, FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType>> const &arg) noexcept {
        os << "gp(";
        detail::write(os, arg);
        os << ")";
        return os;
    }

    template<typename ScalarType, typename MetricSpaceType>
    std::ostream & operator <<(std::ostream &os, FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType>> const &arg) noexcept {
        os << "op(";
        detail::write(os, arg);
        os << ")";
        return os;
    }

}

#endif // __TBGAL_WRITE_HPP__
