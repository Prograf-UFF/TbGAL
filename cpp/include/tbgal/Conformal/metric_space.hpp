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

#ifndef __TBGAL_CONFORMAL_METRIC_SPACE_HPP__
#define __TBGAL_CONFORMAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime_, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime_ = BaseSpaceDimensionsAtCompileTime_>
    class ConformalMetricSpace : public BaseSignedMetricSpace<BaseSpaceDimensionsAtCompileTime_ != Dynamic ? BaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic, 1, MaxBaseSpaceDimensionsAtCompileTime_ != Dynamic ? MaxBaseSpaceDimensionsAtCompileTime_ + 2 : Dynamic> {
    private:

        using Super = BaseSignedMetricSpace<BaseSpaceDimensionsAtCompileTime_ != Dynamic ? BaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic, 1, MaxBaseSpaceDimensionsAtCompileTime_ != Dynamic ? MaxBaseSpaceDimensionsAtCompileTime_ + 2 : Dynamic>;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = typename Super::ScalarType;

        constexpr static IndexType DimensionsAtCompileTime = Super::DimensionsAtCompileTime;
        constexpr static IndexType MaxDimensionsAtCompileTime = Super::MaxDimensionsAtCompileTime;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime_;
        constexpr static IndexType MaxBaseSpaceDimensionsAtCompileTime = MaxBaseSpaceDimensionsAtCompileTime_;

        inline ConformalMetricSpace(ConformalMetricSpace const &) = default;
        inline ConformalMetricSpace(ConformalMetricSpace &&) = default;

        inline ConformalMetricSpace(IndexType base_space_dimensions) :
            Super(base_space_dimensions + 1, 1),
            basis_vectors_str_(),
            from_actual_to_signed_metric_(),
            from_signed_to_actual_metric_() {
            update_basis_vectors_str(base_space_dimensions);
        }

        inline ConformalMetricSpace() :
            ConformalMetricSpace((BaseSpaceDimensionsAtCompileTime != Dynamic) ? BaseSpaceDimensionsAtCompileTime : 0) {
        }

        inline ConformalMetricSpace & operator=(ConformalMetricSpace const &) = default;
        inline ConformalMetricSpace & operator=(ConformalMetricSpace &&) = default;

        inline std::string const & basis_vector_str(IndexType index) const override {
            return basis_vectors_str_[index];
        }
        
        inline IndexType base_space_dimensions() const {
            return Super::dimensions() - 2;
        }

        inline void set_base_space_dimensions(IndexType base_space_dimensions) {
            Super::set_dimensions(base_space_dimensions + 1, 1);
            update_basis_vectors_str(base_space_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType base_space_dimensions) {
            basis_vectors_str_.resize(base_space_dimensions + 2);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "no";
            basis_vectors_str_[base_space_dimensions + 1] = "ni";

            auto aux = 1 / sqrt(ScalarType(2));
            from_actual_to_signed_metric_ = detail::make_identity_matrix<ScalarType, DimensionsAtCompileTime, MaxDimensionsAtCompileTime>(base_space_dimensions + 2);
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions, base_space_dimensions) = aux;
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions + 1, base_space_dimensions) = aux;
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions, base_space_dimensions + 1) = -aux;
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions + 1, base_space_dimensions + 1) = aux;

            from_signed_to_actual_metric_ = detail::transpose(from_actual_to_signed_metric_);
        }

        std::vector<std::string> basis_vectors_str_;
        detail::matrix_type_t<ScalarType, DimensionsAtCompileTime, DimensionsAtCompileTime, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime> from_actual_to_signed_metric_;
        detail::matrix_type_t<ScalarType, DimensionsAtCompileTime, DimensionsAtCompileTime, MaxDimensionsAtCompileTime, MaxDimensionsAtCompileTime> from_signed_to_actual_metric_;

        template<typename SomeMetricSpaceType> friend struct detail::from_actual_to_signed_metric_impl;
        template<typename SomeMetricSpaceType> friend struct detail::from_signed_to_actual_metric_impl;
    };

    namespace detail {

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct from_actual_to_signed_metric_impl<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_actual_metric) {
                return detail::prod(space_ptr->from_actual_to_signed_metric_, std::move(factors_in_actual_metric));
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct from_signed_to_actual_metric_impl<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) {
                return detail::prod(space_ptr->from_signed_to_actual_metric_, std::move(factors_in_signed_metric));
            }
        };

    }

}

#endif // __TBGAL_CONFORMAL_METRIC_SPACE_HPP__
