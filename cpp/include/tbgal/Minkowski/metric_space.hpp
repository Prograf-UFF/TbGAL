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

#ifndef __TBGAL_MINKOWSKI_METRIC_SPACE_HPP__
#define __TBGAL_MINKOWSKI_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime_, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime_ = BaseSpaceDimensionsAtCompileTime_>
    class MinkowskiMetricSpace : public BaseSignedMetricSpace<BaseSpaceDimensionsAtCompileTime_ != Dynamic ? BaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic, 1, MaxBaseSpaceDimensionsAtCompileTime_ != Dynamic ? MaxBaseSpaceDimensionsAtCompileTime_ + 2 : Dynamic> {
    private:

        using Super = BaseSignedMetricSpace<BaseSpaceDimensionsAtCompileTime_ != Dynamic ? BaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic, 1, MaxBaseSpaceDimensionsAtCompileTime_ != Dynamic ? MaxBaseSpaceDimensionsAtCompileTime_ + 2 : Dynamic>;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = typename Super::ScalarType;

        constexpr static IndexType DimensionsAtCompileTime = Super::DimensionsAtCompileTime;
        constexpr static IndexType MaxDimensionsAtCompileTime = Super::MaxDimensionsAtCompileTime;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime_;
        constexpr static IndexType MaxBaseSpaceDimensionsAtCompileTime = MaxBaseSpaceDimensionsAtCompileTime_;

        inline MinkowskiMetricSpace(MinkowskiMetricSpace const &) = default;
        inline MinkowskiMetricSpace(MinkowskiMetricSpace &&) = default;

        inline MinkowskiMetricSpace(IndexType base_space_dimensions) noexcept :
            Super(base_space_dimensions + 1, 1),
            basis_vectors_str_() {
            update_basis_vectors_str(base_space_dimensions);
        }

        inline MinkowskiMetricSpace() noexcept :
            MinkowskiMetricSpace((BaseSpaceDimensionsAtCompileTime != Dynamic) ? BaseSpaceDimensionsAtCompileTime : 0) {
        }

        inline MinkowskiMetricSpace & operator=(MinkowskiMetricSpace const &) = default;
        inline MinkowskiMetricSpace & operator=(MinkowskiMetricSpace &&) = default;

        inline std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        inline IndexType base_space_dimensions() const noexcept {
            return Super::dimensions() - 2;
        }

        inline void set_base_space_dimensions(IndexType base_space_dimensions) noexcept {
            Super::set_dimensions(base_space_dimensions + 1, 1);
            update_basis_vectors_str(base_space_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType base_space_dimensions) noexcept {
            basis_vectors_str_.resize(base_space_dimensions + 2);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "ep";
            basis_vectors_str_[base_space_dimensions + 1] = "em";
        }

        std::vector<std::string> basis_vectors_str_;
    };

    namespace detail {

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct from_actual_to_signed_metric_impl<MinkowskiMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(MinkowskiMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const *, MatrixType &&factors_in_actual_metric) noexcept {
                return std::move(factors_in_actual_metric);
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct from_signed_to_actual_metric_impl<MinkowskiMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(MinkowskiMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) noexcept {
                return std::move(factors_in_signed_metric);
            }
        };

    }

}

#endif // __TBGAL_MINKOWSKI_METRIC_SPACE_HPP__
