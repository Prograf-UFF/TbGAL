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

#ifndef __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
#define __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime_, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime_ = BaseSpaceDimensionsAtCompileTime_>
    class HomogeneousMetricSpace : public BaseSignedMetricSpace<BaseSpaceDimensionsAtCompileTime_ != Dynamic ? BaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic, 0, MaxBaseSpaceDimensionsAtCompileTime_ != Dynamic ? MaxBaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic> {
    private:

        using Super = BaseSignedMetricSpace<BaseSpaceDimensionsAtCompileTime_ != Dynamic ? BaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic, 0, MaxBaseSpaceDimensionsAtCompileTime_ != Dynamic ? MaxBaseSpaceDimensionsAtCompileTime_ + 1 : Dynamic>;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = typename Super::ScalarType;

        constexpr static IndexType DimensionsAtCompileTime = Super::DimensionsAtCompileTime;
        constexpr static IndexType MaxDimensionsAtCompileTime = Super::MaxDimensionsAtCompileTime;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime_;
        constexpr static IndexType MaxBaseSpaceDimensionsAtCompileTime = MaxBaseSpaceDimensionsAtCompileTime_;

        inline HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        inline HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        inline HomogeneousMetricSpace(IndexType base_space_dimensions) :
            Super(base_space_dimensions + 1, 0),
            basis_vectors_str_() {
            update_basis_vectors_str(base_space_dimensions);
        }

        inline HomogeneousMetricSpace() :
            HomogeneousMetricSpace((BaseSpaceDimensionsAtCompileTime != Dynamic) ? BaseSpaceDimensionsAtCompileTime : 0) {
        }

        inline HomogeneousMetricSpace & operator=(HomogeneousMetricSpace const &) = default;
        inline HomogeneousMetricSpace & operator=(HomogeneousMetricSpace &&) = default;

        inline std::string const & basis_vector_str(IndexType index) const override {
            return basis_vectors_str_[index];
        }
        
        inline IndexType base_space_dimensions() const {
            return Super::dimensions() - 1;
        }

        inline void set_base_space_dimensions(IndexType base_space_dimensions) {
            Super::set_dimensions(base_space_dimensions + 1, 0);
            update_basis_vectors_str(base_space_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType base_space_dimensions) {
            basis_vectors_str_.resize(base_space_dimensions + 1);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "ep";
        }

        std::vector<std::string> basis_vectors_str_;
    };

}

#endif // __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
