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

#ifndef __TBGAL_SIGNED_METRIC_SPACE_HPP__
#define __TBGAL_SIGNED_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType PDimensionsAtCompileTime_, DefaultIndexType QDimensionsAtCompileTime_, DefaultIndexType MaxDimensionsAtCompileTime_ = ((PDimensionsAtCompileTime_ != Dynamic && QDimensionsAtCompileTime_ != Dynamic) ? (PDimensionsAtCompileTime_ + QDimensionsAtCompileTime_) : Dynamic)>
    class SignedMetricSpace : public BaseSignedMetricSpace<PDimensionsAtCompileTime_, QDimensionsAtCompileTime_, MaxDimensionsAtCompileTime_> {
    private:

        using Super = BaseSignedMetricSpace<PDimensionsAtCompileTime_, QDimensionsAtCompileTime_, MaxDimensionsAtCompileTime_>;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = typename Super::ScalarType;

        constexpr static IndexType PDimensionsAtCompileTime = Super::PDimensionsAtCompileTime;
        constexpr static IndexType QDimensionsAtCompileTime = Super::QDimensionsAtCompileTime;
        constexpr static IndexType MaxDimensionsAtCompileTime = Super::MaxDimensionsAtCompileTime;

        inline SignedMetricSpace(SignedMetricSpace const &) = default;
        inline SignedMetricSpace(SignedMetricSpace &&) = default;

        inline SignedMetricSpace(IndexType p_dimensions, IndexType q_dimensions) noexcept :
            Super(p_dimensions, q_dimensions),
            basis_vectors_str_(),
            p_dimensions_(p_dimensions),
            q_dimensions_(q_dimensions) {
            update_basis_vectors_str(p_dimensions, q_dimensions);
        }

        inline SignedMetricSpace() noexcept :
            SignedMetricSpace((PDimensionsAtCompileTime != Dynamic && QDimensionsAtCompileTime != Dynamic) ? PDimensionsAtCompileTime : 0, (PDimensionsAtCompileTime != Dynamic && QDimensionsAtCompileTime != Dynamic) ? QDimensionsAtCompileTime : 0) {
        }

        inline SignedMetricSpace & operator=(SignedMetricSpace const &) = default;
        inline SignedMetricSpace & operator=(SignedMetricSpace &&) = default;

        inline IndexType p_dimensions() const noexcept {
            return p_dimensions_;
        }

        inline IndexType q_dimensions() const noexcept {
            return q_dimensions_;
        }

        inline std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        inline void set_dimensions(IndexType p_dimensions, IndexType q_dimensions) noexcept {
            Super::set_dimensions(p_dimensions, q_dimensions);
            p_dimensions_ = p_dimensions;
            q_dimensions_ = q_dimensions;
            update_basis_vectors_str(p_dimensions, q_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType p_dimensions, IndexType q_dimensions) noexcept {
            basis_vectors_str_.resize(base_space_dimensions + 1);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "ep";
        }

        std::vector<std::string> basis_vectors_str_;
        IndexType p_dimensions_;
        IndexType q_dimensions_;
    };

}

#endif // __TBGAL_SIGNED_METRIC_SPACE_HPP__
