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

        constexpr static IndexType DimensionsAtCompileTime = Super::DimensionsAtCompileTime;
        constexpr static IndexType MaxDimensionsAtCompileTime = Super::MaxDimensionsAtCompileTime;

        inline SignedMetricSpace(SignedMetricSpace const &) = default;
        inline SignedMetricSpace(SignedMetricSpace &&) = default;

        inline SignedMetricSpace(IndexType p_dimensions, IndexType q_dimensions) :
            Super(p_dimensions, q_dimensions),
            basis_vectors_str_(),
            p_dimensions_(p_dimensions),
            q_dimensions_(q_dimensions) {
            update_basis_vectors_str(p_dimensions, q_dimensions);
        }

        inline SignedMetricSpace() :
            SignedMetricSpace((PDimensionsAtCompileTime != Dynamic && QDimensionsAtCompileTime != Dynamic) ? PDimensionsAtCompileTime : 0, (PDimensionsAtCompileTime != Dynamic && QDimensionsAtCompileTime != Dynamic) ? QDimensionsAtCompileTime : 0) {
        }

        inline SignedMetricSpace & operator=(SignedMetricSpace const &) = default;
        inline SignedMetricSpace & operator=(SignedMetricSpace &&) = default;

        inline IndexType p_dimensions() const {
            return p_dimensions_;
        }

        inline IndexType q_dimensions() const {
            return q_dimensions_;
        }

        inline std::string const & basis_vector_str(IndexType index) const override {
            return basis_vectors_str_[index];
        }
        
        inline void set_dimensions(IndexType p_dimensions, IndexType q_dimensions) {
            Super::set_dimensions(p_dimensions, q_dimensions);
            p_dimensions_ = p_dimensions;
            q_dimensions_ = q_dimensions;
            update_basis_vectors_str(p_dimensions, q_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType p_dimensions, IndexType q_dimensions) {
            basis_vectors_str_.resize(p_dimensions + q_dimensions);
            for (IndexType ind = 0, ind_end = p_dimensions + q_dimensions; ind != ind_end; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        std::vector<std::string> basis_vectors_str_;
        IndexType p_dimensions_;
        IndexType q_dimensions_;
    };

    namespace detail {

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
        struct from_actual_to_signed_metric_impl<SignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(SignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> const *, MatrixType &&factors_in_actual_metric) {
                return std::move(factors_in_actual_metric);
            }
        };

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
        struct from_signed_to_actual_metric_impl<SignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(SignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) {
                return std::move(factors_in_signed_metric);
            }
        };

    }

}

#endif // __TBGAL_SIGNED_METRIC_SPACE_HPP__
