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

#ifndef __TBGAL_METRIC_SPACE_HPP__
#define __TBGAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
    class BaseSignedMetricSpace;
    
    namespace detail {

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) apply_signed_metric(BaseSignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> const *, MatrixType &&);

    }

    template<DefaultIndexType PDimensionsAtCompileTime_, DefaultIndexType QDimensionsAtCompileTime_, DefaultIndexType MaxDimensionsAtCompileTime_ = ((PDimensionsAtCompileTime_ != Dynamic && QDimensionsAtCompileTime_ != Dynamic) ? (PDimensionsAtCompileTime_ + QDimensionsAtCompileTime_) : Dynamic)>
    class BaseSignedMetricSpace {
    public:

        static_assert(PDimensionsAtCompileTime_ == Dynamic || QDimensionsAtCompileTime_ == Dynamic || (PDimensionsAtCompileTime_ >= 0 && QDimensionsAtCompileTime_ >= 0 && (PDimensionsAtCompileTime_ + QDimensionsAtCompileTime_) > 0 && (PDimensionsAtCompileTime_ + QDimensionsAtCompileTime_) == MaxDimensionsAtCompileTime_), "Invalid number of dimensions.");

        using IndexType = DefaultIndexType;
        using ScalarType = DefaultScalarType;

        constexpr static DefaultIndexType PDimensionsAtCompileTime = PDimensionsAtCompileTime_;
        constexpr static DefaultIndexType QDimensionsAtCompileTime = QDimensionsAtCompileTime_;

        constexpr static IndexType DimensionsAtCompileTime = (PDimensionsAtCompileTime != Dynamic && QDimensionsAtCompileTime != Dynamic) ? (PDimensionsAtCompileTime + QDimensionsAtCompileTime) : Dynamic;
        constexpr static IndexType MaxDimensionsAtCompileTime = MaxDimensionsAtCompileTime_;

        inline IndexType dimensions() const {
            return detail::rows(signed_metric_);
        }

        virtual std::string const & basis_vector_str(IndexType) const = 0;

    protected:

        inline BaseSignedMetricSpace(BaseSignedMetricSpace const &) = default;
        inline BaseSignedMetricSpace(BaseSignedMetricSpace &&) = default;

        inline BaseSignedMetricSpace(IndexType p_dimensions, IndexType q_dimensions) :
            signed_metric_() {
            update_basis(p_dimensions, q_dimensions);
        }

        inline BaseSignedMetricSpace & operator=(BaseSignedMetricSpace const &) = default;
        inline BaseSignedMetricSpace & operator=(BaseSignedMetricSpace &&) = default;

        inline void set_dimensions(IndexType p_dimensions, IndexType q_dimensions) {
            update_basis(p_dimensions, q_dimensions);
        }

    private:

        inline void update_basis(IndexType p_dimensions, IndexType q_dimensions) {
            if (!(p_dimensions >= 0 && q_dimensions >= 0 && (p_dimensions + q_dimensions) >= 0 && (PDimensionsAtCompileTime == Dynamic || PDimensionsAtCompileTime == p_dimensions) && (QDimensionsAtCompileTime == Dynamic || QDimensionsAtCompileTime == q_dimensions) && (MaxDimensionsAtCompileTime == Dynamic || MaxDimensionsAtCompileTime >= (p_dimensions + q_dimensions)))) {
                throw NotSupportedError("Invalid values for the (p, q) signature of the current space.");
            }
            signed_metric_ = detail::make_diagonal_matrix<ScalarType, DimensionsAtCompileTime, MaxDimensionsAtCompileTime>(p_dimensions + q_dimensions);
            for (IndexType ind = 0; ind != p_dimensions; ++ind) {
                detail::coeff(signed_metric_, ind, ind) = 1;
            }
            for (IndexType ind = p_dimensions; ind != (p_dimensions + q_dimensions); ++ind) {
                detail::coeff(signed_metric_, ind, ind) = -1;
            }
        }

        detail::diagonal_matrix_type_t<ScalarType, DimensionsAtCompileTime, MaxDimensionsAtCompileTime> signed_metric_;

        template<DefaultIndexType SomePDimensionsAtCompileTime, DefaultIndexType SomeQDimensionsAtCompileTime, DefaultIndexType SomeMaxDimensionsAtCompileTime, typename SomeMatrixType> friend constexpr decltype(auto) detail::apply_signed_metric(BaseSignedMetricSpace<SomePDimensionsAtCompileTime, SomeQDimensionsAtCompileTime, SomeMaxDimensionsAtCompileTime> const *, SomeMatrixType &&);
    };

    template<DefaultIndexType PDimensionsAtCompileTime_, DefaultIndexType MaxDimensionsAtCompileTime_>
    class BaseSignedMetricSpace<PDimensionsAtCompileTime_, 0, MaxDimensionsAtCompileTime_> {
    public:

        static_assert(PDimensionsAtCompileTime_ == Dynamic || (PDimensionsAtCompileTime_ > 0 && PDimensionsAtCompileTime_ == MaxDimensionsAtCompileTime_), "Invalid number of dimensions.");

        using IndexType = DefaultIndexType;
        using ScalarType = DefaultScalarType;

        constexpr static DefaultIndexType PDimensionsAtCompileTime = PDimensionsAtCompileTime_;
        constexpr static DefaultIndexType QDimensionsAtCompileTime = 0;

        constexpr static IndexType DimensionsAtCompileTime = (PDimensionsAtCompileTime != Dynamic) ? PDimensionsAtCompileTime : Dynamic;
        constexpr static IndexType MaxDimensionsAtCompileTime = MaxDimensionsAtCompileTime_;

        inline IndexType dimensions() const {
            return dimensions_;
        }

        virtual std::string const & basis_vector_str(IndexType) const = 0;

    protected:

        inline BaseSignedMetricSpace(BaseSignedMetricSpace const &) = default;
        inline BaseSignedMetricSpace(BaseSignedMetricSpace &&) = default;

        inline BaseSignedMetricSpace(IndexType p_dimensions, IndexType q_dimensions) :
            dimensions_(p_dimensions) {
            if (!(p_dimensions >= 0 && q_dimensions == 0 && (PDimensionsAtCompileTime == Dynamic || PDimensionsAtCompileTime == p_dimensions) && (MaxDimensionsAtCompileTime == Dynamic || MaxDimensionsAtCompileTime >= p_dimensions))) {
                throw NotSupportedError("Invalid values for the (p, q) signature of the current space.");
            }
        }

        inline BaseSignedMetricSpace & operator=(BaseSignedMetricSpace const &) = default;
        inline BaseSignedMetricSpace & operator=(BaseSignedMetricSpace &&) = default;

        inline void set_dimensions(IndexType p_dimensions, IndexType q_dimensions) {
            if (!(p_dimensions >= 0 && q_dimensions == 0 && (PDimensionsAtCompileTime == Dynamic || PDimensionsAtCompileTime == p_dimensions) && (MaxDimensionsAtCompileTime == Dynamic || MaxDimensionsAtCompileTime >= p_dimensions))) {
                throw NotSupportedError("Invalid value for the dimensionaliry of the current space.");
            }
            dimensions_ = p_dimensions;
        }

    private:

        IndexType dimensions_;
    };

    namespace detail {

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) apply_signed_metric(BaseSignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) {
            return prod(space_ptr->signed_metric_, std::move(factors_in_signed_metric));
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) apply_signed_metric(BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime> const *, MatrixType &&factors_in_signed_metric) {
            return std::move(factors_in_signed_metric);
        }

        template<typename MetricSpaceType, typename MatrixType, typename = std::enable_if_t<MetricSpaceType::QDimensionsAtCompileTime != 0, int> >
        constexpr decltype(auto) from_actual_to_signed_metric(MetricSpaceType const *space_ptr, MatrixType &&factors_in_actual_metric) {
            return from_actual_to_signed_metric_impl<MetricSpaceType>::eval(space_ptr, std::move(factors_in_actual_metric));
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) from_actual_to_signed_metric(BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime> const *, MatrixType &&factors_in_actual_metric) {
            return std::move(factors_in_actual_metric);
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) from_outer_to_geometric_factors(BaseSignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) {
            auto new_factors_in_signed_metric = evaluate(prod(factors_in_signed_metric, es_eigenvectors_matrix(evaluate(prod(transpose(factors_in_signed_metric), apply_signed_metric(space_ptr, factors_in_signed_metric))))));
            return std::make_tuple(determinant(prod(transpose(new_factors_in_signed_metric), factors_in_signed_metric)), new_factors_in_signed_metric);
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) from_outer_to_geometric_factors(BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime> const *, MatrixType &&factors_in_signed_metric) {
            using MetricSpaceType = BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime>;
            using ResultingScalarType = std::common_type_t<scalar_type_t<MatrixType>, typename MetricSpaceType::ScalarType>;
            return std::make_tuple(ResultingScalarType(1), std::move(factors_in_signed_metric));
        }

        template<typename MetricSpaceType, typename MatrixType, typename = std::enable_if_t<MetricSpaceType::QDimensionsAtCompileTime != 0, int> >
        constexpr decltype(auto) from_signed_to_actual_metric(MetricSpaceType const *space_ptr, MatrixType &&factors_in_signed_metric) {
            return from_signed_to_actual_metric_impl<MetricSpaceType>::eval(space_ptr, std::move(factors_in_signed_metric));
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) from_signed_to_actual_metric(BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime> const *, MatrixType &&factors_in_signed_metric) {
            return std::move(factors_in_signed_metric);
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType QDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) metric_factor(BaseSignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) {
            using MetricSpaceType = BaseSignedMetricSpace<PDimensionsAtCompileTime, QDimensionsAtCompileTime, MaxDimensionsAtCompileTime>;
            using ResulingScalarType = std::common_type_t<scalar_type_t<MatrixType>, typename MetricSpaceType::ScalarType>;
            using IndexType = index_type_t<MatrixType>;
            
            constexpr DefaultIndexType DimensionsAtCompileTime = MetricSpaceType::DimensionsAtCompileTime;

            ResulingScalarType result = 1;
            for (IndexType col = 0, col_end = cols(factors_in_signed_metric); col != col_end; ++col) {
                auto column = block<DimensionsAtCompileTime, 1>(factors_in_signed_metric, 0, col, space_ptr->dimensions(), 1);
                result *= coeff(prod(transpose(column), apply_signed_metric(space_ptr, column)), 0, 0);
            }
            return result;
        }

        template<DefaultIndexType PDimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime, typename MatrixType>
        constexpr decltype(auto) metric_factor(BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime> const *space_ptr, MatrixType &&factors_in_signed_metric) {
            using MetricSpaceType = BaseSignedMetricSpace<PDimensionsAtCompileTime, 0, MaxDimensionsAtCompileTime>;
            using ResulingScalarType = std::common_type_t<scalar_type_t<MatrixType>, typename MetricSpaceType::ScalarType>;
            return ResulingScalarType(1);
        }

    }

}

#endif // __TBGAL_METRIC_SPACE_HPP__
