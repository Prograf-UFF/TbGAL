#ifndef __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
#define __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime_, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime_ = BaseSpaceDimensionsAtCompileTime_>
    class HomogeneousMetricSpace : public MetricSpace<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime_, MaxBaseSpaceDimensionsAtCompileTime_> > {
    private:

        static_assert(BaseSpaceDimensionsAtCompileTime_ == Dynamic || (BaseSpaceDimensionsAtCompileTime_ >= 0 && BaseSpaceDimensionsAtCompileTime_ == MaxBaseSpaceDimensionsAtCompileTime_), "Invalid number of base dimensions.");

        using Super = MetricSpace<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime_, MaxBaseSpaceDimensionsAtCompileTime_> >;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = DefaultScalarType;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime_;
        constexpr static IndexType MaxBaseSpaceDimensionsAtCompileTime = MaxBaseSpaceDimensionsAtCompileTime_;

        constexpr static IndexType DimensionsAtCompileTime = (BaseSpaceDimensionsAtCompileTime != Dynamic) ? (BaseSpaceDimensionsAtCompileTime + 1) : Dynamic;
        constexpr static IndexType MaxDimensionsAtCompileTime = (MaxBaseSpaceDimensionsAtCompileTime != Dynamic) ? (MaxBaseSpaceDimensionsAtCompileTime + 1) : Dynamic;

        inline HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        inline HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        inline HomogeneousMetricSpace(IndexType base_space_dimensions) noexcept :
            Super(),
            basis_vectors_str_() {
            update_basis_vectors_str(base_space_dimensions);
        }

        inline HomogeneousMetricSpace() noexcept :
            HomogeneousMetricSpace((BaseSpaceDimensionsAtCompileTime != Dynamic) ? BaseSpaceDimensionsAtCompileTime : 0) {
        }

        inline HomogeneousMetricSpace & operator=(HomogeneousMetricSpace const &) = default;
        inline HomogeneousMetricSpace & operator=(HomogeneousMetricSpace &&) = default;

        inline std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        inline IndexType base_space_dimensions() const noexcept {
            return basis_vectors_str_.size() - 1;
        }

        inline IndexType dimensions() const noexcept override {
            return basis_vectors_str_.size();
        }

        inline void set_base_space_dimensions(IndexType base_space_dimensions) noexcept {
            update_basis_vectors_str(base_space_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType base_space_dimensions) noexcept {
            assert(base_space_dimensions >= 0 && (BaseSpaceDimensionsAtCompileTime == Dynamic || base_space_dimensions == BaseSpaceDimensionsAtCompileTime) && (MaxBaseSpaceDimensionsAtCompileTime == Dynamic || base_space_dimensions <= MaxBaseSpaceDimensionsAtCompileTime));
            basis_vectors_str_.resize(base_space_dimensions + 1);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "ep";
        }

        std::vector<std::string> basis_vectors_str_;
    };

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
    struct is_metric_space<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > :
        std::true_type {
    };

    namespace detail {

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct apply_signed_metric_impl<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static MatrixType const & eval(HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const &, MatrixType const &factors_in_signed_metric) noexcept {
                return factors_in_signed_metric;
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct from_actual_to_signed_metric_impl<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const &, MatrixType &&factors_in_actual_metric) noexcept {
                return std::move(factors_in_actual_metric);
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct from_signed_to_actual_metric_impl<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const &, MatrixType &&factors_in_signed_metric) noexcept {
                return std::move(factors_in_signed_metric);
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime, DefaultIndexType MaxBaseSpaceDimensionsAtCompileTime>
        struct metric_factor_impl<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime> const &space, MatrixType const &) noexcept {
                using ResultingScalarType = std::common_type_t<scalar_type_t<MatrixType>, typename HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime, MaxBaseSpaceDimensionsAtCompileTime>::ScalarType>;
                return ResultingScalarType(1);
            }
        };

    }

}

#endif // __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
