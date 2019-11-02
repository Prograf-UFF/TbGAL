#ifndef __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
#define __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType DimensionsAtCompileTime_, DefaultIndexType MaxDimensionsAtCompileTime_ = DimensionsAtCompileTime_>
    class EuclideanMetricSpace : public MetricSpace<EuclideanMetricSpace<DimensionsAtCompileTime_, MaxDimensionsAtCompileTime_> > {
    private:

        static_assert(DimensionsAtCompileTime_ == Dynamic || (DimensionsAtCompileTime_ > 0 && DimensionsAtCompileTime_ == MaxDimensionsAtCompileTime_), "Invalid number of base dimensions.");

        using Super = MetricSpace<EuclideanMetricSpace<DimensionsAtCompileTime_, MaxDimensionsAtCompileTime_> >;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = DefaultScalarType;

        constexpr static IndexType DimensionsAtCompileTime = DimensionsAtCompileTime_;
        constexpr static IndexType MaxDimensionsAtCompileTime = MaxDimensionsAtCompileTime_;

        inline EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        inline EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        inline EuclideanMetricSpace(IndexType dimensions) noexcept :
            Super(),
            basis_vectors_str_(0) {
            update_basis_vectors_str(dimensions);
        }

        inline EuclideanMetricSpace() noexcept :
            EuclideanMetricSpace((DimensionsAtCompileTime != Dynamic) ? DimensionsAtCompileTime : 0) {
        }

        inline EuclideanMetricSpace & operator=(EuclideanMetricSpace const &) = default;
        inline EuclideanMetricSpace & operator=(EuclideanMetricSpace &&) = default;

        inline std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }

        inline IndexType dimensions() const noexcept override {
            return basis_vectors_str_.size();
        }

        inline void set_dimensions(IndexType dimensions) noexcept {
            update_basis_vectors_str(dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType dimensions) noexcept {
            assert(dimensions > 0 && (DimensionsAtCompileTime == Dynamic || dimensions == DimensionsAtCompileTime) && (MaxDimensionsAtCompileTime == Dynamic || dimensions <= MaxDimensionsAtCompileTime));
            basis_vectors_str_.resize(dimensions);
            for (IndexType ind = 0; ind != dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        std::vector<std::string> basis_vectors_str_;
    };

    template<DefaultIndexType DimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
    struct is_metric_space<EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> > :
        std::true_type {
    };

    namespace detail {

        template<DefaultIndexType DimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
        struct apply_signed_metric_impl<EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static MatrixType const & eval(EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> const &, MatrixType const &factors_in_signed_metric) noexcept {
                return factors_in_signed_metric;
            }
        };

        template<DefaultIndexType DimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
        struct from_actual_to_signed_metric_impl<EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static MatrixType const & eval(EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> const &, MatrixType const &factors_in_actual_metric) noexcept {
                return factors_in_actual_metric;
            }
        };

        template<DefaultIndexType DimensionsAtCompileTime, DefaultIndexType MaxDimensionsAtCompileTime>
        struct from_signed_to_actual_metric_impl<EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static MatrixType const & eval(EuclideanMetricSpace<DimensionsAtCompileTime, MaxDimensionsAtCompileTime> const &, MatrixType const &factors_in_signed_metric) noexcept {
                return factors_in_signed_metric;
            }
        };

    }

}

#endif // __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
