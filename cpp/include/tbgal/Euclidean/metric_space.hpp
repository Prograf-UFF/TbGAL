#ifndef __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
#define __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType DimensionsAtCompileTime_>
    class EuclideanMetricSpace : public MetricSpace<EuclideanMetricSpace<DimensionsAtCompileTime_> > {
    private:

        static_assert(DimensionsAtCompileTime_ > 0, "Invalid number of dimensions.");

        using Super = MetricSpace<EuclideanMetricSpace<DimensionsAtCompileTime_> >;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = DefaultScalarType;

        constexpr static IndexType DimensionsAtCompileTime = DimensionsAtCompileTime_;

        constexpr EuclideanMetricSpace() noexcept :
            Super(),
            basis_vectors_str_() {
            for (IndexType ind = 0; ind != DimensionsAtCompileTime; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        constexpr EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        constexpr EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        constexpr EuclideanMetricSpace & operator=(EuclideanMetricSpace const &) = default;
        constexpr EuclideanMetricSpace & operator=(EuclideanMetricSpace &&) = default;

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        constexpr IndexType dimensions() const noexcept override {
            return DimensionsAtCompileTime;
        }

    private:

        std::array<std::string, DimensionsAtCompileTime> basis_vectors_str_;
    };

    template<>
    class EuclideanMetricSpace<Dynamic> : public MetricSpace<EuclideanMetricSpace<Dynamic> > {
    private:

        using Super = MetricSpace<EuclideanMetricSpace<Dynamic> >;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = DefaultScalarType;

        constexpr static IndexType DimensionsAtCompileTime = Dynamic;

        inline EuclideanMetricSpace() noexcept :
            Super(),
            basis_vectors_str_(0) {
        }

        inline EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        inline EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        inline EuclideanMetricSpace(IndexType dimensions) noexcept :
            Super(),
            basis_vectors_str_(0) {
            update_basis_vectors_str(dimensions);
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
            assert(dimensions > 0);
            basis_vectors_str_.resize(dimensions);
            for (IndexType ind = 0; ind != dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        std::vector<std::string> basis_vectors_str_;
    };

    template<DefaultIndexType DimensionsAtCompileTime>
    struct is_metric_space<EuclideanMetricSpace<DimensionsAtCompileTime> > :
        std::true_type {
    };

    namespace detail {

        template<DefaultIndexType DimensionsAtCompileTime>
        struct from_actual_to_orthogonal_metric_impl<EuclideanMetricSpace<DimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static MatrixType const & eval(EuclideanMetricSpace<DimensionsAtCompileTime> const &, MatrixType const &factors) noexcept {
                return factors;
            }
        };

        template<DefaultIndexType DimensionsAtCompileTime>
        struct from_orthogonal_to_actual_metric_impl<EuclideanMetricSpace<DimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static MatrixType const & eval(EuclideanMetricSpace<DimensionsAtCompileTime> const &, MatrixType const &factors) noexcept {
                return factors;
            }
        };

        template<DefaultIndexType DimensionsAtCompileTime>
        struct orthogonal_metric_factor_impl<EuclideanMetricSpace<DimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(EuclideanMetricSpace<DimensionsAtCompileTime> const &, MatrixType const &, DefaultIndexType) noexcept {
                using ScalarType = typename EuclideanMetricSpace<DimensionsAtCompileTime>::ScalarType;
                return ScalarType(1);
            }
        };

    }

}

#endif // __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
