#ifndef __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
#define __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename ScalarType, DefaultIndexType DimensionsAtCompileTime>
    class EuclideanMetricSpace : public MetricSpace<detail::identity_matrix_type_t<ScalarType, DimensionsAtCompileTime> > {
    private:

        static_assert(DimensionsAtCompileTime > 0, "Invalid number of dimensions.");

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, DimensionsAtCompileTime> >;
    
    public:

        using IndexType = typename Super::IndexType;

        constexpr EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        constexpr EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        constexpr EuclideanMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, DimensionsAtCompileTime>(DimensionsAtCompileTime)),
            basis_vectors_str_() {
            for (DefaultIndexType ind = 0; ind != DimensionsAtCompileTime; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }

    private:

        std::array<std::string, DimensionsAtCompileTime> basis_vectors_str_;
    };

    template<typename ScalarType>
    class EuclideanMetricSpace<ScalarType, Dynamic> : public MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> > {
    private:

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> >;
    
    public:

        using IndexType = typename Super::IndexType;

        constexpr EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        constexpr EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        constexpr EuclideanMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(0)),
            basis_vectors_str_() {
        }

        constexpr EuclideanMetricSpace(IndexType dimensions) noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(dimensions)),
            basis_vectors_str_() {
            assert(dimensions > 0);
            update_basis_vectors_str();
        }

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }

        constexpr void set_dimensions(IndexType dimensions) noexcept {
            assert(dimensions > 0);
            Super::metric_matrix_ = detail::make_identity_matrix<ScalarType, Dynamic>(dimensions);
            update_basis_vectors_str();
        }

    private:

        constexpr void update_basis_vectors_str() noexcept {
            basis_vectors_str_.resize(dimensions());
            for (IndexType ind = 0; ind != dimensions(); ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        std::vector<std::string> basis_vectors_str_;
    };

    namespace detail {

        template<typename ScalarType, DefaultIndexType DimensionsAtCompileTime>
        struct is_metric_space<EuclideanMetricSpace<ScalarType, DimensionsAtCompileTime> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
