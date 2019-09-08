#ifndef __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
#define __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__

namespace tbgal {

    template<typename ScalarType, DefaultIndexType BaseDimensionsAtCompileTime>
    class HomogeneousMetricSpace : public MetricSpace<detail::identity_matrix_type_t<ScalarType, BaseDimensionsAtCompileTime + 1> > {
    private:

        static_assert(BaseDimensionsAtCompileTime > 0, "Invalid number of base dimensions.");

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, BaseDimensionsAtCompileTime + 1> >;
    
    public:

        using IndexType = typename Super::IndexType;

        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        constexpr HomogeneousMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, BaseDimensionsAtCompileTime + 1>(BaseDimensionsAtCompileTime + 1)),
            basis_vectors_str_() {
            for (DefaultIndexType ind = 0; ind != BaseDimensionsAtCompileTime; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[BaseDimensionsAtCompileTime] = "ep";
        }

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }

    private:

        std::array<std::string, BaseDimensionsAtCompileTime + 1> basis_vectors_str_;
    };

    template<typename ScalarType>
    class HomogeneousMetricSpace<ScalarType, Dynamic> : public MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> > {
    private:

        using Super = MetricSpace<detail::identity_matrix_type_t<ScalarType, Dynamic> >;
    
    public:

        using IndexType = typename Super::IndexType;

        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        constexpr HomogeneousMetricSpace() noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(1)),
            basis_vectors_str_() {
        }

        constexpr HomogeneousMetricSpace(IndexType base_dimensions) noexcept :
            Super(detail::make_identity_matrix<ScalarType, Dynamic>(base_dimensions + 1)),
            basis_vectors_str_() {
            assert(base_dimensions > 0);
            update_basis_vectors_str();
        }

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }

        constexpr void set_base_dimensions(IndexType base_dimensions) noexcept {
            assert(base_dimensions > 0);
            Super::metric_matrix_ = detail::make_identity_matrix<ScalarType, Dynamic>(base_dimensions + 1);
        }

    private:

        constexpr void update_basis_vectors_str() noexcept {
            basis_vectors_str_.resize(dimensions());
            for (IndexType ind = 1; ind != dimensions(); ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind);
            }
            basis_vectors_str_[dimensions() - 1] = "ep";
        }

        std::vector<std::string> basis_vectors_str_;
    };

    namespace detail {

        template<typename ScalarType, DefaultIndexType BaseDimensionsAtCompileTime>
        struct is_metric_space<HomogeneousMetricSpace<ScalarType, BaseDimensionsAtCompileTime> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
