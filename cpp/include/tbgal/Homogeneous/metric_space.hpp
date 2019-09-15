#ifndef __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
#define __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime_>
    class HomogeneousMetricSpace : public MetricSpace<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime_> > {
    private:

        static_assert(BaseSpaceDimensionsAtCompileTime_ >= 0, "Invalid number of base dimensions.");

        using Super = MetricSpace<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime_> >;
    
    public:

        using IndexType = typename Super::IndexType;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime_;
        constexpr static IndexType DimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime + 1;

        constexpr HomogeneousMetricSpace() noexcept :
            Super(),
            basis_vectors_str_() {
            for (IndexType ind = 0; ind != BaseSpaceDimensionsAtCompileTime; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[BaseSpaceDimensionsAtCompileTime] = "ep";
        }

        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        constexpr HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        constexpr HomogeneousMetricSpace & operator=(HomogeneousMetricSpace const &) = default;
        constexpr HomogeneousMetricSpace & operator=(HomogeneousMetricSpace &&) = default;

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        constexpr IndexType base_space_dimensions() const noexcept {
            return BaseSpaceDimensionsAtCompileTime;
        }

        constexpr IndexType dimensions() const noexcept override {
            return BaseSpaceDimensionsAtCompileTime + 1;
        }

    private:

        std::array<std::string, BaseSpaceDimensionsAtCompileTime + 1> basis_vectors_str_;
    };

    template<>
    class HomogeneousMetricSpace<Dynamic> : public MetricSpace<HomogeneousMetricSpace<Dynamic> > {
    private:

        using Super = MetricSpace<HomogeneousMetricSpace<Dynamic> >;
    
    public:

        using IndexType = typename Super::IndexType;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = Dynamic;
        constexpr static IndexType DimensionsAtCompileTime = Dynamic;

        inline HomogeneousMetricSpace() noexcept :
            Super(),
            basis_vectors_str_() {
            update_basis_vectors_str(0);
        }

        inline HomogeneousMetricSpace(HomogeneousMetricSpace const &) = default;
        inline HomogeneousMetricSpace(HomogeneousMetricSpace &&) = default;

        inline HomogeneousMetricSpace(IndexType base_space_dimensions) noexcept :
            Super(),
            basis_vectors_str_() {
            update_basis_vectors_str(base_space_dimensions);
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
            assert(base_space_dimensions >= 0);
            basis_vectors_str_.resize(base_space_dimensions + 1);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "ep";
        }

        std::vector<std::string> basis_vectors_str_;
    };

    namespace detail {

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime>
        struct is_metric_space<HomogeneousMetricSpace<BaseSpaceDimensionsAtCompileTime> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_HOMOGENEOUS_METRIC_SPACE_HPP__
