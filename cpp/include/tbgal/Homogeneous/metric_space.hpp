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

        inline HomogeneousMetricSpace(IndexType base_space_dimensions) noexcept :
            Super(base_space_dimensions + 1, 0),
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
            return Super::dimensions() - 1;
        }

        inline void set_base_space_dimensions(IndexType base_space_dimensions) noexcept {
            Super::set_dimensions(base_space_dimensions + 1, 0);
            update_basis_vectors_str(base_space_dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType base_space_dimensions) noexcept {
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
