#ifndef __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
#define __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType DimensionsAtCompileTime_, DefaultIndexType MaxDimensionsAtCompileTime_ = DimensionsAtCompileTime_>
    class EuclideanMetricSpace : public BaseSignedMetricSpace<DimensionsAtCompileTime_, 0, MaxDimensionsAtCompileTime_> {
    private:

        using Super = BaseSignedMetricSpace<DimensionsAtCompileTime_, 0, MaxDimensionsAtCompileTime_>;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = typename Super::ScalarType;

        constexpr static IndexType DimensionsAtCompileTime = Super::DimensionsAtCompileTime;
        constexpr static IndexType MaxDimensionsAtCompileTime = Super::MaxDimensionsAtCompileTime;

        inline EuclideanMetricSpace(EuclideanMetricSpace const &) = default;
        inline EuclideanMetricSpace(EuclideanMetricSpace &&) = default;

        inline EuclideanMetricSpace(IndexType dimensions) noexcept :
            Super(dimensions, 0),
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

        inline void set_dimensions(IndexType dimensions) noexcept {
            Super::set_dimensions(dimensions, 0);
            update_basis_vectors_str(dimensions);
        }

    private:

        inline void update_basis_vectors_str(IndexType dimensions) noexcept {
            basis_vectors_str_.resize(dimensions);
            for (IndexType ind = 0; ind != dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
        }

        std::vector<std::string> basis_vectors_str_;
    };

}

#endif // __TBGAL_EUCLIDEAN_METRIC_SPACE_HPP__
