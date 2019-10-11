#ifndef __TBGAL_CONFORMAL_METRIC_SPACE_HPP__
#define __TBGAL_CONFORMAL_METRIC_SPACE_HPP__

namespace tbgal {

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime_>
    class ConformalMetricSpace : public MetricSpace<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime_> > {
    private:

        static_assert(BaseSpaceDimensionsAtCompileTime_ >= 0, "Invalid number of base dimensions.");

        using Super = MetricSpace<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime_> >;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = DefaultScalarType;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime_;
        constexpr static IndexType DimensionsAtCompileTime = BaseSpaceDimensionsAtCompileTime + 2;

        constexpr ConformalMetricSpace() noexcept :
            Super(),
            basis_vectors_str_(),
            from_actual_to_signed_metric_(detail::make_identity_matrix<ScalarType, DimensionsAtCompileTime>(DimensionsAtCompileTime)),
            from_signed_to_actual_metric_() {
            for (IndexType ind = 0; ind != BaseSpaceDimensionsAtCompileTime; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[BaseSpaceDimensionsAtCompileTime] = "no";
            basis_vectors_str_[BaseSpaceDimensionsAtCompileTime + 1] = "ni";
            auto aux = 1 / sqrt(ScalarType(2));
            detail::coeff(from_actual_to_signed_metric_, BaseSpaceDimensionsAtCompileTime, BaseSpaceDimensionsAtCompileTime) = aux;
            detail::coeff(from_actual_to_signed_metric_, BaseSpaceDimensionsAtCompileTime + 1, BaseSpaceDimensionsAtCompileTime) = aux;
            detail::coeff(from_actual_to_signed_metric_, BaseSpaceDimensionsAtCompileTime, BaseSpaceDimensionsAtCompileTime + 1) = -aux;
            detail::coeff(from_actual_to_signed_metric_, BaseSpaceDimensionsAtCompileTime + 1, BaseSpaceDimensionsAtCompileTime + 1) = aux;
            from_signed_to_actual_metric_ = detail::transpose(from_actual_to_signed_metric_);
        }

        constexpr ConformalMetricSpace(ConformalMetricSpace const &) = default;
        constexpr ConformalMetricSpace(ConformalMetricSpace &&) = default;

        constexpr ConformalMetricSpace & operator=(ConformalMetricSpace const &) = default;
        constexpr ConformalMetricSpace & operator=(ConformalMetricSpace &&) = default;

        constexpr std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        constexpr IndexType base_space_dimensions() const noexcept {
            return BaseSpaceDimensionsAtCompileTime;
        }

        constexpr IndexType dimensions() const noexcept override {
            return BaseSpaceDimensionsAtCompileTime + 2;
        }

    private:

        std::array<std::string, DimensionsAtCompileTime> basis_vectors_str_;
        detail::matrix_type_t<ScalarType, DimensionsAtCompileTime, DimensionsAtCompileTime> from_actual_to_signed_metric_;
        detail::matrix_type_t<ScalarType, DimensionsAtCompileTime, DimensionsAtCompileTime> from_signed_to_actual_metric_;

        template<DefaultIndexType SomeBaseSpaceDimensionsAtCompileTime> friend struct detail::from_actual_to_signed_metric_impl;
        template<DefaultIndexType SomeBaseSpaceDimensionsAtCompileTime> friend struct detail::from_signed_to_actual_metric_impl;
    };

    template<>
    class ConformalMetricSpace<Dynamic> : public MetricSpace<ConformalMetricSpace<Dynamic> > {
    private:

        using Super = MetricSpace<ConformalMetricSpace<Dynamic> >;
    
    public:

        using IndexType = typename Super::IndexType;
        using ScalarType = DefaultScalarType;

        constexpr static IndexType BaseSpaceDimensionsAtCompileTime = Dynamic;
        constexpr static IndexType DimensionsAtCompileTime = Dynamic;

        inline ConformalMetricSpace() noexcept :
            Super(),
            basis_vectors_str_() {
            update_basis_vectors_str(0);
        }

        inline ConformalMetricSpace(ConformalMetricSpace const &) = default;
        inline ConformalMetricSpace(ConformalMetricSpace &&) = default;

        inline ConformalMetricSpace(IndexType base_space_dimensions) noexcept :
            Super(),
            basis_vectors_str_(),
            from_actual_to_signed_metric_(),
            from_signed_to_actual_metric_() {
            update_basis_vectors_str(base_space_dimensions);
        }

        inline ConformalMetricSpace & operator=(ConformalMetricSpace const &) = default;
        inline ConformalMetricSpace & operator=(ConformalMetricSpace &&) = default;

        inline std::string const & basis_vector_str(IndexType index) const noexcept override {
            return basis_vectors_str_[index];
        }
        
        inline IndexType base_space_dimensions() const noexcept {
            return basis_vectors_str_.size() - 2;
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
            basis_vectors_str_.resize(base_space_dimensions + 2);
            for (IndexType ind = 0; ind != base_space_dimensions; ++ind) {
                basis_vectors_str_[ind] = "e" + std::to_string(ind + 1);
            }
            basis_vectors_str_[base_space_dimensions] = "no";
            basis_vectors_str_[base_space_dimensions + 1] = "ni";

            auto aux = 1 / sqrt(ScalarType(2));
            from_actual_to_signed_metric_ = detail::make_identity_matrix<ScalarType, DimensionsAtCompileTime>(base_space_dimensions + 2);
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions, base_space_dimensions) = aux;
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions + 1, base_space_dimensions) = aux;
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions, base_space_dimensions + 1) = -aux;
            detail::coeff(from_actual_to_signed_metric_, base_space_dimensions + 1, base_space_dimensions + 1) = aux;

            from_signed_to_actual_metric_ = detail::transpose(from_actual_to_signed_metric_);
        }

        std::vector<std::string> basis_vectors_str_;
        detail::matrix_type_t<ScalarType, DimensionsAtCompileTime, DimensionsAtCompileTime> from_actual_to_signed_metric_;
        detail::matrix_type_t<ScalarType, DimensionsAtCompileTime, DimensionsAtCompileTime> from_signed_to_actual_metric_;

        template<DefaultIndexType SomeBaseSpaceDimensionsAtCompileTime> friend struct detail::from_actual_to_signed_metric_impl;
        template<DefaultIndexType SomeBaseSpaceDimensionsAtCompileTime> friend struct detail::from_signed_to_actual_metric_impl;
    };

    template<DefaultIndexType BaseSpaceDimensionsAtCompileTime>
    struct is_metric_space<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> > :
        std::true_type {
    };

    namespace detail {

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime>
        struct apply_signed_metric_impl<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> const &, MatrixType const &factors_in_signed_metric) noexcept {
                using IndexType = index_type_t<MatrixType>;
                MatrixType result = factors_in_signed_metric;
                for (IndexType col = 0, col_end = cols(result), last_row = rows(result) - 1; col != col_end; ++col) {
                    coeff(result, last_row, col) *= -1;
                }
                return result;
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime>
        struct from_actual_to_signed_metric_impl<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> const &space, MatrixType const &factors_in_actual_metric) noexcept {
                return detail::prod(space.from_actual_to_signed_metric_, factors_in_actual_metric);
            }
        };

        template<DefaultIndexType BaseSpaceDimensionsAtCompileTime>
        struct from_signed_to_actual_metric_impl<ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> > {
            template<typename MatrixType>
            constexpr static decltype(auto) eval(ConformalMetricSpace<BaseSpaceDimensionsAtCompileTime> const &space, MatrixType const &factors_in_signed_metric) noexcept {
                return detail::prod(space.from_signed_to_actual_metric_, factors_in_signed_metric);
            }
        };

    }

}

#endif // __TBGAL_CONFORMAL_METRIC_SPACE_HPP__
