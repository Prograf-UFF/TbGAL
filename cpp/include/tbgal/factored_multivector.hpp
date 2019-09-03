#ifndef __TBGAL_MULTIVECTOR_HPP__
#define __TBGAL_MULTIVECTOR_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    class FactoredMultivector final {
    public:

        static_assert(detail::is_factoring_product_v<FactoringProductType>, "Invalid FactoringProductType.");

        using ScalarType = detail::scalar_type_t<SquareMatrixType>;
        using IndexType = detail::index_type_t<SquareMatrixType>;

        using SpaceType = typename FactoringProductType::SpaceType;

        constexpr FactoredMultivector(FactoredMultivector const &) = default;
        constexpr FactoredMultivector(FactoredMultivector &&) = default;

        constexpr FactoredMultivector(SpaceType const &space) noexcept :
            space_(space),
            scalar_(0),
            factors_(detail::make_identity_matrix<ScalarType, SpaceType::DimensionsAtCompileTime>(space.dimensions())),
            factors_count_(0) {
        }

        template<typename OtherScalarType>
        constexpr FactoredMultivector(SpaceType const &space, OtherScalarType &&scalar) noexcept :
            space_(space),
            scalar_(std::move(scalar)),
            factors_(detail::make_identity_matrix<ScalarType, SpaceType::DimensionsAtCompileTime>(space.dimensions())),
            factors_count_(0) {
        }

        template<typename OtherFactoringProductType, typename OtherSquareMatrixType>
        constexpr FactoredMultivector(FactoredMultivector<OtherFactoringProductType, OtherSquareMatrixType> const &other) noexcept;
        //TODO [FUTURE] Specialization.

        constexpr SpaceType& space() const noexcept {
            return space_;
        }

        constexpr ScalarType scalar() const noexcept {
            return scalar_;
        }

        constexpr SquareMatrixType factors() const noexcept {
            return factors_;
        }

        constexpr IndexType factors_count() const noexcept {
            return factors_count_;
        }

    private:

        constexpr FactoredMultivector(SpaceType const &space, ScalarType &&scalar, SquareMatrixType &&factors, IndexType &&factors_count) noexcept :
            space_(space),
            scalar_(std::move(scalar)),
            factors_(std::move(factors)),
            factors_count_(std::move(factors_count)) {
        }

    private:

        ScalarType scalar_;

        SquareMatrixType factors_;
        IndexType factors_count_;

        SpaceType &space_;

        template<typename SomeMetricSpaceType, typename SomeScalarType, typename> friend constexpr decltype(auto) scalar(SomeMetricSpaceType const &, SomeScalarType const &);
        template<typename SomeMetricSpaceType, typename FirstScalarType, typename... NextScalarTypes> friend constexpr decltype(auto) vector(SomeMetricSpaceType const &, FirstScalarType const &arg1, NextScalarTypes const &...);

        template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType> friend constexpr decltype(auto) ADD(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &);
        template<typename SomeFactoringProductType, typename SomeSquareMatrixType> friend constexpr decltype(auto) DUAL(FactoredMultivector<SomeFactoringProductType, SomeSquareMatrixType> const &);
        template<typename FirstScalarType, typename SecondFactoringProduct, typename SecondSquareMatrixType, typename> friend constexpr decltype(auto) LCONT(FirstScalarType const &, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &);
        template<typename SomeFactoringProductType, typename SomeSquareMatrixType> friend constexpr decltype(auto) REVERSE(FactoredMultivector<SomeFactoringProductType, SomeSquareMatrixType> const &);
        template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType> friend constexpr decltype(auto) SUB(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &);
        template<typename SomeFactoringProductType, typename SomeSquareMatrixType> friend constexpr FactoredMultivector<SomeFactoringProductType, SomeSquareMatrixType> UNARY_MINUS(FactoredMultivector<SomeFactoringProductType, SomeSquareMatrixType> const &);
    };

    namespace detail {

        template<typename FactoringProductType, typename SquareMatrixType>
        struct is_multivector<FactoredMultivector<FactoringProductType, SquareMatrixType> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_MULTIVECTOR_HPP__
