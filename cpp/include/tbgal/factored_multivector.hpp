#ifndef __TBGAL_MULTIVECTOR_HPP__
#define __TBGAL_MULTIVECTOR_HPP__

namespace tbgal {

    template<typename FactoringProductType, typename SquareMatrixType>
    class FactoredMultivector {
    public:

        static_assert(detail::is_factoring_product_v<FactoringProductType>, "Invalid FactoringProductType.");

        using ScalarType = detail::scalar_type_t<SquareMatrixType>;
        using IndexType = detail::index_type_t<SquareMatrixType>;

        using SpaceType = typename FactoringProductType::SpaceType;

        constexpr FactoredMultivector(FactoredMultivector const &) = default;
        constexpr FactoredMultivector(FactoredMultivector &&) = default;

        constexpr FactoredMultivector(SpaceType const &space) noexcept :
            scalar_(0),
            factors_(detail::make_identity_matrix<ScalarType>(space)),
            factors_count_(0),
            space_(space) {
        }

        template<typename OtherScalarType>
        constexpr FactoredMultivector(OtherScalarType &&scalar, SpaceType const &space) noexcept :
            scalar_(std::move(scalar)),
            factors_(detail::make_identity_matrix<ScalarType>(space)),
            factors_count_(0),
            space_(space) {
        }

        template<typename OtherFactoringProductType, typename OtherSquareMatrixType>
        constexpr FactoredMultivector(FactoredMultivector<OtherFactoringProductType, OtherSquareMatrixType> const &other) noexcept;
        //TODO Especializar

        constexpr ScalarType scalar() const noexcept {
            return scalar_;
        }

        constexpr SquareMatrixType factors() const noexcept {
            return factors_;
        }

        constexpr IndexType factors_count() const noexcept {
            return factors_count_;
        }

        constexpr SpaceType& space() const noexcept {
            return space_;
        }

    private:

        constexpr FactoredMultivector(ScalarType &&scalar, SquareMatrixType &&factors, IndexType &&factors_count, SpaceType const &space) noexcept :
            scalar_(std::move(scalar)),
            factors_(std::move(factors)),
            factors_count_(std::move(factors_count)),
            space_(space) {
        }

    private:

        ScalarType scalar_;

        SquareMatrixType factors_;
        IndexType factors_count_;

        SpaceType &space_;

        template<typename ScalarType, typename MetricSpaceType, typename> friend decltype(auto) scalar(ScalarType const &, MetricSpaceType const &);

        template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType> friend decltype(auto) ADD(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &);
        template<typename SomeMetricSpaceType, typename SomeSquareMatrixType> friend decltype(auto) DUAL(FactoredMultivector<OuterProduct<SomeMetricSpaceType>, SomeSquareMatrixType> const &);
        template<typename FirstType, typename SecondFactoringProduct, typename SecondSquareMatrixType, typename> friend decltype(auto) LCONT(FirstType const &, FactoredMultivector<SecondFactoringProduct, SecondSquareMatrixType> const &);
        template<typename FirstFactoringProductType, typename FirstSquareMatrixType, typename SecondFactoringProductType, typename SecondSquareMatrixType> friend decltype(auto) SUB(FactoredMultivector<FirstFactoringProductType, FirstSquareMatrixType> const &, FactoredMultivector<SecondFactoringProductType, SecondSquareMatrixType> const &);
        template<typename SomeFactoringProductType, typename SomeSquareMatrixType> friend FactoredMultivector<SomeFactoringProductType, SomeSquareMatrixType> UNARY_MINUS(FactoredMultivector<SomeFactoringProductType, SomeSquareMatrixType> const &);
    };

    namespace detail {

        template<typename FactoringProductType, typename SquareMatrixType>
        struct is_multivector<FactoredMultivector<FactoringProductType, SquareMatrixType> > :
            std::true_type {
        };

    }

}

#endif // __TBGAL_MULTIVECTOR_HPP__
