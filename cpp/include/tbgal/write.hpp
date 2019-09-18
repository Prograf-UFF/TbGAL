#ifndef __TBGAL_WRITE_HPP__
#define __TBGAL_WRITE_HPP__

namespace tbgal {

    namespace detail {

        template<typename FactoringProductType, typename SquareMatrixType>
        std::ostream & write(std::ostream &os, FactoredMultivector<FactoringProductType, SquareMatrixType> const &arg) noexcept {
            using IndexType = typename FactoredMultivector<FactoringProductType, SquareMatrixType>::IndexType;
            os << arg.scalar();
            if (arg.factors_count() > 0) {
                auto factors = from_orthogonal_to_actual_metric(arg.space(), arg.factors());
                for (IndexType ind = 0; ind != arg.factors_count(); ++ind) {
                    os << ", (";
                    for (IndexType dim = 0; dim != arg.space().dimensions(); ++dim) {
                        auto value = coeff(factors, dim, ind);
                        if (dim > 0) {
                            os << (value >= 0 ? " + " : " - ") << std::abs(value);
                        }
                        else {
                            os << value;
                        }
                        os << " * " << arg.space().basis_vector_str(dim);
                    }
                    os << ")";
                }
            }
            return os;
        }

    }

    template<typename MetricSpaceType, typename SquareMatrixType>
    std::ostream & operator <<(std::ostream &os, FactoredMultivector<GeometricProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        os << "GP(";
        detail::write(os, arg);
        os << ")";
        return os;
    }

    template<typename MetricSpaceType, typename SquareMatrixType>
    std::ostream & operator <<(std::ostream &os, FactoredMultivector<OuterProduct<MetricSpaceType>, SquareMatrixType> const &arg) noexcept {
        os << "OP(";
        detail::write(os, arg);
        os << ")";
        return os;
    }

}

#endif // __TBGAL_WRITE_HPP__
