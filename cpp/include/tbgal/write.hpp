#ifndef __TBGAL_WRITE_HPP__
#define __TBGAL_WRITE_HPP__

namespace tbgal {

    namespace detail {

        template<typename ScalarType, typename FactoringProductType>
        std::ostream & write(std::ostream &os, FactoredMultivector<ScalarType, FactoringProductType> const &arg) noexcept {
            using IndexType = typename FactoredMultivector<ScalarType, FactoringProductType>::IndexType;
            os << arg.scalar();
            if (arg.factors_count() > 0) {
                auto factors = arg.factors_in_actual_metric();
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

    template<typename ScalarType, typename MetricSpaceType>
    std::ostream & operator <<(std::ostream &os, FactoredMultivector<ScalarType, GeometricProduct<MetricSpaceType>> const &arg) noexcept {
        os << "gp(";
        detail::write(os, arg);
        os << ")";
        return os;
    }

    template<typename ScalarType, typename MetricSpaceType>
    std::ostream & operator <<(std::ostream &os, FactoredMultivector<ScalarType, OuterProduct<MetricSpaceType>> const &arg) noexcept {
        os << "op(";
        detail::write(os, arg);
        os << ")";
        return os;
    }

}

#endif // __TBGAL_WRITE_HPP__
