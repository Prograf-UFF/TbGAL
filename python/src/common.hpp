/* Copyright (C) Eduardo Vera Sousa and Leandro Augusto Frata Fernandes
 *
 * authors    : Sousa, Eduardo V.
 *              Fernandes, Leandro A. F.
 * repository : https://github.com/Prograf-UFF/TbGAL
 *
 * This file is part of the Tensor-based Geometric Algebra Library (TbGAL).
 *
 * TbGAL is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * TbGAL is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU General Public License for more details.
 * You should have received a copy of the GNU General Public License
 * along with TbGAL. If not, see <https://www.gnu.org/licenses/>.
 */

#ifndef __TBGAL_PYTHON_COMMON_HPP__
#define __TBGAL_PYTHON_COMMON_HPP__

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>

namespace py = boost::python;
namespace np = boost::python::numpy;

namespace tbgal {

    namespace python {

        template<typename ValueType>
        class sequence_const_iterator {
        public:

            using difference_type = py::ssize_t;
            using value_type = ValueType;
            using pointer = ValueType*;
            using reference = ValueType&;
            using iterator_category = std::random_access_iterator_tag;

            inline sequence_const_iterator(py::object const *object_ptr, difference_type index) :
                object_ptr_(object_ptr),
                index_(index) {
            }

            inline sequence_const_iterator(sequence_const_iterator const &other) = default;
            inline sequence_const_iterator(sequence_const_iterator &&other) = default;

            inline sequence_const_iterator& operator++() {
                ++index_;
                return *this;
            }

            inline sequence_const_iterator& operator--() {
                --index_;
                return *this;
            }

            inline sequence_const_iterator operator++(int) {
                sequence_const_iterator retval = *this;
                ++(*this);
                return retval;
            }

            inline sequence_const_iterator operator--(int) {
                sequence_const_iterator retval = *this;
                --(*this);
                return retval;
            }

            inline bool operator==(sequence_const_iterator const &other) const {
                assert(object_ptr_ == other.object_ptr_);
                return index_ == other.index_;
            }

            inline bool operator!=(sequence_const_iterator const &other) const {
                assert(object_ptr_ == other.object_ptr_);
                return index_ != other.index_;
            }

            inline value_type operator*() const {
                return py::extract<value_type>(object_ptr_->operator[](index_))();
            }

            inline difference_type operator-(sequence_const_iterator const &other) const {
                return index_ - other.index_;
            }

        private:

            py::object const *object_ptr_;
            difference_type index_;

            friend class boost::iterator_core_access;
        };

        template<typename ValueType>
        sequence_const_iterator<ValueType> begin(py::object const *object_ptr) {
            return sequence_const_iterator<ValueType>(object_ptr, 0);
        }

        template<typename ValueType>
        sequence_const_iterator<ValueType> end(py::object const *object_ptr) {
            return sequence_const_iterator<ValueType>(object_ptr, py::len(*object_ptr));
        }

    }

}

#define PY_TBGAL_INITIALIZE() \
    np::initialize(); \
    py::register_exception_translator<NotSupportedError>(+[](NotSupportedError const &error) -> void { PyErr_SetString(PyExc_RuntimeError, error.what()); })

#define PY_TBGAL_EXPOSE_CORE(METRIC_SPACE_TYPE) \
    _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_CLASSES(METRIC_SPACE_TYPE); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS(METRIC_SPACE_TYPE); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS(METRIC_SPACE_TYPE)

#define _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_CLASSES(METRIC_SPACE_TYPE) \
    _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_CLASS("FactoredMultivector_GeometricProduct", DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>) \
    _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_MULTIVECTOR(DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>); \
    \
    _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_CLASS("FactoredMultivector_OuterProduct", DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>) \
    _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_MULTIVECTOR(DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>)

#define _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_CLASS(NAME, SCALAR_TYPE, FACTORING_PRODUCT_TYPE) \
    py::class_<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >(NAME) \
        .def("scalar", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { \
            return DefaultScalarType(arg.scalar()); \
        }) \
        .def("factors", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { \
            using IndexType = FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>::IndexType; \
            auto factors = arg.factors_in_actual_metric(); \
            IndexType const rows = detail::rows(factors), cols = detail::cols(factors); \
            np::ndarray result = np::empty(py::make_tuple(rows, cols), np::dtype::get_builtin<SCALAR_TYPE>()); \
            SCALAR_TYPE *data = (SCALAR_TYPE*)result.get_data(); \
            for (IndexType row = 0; row != rows; ++row) { \
                for (IndexType col = 0; col != cols; ++col, ++data) { \
                    *data = detail::coeff(factors, row, col); \
                } \
            } \
            return result; \
        }) \
        .def("factors_count", &FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE>::factors_count) \
        .def("__repr__", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { \
            std::ostringstream os; \
            os << arg; \
            return os.str(); \
        }) \
        .def(+py::self) \
        .def(-py::self) \
        .def(~py::self) \
        .def(py::self + py::self) \
        .def(py::self - py::self) \
        .def(py::self ^ py::self) \
        .def(py::self * py::self) \
        .def(py::self / py::self) \
        _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_NATIVE(std::int16_t) \
        _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_NATIVE(std::int32_t) \
        _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_NATIVE(std::int64_t) \
        _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_NATIVE(DefaultScalarType)

#define _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_MULTIVECTOR(SCALAR_TYPE, FACTORING_PRODUCT_TYPE) \
    .def(py::self + py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >()) \
    .def(py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >() + py::self) \
    .def(py::self - py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >()) \
    .def(py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >() - py::self) \
    .def(py::self ^ py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >()) \
    .def(py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >() ^ py::self) \
    .def(py::self * py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >()) \
    .def(py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >() * py::self) \
    .def(py::self / py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >()) \
    .def(py::other<FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> >() / py::self)

#define _PY_TBGAL_EXPOSE_FACTORED_MULTIVECTOR_SELF_OPERATIONS_WITH_NATIVE(SCALAR_TYPE) \
    .def(py::self + py::other<SCALAR_TYPE>()) \
    .def(py::other<SCALAR_TYPE>() + py::self) \
    .def(py::self - py::other<SCALAR_TYPE>()) \
    .def(py::other<SCALAR_TYPE>() - py::self) \
    .def(py::self ^ py::other<SCALAR_TYPE>()) \
    .def(py::other<SCALAR_TYPE>() ^ py::self) \
    .def(py::self * py::other<SCALAR_TYPE>()) \
    .def(py::other<SCALAR_TYPE>() * py::self) \
    .def(py::self / py::other<SCALAR_TYPE>()) \
    .def(py::other<SCALAR_TYPE>() / py::self)

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS(METRIC_SPACE_TYPE) \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_SOMETHING(METRIC_SPACE_TYPE, DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_SOMETHING(METRIC_SPACE_TYPE, DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_SOMETHING(METRIC_SPACE_TYPE, std::int16_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_SOMETHING(METRIC_SPACE_TYPE, std::int32_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_SOMETHING(METRIC_SPACE_TYPE, std::int64_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_SOMETHING(METRIC_SPACE_TYPE, DefaultScalarType)

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_SOMETHING(METRIC_SPACE_TYPE, FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE) \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_NATIVE(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, std::int16_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_NATIVE(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, std::int32_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_NATIVE(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, std::int64_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_NATIVE(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, DefaultScalarType); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>)

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_SOMETHING(METRIC_SPACE_TYPE, FIRST_SCALAR_TYPE) \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_NATIVE(FIRST_SCALAR_TYPE, std::int16_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_NATIVE(FIRST_SCALAR_TYPE, std::int32_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_NATIVE(FIRST_SCALAR_TYPE, std::int64_t); \
    _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_NATIVE(FIRST_SCALAR_TYPE, DefaultScalarType)

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE) \
    py::def("add", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return add(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("addition", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return addition(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("sub", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return sub(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("subtraction", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return subtraction(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("gp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return gp(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("igp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return igp(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("dot", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return dot(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("hip", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return hip(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("lcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return lcont(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("rcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return rcont(arg1, arg2); }), py::args("arg1", "arg2"); \
    py::def("sp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return sp(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("op", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return op(arg1, arg2); }, py::args("arg1", "arg2")); \
    py::def("apply_even_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &versor, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg) { return apply_even_versor(versor, arg); }, py::args("versor", "arg")); \
    py::def("apply_odd_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &versor, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg) { return apply_odd_versor(versor, arg); }, py::args("versor", "arg")); \
    py::def("apply_rotor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &rotor, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg) { return apply_rotor(rotor, arg); }, py::args("rotor", "arg"))

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_NATIVE(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, SECOND_SCALAR_TYPE) \
    py::def("add", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return add(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("addition", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return addition(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("sub", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return sub(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("subtraction", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return subtraction(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("gp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return gp(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("igp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return igp(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("dot", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return dot(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("hip", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return hip(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("lcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return lcont(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("rcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return rcont(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("sp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return sp(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("op", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &arg1, SECOND_SCALAR_TYPE arg2) { return op(arg1, DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("apply_even_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &versor, SECOND_SCALAR_TYPE arg) { return apply_even_versor(versor, DefaultScalarType(arg)); }, py::args("versor", "arg")); \
    py::def("apply_odd_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &versor, SECOND_SCALAR_TYPE arg) { return apply_odd_versor(versor, DefaultScalarType(arg)); }, py::args("versor", "arg")); \
    py::def("apply_rotor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &rotor, SECOND_SCALAR_TYPE arg) { return apply_rotor(rotor, DefaultScalarType(arg)); }, py::args("rotor", "arg"))

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE) \
    py::def("add", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return add(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("addition", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return addition(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("sub", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return sub(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("subtraction", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return subtraction(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("gp", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return gp(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("igp", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return igp(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("dot", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return dot(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("hip", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return hip(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("lcont", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return lcont(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("rcont", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return rcont(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("sp", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return sp(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("op", +[](FIRST_SCALAR_TYPE arg1, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg2) { return op(DefaultScalarType(arg1), arg2); }, py::args("arg1", "arg2")); \
    py::def("apply_even_versor", +[](FIRST_SCALAR_TYPE versor, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg) { return apply_even_versor(DefaultScalarType(versor), arg); }, py::args("versor", "arg")); \
    py::def("apply_odd_versor", +[](FIRST_SCALAR_TYPE versor, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg) { return apply_odd_versor(DefaultScalarType(versor), arg); }, py::args("versor", "arg")); \
    py::def("apply_rotor", +[](FIRST_SCALAR_TYPE rotor, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &arg) { return apply_rotor(DefaultScalarType(rotor), arg); }, py::args("rotor", "arg"))

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_NATIVE(FIRST_SCALAR_TYPE, SECOND_SCALAR_TYPE) \
    py::def("add", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return add(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("addition", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return addition(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("sub", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return sub(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("subtraction", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return subtraction(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("gp", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return gp(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("igp", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return igp(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("dot", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return dot(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("hip", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return hip(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("lcont", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return lcont(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("rcont", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return rcont(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("sp", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return sp(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("op", +[](FIRST_SCALAR_TYPE arg1, SECOND_SCALAR_TYPE arg2) { return op(DefaultScalarType(arg1), DefaultScalarType(arg2)); }, py::args("arg1", "arg2")); \
    py::def("apply_even_versor", +[](FIRST_SCALAR_TYPE versor, SECOND_SCALAR_TYPE arg) { return apply_even_versor(DefaultScalarType(versor), DefaultScalarType(arg)); }, py::args("versor", "arg")); \
    py::def("apply_odd_versor", +[](FIRST_SCALAR_TYPE versor, SECOND_SCALAR_TYPE arg) { return apply_odd_versor(DefaultScalarType(versor), DefaultScalarType(arg)); }, py::args("versor", "arg")); \
    py::def("apply_rotor", +[](FIRST_SCALAR_TYPE rotor, SECOND_SCALAR_TYPE arg) { return apply_rotor(DefaultScalarType(rotor), DefaultScalarType(arg)); }, py::args("rotor", "arg"))

#define _PY_TBGAL_EXPOSE_UNARY_OPERATIONS(METRIC_SPACE_TYPE) \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_MULTIVECTOR(DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_MULTIVECTOR(DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(std::int16_t); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(std::int32_t); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(std::int64_t); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(DefaultScalarType)

#define _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_MULTIVECTOR(SCALAR_TYPE, FACTORING_PRODUCT_TYPE) \
    py::def("is_blade", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return is_blade(arg); }, py::args("arg")); \
    py::def("is_zero", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return is_zero(arg); }, py::args("arg")); \
    py::def("unary_minus", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return unary_minus(arg); }, py::args("arg")); \
    py::def("unary_plus", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return unary_plus(arg); }, py::args("arg")); \
    py::def("conjugate", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return conjugate(arg); }, py::args("arg")); \
    py::def("involute", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return involute(arg); }, py::args("arg")); \
    py::def("reverse", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return reverse(arg); }, py::args("arg")); \
    py::def("dual", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return dual(arg); }, py::args("arg")); \
    py::def("undual", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return undual(arg); }, py::args("arg")); \
    py::def("inverse", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return inverse(arg); }, py::args("arg")); \
    py::def("inv", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return inv(arg); }, py::args("arg")); \
    py::def("unit", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return unit(arg); }, py::args("arg")); \
    py::def("rnorm_sqr", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return rnorm_sqr(arg); }, py::args("arg")); \
    py::def("rnorm", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return rnorm(arg); }, py::args("arg"))

#define _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(SCALAR_TYPE) \
    py::def("is_blade", +[](SCALAR_TYPE arg) { return is_blade(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("is_zero", +[](SCALAR_TYPE arg) { return is_zero(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("unary_minus", +[](SCALAR_TYPE arg) { return unary_minus(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("unary_plus", +[](SCALAR_TYPE arg) { return unary_plus(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("conjugate", +[](SCALAR_TYPE arg) { return conjugate(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("involute", +[](SCALAR_TYPE arg) { return involute(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("reverse", +[](SCALAR_TYPE arg) { return reverse(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("dual", +[](SCALAR_TYPE arg) { return dual(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("undual", +[](SCALAR_TYPE arg) { return undual(DefaultScalarType(arg)); }), py::args("arg"); \
    py::def("inverse", +[](SCALAR_TYPE arg) { return inverse(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("inv", +[](SCALAR_TYPE arg) { return inv(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("unit", +[](SCALAR_TYPE arg) { return unit(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("rnorm_sqr", +[](SCALAR_TYPE arg) { return rnorm_sqr(DefaultScalarType(arg)); }, py::args("arg")); \
    py::def("rnorm", +[](SCALAR_TYPE arg) { return rnorm(DefaultScalarType(arg)); }, py::args("arg"))

#define PY_TBGAL_EXPOSE_UTILS() \
    py::def("e", +[](DefaultIndexType index) { return e(index); }, py::args("index")); \
    py::def("scalar", +[](std::int16_t value) { return scalar(DefaultScalarType(value)); }, py::args("value")); \
    py::def("scalar", +[](std::int32_t value) { return scalar(DefaultScalarType(value)); }, py::args("value")); \
    py::def("scalar", +[](std::int64_t value) { return scalar(DefaultScalarType(value)); }, py::args("value")); \
    py::def("scalar", +[](DefaultScalarType value) { return scalar(value); }, py::args("value")); \
    py::def("vector", py::raw_function(+[](py::tuple const &args, py::dict const &) { return vector(tbgal::python::begin<DefaultScalarType>(&args), tbgal::python::end<DefaultScalarType>(&args)); }));

#define PY_TBGAL_EXPOSE_FUNCTION(NAME, FUNCTION) \
    py::def(NAME, &FUNCTION)

#define PY_TBGAL_EXPOSE_GLOBAL_VARIABLE(NAME, INSTANCE) \
    py::scope().attr(NAME) = py::ptr(&INSTANCE)

#endif // __TBGAL_PYTHON_COMMON_HPP__
