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

#ifndef __TBGAL_PYTHON_MACRO_HPP__
#define __TBGAL_PYTHON_MACRO_HPP__

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
    py::def("add", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return add(lhs, rhs); }); \
    py::def("addition", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return addition(lhs, rhs); }); \
    py::def("sub", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return sub(lhs, rhs); }); \
    py::def("subtraction", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return subtraction(lhs, rhs); }); \
    py::def("gp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return gp(lhs, rhs); }); \
    py::def("igp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return igp(lhs, rhs); }); \
    py::def("dot", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return dot(lhs, rhs); }); \
    py::def("hip", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return hip(lhs, rhs); }); \
    py::def("lcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return lcont(lhs, rhs); }); \
    py::def("rcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return rcont(lhs, rhs); }); \
    py::def("sp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return sp(lhs, rhs); }); \
    py::def("op", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return op(lhs, rhs); }); \
    py::def("apply_even_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return apply_even_versor(lhs, rhs); }); \
    py::def("apply_odd_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return apply_odd_versor(lhs, rhs); }); \
    py::def("apply_rotor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return apply_rotor(lhs, rhs); })

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_MULTIVECTOR_AND_NATIVE(FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE, SECOND_SCALAR_TYPE) \
    py::def("add", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return add(lhs, DefaultScalarType(rhs)); }); \
    py::def("addition", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return addition(lhs, DefaultScalarType(rhs)); }); \
    py::def("sub", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return sub(lhs, DefaultScalarType(rhs)); }); \
    py::def("subtraction", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return subtraction(lhs, DefaultScalarType(rhs)); }); \
    py::def("gp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return gp(lhs, DefaultScalarType(rhs)); }); \
    py::def("igp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return igp(lhs, DefaultScalarType(rhs)); }); \
    py::def("dot", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return dot(lhs, DefaultScalarType(rhs)); }); \
    py::def("hip", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return hip(lhs, DefaultScalarType(rhs)); }); \
    py::def("lcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return lcont(lhs, DefaultScalarType(rhs)); }); \
    py::def("rcont", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return rcont(lhs, DefaultScalarType(rhs)); }); \
    py::def("sp", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return sp(lhs, DefaultScalarType(rhs)); }); \
    py::def("op", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return op(lhs, DefaultScalarType(rhs)); }); \
    py::def("apply_even_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return apply_even_versor(lhs, DefaultScalarType(rhs)); }); \
    py::def("apply_odd_versor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return apply_odd_versor(lhs, DefaultScalarType(rhs)); }); \
    py::def("apply_rotor", +[](FactoredMultivector<FIRST_SCALAR_TYPE, FIRST_FACTORING_PRODUCT_TYPE> const &lhs, SECOND_SCALAR_TYPE rhs) { return apply_rotor(lhs, DefaultScalarType(rhs)); })

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_MULTIVECTOR(FIRST_SCALAR_TYPE, SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE) \
    py::def("add", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return add(DefaultScalarType(lhs), rhs); }); \
    py::def("addition", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return addition(DefaultScalarType(lhs), rhs); }); \
    py::def("sub", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return sub(DefaultScalarType(lhs), rhs); }); \
    py::def("subtraction", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return subtraction(DefaultScalarType(lhs), rhs); }); \
    py::def("gp", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return gp(DefaultScalarType(lhs), rhs); }); \
    py::def("igp", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return igp(DefaultScalarType(lhs), rhs); }); \
    py::def("dot", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return dot(DefaultScalarType(lhs), rhs); }); \
    py::def("hip", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return hip(DefaultScalarType(lhs), rhs); }); \
    py::def("lcont", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return lcont(DefaultScalarType(lhs), rhs); }); \
    py::def("rcont", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return rcont(DefaultScalarType(lhs), rhs); }); \
    py::def("sp", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return sp(DefaultScalarType(lhs), rhs); }); \
    py::def("op", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return op(DefaultScalarType(lhs), rhs); }); \
    py::def("apply_even_versor", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return apply_even_versor(DefaultScalarType(lhs), rhs); }); \
    py::def("apply_odd_versor", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return apply_odd_versor(DefaultScalarType(lhs), rhs); }); \
    py::def("apply_rotor", +[](FIRST_SCALAR_TYPE lhs, FactoredMultivector<SECOND_SCALAR_TYPE, SECOND_FACTORING_PRODUCT_TYPE> const &rhs) { return apply_rotor(DefaultScalarType(lhs), rhs); })

#define _PY_TBGAL_EXPOSE_BINARY_OPERATIONS_WITH_NATIVE_AND_NATIVE(FIRST_SCALAR_TYPE, SECOND_SCALAR_TYPE) \
    py::def("add", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return add(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("addition", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return addition(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("sub", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return sub(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("subtraction", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return subtraction(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("gp", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return gp(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("igp", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return igp(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("dot", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return dot(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("hip", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return hip(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("lcont", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return lcont(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("rcont", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return rcont(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("sp", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return sp(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("op", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return op(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("apply_even_versor", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return apply_even_versor(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("apply_odd_versor", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return apply_odd_versor(DefaultScalarType(lhs), DefaultScalarType(rhs)); }); \
    py::def("apply_rotor", +[](FIRST_SCALAR_TYPE lhs, SECOND_SCALAR_TYPE rhs) { return apply_rotor(DefaultScalarType(lhs), DefaultScalarType(rhs)); })

#define _PY_TBGAL_EXPOSE_UNARY_OPERATIONS(METRIC_SPACE_TYPE) \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_MULTIVECTOR(DefaultScalarType, GeometricProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_MULTIVECTOR(DefaultScalarType, OuterProduct<METRIC_SPACE_TYPE>); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(std::int16_t); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(std::int32_t); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(std::int64_t); \
    _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(DefaultScalarType)

#define _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_MULTIVECTOR(SCALAR_TYPE, FACTORING_PRODUCT_TYPE) \
    py::def("is_blade", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return is_blade(arg); }); \
    py::def("is_zero", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return is_zero(arg); }); \
    py::def("unary_minus", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return unary_minus(arg); }); \
    py::def("unary_plus", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return unary_plus(arg); }); \
    py::def("conjugate", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return conjugate(arg); }); \
    py::def("involute", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return involute(arg); }); \
    py::def("reverse", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return reverse(arg); }); \
    py::def("dual", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return dual(arg); }); \
    py::def("undual", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return undual(arg); }); \
    py::def("inverse", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return inverse(arg); }); \
    py::def("inv", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return inv(arg); }); \
    py::def("unit", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return unit(arg); }); \
    py::def("rnorm_sqr", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return rnorm_sqr(arg); }); \
    py::def("rnorm", +[](FactoredMultivector<SCALAR_TYPE, FACTORING_PRODUCT_TYPE> const &arg) { return rnorm(arg); })

#define _PY_TBGAL_EXPOSE_UNARY_OPERATIONS_WITH_NATIVE(SCALAR_TYPE) \
    py::def("is_blade", +[](SCALAR_TYPE arg) { return is_blade(DefaultScalarType(arg)); }); \
    py::def("is_zero", +[](SCALAR_TYPE arg) { return is_zero(DefaultScalarType(arg)); }); \
    py::def("unary_minus", +[](SCALAR_TYPE arg) { return unary_minus(DefaultScalarType(arg)); }); \
    py::def("unary_plus", +[](SCALAR_TYPE arg) { return unary_plus(DefaultScalarType(arg)); }); \
    py::def("conjugate", +[](SCALAR_TYPE arg) { return conjugate(DefaultScalarType(arg)); }); \
    py::def("involute", +[](SCALAR_TYPE arg) { return involute(DefaultScalarType(arg)); }); \
    py::def("reverse", +[](SCALAR_TYPE arg) { return reverse(DefaultScalarType(arg)); }); \
    py::def("dual", +[](SCALAR_TYPE arg) { return dual(DefaultScalarType(arg)); }); \
    py::def("undual", +[](SCALAR_TYPE arg) { return undual(DefaultScalarType(arg)); }); \
    py::def("inverse", +[](SCALAR_TYPE arg) { return inverse(DefaultScalarType(arg)); }); \
    py::def("inv", +[](SCALAR_TYPE arg) { return inv(DefaultScalarType(arg)); }); \
    py::def("unit", +[](SCALAR_TYPE arg) { return unit(DefaultScalarType(arg)); }); \
    py::def("rnorm_sqr", +[](SCALAR_TYPE arg) { return rnorm_sqr(DefaultScalarType(arg)); }); \
    py::def("rnorm", +[](SCALAR_TYPE arg) { return rnorm(DefaultScalarType(arg)); })

#define PY_TBGAL_EXPOSE_UTILS() \
    py::def("e", +[](DefaultIndexType index) { return e(index); }); \
    py::def("scalar", +[](std::int16_t value) { return scalar(DefaultScalarType(value)); }); \
    py::def("scalar", +[](std::int32_t value) { return scalar(DefaultScalarType(value)); }); \
    py::def("scalar", +[](std::int64_t value) { return scalar(DefaultScalarType(value)); }); \
    py::def("scalar", +[](DefaultScalarType value) { return scalar(value); }); \
    py::def("vector", py::raw_function(+[](py::tuple const &args, py::dict const &) { \
        std::vector<DefaultScalarType> dummy; \
        dummy.assign(py::stl_input_iterator<DefaultScalarType>(args), py::stl_input_iterator<DefaultScalarType>()); \
        return vector(dummy.begin(), dummy.end()); \
    }));

#define PY_TBGAL_EXPOSE_FUNCTION(NAME, FUNCTION) \
    py::def(NAME, &FUNCTION)

#define PY_TBGAL_EXPOSE_VARIABLE(NAME, INSTANCE) \
    py::scope().attr(NAME) = py::ptr(&INSTANCE)

#endif // __TBGAL_PYTHON_MACRO_HPP__
