#include <boost/python/module.hpp>
#include <boost/python/def.hpp>

#include <boost/python.hpp>
#include <boost/python/object.hpp>
#include <boost/python/args.hpp>
#include <boost/python/list.hpp>
#include <boost/python/tuple.hpp>
#include <boost/python/dict.hpp>
#include "boost/python/stl_iterator.hpp"
#include <boost/python/raw_function.hpp>

using namespace tbgal;
namespace python = boost::python;

/*

Hestenes' inner product         - ok
Dot product                     - ok
left contraction                - ok
right contraction               - ok
Reversion                       - ok
Inversion                       - ok
Reverse norm                    - ok
Squared reverse norm            - ok
Dualization                     - ok
Undualization                   - ok
Scalar product                  - ok
Outer product                   - ok
Geometric product               - ok
Addition                        - ok?
Subtraction                     - ok?
Unary plus                      - ok
Unary minus                     - ok

*/

template<typename T>
inline
std::vector< T > py_list_to_std_vector( const boost::python::object& iterable )
{
    return std::vector< T >( boost::python::stl_input_iterator< T >( iterable ),
                             boost::python::stl_input_iterator< T >( ) );
}

template<typename T>
T py_hip(const T &lhs, const T &rhs) {
    return tbgal::hip(lhs, rhs);
}

template<typename T>
T py_dot(const T &lhs, const T &rhs) {
    return tbgal::dot(lhs, rhs);
}

template<typename T>
T py_lcont(const T &lhs, const T &rhs) {
    return tbgal::lcont(lhs, rhs);
}

template<typename T>
T py_rcont(const T &lhs, const T &rhs) {
    return tbgal::rcont(lhs, rhs);
}

template<typename T>
T py_reverse(const T &mv) {
    return tbgal::reverse(mv);
}

template<typename T>
T py_inverse(const T &mv) {
    return tbgal::inv(mv);
}

template<typename T>
double py_rnorm(const T &mv) {
    return tbgal::rnorm(mv);
}

template<typename T>
double py_rnorm_sqr(const T &mv) {
    return tbgal::rnorm_sqr(mv);
}

template<typename T>
T py_dual(const T &mv) {
    return tbgal::dual(mv);
}

template<typename T>
T py_undual(const T &mv) {
    return tbgal::undual(mv);
}

template<typename T>
double py_sp(const T &lhs, const T &rhs) {
    return tbgal::sp(lhs, rhs);
}
