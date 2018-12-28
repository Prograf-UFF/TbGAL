# Copyright 2017 Leandro Augusto Frata Fernandes.
# Distributed under the GNU General Public License, Version 3.
#
# Last update: August 9th, 2017.
#
# This implementation is part of the material accompanying the book:
#
#    Álgebra Geométrica e Aplicações
#    Fernandes, L. A. F.; Lavor, C.; Oliveira, M. M.
#
#    Sociedade Brasileira de Matemática Aplicada e Computacional (SBMAC)
#    Notas em Matemática Aplicada (NoMA), Volume 85, 2017
#
#    ISSN: 2175-3385 | e-ISSN: 2236-5915
#    ISBN: 978-85-8215-081-8 | e-ISBN: 978-85-8215-080-1
#
#    http://www.ic.uff.br/algebrageometrica

import abc
import math
import numbers

DEFAULT_METRIC = None  # Default metric.
DEFAULT_TOLERANCE = 1.0e-8  # Default tolerance value for floating-point comparison.


class Metric(metaclass=abc.ABCMeta):
    """The superclass for metric matrix definition.
    """

    def __init__(self):
        """Default constructor.
        """

    @abc.abstractmethod
    def matrix_entry(self, row, col):
        """Return a matrix entry.
        :rtype: float
        """


class OrthogonalMetric(Metric):
    """The superclass for orhogonal matrix definition.
    """

    def __init__(self):
        """Default constructor.
        """
        super().__init__()

    def matrix_entry(self, row,  col):
        """Return a matrix entry.
        """
        assert isinstance(row, int) and row >= 0
        assert isinstance(col, int) and col >= 0

        if row == col:
            return self.diagonal_entry(row)
        else:
            return 0

    @abc.abstractmethod
    def diagonal_entry(self, index):
        """Return a matrix entry on the main diagonal.
        :rtype: float
        """

    @abc.abstractmethod
    def metric_factor(self, mask):
        """Return the metric factor computed for the given basis blade.
        :rtype: float
        """


class EuclideanMetric(OrthogonalMetric):
    """Euclidean metric matrix.
    """

    def __init__(self):
        """Default constructor.
        """
        super().__init__()

    def diagonal_entry(self, index):
        """Return a matrix entry on the main diagonal.
        """
        assert isinstance(index, int) and index >= 0

        return 1

    def metric_factor(self, mask):
        """Return the metric factor computed for the given basis blade.
        """
        assert isinstance(mask, int)

        return 1


class SignedMetric(OrthogonalMetric):
    """Orhogonal metric matrix having a given p,q signature.
    """

    def __init__(self, p, q):
        """Default constructor.
        """
        assert isinstance(p, int) and p >= 0
        assert isinstance(q, int) and q >= 0

        super().__init__()
        self.__p = p
        self.__q = q

    def diagonal_entry(self, index):
        """Return a matrix entry on the main diagonal.
        """
        assert isinstance(index, int) and index >= 1

        if index <= self.__p:
            return +1
        elif index <= self.__p + self.__q:
            return -1
        else:
            return 0

    def metric_factor(self, mask):
        """Return the metric factor computed for the given basis blade.
        """
        assert isinstance(mask, int)

        index = 0
        result = 1
        while mask != 0:
            if (mask & 1) != 0:
                result *= self.diagonal_entry(index)
            mask >>= 1
            index += 1
        return result


class Multivector:
    """Multivector defines a sparce multivector structure.
    """

    def __init__(self, arg=None):
        """Multivector()
           Multivector(scalar)
           Multivector({mask1: coefficient1, mask2: coefficient2, ...})
        """
        from collections import OrderedDict

        if arg is None or (isinstance(arg, numbers.Number) and arg == 0):
            self.__components = {}
        elif isinstance(arg, numbers.Number):
            self.__components = {0: arg}
        elif isinstance(arg, dict):
            for mask, coef in arg.items():
                assert isinstance(mask, int) and isinstance(coef, numbers.Number)
            self.__components = OrderedDict(sorted(arg.copy().items(), key=lambda comp: comp[0]))
        else:
            raise ValueError

    def __mul__(self, other):
        """self * other -> gp_em(self, other), where self or other is a scalar value.
        """
        assert is_scalar(self) or is_scalar(other)

        return gp_em(self, other)

    def __rmul__(self, other):
        """other * self -> gp_em(other, self), where self or other is a scalar value.
        """
        assert is_scalar(self) or is_scalar(other)

        return gp_em(other, self)

    def __pos__(self):
        """+self
        """
        return self

    def __add__(self, other):
        """self + other
        """
        return merge(self, other, lambda lhs_coef, rhs_coef: lhs_coef + rhs_coef)

    def __radd__(self, other):
        """other + self
        """
        return merge(other, self, lambda lhs_coef, rhs_coef: lhs_coef + rhs_coef)

    def __neg__(self):
        """-self
        """
        return Multivector({mask: -coef for mask, coef in self.components()})

    def __sub__(self, other):
        """self - other
        """
        return merge(self, other, lambda lhs_coef, rhs_coef: lhs_coef - rhs_coef)

    def __rsub__(self, other):
        """other - self
        """
        return merge(other, self, lambda lhs_coef, rhs_coef: lhs_coef - rhs_coef)

    def __xor__(self, other):
        """self ^ other -> op(self, other)
        """
        return op(self, other)

    def __rxor__(self, other):
        """other ^ self -> op(other, self)
        """
        return op(other, self)

    def tostring(self):
        """str(self)
        """
        result = ''
        is_first_comp = True

        for mask, coef in self.components():
            if coef != 0:
                if is_first_comp:
                    if coef < 0:
                        result += '-'
                    is_first_comp = False
                else:
                    if coef < 0:
                        result += ' - '
                    else:
                        result += ' + '

                if mask == 0 or abs(coef) != 1:
                    result += str(abs(coef))

                if mask != 0:
                    if abs(coef) != 1:
                        result += ' * '

                    index = 0
                    is_first_factor = True
                    while mask != 0:
                        if (mask & 1) != 0:
                            if not is_first_factor:
                                result += '^'
                            result += 'e%d' % index
                            is_first_factor = False
                        index += 1
                        mask >>= 1

        if is_first_comp:
            return str(0)
        else:
            return result

    
    def __repr__(self):
        return self.tostring()

    def __str__(self):
        return self.tostring()
   
    def components(self):
        """Return a set like object providing a view on multivector's components.
        """
        return self.__components.items()


def __population_count(mask):
    """Population (number of 1 bits) count.
    """
    mask = mask - ((mask >> 1) & 0x55555555)
    mask = (mask & 0x33333333) + ((mask >> 2) & 0x33333333)
    return (((mask + (mask >> 4) & 0xF0F0F0F) * 0x1010101) & 0xffffffff) >> 24


def __reordering_sign(lhs_mask, rhs_mask):
    """Sign change due to canonical reordering.
    """
    lhs_mask >>= 1
    changes = 0
    while lhs_mask != 0:
        changes += __population_count(lhs_mask & rhs_mask)
        lhs_mask >>= 1
    if (changes & 1) == 0:
        return +1
    else:
        return -1


def __used_basis_vectors(arg):
    """Extract the basis vectors used by the given mutivector.
    """
    bits = 0
    for mask, coef in arg.components():
        bits ^= mask

    i = 0
    result = []
    while bits != 0:
        if (bits & 1) == 1:
            result.append(e(i))
        bits >>= 1
        i += 1

    return result


def conjugation(arg):
    """Clifford conjugation operation.
    """
    def change(mask):
        k = __population_count(mask)
        return (((k * (k + 1)) >> 1) & 1) != 0

    return graded_uminus(arg, change)


def delta(lhs, rhs, metric=DEFAULT_METRIC, tol=DEFAULT_TOLERANCE):
    """Delta product.
    """
    m = gp(lhs, rhs, metric)

    max_grade = -1
    for mask, coef in m.components():
        if abs(coef) > tol:
            curr_grade = __population_count(mask)
            if max_grade < curr_grade:
                max_grade = curr_grade

    return take_grade(m, max_grade)


def delta_em(lhs, rhs, tol=DEFAULT_TOLERANCE):
    """Delta product under Euclidean metric.
    """
    return delta(lhs, rhs, EuclideanMetric(), tol)


def dual(arg, pseudoscalar, metric=DEFAULT_METRIC):
    """Dualization.
    """
    return gp(arg, inv(pseudoscalar, metric), metric)


def dual_em(arg, pseudoscalar):
    """Dualization under Euclidean metric.
    """
    return dual(arg, pseudoscalar, EuclideanMetric())


def e(index):
    """Return the i-th basis vector.
    """
    assert isinstance(index, int) and index >= 0

    return Multivector({1 << index: 1})


def equal(lhs, rhs, tol=DEFAULT_TOLERANCE):
    """Equality check function.
    """
    if not isinstance(lhs, Multivector):
        lhs = Multivector(lhs)
    if not isinstance(rhs, Multivector):
        rhs = Multivector(rhs)

    lhs_itr = iter(lhs.components())
    rhs_itr = iter(rhs.components())
    lhs_comp = next(lhs_itr, None)
    rhs_comp = next(rhs_itr, None)
    while lhs_comp is not None and rhs_comp is not None:
        if lhs_comp[0] < rhs_comp[0]:
            if abs(lhs_comp[1]) > tol:
                return False
            lhs_comp = next(lhs_itr, None)
        elif lhs_comp[0] > rhs_comp[0]:
            if abs(rhs_comp[1]) > tol:
                return False
            rhs_comp = next(rhs_itr, None)
        else:
            if abs(lhs_comp[1] - rhs_comp[1]) > tol:
                return False
            lhs_comp = next(lhs_itr, None)
            rhs_comp = next(rhs_itr, None)

    while lhs_comp is not None:
        if abs(lhs_comp[1]) > tol:
            return False
        lhs_comp = next(lhs_itr, None)

    while rhs_comp is not None:
        if abs(rhs_comp[1]) > tol:
            return False
        rhs_comp = next(rhs_itr, None)

    return True


def gp(lhs, rhs, metric=DEFAULT_METRIC):
    """Geometric product.
    """
    return graded_gp(lhs, rhs, metric,
                     lambda lhs_grade, rhs_grade, result_grade: True)


def gp_em(lhs, rhs):
    """Geometric product under Euclidean metric.
    """
    return gp(lhs, rhs, EuclideanMetric())


def grade(arg, tol=DEFAULT_TOLERANCE):
    """Return the grade of the multivector, -1 for mixed grade, or None for zero multivector.
    """
    if not isinstance(arg, Multivector):
        arg = Multivector(arg)

    result = None
    for mask, coef in arg.components():
        if abs(coef) > tol:
            curr = __population_count(mask)
            if result is not None:
                if result != curr:
                    return -1
            else:
                result = curr
    return result


def graded_gp(lhs, rhs, metric, keep):
    """General implementation of bilinear products.
    """
    if not isinstance(lhs, Multivector):
        lhs = Multivector(lhs)
    if not isinstance(rhs, Multivector):
        rhs = Multivector(rhs)
    assert isinstance(metric, Metric)

    if isinstance(metric, OrthogonalMetric):
        result = {}

        for lhs_mask, lhs_coef in lhs.components():
            for rhs_mask, rhs_coef in rhs.components():
                mask = lhs_mask ^ rhs_mask
                if keep(lhs_mask, rhs_mask, mask):
                    coef = lhs_coef * rhs_coef * __reordering_sign(lhs_mask, rhs_mask)\
                           * metric.metric_factor(lhs_mask & rhs_mask)
                    if mask in result:
                        result[mask] += coef
                    else:
                        result[mask] = coef

        return Multivector(result)
    else:
        raise NotImplemented  # Non-orthogonal metric


def graded_uminus(rhs, change):
    """General implementation of grade-based sign-change operations.
    """
    if not isinstance(rhs, Multivector):
        rhs = Multivector(rhs)

    return Multivector({mask: (-coef if change(mask) else coef) for mask, coef in rhs.components()})


def igp(lhs, rhs, metric=DEFAULT_METRIC, tol=DEFAULT_TOLERANCE):
    """Inverse geometric product.
    """
    return gp(lhs, inv(rhs, metric, tol), metric)


def igp_em(lhs, rhs, tol=DEFAULT_TOLERANCE):
    """Inverse geometric product under Euclidean metric.
    """
    return igp(lhs, rhs, EuclideanMetric(), tol)


def inv(arg, metric=DEFAULT_METRIC, tol=DEFAULT_TOLERANCE):
    """Blade inverse.
    """
    rev = reversion(arg)
    return gp(rev, 1 / native(scp(arg, rev, metric), tol), metric)


def inv_em(arg, tol=DEFAULT_TOLERANCE):
    """Blade inverse under Euclidean metric.
    """
    return inv(arg, EuclideanMetric(), tol)


def involution(arg):
    """Grade involution operation.
    """
    return graded_uminus(arg, lambda mask: (__population_count(mask) & 1) != 0)


def is_blade(arg, tol=DEFAULT_TOLERANCE):
    """Returns whether the given multivetor is a blade.
    """
    return grade(arg, tol) != -1 and is_versor_em(arg, tol)


def is_scalar(arg, tol=DEFAULT_TOLERANCE):
    """Return whether the given argument is a scalar value.
    """
    assert isinstance(tol, numbers.Number)

    if not isinstance(arg, Multivector):
        arg = Multivector(arg)

    for mask, coef in arg.components():
        if mask != 0 and abs(coef) > tol:
            return False
    return True


def is_versor(arg, metric=DEFAULT_METRIC, tol=DEFAULT_TOLERANCE):
    """Returns whether the given multivetor is a versor.
    """
    if not isinstance(arg, Multivector):
        arg = Multivector(arg)

    arg_hat = involution(arg)
    inv_arg = inv(arg, metric, tol)

    tmp1 = gp(arg_hat, inv_arg, metric)
    if grade(tmp1, tol) != 0:
        return False

    tmp2 = gp(inv_arg, arg_hat, metric)
    if not equal(tmp1, tmp2, tol):
        return False

    arg_tild = reversion(arg)
    if isinstance(metric, OrthogonalMetric):
        basis_vectors = __used_basis_vectors(arg)
    else:
        basis_vectors = metric.basis_vectors()
    for v in basis_vectors:
        if grade(gp(gp(arg_hat, v, metric), arg_tild, metric), tol) != 1:
            return False

    return True


def is_versor_em(arg, tol=DEFAULT_TOLERANCE):
    """Returns whether the given multivetor is a versor under Euclidean metric.
    """
    return is_versor(arg, EuclideanMetric(), tol)


def is_zero(arg, tol=DEFAULT_TOLERANCE):
    """Return whether the given argument is zero.
    """
    assert isinstance(tol, numbers.Number)

    if not isinstance(arg, Multivector):
        arg = Multivector(arg)

    for mask, coef in arg.components():
        if abs(coef) > tol:
            return False
    return True


def lcont(lhs, rhs, metric=DEFAULT_METRIC):
    """Left contraction.
    """
    return graded_gp(lhs, rhs, metric,
                     lambda lhs_grade, rhs_grade, result_grade: result_grade == (rhs_grade - lhs_grade))


def lcont_em(lhs, rhs):
    """Left contraction under Euclidean metric.
    """
    return lcont(lhs, rhs, EuclideanMetric())


def merge(lhs, rhs, func):
    """General implementation of add and sub operations.
    """
    if not isinstance(lhs, Multivector):
        lhs = Multivector(lhs)
    if not isinstance(rhs, Multivector):
        rhs = Multivector(rhs)

    result = {}

    lhs_itr = iter(lhs.components())
    rhs_itr = iter(rhs.components())
    lhs_comp = next(lhs_itr, None)
    rhs_comp = next(rhs_itr, None)
    while lhs_comp is not None and rhs_comp is not None:
        if lhs_comp[0] < rhs_comp[0]:
            result[lhs_comp[0]] = func(lhs_comp[1], 0)
            lhs_comp = next(lhs_itr, None)
        elif lhs_comp[0] > rhs_comp[0]:
            result[rhs_comp[0]] = func(0, rhs_comp[1])
            rhs_comp = next(rhs_itr, None)
        else:
            result[lhs_comp[0]] = func(lhs_comp[1], rhs_comp[1])
            lhs_comp = next(lhs_itr, None)
            rhs_comp = next(rhs_itr, None)

    while lhs_comp is not None:
        result[lhs_comp[0]] = func(lhs_comp[1], 0)
        lhs_comp = next(lhs_itr, None)

    while rhs_comp is not None:
        result[rhs_comp[0]] = func(0, rhs_comp[1])
        rhs_comp = next(rhs_itr, None)

    return Multivector(result)


def native(arg, tol=DEFAULT_TOLERANCE):
    """Scalar multivector to Number conversion.
    """
    assert isinstance(arg, Multivector) and is_scalar(arg, tol)

    comp = next(iter(arg.components()), None)
    if comp is not None and comp[0] == 0:
        return comp[1]
    else:
        return 0


def op(lhs, rhs):
    """Outer (wedge) product.
    """
    return graded_gp(lhs, rhs, EuclideanMetric(),
                     lambda lhs_grade, rhs_grade, result_grade: result_grade == (lhs_grade + rhs_grade))


def reversion(arg):
    """Reversion operation.
    """
    def change(mask):
        k = __population_count(mask)
        return (((k * (k - 1)) >> 1) & 1) != 0

    return graded_uminus(arg, change)


def rcont(lhs, rhs, metric=DEFAULT_METRIC):
    """Right contraction.
    """
    return graded_gp(lhs, rhs, metric,
                     lambda lhs_grade, rhs_grade, result_grade: result_grade == (lhs_grade - rhs_grade))


def rcont_em(lhs, rhs):
    """Right contraction under Euclidean metric.
    """
    return rcont(lhs, rhs, EuclideanMetric())


def rnorm(arg, metric=DEFAULT_METRIC, tol=DEFAULT_TOLERANCE):
    """Reverse norm.
    """
    return math.sqrt(sqr_rnorm(arg, metric, tol))


def rnorm_em(arg, tol=DEFAULT_TOLERANCE):
    """Reverse norm under Euclidean metric.
    """
    return rnorm(arg, EuclideanMetric(), tol)


def scalar(arg):
    """Number to scalar multivector conversion.
    """
    assert isinstance(arg, numbers.Number)

    return Multivector(arg)


def scp(lhs, rhs, metric=DEFAULT_METRIC):
    """Scalar product.
    """
    return graded_gp(lhs, rhs, metric,
                     lambda lhs_grade, rhs_grade, result_grade: result_grade == 0)


def scp_em(lhs, rhs):
    """Scalar product under Euclidean metric.
    """
    return scp(lhs, rhs, EuclideanMetric())


def sqr_rnorm(arg, metric=DEFAULT_METRIC, tol=DEFAULT_TOLERANCE):
    """Squared reverse norm.
    """
    return native(scp(arg, reversion(arg), metric), tol)


def sqr_rnorm_em(arg, tol=DEFAULT_TOLERANCE):
    """Squared reverse norm under Euclidean metric.
    """
    return sqr_rnorm(arg, EuclideanMetric(), tol)


def take_grade(arg, k):
    """Grade extraction operation.
    """
    assert isinstance(k, int)

    if not isinstance(arg, Multivector):
        arg = Multivector(arg)

    return Multivector({mask: coef for mask, coef in arg.components() if __population_count(mask) == k})


def undual(arg, pseudoscalar, metric=DEFAULT_METRIC):
    """Undualization.
    """
    return gp(arg, pseudoscalar, metric)


def undual_em(arg, pseudoscalar):
    """Undualization under Euclidean metric.
    """
    return undual(arg, pseudoscalar, EuclideanMetric())
