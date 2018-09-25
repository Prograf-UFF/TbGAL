import os
import sys
sys.path.append(os.path.join("..", "..", "..", "build"))
from pycliff.metric import EuclideanMetric
from pycliff.multivector import *
from math import sin, cos, atan2, sqrt, pi

def setup():
    DIMS = 5
    Multivector.set_N(DIMS)
    generate_T(EuclideanMetric())


def versor_factorization_test():
    B = (e(1)^e(2))+(e(1)^e(3)+(e(2)^e(3)))
    B = B * (1.0 / NORM(B))

    V = cos(pi/4.0) - (sin(pi/4.0)*B)
    vectors = fact_versor(V)

    V_new = vectors[0]
    for v in vectors[1:]:
        V_new = GP(V_new, v)

    assert V == V_new


def versor_factorization_apply_test():
    B = (e(1)^e(2))+(e(1)^e(3)+(e(2)^e(3)))
    B = B * (1.0 / NORM(B))

    V = cos(pi/4.0) - (sin(pi/4.0)*B)

    k1 = (IGP(GP(V, e(1)), V))

    vectors = fact_versor(V)
    v1 = vectors[0]
    v2 = vectors[1]

    k2  = (IGP(GP(v1, IGP(GP(v2, e(1)), v2)), v1))

    assert k1 == k2


def teardown():
    pass
