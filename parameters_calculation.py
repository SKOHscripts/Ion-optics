from math import *


def gamma(a, b, eps):
    return (eps**2 + a**2 * eps**2) / (b * eps**2)


def f(L, V0, V, R):
    return (R**2 * V) / (L * (V0 / 2))
