from math import *


def gamma(a, b, eps):
    return (eps**2 + a**2 * eps**2) / (b * eps**2)


def k(V0, V, R):
    return V0 / (2 * V * R**2)
