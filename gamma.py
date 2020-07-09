from math import *


def gamma(a, b, eps):
    return (eps**2 + a**2 * eps**2) / (b * eps**2)


print(gamma(1.4, 2.3, 30))
