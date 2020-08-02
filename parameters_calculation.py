'''
@file : parameters_calculation.py
@brief :  Functions to calculate some other parameters. 

@author : Corentin MICHEL
creation : 20/07/2020
'''

from math import *


def gamma(a, b, eps):
    return (eps**2 + a**2 * eps**2) / (b * eps**2)


def k(V0, V, R):
    return sqrt(V0 / (2 * V * R**2))


def f(L, V0, V, R):  # L en m, R en m
    return (R**2 * V) / (L * (V0 / 2))
