from pylab import *
import numpy as np
from math import *

import parameters_calculation
import parameters

#############################################################################
#
# Input signal
#
#

def Input(a, b, g, eps):

    sigma = np.array([[b, -a],
                      [-a, g]])

    return sigma


print("Input signal: \n", Input(parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))

#############################################################################
#
# Drift
#
#


def Drift(L, a, b, g, eps):
    T = np.array([[1, L], [0, 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])
    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Drift: \n", Drift(parameters.drift_L, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon))


#############################################################################
#
# Thin lens
#

def Lens(f, a, b, g, eps):
    T = np.array([[1, 0], [-1 / f, 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Lens: \n", Lens(parameters.lens_f, parameters.lens_alpha, parameters.lens_beta, parameters.lens_gamma, parameters.epsilon))


#############################################################################
#
# Thin lens + drift
#
#
# f=focale de la lentille, L=longueur du drift, a=alpha, b=beta, g=gamma, eps=epsilon

def LensDrift(f, L, a, b, g, eps):
    T = np.array([[1 - (L / f), L], [-1 / f, 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])
    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Lens + Drift: \n", LensDrift(parameters.LensDrift_f, parameters.LensDrift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))


#############################################################################
#
# Magnetic dipole for (x,x') plan
#
#
# phi=angle du dipôle, p=longueur, a=alpha, b=beta, g=gamma, eps=epsilon

def DipoleMag_x(phi, p, a, b, g, eps):
    T = np.array([[cos(phi), p * sin(phi)],
                  [(1 / p) * sin(phi), cos(phi)]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Magnetic dipole for x plan: \n", DipoleMag_x(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon))

#############################################################################
#
# Magnetic dipole for (y,y') plan
#
#
# phi=angle du dipôle, p=longueur, a=alpha, b=beta, g=gamma, eps=epsilon


def DipoleMag_y(phi, p, a, b, g, eps):
    T = np.array([[1, p * phi], [0, 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Magnetic dipole for y plan: \n", DipoleMag_y(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon))


#############################################################################
#
# Einzel lens
#
#
def Einzel(L1, f, L2, a, b, g, eps):
    T = np.array([[1 - (L2 / f), L1 * (1 - (L2 / f)) + L2],
                  [-1 / f, -(L1 / f) + 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = T.dot(sigma_in)

    return np.dot(A, Ttr)


print("")
print("Einzel: \n", Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon))

#############################################################################
#
# Electrostatic quadrupole (WITH A DRIFT)
#
#


def Quadru_conv(k, L1, L2, a, b, g, eps):
    T = np.array([[np.cos(k * L1), (1 / k) * np.sin(k * L1)],
                  [-k * np.sin(k * L1), np.cos(k * L1)]])
    Ttr = np.transpose(T)
    Drifttr = np.transpose(Drift(L2, a, b, g, eps))

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    return Drift(L2, a, b, g, eps).dot(T).dot(sigma_in).dot(Ttr).dot(Drifttr)


print("")
print("Convergent electrostatic quadrupole: \n", Quadru_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


def Quadru_div(k, L1, L2, a, b, g, eps):
    T = np.array([[np.cosh(k * L1), (1 / k) * np.sinh(k * L1)],
                  [k * np.sinh(k * L1), np.cosh(k * L1)]])
    Ttr = T.transpose()
    Drifttr = Drift(L2, a, b, g, eps).transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    return Drift(L2, a, b, g, eps).dot(T).dot(sigma_in).dot(Ttr).dot(Drifttr)


print("")
print("Divergent electrostatic quadrupole: \n", Quadru_div(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


'''
If the incoming beam converges, the quadrupole lens behaves as a convergent lens in x and a divergent lens in y, and vice versa if the beam is divergent.
'''

if parameters.quadru_alpha >= 0:
    QUADRU_X = Quadru_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)
    QUADRU_Y = Quadru_div(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)
else:
    QUADRU_Y = Quadru_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)
    QUADRU_X = Quadru_div(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)


def Quadru_doubletPM(k1, k2, L1, L2, L3, L4, a, b, g, eps):
    if parameters.quadru_alpha >= 0:
        QUADRU_X = Quadru_conv(k1, L1, L2, a, b, g, eps)
        QUADRU_Y = Quadru_div(k2, L3, L4, a, b, g, eps)
    else:
        QUADRU_Y = Quadru_conv(k1, L3, L4, a, b, g, eps)
        QUADRU_X = Quadru_div(k2, L1, L2, a, b, g, eps)

    T1 = QUADRU_X
    T2 = QUADRU_Y
    T1tr = T1.transpose()
    T2tr = T2.transpose()

    T = np.dot(T1, T2)

    sigma_in = np.array([[b, -a],
                         [-a, g]])
    Ttr = np.dot(T1tr, T2tr)

    return np.dot(T, sigma_in, Ttr)


print("")
print("Doublet of quadrupoles for x plan: \n", Quadru_doubletPM(parameters.doublet_k1, parameters.doublet_k2, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon))


def Quadru_doubletMP(k1, k2, L1, L2, L3, L4, a, b, g, eps):
    if parameters.quadru_alpha >= 0:
        QUADRU_X = Quadru_conv(k1, L1, L2, a, b, g, eps)
        QUADRU_Y = Quadru_div(k2, L3, L4, a, b, g, eps)
    else:
        QUADRU_Y = Quadru_conv(k1, L3, L4, a, b, g, eps)
        QUADRU_X = Quadru_div(k2, L1, L2, a, b, g, eps)

    T1 = QUADRU_X
    T2 = QUADRU_Y
    T1tr = T1.transpose()
    T2tr = T2.transpose()

    T = np.dot(T1, T2)

    sigma_in = np.array([[b, -a],
                         [-a, g]])
    Ttr = np.dot(T1tr, T2tr)

    return np.dot(T, sigma_in, Ttr)


print("")
print("Doublet of quadrupoles for y plan: \n", Quadru_doubletMP(parameters.doublet_k1, parameters.doublet_k2, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon))


def Quadru_doubletPM_approx(f1, f2, L1, L2, V, R, a, b, g, eps):
    return Drift(L1 / 2, a, b, g, eps).dot(LensDrift(f1, L2, a, b, g, eps)).dot(LensDrift(-f2, L2 / 2, a, b, g, eps))


print("")
print("Doublet_approxPM: \n", Quadru_doubletPM_approx(parameters_calculation.f(parameters.quadru_L, 900, parameters.V, parameters.R), parameters_calculation.f(parameters.quadru_L, 2000, parameters.V, parameters.R), parameters.quadru_L, parameters.quadru_drift_L, parameters.V, parameters.R, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


def Quadru_doubletMP_approx(f, L1, L2, V, R, a, b, g, eps):
    return Drift(L1 / 2, a, b, g, eps).dot(LensDrift(-f, L2, a, b, g, eps)).dot(LensDrift(f, L2, a, b, g, eps))


print("")
print("Doublet_approxMP: \n", Quadru_doubletMP_approx(-parameters_calculation.f(parameters.quadru_L, 900, parameters.V, parameters.R), parameters.quadru_L, parameters.quadru_drift_L, parameters.V, parameters.R, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))
