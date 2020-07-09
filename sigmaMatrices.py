from pylab import *
import numpy as np
from math import *

import gamma
import parameters

#############################################################################
#
# Input signal
#
#


def Input(a, b, g, eps):

    sigma = np.array([[b, -a],
                      [-a, g]])

    return eps / math.pi * sigma


print("Input signal: \n", Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon))

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
    C = np.dot(T, sigma_in)

    return eps / math.pi * np.dot(C, Ttr)


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

    return eps / math.pi * np.dot(np.dot(T, sigma_in), Ttr)


print("")
print("Lens: \n", Lens(parameters.lens_f, parameters.lens_alpha, parameters.lens_beta, parameters.lens_gamma, parameters.epsilon))


#############################################################################
#
# Thin lens + drift
#
#
# f=focale de la lentille, L=longueur du drift, a=alpha, b=beta, g=gamma, eps=epsilon

def LensDrift(f, L, a, b, g, eps):
    T = np.array([[1 - L / f, L], [-1 / f, 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    return eps / math.pi * np.dot(np.dot(T, sigma_in), Ttr)


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

    return eps / math.pi * np.dot(np.dot(T, sigma_in), Ttr)


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

    return eps / math.pi * np.dot(np.dot(T, sigma_in), Ttr)


print("")
print("Magnetic dipole for y plan: \n", DipoleMag_y(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon))


#############################################################################
#
# Einzel lens
#
#
def Einzel(L1, f, L2, a, b, g, eps):
    T = np.array([[1 - L2 / f, L1 * (1 - L2 / f) + L2],
                  [-1 / f, -L1 / f + 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    return eps / math.pi * np.dot(np.dot(T, sigma_in), Ttr)


print("")
print("Einzel: \n", Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon))

#############################################################################
#
# Einzel lens
#
#


def Quadri_conv(, a, b, g, eps):
    T = np.array([[1 - L2 / f, L1 * (1 - L2 / f) + L2],
                  [-1 / f, -L1 / f + 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    return eps / math.pi * np.dot(np.dot(T, sigma_in), Ttr)


print("")
print("Quadripôle électrostatique convergent: \n", Quadri_conv(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon))
