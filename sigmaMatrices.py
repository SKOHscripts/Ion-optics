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
    T = np.array([[1 - L2 / f, L1 * (1 - L2 / f) + L2],
                  [-1 / f, -L1 / f + 1]])

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Einzel: \n", Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon))

#############################################################################
#
# Electrostatic quadrupole (WITH A DRIFT)
#
#


def Quadri_conv(k, L1, L2, a, b, g, eps):
    T = np.array([[cos(sqrt(k) * L1), (1 / sqrt(k)) * sin(sqrt(k) * L1)],
                  [-sqrt(k) * sin(sqrt(k) * L1), cos(sqrt(k) * L1)]])

    # print(f"T={T}")
    T = np.dot(Drift(L2, a, b, g, eps), T)

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Convergent electrostatic quadripole: \n", Quadri_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


def Quadri_div(k, L1, L2, a, b, g, eps):
    T = np.array([[cosh(sqrt(k) * L1), (1 / sqrt(k)) * sinh(sqrt(k) * L1)],
                  [sqrt(k) * sinh(sqrt(k) * L1), cosh(sqrt(k) * L1)]])

    # print(f"T={T}")
    T = np.dot(Drift(L2, a, b, g, eps), T)

    Ttr = T.transpose()

    sigma_in = np.array([[b, -a],
                         [-a, g]])

    A = np.dot(T, sigma_in)

    return np.dot(A, Ttr)


print("")
print("Divergent electrostatic quadripole: \n", Quadri_div(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


'''
If the incoming beam converges, the quadrupole lens behaves as a convergent lens in x and a divergent lens in y, and vice versa if the beam is divergent.
'''
if parameters.quadru_alpha >= 0:
    QUADRU_X = Quadri_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)
    QUADRU_Y = Quadri_div(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)
else:
    QUADRU_Y = Quadri_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)
    QUADRU_X = Quadri_div(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)


# def Quadru_doubletPM(k, L1, L2, L3, L4, a, b, g, eps):
#     T1 = QUADRU_X
#     T2 = QUADRU_Y

#     # print(f"T={T}")

#     T = np.dot(T2, T1)

#     print(f"T={T}")

#     Ttr = T.transpose()
#     print(f"Ttr={Ttr}")

#     sigma_in = np.array([[b, -a],
#                          [-a, g]])

#     return T.dot(sigma_in).dot(Ttr)


# print("")
# print("Doublet of quadripoles : \n", Quadru_doubletPM(parameters.doublet_k, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon))


# def Quadru_doubletMP(k, L1, L2, L3, L4, a, b, g, eps):
#     T1 = QUADRU_X
#     T2 = QUADRU_Y

#     # print(f"T={T}")

#     T = np.dot(T2, T1)

#     print(f"T={T}")

#     Ttr = T.transpose()
#     print(f"Ttr={Ttr}")

#     sigma_in = np.array([[b, -a],
#                          [-a, g]])

#     return T.dot(sigma_in).dot(Ttr)


# print("")
# print("Doublet of quadripoles : \n", Quadru_doubletMP(parameters.doublet_k, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon))
