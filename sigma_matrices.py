'''
@file : sigma_matrices.py
@brief : Here we calculate the sigma matrices of every optical element of the pipeline thanks to their transfert matrices. 

@author : Corentin MICHEL
creation : 20/07/2020
'''

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
  D = np.array([[1, L2],
                [0, 1]])

  sigma_in = np.array([[b, -a],
                       [-a, g]])

  return D.dot(T).dot(sigma_in).dot(T.transpose()).dot(D.transpose())


print("")
print("Convergent electrostatic quadrupole: \n", Quadru_conv(parameters.quadru_k, parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


def Quadru_div(k, L1, L2, a, b, g, eps):
  T = np.array([[np.cosh(k * L1), (1 / k) * np.sinh(k * L1)],
                [k * np.sinh(k * L1), np.cosh(k * L1)]])
  D = np.array([[1, L2],
                [0, 1]])

  sigma_in = np.array([[b, -a],
                       [-a, g]])

  return D.dot(T).dot(sigma_in).dot(T.transpose()).dot(D.transpose())


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


#####################


def Quadru_doubletPM(k1, k2, Lq1, Lq2, Ld1, Ld2, a, b, g, eps):
  L1 = np.array([[1, Ld1],
                 [0, 1]])
  L2 = np.array([[1, Ld2],
                 [0, 1]])

  T1 = np.array([[cos(k1 * Lq1), 1 / k1 * sin(k1 * Lq1)],
                 [-k1 * sin(k1 * Lq1), cos(k1 * Lq1)]])

  T2 = np.array([[cosh(k2 * Lq2), 1 / k2 * sinh(k2 * Lq2)],
                 [k2 * sinh(k2 * Lq2), cosh(k2 * Lq2)]])
  sigma_in = np.array([[b, -a],
                       [-a, g]])

  return L2.dot(T2).dot(L1).dot(T1).dot(sigma_in).dot(T1.transpose()).dot(L1.transpose()).dot(T2.transpose()).dot(L2.transpose())


print("")
print("Doublet of quadrupoles for x plan: \n", Quadru_doubletPM(parameters.doublet_k1, parameters.doublet_k2, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon))


def Quadru_doubletMP(k1, k2, Lq1, Lq2, Ld1, Ld2, a, b, g, eps):
  L1 = np.array([[1, Ld1],
                 [0, 1]])
  L2 = np.array([[1, Ld2],
                 [0, 1]])

  T1 = np.array([[cosh(k1 * Lq1), 1 / k1 * sinh(k1 * Lq1)],
                 [k1 * sinh(k1 * Lq1), cosh(k1 * Lq1)]])
  T2 = np.array([[cos(k2 * Lq2), 1 / k2 * sin(k2 * Lq2)],
                 [-k2 * sin(k2 * Lq2), cos(k2 * Lq2)]])

  sigma_in = np.array([[b, -a],
                       [-a, g]])

  return L2.dot(T2).dot(L1).dot(T1).dot(sigma_in).dot(T1.transpose()).dot(L1.transpose()).dot(T2.transpose()).dot(L2.transpose())


print("")
print("Doublet of quadrupoles for y plan: \n", Quadru_doubletMP(parameters.doublet_k1, parameters.doublet_k2, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon))


def Quadru_doubletPM_approx(f1, f2, L1, L2, a, b, g, eps):
  Td = np.array([[1, L1 / 2],
                 [0, 1]])
  Tdl1 = np.array([[1 - (L2 / 2 * f1), L2 / 2],
                   [-1 / f1, 1]])

  Tdl2 = np.array([[1 + (L2 / (2 * f2)), L2 / 2],
                   [1 / f2, 1]])
  T = Td.dot(Tdl1).dot(Tdl2)
  sigma_in = np.array([[b, -a],
                       [-a, g]])
  return T.dot(sigma_in).dot(T.transpose())


print("")
print("Doublet_approxPM: \n", Quadru_doubletPM_approx(parameters_calculation.f(parameters.quadru_L, 900, 20000, 0.035), parameters_calculation.f(parameters.quadru_L, 900, 20000, 0.035), parameters.quadru_drift_L, parameters.quadru_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


def Quadru_doubletMP_approx(f1, f2, L1, L2, a, b, g, eps):
  Td = np.array([[1, L1 / 2],
                 [0, 1]])
  Tdl1 = np.array([[1 + (L2 / 2 * f1), L2 / 2],
                   [1 / f1, 1]])
  Tdl2 = np.array([[1 - (L2 / (2 * f2)), L2 / 2],
                   [-1 / f2, 1]])
  T = Td.dot(Tdl1).dot(Tdl2)
  sigma_in = np.array([[b, -a],
                       [-a, g]])
  return T.dot(sigma_in).dot(T.transpose())


print("")
print("Doublet_approxMP: \n", Quadru_doubletMP_approx(parameters_calculation.f(parameters.quadru_L, 900, 20000, 0.035), parameters_calculation.f(parameters.quadru_L, 900, 20000, 0.035), parameters.quadru_drift_L, parameters.quadru_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon))


#####################################################################
#####################################################################

def Quadru_tripletPMP(k1, k2, k3, Lq1, Lq2, Lq3, Ld1, Ld2, Ld3, a, b, g, eps):
  L1 = np.array([[1, Ld1],
                 [0, 1]])
  L2 = np.array([[1, Ld2],
                 [0, 1]])
  L3 = np.array([[1, Ld3],
                 [0, 1]])

  T1 = np.array([[cos(k1 * Lq1), 1 / k1 * sin(k1 * Lq1)],
                 [-k1 * sin(k1 * Lq1), cos(k1 * Lq1)]])
  T2 = np.array([[cosh(k2 * Lq2), 1 / k2 * sinh(k2 * Lq2)],
                 [k2 * sinh(k2 * Lq2), cosh(k2 * Lq2)]])
  T3 = np.array([[cos(k3 * Lq3), 1 / k3 * sin(k3 * Lq3)],
                 [-k3 * sin(k3 * Lq3), cos(k3 * Lq3)]])

  sigma_in = np.array([[b, -a], [-a, g]])

  return L3.dot(T3).dot(L2).dot(T2).dot(L1).dot(T1).dot(sigma_in).dot(T1.transpose()).dot(L1.transpose()).dot(T2.transpose()).dot(L2.transpose()).dot(T3.transpose()).dot(L3.transpose())


print("")
print("Triplet of quadrupoles for x plan: \n", Quadru_tripletPMP(parameters.triplet_k1, parameters.triplet_k2, parameters.triplet_k3, parameters.triplet_Lq1, parameters.triplet_Lq2, parameters.triplet_Lq3, parameters.triplet_Ld1, parameters.triplet_Ld2, parameters.triplet_Ld3, parameters.triplet_alpha, parameters.triplet_beta, parameters.triplet_gamma, parameters.epsilon))


def Quadru_tripletMPM(k1, k2, k3, Lq1, Lq2, Lq3, Ld1, Ld2, Ld3, a, b, g, eps):
  L1 = np.array([[1, Ld1],
                 [0, 1]])
  L2 = np.array([[1, Ld2],
                 [0, 1]])
  L3 = np.array([[1, Ld3],
                 [0, 1]])
  T1 = np.array([[cosh(k1 * Lq1), 1 / k1 * sinh(k1 * Lq1)],
                 [k1 * sinh(k1 * Lq1), cosh(k1 * Lq1)]])
  T2 = np.array([[cos(k2 * Lq2), 1 / k2 * sin(k2 * Lq2)],
                 [-k2 * sin(k2 * Lq2), cos(k2 * Lq2)]])
  T3 = np.array([[cosh(k3 * Lq3), 1 / k3 * sinh(k3 * Lq3)],
                 [k3 * sinh(k3 * Lq3), cosh(k3 * Lq3)]])

  sigma_in = np.array([[b, -a], [-a, g]])

  return L3.dot(T3).dot(L2).dot(T2).dot(L1).dot(T1).dot(sigma_in).dot(T1.transpose()).dot(L1.transpose()).dot(T2.transpose()).dot(L2.transpose()).dot(T3.transpose()).dot(L3.transpose())


print("")
print("Triplet of quadrupoles for y plan: \n", Quadru_tripletMPM(parameters.triplet_k1, parameters.triplet_k2, parameters.triplet_k3, parameters.triplet_Lq1, parameters.triplet_Lq2, parameters.triplet_Lq3, parameters.triplet_Ld1, parameters.triplet_Ld2, parameters.triplet_Ld3, parameters.triplet_alpha, parameters.triplet_beta, parameters.triplet_gamma, parameters.epsilon))

zx = abs(Drift(parameters.optim_drift_L, parameters.optim_drift_alpha, parameters.optim_drift_beta, parameters.optim_drift_gamma, parameters.epsilon).dot(LensDrift(parameters_calculation.f(parameters.optim_LensDrift_L, 733, parameters.optim_V, parameters.optim_R), parameters.optim_drift_L + parameters.optim_dist, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + LensDrift(-parameters_calculation.f(parameters.optim_LensDrift_L, 1250, parameters.optim_V, parameters.optim_R), parameters.optim_drift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])

zy = abs(Drift(parameters.optim_drift_L, parameters.optim_drift_alpha, parameters.optim_drift_beta, parameters.optim_drift_gamma, parameters.epsilon).dot(LensDrift(-parameters_calculation.f(parameters.optim_LensDrift_L, 733, parameters.optim_V, parameters.optim_R), 2 * parameters.optim_drift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + LensDrift(parameters_calculation.f(parameters.optim_LensDrift_L, 1250, parameters.optim_V, parameters.optim_R), 2 * parameters.optim_drift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])

print(zx * zy)
