import parameters_calculation
from math import *

epsilon = 30  # Input emittance for the ion beam.
n_std = 2.2977  # n_std : float 	The number of standard deviations to determine the ellipse's radiuses. We take chi2=2.2977 for 68.3% (epsilon rms)
mu = 0, 0  # Two_dimensional mean (we take 0,0 by default, because that means that the beam is in average in the center of the ion gun).
scale = 1, 1  # Dimension we give to the beam.

#############################################################################
#
# Input signal
#

input_alpha = 0
input_beta = 2.3
input_gamma = parameters_calculation.gamma(input_alpha, input_beta, epsilon)

#############################################################################
#
# Drift parameters
#

drift_L = 0.424
drift_alpha = input_alpha
drift_beta = input_beta
drift_gamma = parameters_calculation.gamma(drift_alpha, drift_beta, epsilon)

#############################################################################
#
# Thin lens parameters
#

lens_f = 0.520
lens_alpha = input_alpha
lens_beta = input_beta
lens_gamma = parameters_calculation.gamma(lens_alpha, lens_beta, epsilon)

#############################################################################
#
# Thin lens + drift parameters
#

LensDrift_f = 0.520
LensDrift_L = 0.100
LensDrift_alpha = input_alpha
LensDrift_beta = input_beta
LensDrift_gamma = parameters_calculation.gamma(LensDrift_alpha, LensDrift_beta, epsilon)

#############################################################################
#
# Magnetic dipole parameters
#

dipoleMag_phi = 90
dipoleMag_p = 0.0260
dipoleMag_alpha = input_alpha
dipoleMag_beta = input_beta
dipoleMag_gamma = parameters_calculation.gamma(dipoleMag_alpha, dipoleMag_beta, epsilon)

#############################################################################
#
# Einzel lens parameters
#

einzel_L1 = 0.424
einzel_f = 0.520
einzel_L2 = 0.424
einzel_alpha = input_alpha
einzel_beta = input_beta
einzel_gamma = parameters_calculation.gamma(einzel_alpha, einzel_beta, epsilon)

#############################################################################
#
# Electrostatic quadrupole parameters
#
V0 = 1800
V = 20000
R = 0.035
quadru_L = 0.100
quadru_drift_L = 0.100
quadru_alpha = input_alpha
quadru_beta = input_beta
quadru_gamma = parameters_calculation.gamma(quadru_alpha, quadru_beta, epsilon)
quadru_k = sqrt(parameters_calculation.k(V0, V, R))  # k [m(-1/2)]
print(f"k= {quadru_k} m(-1/2)")

#############################################################################
#
# Doublet parameters
#

V0 = 1800
V = 20000
R = 0.035
doublet_L1 = 0.100
doublet_L2 = 0.040
doublet_L3 = 0.100
doublet_L4 = 0.100
doublet_alpha = input_alpha
doublet_beta = input_beta
doublet_gamma = parameters_calculation.gamma(doublet_alpha, doublet_beta, epsilon)
doublet_k1 = sqrt(parameters_calculation.k(V0, V, R))
doublet_k2 = sqrt(parameters_calculation.k(V0, V, R))
