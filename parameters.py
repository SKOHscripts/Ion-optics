import gamma
import math

epsilon = 30  # Input emittance for the ion beam.
n_std = 2  # n_std : float 	The number of standard deviations to determine the ellipse's radiuses.
mu = 0, 0  # Two_dimensional mean (we take 0,0 by default, because that means that the beam is in average in the center of the ion gun).
scale = 1, 1  # Dimension we give to the beam.

#############################################################################
#
# Drift parameters
#

drift_L = 424
drift_alpha = 1.4
drift_beta = 2.3
drift_gamma = gamma.gamma(drift_alpha, drift_beta, epsilon)

#############################################################################
#
# Thin lens parameters
#

lens_f = 520
lens_alpha = 1.4
lens_beta = 2.3
lens_gamma = gamma.gamma(lens_alpha, lens_beta, epsilon)

#############################################################################
#
# Thin lens + drift parameters
#

LensDrift_f = 0.520
LensDrift_L = 0.100
LensDrift_alpha = 1.4
LensDrift_beta = 2.3
LensDrift_gamma = gamma.gamma(LensDrift_alpha, LensDrift_beta, epsilon)

#############################################################################
#
# Magnetic dipole parameters
#

dipoleMag_phi = 45
dipoleMag_p = 0.020
dipoleMag_alpha = 1.4
dipoleMag_beta = 2.3
dipoleMag_gamma = gamma.gamma(dipoleMag_alpha, dipoleMag_beta, epsilon)

#############################################################################
#
# Einzel lens parameters
#

einzel_L1 = 0.424
einzel_f = 0.520
einzel_L2 = 0.424
einzel_alpha = 1.4
einzel_beta = 2.3
einzel_gamma = gamma.gamma(einzel_alpha, einzel_beta, epsilon)
