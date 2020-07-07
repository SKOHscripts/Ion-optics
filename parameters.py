import calcul_gamma

epsilon = 30

#############################################################################
#
# Drift parameters
#

drift_L = 424
drift_alpha = 1.4
drift_beta = 2.3
drift_gamma = calcul_gamma.gamma(drift_alpha, drift_beta, epsilon)

#############################################################################
#
# Thin lens parameters
#

lens_f = 520
lens_alpha = 1.4
lens_beta = 2.3
lens_gamma = calcul_gamma.gamma(lens_alpha, lens_beta, epsilon)

#############################################################################
#
# Thin lens + drift parameters
#

LensDrift_f = 520
LensDrift_L = 100
LensDrift_alpha = 1.4
LensDrift_beta = 2.3
LensDrift_gamma = calcul_gamma.gamma(LensDrift_alpha, LensDrift_beta, epsilon)

#############################################################################
#
# Magnetic dipole parameters
#

dipoleMag_phi = 45
dipoleMag_p = 20
dipoleMag_alpha = 1.4
dipoleMag_beta = 2.3
dipoleMag_gamma = calcul_gamma.gamma(dipoleMag_alpha, dipoleMag_beta, epsilon)

#############################################################################
#
# Einzel lens parameters
#

einzel_L1 = 424
einzel_f = 520
einzel_L2 = 424
einzel_alpha = 1.4
einzel_beta = 2.3
einzel_gamma = calcul_gamma.gamma(einzel_alpha, einzel_beta, epsilon)
