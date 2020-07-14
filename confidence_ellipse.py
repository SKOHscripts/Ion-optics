import numpy as np
import pylab
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

import sigma_matrices
import parameters_calculation
import parameters

###############################################################################
#
# Inspired from:
#
# https://matplotlib.org/3.1.0/gallery/statistics/confidence_ellipse.html


###############################################################################
#
# The plotting function itself
# """"""""""""""""""""""""""""
#
# This function plots the confidence ellipse of the covariance of the given
# array-like variables x and y. The ellipse is plotted into the given
# axes-object ax.
#
# The radiuses of the ellipse can be controlled by n_std which is the number
# of standard deviations. The default value is 3 which makes the ellipse
# enclose 99.7% of the points (given the data is normally distributed
# like in these examples).


def confidence_ellipse(x, y, ax, n_std=parameters.n_std, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of `x` and `y`

    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1] / np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
                      width=ell_radius_x * 2,
                      height=ell_radius_y * 2,
                      facecolor=facecolor,
                      **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * np.sqrt(n_std)
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * np.sqrt(n_std)
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)

    return ax.add_patch(ellipse)


###############################################################################
#
# A helper function to create a correlated dataset
# """"""""""""""""""""""""""""""""""""""""""""""""
#
# Creates a random two-dimesional dataset with the specified
# two-dimensional mean (mu) and dimensions (scale).
# The correlation can be controlled by the param 'dependency',
# a 2x2 matrix.

def get_correlated_dataset(n, dependency, mu, scale, **kwargs):
    latent = np.random.randn(n, 2)
    dependent = latent.dot(dependency)
    scaled = dependent * scale
    scaled_with_offset = scaled + mu
    # return x and y of the new, correlated dataset
    return scaled_with_offset[:, 0], scaled_with_offset[:, 1]


###############################################################################
#
# Visualization of output beams of different optical elements from calculated "sigma" matrices.
#

np.random.seed(0)

LENS = {

    'Input signal': sigma_matrices.Input(parameters.lens_alpha, parameters.lens_beta, parameters.lens_gamma, parameters.epsilon),

    'Lens': sigma_matrices.Lens(parameters.lens_f, parameters.lens_alpha, parameters.lens_beta, parameters.lens_gamma, parameters.epsilon),

}

DRIFT = {

    'Input signal': sigma_matrices.Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),

    'Drift': sigma_matrices.Drift(parameters.drift_L, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),
}

LENSDRIFT = {

    'Input signal': sigma_matrices.Input(parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon),

    'Lens + drift': sigma_matrices.LensDrift(parameters.LensDrift_f, parameters.LensDrift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon),
}

EINZEL = {

    'Input signal': sigma_matrices.Input(parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon),

    'Einzel lens': sigma_matrices.Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon),
}

DIPOLEMAG = {

    'Input signal': sigma_matrices.Input(parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon),

    'Magnetic dipole for x plan': sigma_matrices.DipoleMag_x(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon),

    'Magnetic dipole for y plan': sigma_matrices.DipoleMag_y(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon),
}

QUADRUPOLE = {

    'Input signal': sigma_matrices.Input(parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon),

    'Electrostatic quadrupole for x plan': sigma_matrices.QUADRU_X,

    'Electrostatic quadrupole for y plan': sigma_matrices.QUADRU_Y,

}

QUADRU_DOUBLET = {

    'Input signal': sigma_matrices.Input(parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon),

    'Doublet for x plan': sigma_matrices.Quadru_doubletPM(parameters.doublet_k1, parameters.doublet_k2, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon),

    'Doublet for y plan': sigma_matrices.Quadru_doubletMP(parameters.doublet_k1, parameters.doublet_k2, parameters.doublet_L1, parameters.doublet_L2, parameters.doublet_L3, parameters.doublet_L4, parameters.doublet_alpha, parameters.doublet_beta, parameters.doublet_gamma, parameters.epsilon),

}


mu = parameters.mu
scale = parameters.scale


fig, axs = plt.subplots(1, 3)
for ax, (title, dependency) in zip(axs, QUADRU_DOUBLET.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    ax.hexbin(x, y, gridsize=50, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-30, 30))
    ax.set_ylim((-10, 10))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()

fig, axs = plt.subplots(1, 3)
for ax, (title, dependency) in zip(axs, QUADRUPOLE.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale, label='Dataset')
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    ax.hexbin(x, y, gridsize=50, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-30, 30))
    ax.set_ylim((-10, 10))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()


fig, axs = plt.subplots(1, 3)
for ax, (title, dependency) in zip(axs, DIPOLEMAG.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    ax.hexbin(x, y, gridsize=50, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-20, 20))
    ax.set_ylim((-5, 5))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()

fig, axs = plt.subplots(1, 2)
for ax, (title, dependency) in zip(axs, EINZEL.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    ax.hexbin(x, y, gridsize=50, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-10, 10))
    ax.set_ylim((-20, 20))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()

fig, axs = plt.subplots(1, 2)
for ax, (title, dependency) in zip(axs, LENSDRIFT.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    ax.hexbin(x, y, gridsize=50, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-10, 10))
    ax.set_ylim((-20, 20))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()

fig, axs = plt.subplots(1, 2)
for ax, (title, dependency) in zip(axs, LENS.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    ax.hexbin(x, y, gridsize=50, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-10, 10))
    ax.set_ylim((-20, 20))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()

fig, axs = plt.subplots(1, 2)
for ax, (title, dependency) in zip(axs, DRIFT.items()):
    x, y = get_correlated_dataset(50000, dependency, mu, scale, label='Dataset')
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)
    ax.hexbin(x, y, gridsize=100, cmap='inferno')
    confidence_ellipse(x, y, ax, label=r'$\epsilon_{rms}=68.3\%$', edgecolor='red', linewidth=2)

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    # ax.set_aspect('equal')
    ax.set_xlim((-10, 10))
    ax.set_ylim((-5, 5))
    ax.set_xlabel('X/Y (mm)')
    ax.set_ylabel("X'/Y' (mrad)")
    ax.legend()

# plt.show()


###############################################################################
#
# Different number of standard deviations
# """""""""""""""""""""""""""""""""""""""
#
#


fig, ax_nstd = plt.subplots(figsize=(6, 6))

dependency_nstd = sigma_matrices.Drift(parameters.drift_L, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon)

mu = parameters.mu
scale = parameters.scale

ax_nstd.axvline(c='grey', lw=1)
ax_nstd.axhline(c='grey', lw=1)
# ax_nstd.set_aspect('equal')
ax_nstd.set_xlim((-10, 10))
ax_nstd.set_ylim((-5, 5))

# ax_nstd.hexbin(x, y, gridsize=50, cmap='inferno')
x, y = get_correlated_dataset(10000, dependency_nstd, mu, scale, label='Dataset')
ax_nstd.scatter(x, y, s=0.5)

'''
https://www.visiondummy.com/2014/04/draw-error-ellipse-representing-covariance-matrix/
'''
'''
chi^2 calculated from 
https://www.fourmilab.ch/rpkp/experiments/analysis/chiCalc.html
'''
confidence_ellipse(x, y, ax_nstd, n_std=2.2977,
                   label=r'$\epsilon_{rms}=68.3\%$', edgecolor='firebrick', linewidth=2)
confidence_ellipse(x, y, ax_nstd, n_std=4.6051,
                   label=r'$2.\epsilon_{rms}=90\%$', edgecolor='firebrick', linestyle='--', linewidth=2)
confidence_ellipse(x, y, ax_nstd, n_std=6.2021,
                   label=r'$4.\epsilon_{rms}=95.5\%$', edgecolor='firebrick', linestyle=':', linewidth=2)

ax_nstd.scatter(mu[0], mu[1], c='red', s=3)
ax_nstd.set_title('Different standard deviations for a drift')
ax_nstd.set_xlabel('X/Y (mm)')
ax_nstd.set_ylabel("X'/Y' (mrad)")
ax_nstd.legend()
plt.show()


###############################################################################
#
# Using the keyword arguments
# """""""""""""""""""""""""""
#
# Use the kwargs specified for matplotlib.patches.Patch in order
# to have the ellipse rendered in different ways.


# fig, ax_kwargs = plt.subplots(figsize=(6, 6))
# dependency_kwargs = sigma_matrices.Drift(parameters.drift_L, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon)

# mu = parameters.mu
# scale = parameters.scale

# ax_kwargs.axvline(c='grey', lw=1)
# ax_kwargs.axhline(c='grey', lw=1)

# x, y = get_correlated_dataset(2000, dependency_kwargs, mu, scale)
# # Plot the ellipse with zorder=0 in order to demonstrate
# # its transparency (caused by the use of alpha).
# confidence_ellipse(x, y, ax_kwargs,
#                    alpha=0.5, label=r'$4.\epsilon_{rms}=95.5\%$', facecolor='pink', edgecolor='pink', zorder=0)

# ax_kwargs.scatter(x, y, s=0.5)
# ax_kwargs.scatter(mu[0], mu[1], c='red', s=3)
# ax_kwargs.set_title(f'Thin lens')
# ax_kwargs.set_aspect('equal')
# # ax_kwargs.set_xlim((-60, 60))
# # ax_kwargs.set_ylim((-10, 10))
# ax_kwargs.set_xlabel('X/Y (mm)')
# ax_kwargs.set_ylabel("X'/Y' (mrad)")
# ax_kwargs.legend()
# fig.subplots_adjust(hspace=0.25)
# plt.show()
