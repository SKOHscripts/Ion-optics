import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse
import matplotlib.transforms as transforms

import sigmaMatrices
import gamma
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
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
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

def get_correlated_dataset(n, dependency, mu, scale):
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

    'Input signal': sigmaMatrices.Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),

    'Lens': sigmaMatrices.Lens(parameters.lens_f, parameters.lens_alpha, parameters.lens_beta, parameters.lens_gamma, parameters.epsilon),

}

DRIFT = {

    'Input signal': sigmaMatrices.Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),

    'Drift': sigmaMatrices.Drift(parameters.drift_L, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),
}

LENSDRIFT = {

    'Input signal': sigmaMatrices.Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),

    'Lens + drift': sigmaMatrices.LensDrift(parameters.LensDrift_f, parameters.LensDrift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon),
}

EINZEL = {

    'Input signal': sigmaMatrices.Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),

    'Einzel lens': sigmaMatrices.Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon),
}

DIPOLEMAG = {

    'Input signal': sigmaMatrices.Input(parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon),

    'Magnetic dipole for x plan': sigmaMatrices.DipoleMag_x(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon),

    'Magnetic dipole for y plan': sigmaMatrices.DipoleMag_y(parameters.dipoleMag_phi, parameters.dipoleMag_p, parameters.dipoleMag_alpha, parameters.dipoleMag_beta, parameters.dipoleMag_gamma, parameters.epsilon),
}

mu = parameters.mu
scale = parameters.scale

fig, axs = plt.subplots(1, 2, figsize=(6, 6))
for ax, (title, dependency) in zip(axs, LENS.items()):
    x, y = get_correlated_dataset(2000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    confidence_ellipse(x, y, ax, edgecolor='red')

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    ax.axis('equal')

fig, axs = plt.subplots(1, 2, figsize=(6, 6))
for ax, (title, dependency) in zip(axs, DRIFT.items()):
    x, y = get_correlated_dataset(2000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    confidence_ellipse(x, y, ax, edgecolor='red')

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    ax.axis('equal')

fig, axs = plt.subplots(1, 2, figsize=(6, 6))
for ax, (title, dependency) in zip(axs, LENSDRIFT.items()):
    x, y = get_correlated_dataset(2000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    confidence_ellipse(x, y, ax, edgecolor='red')

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    ax.axis('equal')

fig, axs = plt.subplots(1, 2, figsize=(6, 6))
for ax, (title, dependency) in zip(axs, EINZEL.items()):
    x, y = get_correlated_dataset(2000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    confidence_ellipse(x, y, ax, edgecolor='red')

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    ax.axis('equal')

fig, axs = plt.subplots(1, 3, figsize=(6, 6))
for ax, (title, dependency) in zip(axs, DIPOLEMAG.items()):
    x, y = get_correlated_dataset(2000, dependency, mu, scale)
    ax.scatter(x, y, s=0.5)

    ax.axvline(c='grey', lw=1)
    ax.axhline(c='grey', lw=1)

    confidence_ellipse(x, y, ax, edgecolor='red')

    ax.scatter(mu[0], mu[1], c='red', s=3)
    ax.set_title(title)
    ax.axis('equal')

plt.show()


###############################################################################
#
# Different number of standard deviations
# """""""""""""""""""""""""""""""""""""""
#
# A plot with n_std = 3 (blue), 2 (purple) and 1 (red)


fig, ax_nstd = plt.subplots(figsize=(6, 6))

dependency_nstd = sigmaMatrices.Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon)

mu = parameters.mu
scale = parameters.scale

ax_nstd.axvline(c='grey', lw=1)
ax_nstd.axhline(c='grey', lw=1)
ax.axis('equal')

x, y = get_correlated_dataset(3000, dependency_nstd, mu, scale)
ax_nstd.scatter(x, y, s=0.5)

confidence_ellipse(x, y, ax_nstd, n_std=1,
                   label=r'$1\sigma=(\beta*\epsilon)^{0.5}$', edgecolor='firebrick', linewidth=2)
confidence_ellipse(x, y, ax_nstd, n_std=2,
                   label=r'$2\sigma=2(\beta*\epsilon)^{0.5}$', edgecolor='firebrick', linestyle='--', linewidth=2)
confidence_ellipse(x, y, ax_nstd, n_std=3,
                   label=r'$3\sigma=3(\beta*\epsilon)^{0.5}$', edgecolor='firebrick', linestyle=':', linewidth=2)

ax_nstd.scatter(mu[0], mu[1], c='red', s=3)
ax_nstd.set_title('Different standard deviations for a thin lens')
ax_nstd.legend()
plt.show()


###############################################################################
#
# Using the keyword arguments
# """""""""""""""""""""""""""
#
# Use the kwargs specified for matplotlib.patches.Patch in order
# to have the ellipse rendered in different ways.


fig, ax_kwargs = plt.subplots(figsize=(6, 6))
dependency_kwargs = sigmaMatrices.Einzel(parameters.einzel_L1, parameters.einzel_f, parameters.einzel_L2, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon)

mu = parameters.mu
scale = parameters.scale

ax_kwargs.axvline(c='grey', lw=1)
ax_kwargs.axhline(c='grey', lw=1)

x, y = get_correlated_dataset(2000, dependency_kwargs, mu, scale)
# Plot the ellipse with zorder=0 in order to demonstrate
# its transparency (caused by the use of alpha).
confidence_ellipse(x, y, ax_kwargs,
                   alpha=0.5, facecolor='pink', edgecolor='pink', zorder=0)

ax_kwargs.scatter(x, y, s=0.5)
ax_kwargs.scatter(mu[0], mu[1], c='red', s=3)
ax_kwargs.set_title(f'Thin lens')
ax.axis('equal')

fig.subplots_adjust(hspace=0.25)
plt.show()
