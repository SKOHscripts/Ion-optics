from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np
from math import *

import parameters
import parameters_calculation
import sigma_matrices

###############################################################################
# Figure
#
fig = plt.figure(figsize=(18, 18), dpi=100)
fig.suptitle('Angle variation according to the voltage of the quadrupoles for x and y plans together', fontsize=20)
ax = fig.gca(projection='3d')

###############################################################################
# X and Y
#
x = np.linspace(-2000, 2000, 2000)
y = np.linspace(-2000, 2000, 2000)
x, y = np.meshgrid(x, y)

###############################################################################
# Z calculation
#

z = abs(sigma_matrices.Drift(0.050, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon).dot(sigma_matrices.LensDrift(parameters_calculation.f(0.100, x, 20000, 0.035), 0.090, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + sigma_matrices.LensDrift(-parameters_calculation.f(0.100, y, 20000, 0.035), 0.050, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1]) * abs(sigma_matrices.Drift(0.050, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon).dot(sigma_matrices.LensDrift(-parameters_calculation.f(0.100, x, 20000, 0.035), 0.100, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + sigma_matrices.LensDrift(parameters_calculation.f(0.100, y, 20000, 0.035), 0.050, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])


# z = sigma_matrices.Einzel(0.050, parameters_calculation.f(0.100, x, 20000, 0.035), 0.050, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon)[0][1] - sigma_matrices.Einzel(0.050, parameters_calculation.f(0.100, y, 20000, 0.035), 0.050, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon)[0][1]


###############################################################################
# Axes parameters
#

# ax.set_xlabel('Quadrupole voltage (1st or 2nd) (V)', fontsize=15)
ax.set_xlim(-3000, 3000)
# ax.set_ylabel('Quadrupole voltage (1st or 2nd) (V)', fontsize=15)
ax.set_ylim(-3000, 3000)
ax.set_zlabel(r'$\Delta \alpha$', fontsize=20)

zmin = np.where(z == np.amin(z))
xmin = x[zmin][0]
ymin = y[zmin][0]


ax.set_xticks([-2000, ymin, abs(xmin), 2000])
ax.set_xticklabels(['-2000 V', f'1st combination= {ymin:.2f} V', f'2nd combination= {abs(xmin):.2f} V', '2000 V'])


ax.set_yticks([-2000, xmin, abs(ymin), 2000])
ax.set_yticklabels(['-2000 V', f'1st combination= {xmin:.2f} V', f'2nd combination= {abs(ymin):.2f} V', '2000 V'])

ax.set_zticks([np.amin(z)])
ax.set_zticklabels([f'min= {np.amin(z):.2E} rad'])


surf = ax.plot_surface(x, y, z, alpha=0.8, cmap=cm.gist_earth, antialiased=False)
cset = ax.contour(x, y, z, zdir='z', alpha=0.7, offset=np.min(z), cmap=cm.gist_earth)
cset = ax.contour(x, y, z, zdir='y', offset=3000, cmap=cm.gist_earth)
cset = ax.contour(x, y, z, zdir='x', offset=-3000, cmap=cm.gist_earth)

cbar = fig.colorbar(surf, cmap=cm.gist_earth)
cbar.set_label('For x and y plans together')

plt.show()
