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

precision = []
for i in range(10, parameters.optim_nb_points):
    x = np.linspace(parameters.optim_min, parameters.optim_max, i)
    y = np.linspace(parameters.optim_min, parameters.optim_max, i)
    x, y = np.meshgrid(x, y)

    ###############################################################################
    # Z calculation
    #

    zx = abs(sigma_matrices.Drift(parameters.optim_drift_L, parameters.optim_drift_alpha, parameters.optim_drift_beta, parameters.optim_drift_gamma, parameters.epsilon).dot(sigma_matrices.LensDrift(parameters_calculation.f(parameters.optim_LensDrift_L, x, parameters.optim_V, parameters.optim_R), parameters.optim_drift_L + parameters.optim_dist, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + sigma_matrices.LensDrift(-parameters_calculation.f(parameters.optim_LensDrift_L, y, parameters.optim_V, parameters.optim_R), parameters.optim_drift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])

    zy = abs(sigma_matrices.Drift(parameters.optim_drift_L, parameters.optim_drift_alpha, parameters.optim_drift_beta, parameters.optim_drift_gamma, parameters.epsilon).dot(sigma_matrices.LensDrift(-parameters_calculation.f(parameters.optim_LensDrift_L, x, parameters.optim_V, parameters.optim_R), 2 * parameters.optim_drift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + sigma_matrices.LensDrift(parameters_calculation.f(parameters.optim_LensDrift_L, y, parameters.optim_V, parameters.optim_R), 2 * parameters.optim_drift_L, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])

    z = zx * zy

    ###############################################################################
    # Axes parameters
    #

    # ax.set_xlabel('Quadrupole voltage (1st or 2nd) (V)', fontsize=15)
    ax.set_xlim(parameters.optim_min, parameters.optim_max)
    # ax.set_ylabel('Quadrupole voltage (1st or 2nd) (V)', fontsize=15)
    ax.set_ylim(parameters.optim_min, parameters.optim_max)
    ax.set_zlabel(r'$\Delta \alpha$', fontsize=20)

    zmin = np.where(z == np.amin(z))
    xmin = x[zmin][0]
    ymin = y[zmin][0]
    precision.append(xmin)
    print(f'{xmin:.2f} ({i}/{parameters.optim_nb_points})')


ax.set_xticks([parameters.optim_min, ymin, abs(xmin), parameters.optim_max])
ax.set_xticklabels([f'{parameters.optim_min} V', f'1st combination= {ymin:.2f} V', f'2nd combination= {abs(xmin):.2f} V', f'{parameters.optim_max} V'])


ax.set_yticks([parameters.optim_min, xmin, abs(ymin), parameters.optim_max])
ax.set_yticklabels([f'{parameters.optim_min} V', f'1st combination= {xmin:.2f} V', f'2nd combination= {abs(ymin):.2f} V', f'{parameters.optim_max} V'])

ax.set_zticks([np.amin(z)])
ax.set_zticklabels([f'min= {np.amin(z):.2E} rad'])


surf = ax.plot_surface(x, y, z, alpha=0.8, cmap=cm.gist_earth)
cset = ax.contour(x, y, z, zdir='z', alpha=0.7, offset=np.min(z), cmap=cm.gist_earth)
cset = ax.contour(x, y, z, zdir='y', offset=parameters.optim_max + 100, cmap=cm.gist_earth)
cset = ax.contour(x, y, z, zdir='x', offset=parameters.optim_min - 100, cmap=cm.gist_earth)

cbar = fig.colorbar(surf, cmap=cm.gist_earth)
cbar.set_label('For x and y plans together')

fig.savefig('pics/optimisation.png', bbox_inches='tight', dpi=100)

xpre = precision
plt.subplots()
print(xpre)
print(len(xpre))
plt.plot(xpre)
plt.suptitle('Precision', fontsize=20)


plt.show()
