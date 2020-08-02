'''
@file : voltage_optimisation.py
@brief : This program uses the approximation of an electrostatic quadrupole by a thin lens edged with drifts to find the voltage couple(s) that allow to have a minimum variation of the emissivity angle at the output with respect to the input. 

@author : Corentin MICHEL
creation : 20/07/2020
'''

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
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

    # zx = abs(sigma_matrices.Quadru_doubletPM_approx(parameters_calculation.f(parameters.quadru_L, x, 20000, 0.035), parameters_calculation.f(parameters.quadru_L, y, 20000, 0.035), parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)[0][1])

    # zy = abs(sigma_matrices.Quadru_doubletMP_approx(parameters_calculation.f(parameters.quadru_L, x, 20000, 0.035), parameters_calculation.f(parameters.quadru_L, y, 20000, 0.035), parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)[0][1])

    # z = zx * zy
    # print(z)
    ###############################################################################
    # Axes parameters
    #

    # ax.set_xlabel('Quadrupole voltage (1st or 2nd) (V)', fontsize=15)
    ax.set_xlim(parameters.optim_min, parameters.optim_max)
    # ax.set_ylabel('Quadrupole voltage (1st or 2nd) (V)', fontsize=15)
    ax.set_ylim(parameters.optim_min, parameters.optim_max)
    ax.set_zlabel(r'$\Delta \alpha$', fontsize=20)
    ax.set_xlabel(r'$1^{st}$ quadrupole (V)', fontsize=20)
    ax.set_ylabel(r'$2^{nd}$ quadrupole (V)', fontsize=20)

    plt.grid(which='both')
    plt.grid(which='minor', alpha=0.3, linestyle='--')
    grid_x_ticks = np.arange(parameters.optim_min, parameters.optim_max, 50)
    grid_y_ticks = np.arange(parameters.optim_min, parameters.optim_max, 50)

    ax.set_xticks(grid_x_ticks, minor=True)
    ax.set_yticks(grid_y_ticks, minor=True)

    zmin = np.where(z == np.amin(z))
    print("zmin=", zmin)
    xmin = x[zmin][0]
    ymin = y[zmin][0]
    print(f'{xmin:.2f} ({i}/{parameters.optim_nb_points})')


extra1 = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
extra2 = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
extra3 = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

plt.legend([extra1, extra2, extra3], (fr'$\Delta \alpha min$= {np.amin(z):.2E} rad', fr'$V_1= {abs(xmin):.2f}$ V', fr'$V_2= {abs(ymin):.2f}$ V'))


surf = ax.plot_surface(x, y, z, alpha=0.6, cmap=cm.gist_earth, antialiased=True)
#cset = ax.contour(x, y, z, zdir='z', alpha=0.7, offset=np.min(z), cmap=cm.gist_earth)
#cset = ax.contour(x, y, z, zdir='y', offset=parameters.optim_max + 100, cmap=cm.gist_earth)
#cset = ax.contour(x, y, z, zdir='x', offset=parameters.optim_min - 100, cmap=cm.gist_earth)

cbar = fig.colorbar(surf, cmap=cm.gist_earth)
cbar.set_label('For x and y plans together')

# fig.savefig('pics/optimisation.png', bbox_inches='tight', dpi=100)


plt.show()
