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
fig = plt.figure(figsize=(18, 10), dpi=100)
fig.suptitle('Angle variation according to the voltage of the quadrupoles', fontsize=20)
ax = fig.gca(projection='3d')

###############################################################################
# X and Y
#
x = np.linspace(-2000, 2000, 100)
y = np.linspace(-2000, 2000, 100)
x, y = np.meshgrid(x, y)

###############################################################################
# Z calculation
#


# z = abs(sigma_matrices.LensDrift(parameters_calculation.f(0.100, x, 20000, 0.035), 0.100, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1] + sigma_matrices.LensDrift(parameters_calculation.f(0.100, y, 20000, 0.035), 0.050, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])

z = abs(sigma_matrices.Drift(0.050, parameters.drift_alpha, parameters.drift_beta, parameters.drift_gamma, parameters.epsilon).dot(sigma_matrices.LensDrift(parameters_calculation.f(0.100, x, 20000, 0.035), 0.100, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon))[0][1] + sigma_matrices.LensDrift(-parameters_calculation.f(0.100, y, 20000, 0.035), 0.050, parameters.LensDrift_alpha, parameters.LensDrift_beta, parameters.LensDrift_gamma, parameters.epsilon)[0][1])

# z = sigma_matrices.Einzel(0.050, parameters_calculation.f(0.100, x, 20000, 0.035), 0.050, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon)[0][1] - sigma_matrices.Einzel(0.050, parameters_calculation.f(0.100, y, 20000, 0.035), 0.050, parameters.einzel_alpha, parameters.einzel_beta, parameters.einzel_gamma, parameters.epsilon)[0][1]


# print(len(z))
print(z)
print(type(z))


###############################################################################
# Axes parameters
#

# ax.zaxis.set_major_locator(LinearLocator(10))
# ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f')

ax.set_xlabel('First quadrupole voltage(V)')
ax.set_xlim(-3000, 3000)
ax.set_ylabel('Second quadrupole voltage (V)')
ax.set_ylim(-3000, 3000)
ax.set_zlabel(r'$\Delta \alpha$')
# ax.set_zlim(np.min(z), np.max(z))

# np.where(np.isin(z[:, 1], min(z)))

# ax.set_xticks([-2000, , 2000])
# ax.set_xticklabels(['-2000', f'MIN={np.min(z)}', '2000'])
# ax.set_yticks([-2000, np.min(z[:, :]), 2000])
# ax.set_yticklabels(['-2000', f'MIN={np.min(z)}', '2000'])


surf = ax.plot_surface(x, y, z, alpha=0.7, cmap=cm.coolwarm, antialiased=True)
# cset = ax.contour(x, y, z, zdir='z', offset=np.min(z), cmap=cm.coolwarm)
cset = ax.contourf(x, y, z, zdir='z', alpha=0.7, offset=np.min(z), cmap=cm.coolwarm)
cset = ax.contour(x, y, z, zdir='y', offset=3000, cmap=cm.coolwarm)
cset = ax.contour(x, y, z, zdir='x', offset=-3000, cmap=cm.coolwarm)

fig.colorbar(surf, cmap=cm.coolwarm)


plt.show()
