from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

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
kx = parameters_calculation.k(x, 20000, 0.035)
ky = parameters_calculation.k(y, 20000, 0.035)
print(kx)

# for i in kx:
#     print(kx[i])
#     z = []
#     z = z.extend(sigma_matrices.Quadri_conv(kx[i], parameters.quadru_L, parameters.quadru_drift_L, parameters.quadru_alpha, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)[0][1])

z = sigma_matrices.Input(kx, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)[0][1] - sigma_matrices.Input(ky, parameters.quadru_beta, parameters.quadru_gamma, parameters.epsilon)[0][1]

# print(len(z))
print(z)


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

ax.set_xticks([-2000, np.min(z[:, 0]), 2000])
ax.set_xticklabels(['-2000', f'MIN={np.min(z)}', '2000'])
ax.set_yticks([-2000, np.min(z[0, :]), 2000])
ax.set_yticklabels(['-2000', f'MIN={np.min(z)}', '2000'])


surf = ax.plot_surface(x, y, z, alpha=0.8, cmap=cm.coolwarm, antialiased=True)
# cset = ax.contour(x, y, z, zdir='z', offset=np.min(z), cmap=cm.coolwarm)
cset = ax.contourf(x, y, z, zdir='z', offset=np.min(z), cmap=cm.coolwarm)
cset = ax.contour(x, y, z, zdir='y', offset=3000, cmap=cm.coolwarm)
cset = ax.contour(x, y, z, zdir='x', offset=-3000, cmap=cm.coolwarm)

fig.colorbar(surf, cmap=cm.coolwarm_r)

plt.show()
