from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np

import parameters
import parameters_calculation
import sigma_matrices


fig = plt.figure()
ax = fig.gca(projection='3d')

x = np.linspace(-2000, 2000, 10)
y = np.linspace(-2000, 2000, 10)
x, y = np.meshgrid(x, y)


R = np.sqrt(x**2 + y**2)
z = np.sin(R)


# print(x)
# print(y)
# print(z)

surf = ax.plot_surface(x, y, z, cmap=cm.coolwarm,
                       linewidth=0, antialiased=True)

# ax.set_zlim(-1.01, 1.01)
ax.zaxis.set_major_locator(LinearLocator(10))
ax.zaxis.set_major_formatter(FormatStrFormatter('%.02f'))
ax.set_xlabel('First quadrupole voltage(V)')
ax.set_ylabel('Second quadrupole voltage (V)')
ax.set_zlabel('$\Delta°$')
ax.set_title('Angle variation according to the voltage of the quadrupoles', fontweight="bold")

fig.colorbar(surf, shrink=1, aspect=15)

plt.show()
