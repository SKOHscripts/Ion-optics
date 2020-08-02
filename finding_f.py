'''
@file : finding_V.py
@brief : Einzel lenses are approximated by thin lenses with drifts. For each lens, we look for the focal length to be applied to obtain a parallel beam at the output, knowing that we emit a parallel beam at the input. We will therefore pass through each lens successively and note the values in a table.  (x,x') plane and (y,y') plane are the same.  

@author : Corentin MICHEL
creation : 20/07/2020
'''

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import parameters
import parameters_calculation
import sigma_matrices
import confidence_ellipse


def find_f(L1, L2, L3, L4, a, b, g, eps, V, R):
    F = [0, 0, 0, 0]
    to_check = parameters.findingf_tocheck

    print('')

    ###########################################################################
    # First Einzel Lens
    #
    alpha1 = [100]
    for i in range(len(to_check)):
        f1 = parameters_calculation.f(L1, to_check[i], parameters.findingf_V, parameters.findingf_R)
        alpha1.append(sigma_matrices.Einzel(L1, f1, L1, a, b, g, eps)[0][1])

    for i in range(len(alpha1)):
        if abs(alpha1[i]) == np.amin([abs(ele) for ele in alpha1]):

            min1 = i
            F[0] = parameters_calculation.f(L1, to_check[min1], parameters.findingf_V, parameters.findingf_R)

            print(f'1er : a= {alpha1[i]} : f= {F[0]}')

    print(F)

    ###########################################################################
    # Second Einzel Lens
    #
    alpha2 = [100]
    for i in range(len(to_check)):
        f2 = parameters_calculation.f(L2, to_check[i], parameters.findingf_V, parameters.findingf_R)
        alpha2.append(sigma_matrices.Einzel(L2, f2, L2, alpha1[min1], b, g, eps)[0][1])

    for i in range(len(alpha2)):
        if abs(alpha2[i]) == np.amin([abs(ele) for ele in alpha2]):

            min2 = i
            F[1] = parameters_calculation.f(L2, to_check[min2], parameters.findingf_V, parameters.findingf_R)

            print(f'2e : a= {alpha2[i]} : f= {F[1]}')

    print(F)

    ###########################################################################
    # Third Einzel Lens
    #

    alpha3 = [100]
    for i in range(len(to_check)):
        f3 = parameters_calculation.f(L3, to_check[i], parameters.findingf_V, parameters.findingf_R)
        alpha3.append(sigma_matrices.Einzel(L3, f3, L3, alpha2[min2], b, g, eps)[0][1])

    for i in range(len(alpha3)):
        if abs(alpha3[i]) == np.amin([abs(ele) for ele in alpha3]):

            min3 = i
            F[2] = parameters_calculation.f(L3, to_check[min3], parameters.findingf_V, parameters.findingf_R)

            print(f'3e : a= {alpha3[i]} : f= {F[2]}')

    print(F)

    ###########################################################################
    # Fourth and last Einzel Lens
    #

    alpha4 = [100]
    for i in range(len(to_check)):
        f4 = parameters_calculation.f(L4, to_check[i], parameters.findingf_V, parameters.findingf_R)
        alpha4.append(sigma_matrices.Einzel(L4, f4, L4, alpha3[min3], b, g, eps)[0][1])

    for i in range(len(alpha4)):
        if abs(alpha4[i]) == np.amin([abs(ele) for ele in alpha4]):

            min4 = i
            F[3] = parameters_calculation.f(L4, to_check[min4], parameters.findingf_V, parameters.findingf_R)

            print(f'4e : a= {alpha4[i]} : f= {F[3]}')

    print(f'f1={F[0]}\nf2={F[1]}\nf3={F[2]}\nf4={F[3]}')

    ###########################################################################
    # Ploting the result
    #
    fig, ax = plt.subplots(1, 1, figsize=(10, 5))

    x = [parameters.findingf_L1, parameters.findingf_L1 + parameters.findingf_L2, parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4]

    plt.plot(x, F, 'o')

    # plt.suptitle('Focal lengths for each lens', fontsize=10)

    ax.set_xlim(-0.2, 3.5)
    ax.set_ylim(-np.amax([abs(ele) for ele in F]) - 0.2, np.amax([abs(ele) for ele in F]) + 0.2)

    grid_x_ticks = np.arange(-0.2, 3.5, 0.1)
    grid_y_ticks = np.arange(-np.amax([abs(ele) for ele in F]) - 0.2, np.amax([abs(ele) for ele in F]) + 0.2, 0.5)

    ax.set_xticks(grid_x_ticks, minor=True)
    ax.set_yticks(grid_y_ticks, minor=True)

    ax.set_ylabel("Focal lengths (m)")
    ax.set_xlabel("Distance between the lens and the entrance (m)")
    extra1 = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)
    extra2 = Rectangle((0, 0), 1, 1, fc="w", fill=False, edgecolor='none', linewidth=0)

    plt.legend([extra1, extra2], (fr'$\alpha_i= {parameters.findingf_alpha:.4f}$', fr'$\alpha_o= {alpha4[min4]:.4f}$'))

    plt.grid(which='both')
    plt.grid(which='minor', alpha=0.3, linestyle='--')

    ###########################################################################
    # Annotate
    #
    bbox = dict(boxstyle="round", fc="0.8")

    ax.annotate('$1^{st}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1, F[0]), (parameters.findingf_L1, F[0]), xytext=(20, 20), textcoords='offset points', bbox=bbox)
    ax.annotate('$2^{nd}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1 + parameters.findingf_L2, F[1]), (parameters.findingf_L1 + parameters.findingf_L2, F[1]), xytext=(20, 20), textcoords='offset points', bbox=bbox)
    ax.annotate('$3^{rd}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, F[2]), (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, F[2]), xytext=(20, 20), textcoords='offset points', bbox=bbox)
    ax.annotate('$4^{th}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, F[3]), (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, F[3]), xytext=(20, 20), textcoords='offset points', bbox=bbox)

    # ax.axvline(x=parameters.findingf_L1, c='k', lw=1)
    # ax.axvline(x=parameters.findingf_L1 + parameters.findingf_L2, c='k', lw=1)
    ax.axvline(x=0, c='red', lw=2)
    ax.axvline(x=1.7, c='red', lw=2)
    ax.axvline(x=3.4, c='red', lw=2)
    # ax.axvline(x=parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, c='k', lw=1)
    # ax.axvline(x=parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, c='k', lw=1)

    ax.annotate("", xy=(parameters.findingf_L1, F[0]), xytext=(parameters.findingf_L1, -F[0]), arrowprops=dict(arrowstyle="<->"))
    ax.annotate("", xy=(parameters.findingf_L1 + parameters.findingf_L2, F[1]), xytext=(parameters.findingf_L1 + parameters.findingf_L2, -F[1]), arrowprops=dict(arrowstyle="<->"))
    ax.annotate("", xy=(parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, F[2]), xytext=(parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, -F[2]), arrowprops=dict(arrowstyle="<->"))
    ax.annotate("", xy=(parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, F[3]), xytext=(parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, -F[3]), arrowprops=dict(arrowstyle="<->"))

    fig.savefig('pics/finding_f.png', bbox_inches='tight', dpi=100)
    plt.show()

    ####################################################################
    # Visualisation

    print(alpha1[min1], alpha2[min2], alpha3[min3])

    ELLIPSES = {

        'Input signal': sigma_matrices.Input(parameters.findingf_alpha, parameters.findingf_beta, parameters_calculation.gamma(parameters.findingf_alpha, parameters.findingf_beta, parameters.findingf_eps), parameters.findingf_eps),

        'In the center': sigma_matrices.Einzel(L2, parameters_calculation.f(L2, to_check[min2], parameters.findingf_V, parameters.findingf_R), L2, alpha1[min1], b, g, eps),

        'After fourth lens': sigma_matrices.Einzel(L4, parameters_calculation.f(L4, to_check[min4], parameters.findingf_V, parameters.findingf_R), L4, alpha3[min3], b, g, eps),

    }

    fig, axs = plt.subplots(1, 3, figsize=(18, 18))
    for ax, (title, dependency) in zip(axs, ELLIPSES.items()):
        x, y = confidence_ellipse.get_correlated_dataset(100000, dependency, parameters.mu, parameters.scale, label='Dataset')
        ax.scatter(x, y, s=0.5)

        ax.axvline(c='grey', lw=1)
        ax.axhline(c='grey', lw=1)
        ax.hexbin(x, y, gridsize=300, cmap='inferno')
        confidence_ellipse.confidence_ellipse(x, y, ax, edgecolor='red', linewidth=2)

        ax.scatter(parameters.mu[0], parameters.mu[1], c='red', s=3)
        ax.set_title(title)
        ax.set_xlim((-3, 3))
        ax.set_ylim((-1, 1))
        ax.set_xlabel('X/Y (mm)')
        ax.set_ylabel("X'/Y' (mrad)")

    fig.savefig('pics/finding_f_beams.png', dpi=100)

    plt.show()


find_f(parameters.findingf_L1, parameters.findingf_L2, parameters.findingf_L3, parameters.findingf_L4, parameters.findingf_alpha, parameters.findingf_beta, parameters_calculation.gamma(parameters.findingf_alpha, parameters.findingf_beta, parameters.findingf_eps), parameters.findingf_eps, parameters.findingf_V, parameters.findingf_R)
