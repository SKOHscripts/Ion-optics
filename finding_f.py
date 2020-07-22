'''
@file : finding_V.py
@brief : Einzel lenses are approximated by thin lenses with drifts. For each lens, we look for the focal length to be applied to obtain a parallel beam at the output, knowing that we emit a parallel beam at the input. We will therefore pass through each lens successively and note the values in a table.  (x,x') plane and (y,y') plane are the same.  

@author : Corentin MICHEL
creation : 20/07/2020
'''

import numpy as np
import matplotlib.pyplot as plt

import parameters
import parameters_calculation
import sigma_matrices


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
        alpha3.append(sigma_matrices.Einzel(L1, f3, L1, alpha2[min2], b, g, eps)[0][1])

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
        alpha4.append(sigma_matrices.Einzel(L1, f4, L1, alpha3[min3], b, g, eps)[0][1])

    for i in range(len(alpha4)):
        if abs(alpha4[i]) == np.amin([abs(ele) for ele in alpha4]):

            min4 = i
            F[3] = parameters_calculation.f(L4, to_check[min4], parameters.findingf_V, parameters.findingf_R)

            print(f'4e : a= {alpha4[i]} : f= {F[3]}')

    print(f'f1={F[0]}\nf2={F[1]}\nf3={F[2]}\nf4={F[3]}')

    ###########################################################################
    # Ploting the result
    #
    fig, ax = plt.subplots(1, 1, figsize=(12, 14))

    y = [parameters.findingf_L1, parameters.findingf_L1 + parameters.findingf_L2, parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4]

    plt.plot(y, F, 'o')

    plt.suptitle('Focal lengths for each lens', fontsize=20)

    grid_x_ticks = np.arange(0, 3.5, 0.1)
    grid_y_ticks = np.arange(0, 3.5, 0.1)

    ax.set_xticks(grid_x_ticks, minor=True)
    ax.set_yticks(grid_y_ticks, minor=True)

    plt.grid(which='both')
    plt.grid(which='minor', alpha=0.3, linestyle='--')

    ###########################################################################
    # Annotate
    #
    bbox = dict(boxstyle="round", fc="0.8")
    arrowprops = dict(arrowstyle="->", connectionstyle="angle,angleA=0,angleB=90,rad=10")

    ax.annotate('$1^{st}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1, F[0]), (parameters.findingf_L1, F[0]), xytext=(20, 20), textcoords='offset points', bbox=bbox, arrowprops=arrowprops)
    ax.annotate('$2^{nd}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1 + parameters.findingf_L2, F[1]), (parameters.findingf_L1 + parameters.findingf_L2, F[1]), xytext=(20, 20), textcoords='offset points', bbox=bbox, arrowprops=arrowprops)
    ax.annotate('$3^{rd}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, F[2]), (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3, F[2]), xytext=(20, 20), textcoords='offset points', bbox=bbox, arrowprops=arrowprops)
    ax.annotate('$4^{th}$ Einzel lens :\n(D= %.2f m, f= %.4f m)' % (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, F[3]), (parameters.findingf_L1 + parameters.findingf_L2 + parameters.findingf_L3 + parameters.findingf_L4, F[3]), xytext=(-7 * 20, 20), textcoords='offset points', bbox=bbox, arrowprops=arrowprops)

    fig.savefig('pics/finding_f.png', bbox_inches='tight', dpi=100)
    plt.show()

    return F


find_f(parameters.findingf_L1, parameters.findingf_L2, parameters.findingf_L3, parameters.findingf_L4, parameters.findingf_alpha, parameters.findingf_beta, parameters_calculation.gamma(parameters.findingf_alpha, parameters.findingf_beta, parameters.findingf_eps), parameters.findingf_eps, parameters.findingf_V, parameters.findingf_R)
