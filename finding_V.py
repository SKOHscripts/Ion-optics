import numpy as np
import parameters
import parameters_calculation
import sigma_matrices


def find_V(L1, L2, L3, L4, a, b, g, eps, V, R):
    V = np.zeros((4, 2))
    to_check = np.linspace(0, 1800, 900)
    alpha1_x = []
    print(alpha1_x)
    for i in range(len(to_check)):
        f1 = parameters_calculation.f(L1, to_check[i], parameters.findingV_V, parameters.findingV_R)
        alpha1_x.append(sigma_matrices.Einzel(L1, f1, L1, a, b, g, eps)[0][1])
        if abs(alpha1_x[i]) <= parameters.findingV_accuracy:
            V[0][0] = to_check[i]

    alpha2_x = []
    for i in range(len(to_check)):
        f2 = parameters_calculation.f(L2, to_check[i], parameters.findingV_V, parameters.findingV_R)
        alpha2_x[i] = sigma_matrices.Einzel(L2, f2, L2, alpha1_x[-1], b, g, eps)[0][1]
        if abs(alpha2_x[i]) <= parameters.findingV_accuracy:
            V[1][0] = to_check[i]

    aplha3_x = []
    for i in range(len(to_check)):
        f3 = parameters_calculation.f(L3, to_check[i], parameters.findingV_V, parameters.findingV_R)
        aplha3_x[i] = sigma_matrices.Einzel(L3, f3, L3, alpha2_x[-1], b, g, eps)[0][1]
        if abs(aplha3_x[i]) <= parameters.findingV_accuracy:
            V[3][0] = to_check[i]

    aplha4_x = []
    for i in range(len(to_check)):
        f4 = parameters_calculation.f(L4, to_check[i], parameters.findingV_V, parameters.findingV_R)
        aplha4_x[i] = sigma_matrices.Einzel(L4, f4, L4, alpha3_x[-1], b, g, eps)[0][1]
        if abs(aplha4_x[i]) <= parameters.findingV_accuracy:
            V[4][0] = to_check[i]

    return V


find_V(0.850, 0.850, 0.850, 0.850, parameters.input_alpha, parameters.input_beta, parameters.input_gamma, parameters.epsilon, parameters.findingV_V, parameters.findingV_R)
