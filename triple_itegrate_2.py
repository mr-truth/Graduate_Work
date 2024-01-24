import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad

W_max = 1
Ed = 1
pump = 0.0625
Nd_array = [10**-2, 10**-3, 10**-4]
t_array = [0, 10, 10**2, 10**3]
E_array = [i for i in np.arange(-2, 6.1, 0.1)]
PL = np.zeros(len(E_array))


for k in range(0, len(Nd_array)):
    Nd = Nd_array[k]
    fig = plt.figure(figsize=(9, 9))
    ax = fig.add_subplot()
    ax.grid()
    for j in range(0, len(t_array)):
        for i in range(0, len(E_array)):
            E = E_array[i]
            t = t_array[j]
            mu = 2 * 0.99 * (Nd ** (1 / 3))

            def integral_1(energy_1):

                def integral(energy_2):
                    function_1 = lambda E: np.exp(-(W_max * np.exp(-4 * Ed / E)) * t - 1) * (E ** (-4))
                    res = np.exp(4 * np.pi * Nd ** 3 * quad(function_1, 0, energy_2)[0])
                    return res

                if ((-3 * energy_1 ** 2) + 2 * energy_1 * E + E ** 2) > 0:
                    e_min = max(energy_1, 0.5 * (((-3 * energy_1 ** 2) + 2 * energy_1 * E + E ** 2) ** 0.5 - energy_1 + E),
                                E - mu)
                else:
                    e_min = energy_1

                e_max = 0.5 * ((5 * energy_1 ** 2 - 2 * energy_1 * E + E ** 2) ** 0.5 + E + energy_1)

                function_2 = lambda energy_2: energy_2 ** (-3) * (energy_1 + energy_2 - E) ** (-3) * np.exp(
                    -32 * np.pi / 3 * Nd * energy_2 ** (-3))

                if e_min < e_max:
                    result = quad(function_2, e_min, e_max)[0] * integral(1)
                else:
                    result = 0.0

                result = result * energy_1 ** (-3) * np.exp(-4 / energy_1 - t * np.exp(-4 / energy_1))

                return result


            a = quad(integral_1, 0, np.inf)[0]

            PL[i] = 512* (1 - pump) * (np.pi ** 2) * (Nd ** 2)/(32 * np.pi * Nd) * a
            print(E)
            print(PL[i])
            print(t)
            print(Nd)


        plt.title(f'Nd = {Nd_array[k]:.2e}', fontsize=34)
        plt.xlabel('E, Ed', fontsize=28, color='black')
        plt.ylabel('P(E, t)/P0', fontsize=28, color='black')
        plt.xlim(-2,6)
        plt.ylim(10**-6, 0.01)
        plt.semilogy(E_array, PL, label=f't={t_array[j]}')
        plt.tick_params(axis='both', which='major', labelsize=34)
        ax.legend(fontsize=34)


plt.show()
