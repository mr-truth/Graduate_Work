import numpy as np
from matplotlib import pyplot as plt
from scipy.integrate import quad
import gmpy2
from gmpy2 import mpfr

gmpy2.get_context().precision = 50



Nd = mpfr('0.01')
compensation = mpfr('0.01')
pump = 0.0625
E_md = 2 * compensation * (Nd ** (2 / 3))
E_array = [i for i in np.arange(-4, 8.1, 0.1)]
t_array = [0, 10, 10**2, 10**3, 10**4, 10**5]
PL = np.zeros(len(E_array))

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot()

for j in range(0, len(t_array)):
    for i in range(0, len(E_array)):

        E = round(E_array[i], 2)
        t = t_array[j]

        def integrate():
            func = lambda energy_1:(energy_1**-5) * (E_md/((4*(E_md/energy_1)**2) + (E - energy_1)**2)) * gmpy2.exp(-4 / energy_1 - 32 * np.pi / 3 * Nd / energy_1 ** 3 - t*gmpy2.exp(-4/energy_1))

            result = quad(func, 0, np.inf, epsabs=1.5e-100)[0]

            result = result * (64 * Nd * pump)/(32 * np.pi * Nd)

            return result


        def integrate_q():
            func = lambda energy: (gmpy2.exp(-t * gmpy2.exp(-4/energy)) - 1)/energy**4

            result =  gmpy2.exp(4 * np.pi * Nd * quad(func, 0, np.inf, epsabs=1.5e-100)[0])

            return result

        a = integrate() * integrate_q()

        PL[i] = a


        #print(E)
        print(PL[i])
    print('-----------------')

    plt.semilogy(E_array, PL)

ax.grid()
plt.show()




