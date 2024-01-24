import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

a = 1
Nd = 0.01
W_max = 1
Ed = 1
t_array = np.geomspace(1e0, 1e5, 15)
E_array = [i for i in np.arange(0.1, 5, 0.01)]
PL = np.zeros(len(E_array))

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot()

with open('outputE, Pl, t.txt', 'w') as file:
    pass

for j in range(len(t_array)):
    for i in range(len(E_array)):
        E = E_array[i]
        t = t_array[j]
        def function_1(r):
            return (np.exp(-W_max * np.exp(-2 * r) * t) - 1) * (r ** 2)

        def integrate_q(n):
            return np.exp(4 * np.pi * n * (a ** 3) * quad(function_1, 0, np.inf)[0])


        def function_2():
            return 32 * np.pi * W_max * Nd * (a ** 3) * ((Ed / E) ** 4) * np.exp(
                (-4 * Ed / E) - (W_max * t * np.exp(-4 * Ed / E)))


        def intensity():
            return function_2() * integrate_q(Nd)

        PL[i] = intensity()

    count = 0
    for x, y in zip(E_array, PL):
        print(f'E{count} = {x} | PL{count} = {y}')
        count += 1
    print("-----------------")

    with open('outputE, Pl, t.txt', 'a') as file:
        file.write('----------------------- \n')
        file.write(f'Number of curve: {j + 1}\n')
        file.write('X[i] values: \n')
        file.write('----------------------- \n')

        for element in E_array:
            file.write(str(element))
            file.write('\n')

        file.write('----------------------- \n')
        file.write('Y[i] values: \n')
        file.write('----------------------- \n')

        for element in PL:
            file.write(str(element))
            file.write('\n')

    plt.xlim(0,4)
    plt.ylim(1e-8, 1e-1)
    plt.xlabel('Energy, Ed', fontsize=15, color='black')
    plt.ylabel('Intensity, p0', fontsize=15, color='black')
    plt.semilogy(E_array, PL, label=f't={t_array[j]:.1e}')
    plt.tick_params(axis='both', which='major', labelsize=16)
    ax.legend(loc='upper right')

ax.grid()
plt.show()
