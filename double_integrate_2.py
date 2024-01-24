import numpy as np
from matplotlib import pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.integrate import quad


t_array = [(10 ** i) for i in np.arange(-1, 10, 0.1)]
array_ = [3e19, 1e19, 3e18, 1e18, 3e17, 1e17, 3e16, 1e16, 3e15, 1e15]
na_array = [i * ((2 * 1.2E-7) ** 3) for i in array_]
PL = np.zeros(len(t_array))
pdf = PdfPages('Figure_2.pdf')

fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot()

with open('outputXY2.txt', 'w') as file:
    pass

for j in range(0, len(na_array)):
    for i in range(0, len(t_array)):
        na = na_array[j]
        t = t_array[i]

        def func(energy):
            return np.exp(-4 / energy) * np.exp(-t * np.exp(-4 / energy)) / energy ** 4

        def integrate():
            return quad(func, 0, np.inf)[0]

        def func_2(energy):
            return (np.exp(-t*np.exp(-4/energy))-1)/energy**4

        def integrate_2():
            return quad(func_2, 0, np.inf)[0]


        PL[i] = (4 * np.pi * na) * integrate() * np.exp(4 * np.pi * na * integrate_2())

    count = 0
    for x, y in zip(t_array, PL):
        print(f't{count} = {x} | PL{count} = {y}')
        count += 1
    print("-----------------")

    with open('outputXY2.txt', 'a') as file:
        file.write('----------------------- \n')
        file.write(f'Number of curve: {j + 1}\n')
        file.write('X[i] values: \n')
        file.write('----------------------- \n')

        for element in t_array:
            file.write(str(element))
            file.write('\n')

        file.write('----------------------- \n')
        file.write('Y[i] values: \n')
        file.write('----------------------- \n')

        for element in PL:
            file.write(str(element))
            file.write('\n')

    plt.ylim(1e-4, 1e6)
    plt.loglog([i / 1e6 for i in t_array], [i*1e6 for i in PL], label=f'Nd={na_array[j]:.1e}')
    plt.xlabel('t, 1/p0', fontsize=15, color='black')
    plt.ylabel('Intensity, p0', fontsize=15, color='black')
    ax.legend(loc='upper right')


pdf.savefig()
pdf.close()
ax.grid()
plt.show()
