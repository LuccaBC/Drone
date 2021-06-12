################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 12.02.2021
#
#   Info: Function to show how Quantization works.
#   RoadMap: Use PyQT to create a GUI interface
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

RMSE = []
bits = [2, 4, 8, 10, 12, 16]
freqs = [1, 2, 10, 30]
graph = ['m--', 'b.-', 'r:', 'g']
label = [' 1 Hz', ' 3 Hz', '10 Hz', '30 Hz']
fig, ax = plt.subplots()
for freq_disc in freqs:
    for N in bits:
        t_amos = linspace(0.0, 1, num=freq_disc + 1)
        y_cont = 2.5 * np.sin(2 * np.pi * t_amos) + 2.5
        dv = 5 / (2 ** N - 1)
        y_amos = np.round((2.5 * np.sin(2 * np.pi * t_amos) + 2.5) / dv) * dv
        error = y_amos - y_cont
        quadError = error ** 2
        sqSum = np.sqrt(np.sum(quadError)/(freq_disc+1))
        RMSE.append(sqSum)
        print("F: %d, N: %d, RMSE = %f" % (freq_disc, N, sqSum))
        # plots
        # fig2, ax2 = plt.subplots()
        # ax2.step(t_amos, y_amos, 'o:', color=(0, 0.6, 0, 1), where='post', linewidth=2.0, markersize=12)
        # ax2.plot(t_amos, y_cont, 'b-', linewidth=3.0)
        # plt.xlabel('Tempo [s]')
        # plt.ylabel('Tensão [V]')
        # plt.show()
    plt.semilogy(bits, RMSE[6*freqs.index(freq_disc):6*(freqs.index(freq_disc)+1)],
             graph[freqs.index(freq_disc)], label=label[freqs.index(freq_disc)])

legend = ax.legend(loc='best', shadow=True, fontsize='x-large')
plt.xlabel('N [bits]')
plt.ylabel('RMSE')
plt.show()







