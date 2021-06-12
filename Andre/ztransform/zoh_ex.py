###############################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 31.03.2021
#
#   Info: Function to show Zero-Order Holder
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np


Tsamples = [0.01, 0.16, 0.32, 1]

time_final = 8
g_s = tf(1, [1, 1])
fig, axs = plt.subplots(2, 2, sharex='col', sharey='row')
n_cont = round(time_final / 0.005) + 1
t_cont = linspace(0, time_final, n_cont)
y_cont, t_c = step(g_s, t_cont)

for Ts in Tsamples:
    idx = Tsamples.index(Ts)
    n_points = round(time_final / Ts) + 1
    time = linspace(0, time_final, n_points)
    g_z = c2d(g_s, Ts, method='zoh')
    print(g_z)
    y_dis, t_dis = step(g_z, time)
    axs[int(np.floor(idx / 2)), idx % 2].plot(t_cont, y_cont, color=[0.5, 0.5, 0.5, 1], label="Contínuo", linewidth=3.0)
    axs[int(np.floor(idx / 2)), idx % 2].step(time, y_dis, color=(0, 0.8, 0, 0.5), label=Ts.__str__(), linewidth=3.0, where='post')
    legend = axs[int(np.floor(idx / 2)), idx % 2].legend(loc='best', shadow=True, fontsize='x-large')
fig.text(0.5, 0.04, 'Tempo [s]', ha='center', va='center')
fig.text(0.06, 0.5, 'Amplitude', ha='center', va='center', rotation='vertical')
plt.show()

