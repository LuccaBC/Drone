################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 18.02.2021
#
#   Info: Function to show magnitude and phase from Zero-Order Holder
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

T = 0.01
w = np.pi*2*np.linspace(0.00001, 250, num=3000)

magnitude = np.abs(np.sin(0.5*T*w)/(0.5*T*w))
magdb = mag2db(magnitude)
cross_wt = np.where(np.abs(magdb+3) <= 0.007)

plt.plot(w*T, magdb, 'b-', linewidth=3.0)
# plt.plot([T*w[0], T*w[-1]], [-3, -3], ':', color=(0.5, 0.5, 0.5, 0.5))
# plt.plot([T*w[cross_wt[0]], T*w[cross_wt[0]]], [np.max(magdb), np.min(magdb)], ':', color=(0.5, 0.5, 0.5, 0.5))
plt.xlabel('wT [rad/s]')
plt.ylabel('Amplitude [dB]')
plt.show()

phase = -0.5*w*T
warp = np.where(phase < -np.pi)
warp2 = np.where(phase < -2*np.pi)
phase[warp2] = phase[warp2] + np.pi
phase[warp] = phase[warp] + np.pi
plt.plot(w*T, phase, 'b-', linewidth=3.0)
# plt.plot([T*w[0], T*w[-1]], [-3, -3], ':', color=(0.5, 0.5, 0.5, 0.5))
# plt.plot([T*w[cross_wt[0]], T*w[cross_wt[0]]], [np.max(magdb), np.min(magdb)], ':', color=(0.5, 0.5, 0.5, 0.5))
plt.xlabel('wT [rad/s]')
plt.ylabel('Phase [rad]')
plt.show()
#####################################
#sinal 3sin(20pit)+1.5
#T=0.01
tempo = np.linspace(0, 0.2, num=4000)
y = 3*np.sin(20*np.pi*tempo) + 3

t_amos = linspace(0.0, 0.2, num=21)
y_amos = 3*np.sin(20*np.pi*t_amos) + 3
y_1_harm = 3*(np.sin(20*np.pi*T/2)/(20*np.pi*T/2))*np.sin(20*np.pi*(tempo-T/2)) + 3
# plots
fig, ax = plt.subplots()
ax.step(t_amos, y_amos, 'o:', color=(0, 0.6, 0, 1), where='post', linewidth=2.0, markersize=12, label='ZOH')
ax.plot(tempo, y, 'b-', linewidth=3.0, label='sinal')
ax.plot(tempo, y_1_harm, 'r:', linewidth=3.0, label='Primeiro Harmônico')
plt.xlabel('Tempo [s]')
plt.ylabel('Amplitude')

legend = ax.legend(loc='best', shadow=True, fontsize='x-large')
plt.show()