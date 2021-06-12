################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 24.02.2021
#
#   Info: Function to show comparison between Forward s approximations
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

# Compensador D(s) = Y/U, s-1/s^2+3s+2
Ts = 0.1
time_final = 8
Npoints = round(time_final / Ts) + 1
sys = tf([1, -1], [1, 3, 2])

s5 = c2d(sys, Ts, method='tustin')
s7 = c2d(sys, Ts, method='forward_diff')
s8 = c2d(sys, Ts, method='backward_diff')

time = np.linspace(0, time_final, Npoints)
y_cont, time = step(sys, time, 0)

control_action = np.ones(Npoints)
y_forward = np.zeros(Npoints)
y_backward = np.zeros(Npoints)

y5, t5 = step(s5, time, 0)
y7, t7 = step(s7, time, 0)
y8, t8 = step(s8, time, 0)

counter = 0
for t in time:
    if counter == 0:
        y_forward[counter] = 0 #(1 / (1 + 2 * Ts)) * (Ts * control_action[counter])
        y_backward[counter] = (1 / (1 + 3 * Ts + 2 * Ts ** 2)) * ((Ts - Ts ** 2) * control_action[counter])
    elif counter == 1:
        y_forward[counter] = Ts * control_action[counter-1] + (2 - 3*Ts)*y_forward[counter - 1]
        y_backward[counter] = (1 / (1 + 3 * Ts + 2 * Ts ** 2)) * (
                (Ts - Ts ** 2) * control_action[counter] - Ts * control_action[counter - 1]
                + (2 + 3 * Ts) * y_backward[counter - 1])
    else:
        y_forward[counter] = (1 / (1 + 2 * Ts)) * (
                Ts * control_action[counter] - (Ts + Ts ** 2) * control_action[counter - 1]
                - (3 * Ts ** 2 - 2 * Ts - 2) * y_forward[counter - 1] - y_forward[counter - 2])
        y_backward[counter] = (1 / (1 + 3 * Ts + 2 * Ts ** 2)) * (
                (Ts - Ts ** 2) * control_action[counter] - Ts * control_action[counter - 1]
                + (2 + 3 * Ts) * y_backward[counter - 1] - y_backward[counter - 2])
        # y_centered[counter] = 2 * Ts * control_action[counter - 1] - 4 * Ts ** 2 * control_action[counter - 2]
        # - 2 * Ts * control_action[counter - 3] - 4 * Ts * y_centered[counter - 1] - (12 * Ts ** 2 - 2) * y_centered[
        #     counter - 2] + 4*Ts*y_centered[counter-3] - y_centered[counter-4]

    counter = counter + 1

fig, ax = plt.subplots()
ax.plot(time, y_cont, color=[0.5, 0.5, 0.5, 1], label="Contínuo", linewidth=3.0)
ax.step(time, y_forward, 'r', label="Para Frente", linewidth=3.0, where='post')
ax.step(time, y_backward+0.01, 'g', label="Para Trás", linewidth=3.0, where='post')
ax.step(time, y7, color=(0, 0, 1, 0.5), label="Para frente s7", linewidth=3.0, where='post')
ax.step(time, y8, color=(1, 0, 0, 0.5), label="Para Trás s8", linewidth=3.0, where='post')
# ax.step(time, y5, color=(0, 1, 0, 0.5), label="Tustin s5", linewidth=3.0, where='post')
legend = ax.legend(loc='best', shadow=True, fontsize='x-large')
plt.xlabel('Tempo [s]')
plt.ylabel('Amplitude')
plt.show(block=True)
