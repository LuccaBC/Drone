################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 07.03.2021
#
#   Info: Exercise at the end od class 4
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np


Ts = 1
time_final = 1000
n_points = round(time_final / Ts) + 1

time = linspace(0, time_final, n_points)

g_z = tf([0.25, 0, 0], [1, 0, -1.25, 0.25], Ts)
y, t = impulse(g_z, time, 0)

print(g_z)

fig, ax = plt.subplots()
ax.stem(time, y, 'rs')

y_med = np.median(y)


# xdata = [time[0], time[1]]
# ydata = [y[0], y[1]]
#
# bbox = dict(boxstyle="round", fc="0.8")
# arrowprops = dict(
#     arrowstyle="->",
#     connectionstyle="angle,angleA=0,angleB=-90,rad=10")
#
# offset = 10
# ax.annotate(
#     f'data = ({xdata[1]:.2f}, {ydata[1]:.2f})',
#     (xdata[1], ydata[1]),
#     xytext=(2*offset, 3*offset), textcoords='offset points',
#     bbox=bbox, arrowprops=arrowprops)
# ax.annotate(
#     f'data = ({xdata[0]:.2f}, {ydata[0]:.2f})',
#     (xdata[0], ydata[0]),
#     xytext=(2*offset, -3*offset), textcoords='offset points',
#     bbox=bbox, arrowprops=arrowprops)
plt.show()

