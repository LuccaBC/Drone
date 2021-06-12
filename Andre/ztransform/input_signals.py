################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 07.03.2021
#
#   Info: Function to show the basic input signals_tf
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np


Ts = 0.1
time_final = 2
n_points = round(time_final / Ts) + 1

time = linspace(0, time_final, n_points)

delay = tf('z')
impulse_tf = (delay**-1)*tf(1, 1)
step_tf = tf(1, [1, 0])
ramp_tf = tf(1, [1, 0, 0])
exp_tf = tf(1, [1, 1])

signals_tf = [impulse_tf, step_tf, ramp_tf, exp_tf]
fig, axs = plt.subplots(2, 2)

for signal in signals_tf:
    idx = signals_tf.index(signal)
    y_c, t_c = impulse(signal, time, 0)
    axs[int(np.floor(idx/2)), idx % 2].stem(t_c, y_c)
plt.show()

