################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 07.03.2021
#
#   Info: Function to show the start value of a advanced signal.
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

T_big = 0.1 #what I will say it is de discretization time

time_final = 2

n_big = round(time_final / T_big) + 1

time_big = linspace(0, time_final, n_big)

step_tf = tf(1, [1, 0])
exp_tf = tf(1, [1, 1])

z_step = c2d(step_tf*1/T_big, T_big, method='impulse')

z_exp = c2d(exp_tf*1/T_big, T_big, method='impulse')

fig, axs = plt.subplots()

y_c, t_c = impulse(z_step, time_big)
axs.step(time_big-0.05, y_c, color=(0.6, 0, 0, 0.5), label="avanço continuo", linewidth=3.0, where='post')
axs.stem(time_big, y_c, 'bo')
plt.xlim(left=-0.1)
plt.xticks([-0.1, 0, 0.5, 1, 1.5, 2])
plt.show()

fig2, axs = plt.subplots()

y_c, t_c = impulse(z_exp, time_big)
axs.plot(time_big-0.05, y_c, color=(0.6, 0, 0, 0.5), label="avanço continuo", linewidth=3.0)
axs.stem(time_big, y_c, 'bo')
plt.xlim(left=-0.1)
plt.xticks([-0.1, 0, 0.5, 1, 1.5, 2])
plt.show()


##########################
T_disc = 0.1

Cz = tf([8.1, 8.1*1.8138], [1, 8.707])
c_f = c2d(Cz, T_disc, method='forward_diff')
c_b = c2d(Cz, T_disc, method='backward_diff')
c_t = c2d(Cz, T_disc, method='tustin')
c_w = c2d(Cz, T_disc, method='tustin', prewarp_frequency=50)
c_m = c2d(Cz, T_disc, method='matched')

print(c_f)
print(c_b)
print(c_t)
print(c_m)
controllers = [c_f, c_b, c_t, c_m]
colors = [(1, 0, 0, 0.5), (0, 1, 0, 0.5), (0, 0, 1, 0.5), (0, 0.5, 0.3, 0.5)]
labels = ['forward', 'backward', 'tustin', 'matched']
plant = tf([1], [1, 1, 0])
plant_d = c2d(plant, T_disc, method='zoh')
time = linspace(0, 5, 51)
fig3, axs = plt.subplots()
for i, c_ds in enumerate(controllers):
    closed = feedback(c_ds*plant_d, 1)
    y, t = step(closed, time)
    plt.step(t, y, color=colors[i], label=labels[i], linewidth=3.0, where='post')
cont_closed = feedback(Cz*plant, 1)
y, t = step(cont_closed, time)
plt.plot(t, y, color=[0.5, 0.5, 0.5, 0.5], label='continuo', linewidth=3.0)
plt.legend(loc='best', shadow=True, fontsize='x-large')
plt.show()
