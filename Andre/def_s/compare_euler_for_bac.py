################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 24.02.2021
#
#   Info: Function to show comparison between Forward Euler and Backward Euler.
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

#Compensador D(s) = U/E = ko(s+a)/(s+b)
#avanço de fase
ko = 3
a = 1
b = 2
T = 0.1
time_final = 3
Npoints = round(time_final/T) + 1
sys = tf([ko, a*ko], [1, b])

time = np.linspace(0, time_final, Npoints)

y_cont, time = step(sys, time, 0)
u_forward = np.zeros(Npoints)
u_backward = np.zeros(Npoints)
u_centered = np.zeros(Npoints)
first = True

for t in time:
    if first:
        first = False
        u_forward[0] = ko
        u_backward[0] = (1/(1+b*T))*(ko*(1+a*T))
        counter = 1
    else:
        u_forward[counter] = (1-b*T)*u_forward[counter-1] + ko*a*T
        u_backward[counter] = (1/(1+b*T))*(u_backward[counter-1] + ko*a*T)
        counter = counter + 1


fig, ax = plt.subplots()
plt.plot(time, y_cont, color=[0.5, 0.5, 0.5, 1], label="Contínuo", linewidth=3.0)
plt.plot(time, u_forward, 'r', label="Para Frente", linewidth=3.0)
plt.plot(time, u_backward, 'g', label="Para Trás", linewidth=3.0)

legend = ax.legend(loc='best', shadow=True, fontsize='x-large')
plt.xlabel('Tempo [s]')
plt.ylabel('Amplitude')
plt.show(block=True)
