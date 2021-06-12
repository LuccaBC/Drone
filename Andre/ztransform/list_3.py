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
import control
from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

############
# montando o tf
##########
num = [3, 3]
den = [1, 4, 6, 4]
Gcont = tf(num, den)

print(Gcont)

zpk_cont = tf2zpk(num, den)

print(zpk_cont)

poles = pole(Gcont)
zeros = zero(Gcont)
k = control.dcgain(Gcont)

print(poles, zeros, k)

pp_cont, zz = pzmap(Gcont, plot=True, grid=True)
zeta = damp(Gcont, doprint=True)
plt.show()

Ts = 0.1
a_back = (1 + 4*Ts + 6*Ts**2 + 4*Ts**3)
num_back = [(3*Ts ** 2) * (1+Ts)/a_back, -(3*Ts ** 2)/a_back, 0, 0]
den_back = [1, -(3+8*Ts+6*Ts**2)/a_back, (3+4*Ts)/a_back, -1/a_back]

G_back = tf(num_back, den_back, Ts)
G_back_d = c2d(Gcont, Ts, method='backward_diff')

print(G_back, G_back_d)
pp, zz = pzmap(G_back_d, plot=True, grid=True)
print(pp)
compare_pp = np.exp(pp_cont*Ts)
print(compare_pp)

num_forw = [(3*Ts ** 2), (3*Ts ** 2)*(Ts-1)]
den_forw = [1, (4*Ts-3), (3-8*Ts+6*Ts**2), 4*Ts**3-1+4*Ts-6*Ts**2]

G_forw = tf(num_forw, den_forw, Ts)
G_forw_d = c2d(Gcont, Ts, method='forward_diff')

print(G_forw, G_forw_d)
pp, zz = pzmap(G_forw_d, plot=True, grid=True)
print(pp)
compare_pp = np.exp(pp_cont*Ts)
print(compare_pp)


a_tust = 8+16*Ts+12*Ts**2+4*Ts**3
Kt = 3*Ts**2/a_tust
a2 = 12*Ts**2+12*Ts**3-16*Ts-24
a1 = 24-16*Ts-12*Ts**2+12*Ts**3
a0 = 16*Ts-8-12*Ts**2+4*Ts**3

num_tust = [Kt*(2+Ts), Kt*(2+3*Ts), Kt*(3*Ts-2), Kt*(Ts-2)]
den_tust = [1, a2/a_tust, a1/a_tust, a0/a_tust]

G_tust = tf(num_tust, den_tust, Ts)
G_tust_d = c2d(Gcont, Ts, method='tustin')

print(G_tust, G_tust_d)
pp, zz = pzmap(G_tust_d, plot=True, grid=True)
print(pp)
compare_pp = np.exp(pp_cont*Ts)
print(compare_pp)




# zpk_back = tf2zpk(num_back, den_back)
pp, zz = pzmap(G_back, plot=True, grid=True)

# print(zpk_back)




#
# Ts = 0.1
#
# impulse_cont = tf(1, 1)
# step_tf = tf(1, [1, 0])
# ramp_tf = tf(1, [1, 0, 0])
# exp_tf = tf(1, [1, 1])
# sin_tf = tf(1, [1, 0, 1])
#
# time_final = 2
# n_points = round(time_final / Ts) + 1
#
# time = linspace(0, time_final, n_points)
#
#
#
# signals_tf = [impulse_tf, step_tf, ramp_tf, exp_tf]
# fig, axs = plt.subplots(2, 2)
#
# for signal in signals_tf:
#     idx = signals_tf.index(signal)
#     y_c, t_c = impulse(signal, time, 0)
#     axs[int(np.floor(idx/2)), idx % 2].stem(t_c, y_c)
# plt.show()

