################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 11.05.2021
#
#   Info: State space model 2
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np

G = tf([57.852, 52.07], [1, -0.81058, 0], 0.01)
print(G)
sys_ss = tf2ss(G)
print(sys_ss)
test_ss, T_transform = ctrl.reachable_form(sys_ss)
Ac = np.round(test_ss.A, decimals=2)
Bc = np.round(test_ss.B, decimals=2)
Cc = np.round(test_ss.C, decimals=2)
Dc = np.round(test_ss.D, decimals=2)

print(Ac, Bc, Cc, Dc)

######################
#testar o step
tfinal = 1
Ts = 0.01
tempo = linspace(0, tfinal, int(tfinal/Ts)+1)
y, t = step(test_ss, tempo)
# plt.plot(t, y)
# plt.show()

#test lsim
X0 = [0, 0]
U = 6*np.ones_like(tempo)

yout, T, xout = lsim(test_ss, U, tempo, X0)

fig, axs = plt.subplots(1, 2)
axs[0].plot(T, yout)
axs[1].step(T, xout[:, 0], 'b', where='post', label='x1')
axs[1].step(T, xout[:, 1], 'r', where='post', label='x2')
axs[1].legend(loc='best', shadow=True, fontsize='x-large')
plt.show()
######################
