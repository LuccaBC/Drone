################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 27.05.2021
#
#   Info: LQR_example based on Digital control of dynamic systems
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np

A = [[-0.2, 0.1, 1], [-0.05, 0, 0], [0, 0, -1]]
B = [[0, 1], [0, 0.7], [1, 0]]
C = [[1, 0, 0], [0, 1, 0]]
D = [[0, 0], [0, 0]]

cont_ss = ss(A, B, C, D)
print(cont_ss)
Ts = 0.2
disc_ss = c2d(cont_ss, Ts)
rounded_ss = ss(np.round(disc_ss.A, 5), np.round(disc_ss.B, 5), disc_ss.C, disc_ss.D, Ts)
print(rounded_ss)
#############
ts1 = 2.4
alfa = np.round(100**(Ts/ts1), 2)
print(alfa)
newPhi = np.round(alfa*rounded_ss.A, 4)
print(newPhi)
newGamma = np.round(alfa*rounded_ss.B, 4)
print(newGamma)
Q = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]
R = [[2, 0], [0, 0.1]]


X, L, G = dare(newPhi, newGamma, Q, R)

#X - Riccati solution
#L - closed loop eigen values
#G - gain matrix

print('Solution')
print('K=', G)
print('autovalores=', L)
print('riccati sol=', X)

#
X0 = [0, 1, 0]
final_time = 5
PHI = disc_ss.A - np.matmul(disc_ss.B, G)
cl_ss = ss(PHI, [[0, 0], [0, 0], [0, 0]], disc_ss.C, disc_ss.D, Ts)
print('A=', np.round(cl_ss.A, 2))
time_d = linspace(0, int(final_time/Ts)*Ts, int(final_time/Ts)+1)
yout, tout, xout = initial(cl_ss, time_d, np.matrix(X0).T, return_x=True)
ctl_action = -G*xout.T
fig, axs = plt.subplots(1, 3)
axs[0].step(tout, xout[:, 0], 'r', where='post', label='x1')
axs[1].step(tout, xout[:, 1], 'b', where='post', label='x2')
axs[2].step(tout, ctl_action[0, :].T, 'r', where='post', label='Uc')
axs[2].step(tout, ctl_action[1, :].T, 'b', where='post', label='Us')
axs[2].legend(loc='best', shadow=True, fontsize='x-large')
axs[0].set_title('x1')
axs[1].set_title('x2')
axs[2].set_title('U')
plt.show()