################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 11.05.2021
#
#   Info: Linearized pendulum example
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np

Ts = [0.07, 0.028]
final_time = 1
l = 0.5
g = 9.81
m = 4
wo2 = g/l
#states x1=\theta, x2 = \theta_dot
theta0 = np.deg2rad(10)
Ac = [[0, 1], [-wo2, 0]]
Bc = [[0], [1/(m*l**2)]]
Cc = [1, 0]
Dc = [0]
cont_ss = ss(Ac, Bc, Cc, Dc)
## Malha aberta
tempo = linspace(0, 4, 4001)
y_c, t_c, x_c = initial(cont_ss, tempo, [theta0, 0], return_x=True)
fig, axs_c = plt.subplots(2, 1)
axs_c[0].plot(t_c, y_c, 'r')
axs_c[1].plot(t_c, x_c[:, 1], 'b')
axs_c[0].set(xlabel='time', ylabel='theta [rad]')
axs_c[1].set(xlabel='time', ylabel='theta_dot [rad/s]')
plt.show()
zeta = 0.83
wd = 2*np.sqrt(wo2)
wn = wd/(np.sqrt(1-zeta**2))
p_c1 = np.complex(-zeta*wn, wd)
p_c2 = np.complex(-zeta*wn, -wd)
for i, _ts in enumerate(Ts):
    fig, axs = plt.subplots(1, 3)
    ss_d = c2d(cont_ss, _ts)
    p_d1 = np.exp(p_c1*_ts)
    p_d2 = np.exp(p_c2*_ts)
    p = pole(ss_d)
    X0 = [[theta0], [0]]
    # print(_ts)
    # print(np.round(ss_d.A, 4))
    # print(np.round(ss_d.B, 4))
    # print(p)

    for j in range(3):
        if j == 0:
            new_ss = ss_d
            init_state = X0
            mT = np.eye(2)
            color = 'r'
            label_t = 'modal'
            offset = -0.1
        elif j == 1:
            new_ss, mT = ctrl.reachable_form(ss_d)
            label_t = 'reachable'
            init_state = np.matmul(mT, X0)
            color = 'g'
            offset = 0
        elif j == 2:
            new_ss, mT = ctrl.observable_form(ss_d)
            label_t = 'observable'
            init_state = np.matmul(mT, X0)
            color = 'b'
            offset = 0.1

        # pzmap(new_ss, grid=True)
        # plt.show()
        K = place(new_ss.A, new_ss.B, [p_d1, p_d2])
        PHI = new_ss.A - np.matmul(new_ss.B, K)
        print('K=', np.round(K, 2))
        cl_ss = ss(PHI, [[0], [0]], new_ss.C, new_ss.D, _ts)
        print('A=', np.round(cl_ss.A, 2))
        # print(zero(cl_ss))
        # print(pole(cl_ss))
        time_d = linspace(0, int(final_time/_ts)*_ts, int(final_time/_ts)+1)
        yout, tout, xout = initial(cl_ss, time_d, init_state, return_x=True)
        xx = np.matmul(np.linalg.inv(mT), xout.T)
        axs[0].step(tout, np.rad2deg(xx[0, :].T)+offset, color, where='post', label=label_t)
        axs[1].step(tout, np.rad2deg(xx[1, :].T)+offset, color, where='post', label=label_t)
        axs[2].step(tout, K.A1[0]*xout[:, 0]+K.A1[1]*xout[:, 1]+offset, color, where='post', label=label_t)
        axs[0].set(xlabel='time', ylabel='theta [deg]')
        axs[1].set(xlabel='time', ylabel='theta_dot [deg/s]')
        axs[2].set(xlabel='time', ylabel='U')
        axs[0].set_title('theta')
        axs[1].set_title('theta_dot')
        axs[2].set_title('U')
        axs[0].legend(loc='best', shadow=True, fontsize='x-large')
        axs[1].legend(loc='best', shadow=True, fontsize='x-large')
        axs[2].legend(loc='best', shadow=True, fontsize='x-large')
    plt.show()

### observador
tt_s = 0.07
ss_d = c2d(cont_ss, tt_s)
p_o1 = np.exp(10*p_c1*tt_s)
p_o2 = np.exp(10*p_c2*tt_s)
X0 = [[theta0], [0]]
Ltr = place(ss_d.A.T, ss_d.C.T, [p_o1, p_o2])
L = Ltr.T
print(L)
# compare
newA = ss_d.A - np.matmul(L, ss_d.C)
ol_ss = ss(newA, [[0], [0]], ss_d.C, ss_d.D, tt_s)
time_o = linspace(0, int(final_time / tt_s) * tt_s, int(final_time / tt_s) + 1)
yout, tout, xout = initial(ol_ss, time_o, init_state, return_x=True)
fig, axs = plt.subplots(2, 1)
axs[0].step(tout, np.rad2deg(xout[:, 0].T), 'r', where='post', label='x1_error')
axs[1].step(tout, np.rad2deg(xout[:, 1].T), 'b', where='post', label='x2_error')
axs[0].set(xlabel='time', ylabel='theta [deg]')
axs[1].set(xlabel='time', ylabel='theta_dot [deg/s]')
axs[0].set_title('theta')
axs[1].set_title('theta_dot')
axs[0].legend(loc='best', shadow=True, fontsize='x-large')
axs[1].legend(loc='best', shadow=True, fontsize='x-large')
plt.show()

#putting all together
#To simulate all together, we shall expande the matrix representation
p_d1 = np.exp(p_c1 * tt_s)
p_d2 = np.exp(p_c2 * tt_s)
K = place(ss_d.A, ss_d.B, [p_d1, p_d2])
all_A = np.concatenate((ss_d.A.A, -np.matmul(ss_d.B, K)), axis=1)
aux_A = np.concatenate((np.matmul(L, ss_d.C), newA-np.matmul(ss_d.B, K)), axis=1)
all_A = np.concatenate((all_A, aux_A), axis=0)

print(all_A)
aug_ss = ss(all_A, [[0], [0], [0], [0]], [1, 0, 0, 0], [0], tt_s)
initial_aug = [[theta0], [0], [0], [0]]
yout, tout, xout = initial(aug_ss, time_o, initial_aug, return_x=True)
fig, axs = plt.subplots(1, 2)
axs[0].step(tout, np.rad2deg(xout[:, 0].T), 'r', where='post', label='x1')
axs[1].step(tout, np.rad2deg(xout[:, 1].T), 'r', where='post', label='x2')
axs[0].step(tout, np.rad2deg(xout[:, 2].T)+1, 'b', where='post', label='x1_bar+1')
axs[1].step(tout, np.rad2deg(xout[:, 3].T)+1, 'b', where='post', label='x2_bar+1')
axs[0].set(xlabel='time', ylabel='theta [deg]')
axs[1].set(xlabel='time', ylabel='theta_dot [deg/s]')
axs[0].set_title('theta')
axs[1].set_title('theta_dot')
axs[0].legend(loc='best', shadow=True, fontsize='x-large')
axs[1].legend(loc='best', shadow=True, fontsize='x-large')
plt.show()

#compensator
d_A = newA-np.matmul(ss_d.B, K)
d_B = L
d_C = K
d_D = [0]
G_d = ss2tf(d_A, d_B, d_C, d_D, tt_s)
print(G_d)

#####
ref = np.deg2rad(2)
all_up = np.concatenate((ss_d.A.A, ss_d.B.A), axis=1)
all_dw = np.concatenate((ss_d.C.A, ss_d.D.A), axis=1)
all_2 = np.concatenate((all_up, all_dw), axis=0)
Nrf = np.matmul(np.linalg.inv(all_2), [[0], [0], [1]])
print(Nrf)

Nx = Nrf[0:2]
Nu = Nrf[2]
Nbar = Nu - np.matmul(K, Nx)
auxB = np.matmul(ss_d.B, Nbar)
newB = np.concatenate((auxB, auxB), axis=0)
aug_ss2 = ss(all_A, newB, [1, 0, 0, 0], [0], tt_s)
print(aug_ss2)


