################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 06.05.2021
#
#   Info: Frequency controller  Example
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

t_cont = 0.001
Ts = 0.2
t_final = 6
vt_cont = linspace(0, t_final, int(t_final/t_cont)+1)
fig, axs_plot = plt.subplots(1, 2)

G_cont = tf([1], [1, 1, 0])
Gd = c2d(G_cont, Ts, method='zoh')
print(Gd)

AA = Gd.num[0][0][0]
BB = Gd.num[0][0][1]
CC = Gd.den[0][0][1]
DD = Gd.den[0][0][2]

num_w = [(BB-AA)*(Ts**2), -4*BB*Ts, 4*(AA+BB)]
den_w = [(1-CC+DD)*(Ts**2), (4-4*DD)*Ts, 4*(1+CC+DD)]

G_w = tf(num_w, den_w)

print(G_w.num[0][0], G_w.den[0][0])

G_w_simp = tf([G_w.num[0][0][1], G_w.num[0][0][2]], [G_w.den[0][0][0], G_w.den[0][0][1], 0])
print(G_w_simp)

ess = 0.3
Kv = 1/ess
bbox = dict(boxstyle="round", fc="0.8")
arrowprops = dict(
    arrowstyle="->",
    connectionstyle="angle,angleA=0,angleB=-90,rad=10")
omega = logspace(start=-1, stop=2, num=10000, base=10)
mag, phase, omega = bode(G_w_simp*Kv, omega)

GM, PM, wg, wp = margin(G_w_simp*Kv)

fig_bode = plt.gcf()
fig_bode.axes[0].plot(wg, -20*np.log10(GM), 'rs')
fig_bode.axes[0].annotate(
    f'data = ({wg:.2f}, {-20*np.log10(GM):.2f})',
    (wg, -20*np.log10(GM)),
    xytext=(20, 30), textcoords='offset points',
    bbox=bbox, arrowprops=arrowprops)

fig_bode.axes[1].plot(wp, -180+PM, 'ro')
fig_bode.axes[1].annotate(
    f'data = ({wp:.2f}, {PM:.2f})',
    (wp, -180+PM),
    xytext=(20, -30), textcoords='offset points',
    bbox=bbox, arrowprops=arrowprops)
####
PM_des = 50
phi_max = PM_des - (PM-np.rad2deg(0.5*wp*Ts))

alfa = (1-np.sin(np.deg2rad(phi_max)))/(1+np.sin(np.deg2rad(phi_max)))

new_mag = -10*np.log10(1/alfa)
print(new_mag)

idx = np.where(abs(20*np.log10(mag)-new_mag) < 0.005)
new_wg = omega[idx]
# offset = 10
plt.show()
tau = 1/(new_wg*np.sqrt(alfa))

Cw = tf([Kv*tau[0], Kv], [alfa*tau[0], 1])
print(Cw)

bode(Cw*G_w_simp)
GM, PM, wg, wp = margin(G_w_simp*Cw)

fig_bode = plt.gcf()
fig_bode.axes[0].plot(wg, -20*np.log10(GM), 'rs')
fig_bode.axes[0].annotate(
    f'data = ({wg:.2f}, {-20*np.log10(GM):.2f})',
    (wg, -20*np.log10(GM)),
    xytext=(20, 30), textcoords='offset points',
    bbox=bbox, arrowprops=arrowprops)

fig_bode.axes[1].plot(wp, -180+PM, 'ro')
fig_bode.axes[1].annotate(
    f'data = ({wp:.2f}, {PM:.2f})',
    (wp, -180+PM),
    xytext=(20, -30), textcoords='offset points',
    bbox=bbox, arrowprops=arrowprops)
plt.show()

Cz = tf([10.66, -8.311], [1, -0.2943], Ts)
bode(Cz*Gd)
plt.show()

R = tf(1, [1, 0])
Rz = c2d(R, Ts)

Tz_mf = feedback(Cz*Gd, 1)

vt_d = linspace(0, int(t_final / Ts) * Ts, int(t_final / Ts) + 1)
ydd, tdd = step(Tz_mf*Rz, vt_d)
y_2, tdd = step(Rz, vt_d)
plt.plot(tdd, y_2, color='b', label=Ts, linewidth=4)
plt.step(tdd, ydd, color='r', where='post', label=Ts, linewidth=4)
plt.show()


