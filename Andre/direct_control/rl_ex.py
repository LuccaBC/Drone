################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 06.05.2021
#
#   Info: Direct Root Locus Example
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

t_cont = 0.001
Ts = [0.64, 0.2, 0.067]
t_final = 6
vt_cont = linspace(0, t_final, int(t_final/t_cont)+1)
fig, axs_plot = plt.subplots(1, 2)

G_cont = tf([1], [1, 1, 0])
T_MA_cont = feedback(G_cont, 1)
U_MA_cont = 1/(G_cont+1)
yc, tc = step(T_MA_cont, vt_cont)
uc, tc = step(U_MA_cont, vt_cont)
axs_plot[0].plot(tc, yc, color=[0.5, 0.5, 0.5, 0.5], label='cont_K=1', linewidth=4)
axs_plot[1].plot(tc, uc, color=[0.5, 0.5, 0.5, 0.5], label='cont_K=1', linewidth=4)

zeta = 0.5
wd = np.pi
wn = wd/np.sqrt(1-zeta**2)
des_pole = complex(-zeta*wn, wd)
color_format = [[1, 0, 0, 0.5], [0, 1, 0, 0.5], [0, 0, 1, 0.5]]
for i, tdisc in enumerate(Ts):
    Gd = c2d(G_cont, tdisc, method='zoh')
    print(Gd)
    zz, pp, kk = tf2zpk(Gd.num[0][0], Gd.den[0][0])
    print("z, p, k:")
    print(zz, pp, kk)
    #Evaluating at gain K=1, closed loop
    Td = feedback(Gd, 1)
    # vt_d = linspace(0, int(t_final/tdisc)*tdisc, int(t_final / tdisc) + 1)
    # ydd, tdd = step(Td, vt_d)
    # plt.step(tdd, ydd, where='post')
    # axs_plot.step(tdd, ydd, color=color_format[i], where='post', label=tdisc, linewidth=4)
    # plt.show()
    print("desired pole:")
    dd_pole = np.exp(des_pole*tdisc)
    print(dd_pole)
    # fig, axs = plt.subplots()
    # rlocus(Gd)
    # plt.plot(dd_pole.real, dd_pole.imag, 'rs')
    # plt.show(block=False)
    #Choosing arbitralery zc = closer to the smaller positive real root
    zc = np.min(pp)-0.1
    phase = np.arctan2(dd_pole.imag, dd_pole.real-zc)
    gain_c = 1/np.abs(dd_pole-zc)
    for zeros in zz:
        phase = phase + np.arctan2(dd_pole.imag-zeros.imag, dd_pole.real-zeros.real)
        gain_c = gain_c/np.abs(dd_pole-zeros)
    for poles in pp:
        phase = phase - np.arctan2(dd_pole.imag-poles.imag, dd_pole.real-poles.real)
        gain_c = gain_c*np.abs(dd_pole-poles)
    if phase < 0:
        phase += np.pi
    else:
        phase -= np.pi
    pc = dd_pole.real - dd_pole.imag/(np.tan(phase))
    gain_c = gain_c*np.abs(dd_pole - pc)
    print(zc, pc, gain_c)
    Cd = tf([gain_c, -gain_c*zc], [1, -pc], tdisc)
    print(Cd)
    # fig, axs = plt.subplots()
    # rlocus(Cd*Gd)
    # plt.plot(dd_pole.real, dd_pole.imag, 'rs')
    # plt.show(block=False)
    Hs = feedback(Cd*0.6*Gd/kk, 1)
    Ud = Cd*0.6/(kk*(1+Cd*Gd))
    vt_d = linspace(0, int(t_final/tdisc)*tdisc, int(t_final / tdisc) + 1)
    ydd, tdd = step(Hs, vt_d)
    udd, tdd = step(Ud, vt_d)
    axs_plot[0].step(tdd, ydd, color=color_format[i], where='post', label=tdisc, linewidth=4)
    axs_plot[1].step(tdd, udd, color=color_format[i], where='post', label=tdisc, linewidth=4)

axs_plot[0].legend(loc='best', shadow=True, fontsize='x-large')
# axs_plot[1].legend(loc='best', shadow=True, fontsize='x-large')
plt.show()
