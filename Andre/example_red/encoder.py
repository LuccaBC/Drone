################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 07.03.2021
#
#   Info: Understanting a encoder measurement
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np
from scipy import integrate


def estimate_vel(w_enc, t_s, window, Np):
    counter = 0
    k = 0
    w_aux = w_enc[0]
    w_e = [0]
    for idx, w in enumerate(w_enc):
        if w < w_aux: #borda de descida
            counter = counter + 1
        if abs(t_s[idx] - window * (k + 1)) < 0.000001:
            ang_desloc = counter / Np
            w_e.append(ang_desloc / window)
            counter = 0
            k = k + 1
        w_aux = w

    return w_e


Ts = 0.05  # 20 ms window
T_encoder = 0.000005  # for emulating encoder signal
t_a = 1  # time ending speed up
t_d = 6  # time starting speed down
a_f = [1, 20, 50, 100]  # acceleration tests
final_time = 9
n = round(final_time / Ts) + 1
time_s = linspace(0, final_time, n)
n_e = round(final_time / T_encoder) + 1
time_e = linspace(0, final_time, n_e)
Np = 12  # 3 pulses per revolution * 4 of quadrature

enc_dividers = np.linspace(0, 1, 2*Np + 1)

accel_signal = np.zeros(n_e)
a_e = np.zeros(n)
DEBUG = True
for a in a_f:
    f_encoder = np.zeros(n_e)
    accel_signal[time_e < t_a] = a
    accel_signal[time_e > t_d] = -(a*t_a)/(final_time-t_d)
    a_e[time_s < t_a] = a
    a_e[time_s > t_d] = -(a*t_a)/(final_time-t_d)
    # vel_signal = integrate.cumtrapz(a_e, time_s, Ts, initial=0)
    f_test = integrate.cumtrapz(accel_signal, time_e, T_encoder, initial=0)
    ang_test = (integrate.cumtrapz(f_test, time_e, T_encoder, initial=0)) % 1
    for i in np.arange(1, 2*Np, 2):
        upperBound = ang_test < enc_dividers[i+1]
        lowerBound = ang_test > enc_dividers[i]
        condition = np.logical_and(lowerBound, upperBound)
        f_encoder[np.where(condition)] = 1

    # if DEBUG:
    #     fig, axs = plt.subplots()
    #     plt.plot(time_e, f_encoder)
    #     plt.plot(time_e, ang_test)
    #     plt.show()
    w_est = estimate_vel(f_encoder, time_e, Ts, Np)
    # w_est.append(w_est[-1])
    if DEBUG:
        fig, axs = plt.subplots()
        plt.plot(time_e, f_test, label='real')
        plt.plot(time_s, w_est, 'r:', label='estimado')
        plt.legend()
        plt.show()
