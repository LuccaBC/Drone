################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 06.05.2021
#
#   Info: non-linear solving
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np

##
# simulated equation
# wdot: -a*w-b*w^2+c*u
# a=3, b= 0.5, c=5

timestep = 0.001
tf = 5
a = 3
b = 0.5
c = 5


def motor_func(omega, voltage):
    w_dot = -a*omega-b*omega**2+c*voltage
    return w_dot

## main
vec_time = np.linspace(0, tf, int(tf/timestep)+1)
omega = 0
dw_k1 = 0
dw_k2 = 0
per_dc = 0
vbat = 11.1
vec_w = [omega]
vec_dw = [0]
omega_ref = 5
sumErr = 0
for count_t, timestamp in enumerate(vec_time[1:]):
    #control part
    err = omega_ref-omega
    sumErr = sumErr + err
    control_action = 1/30*err + 1/1000*sumErr
    if control_action < 0:
        per_dc = 0
    elif control_action > 1:
        per_dc = 1
    else:
        per_dc = control_action
    volt = per_dc*vbat
    #simulation
    dw = motor_func(omega, volt)
    vec_dw.append(dw)
    # first attempt
    gain = timestep / 12
    aux = 5 * dw + 8 * dw_k1 - dw_k2
    omega = omega + gain * aux
    vec_w.append(omega)

plt.plot(vec_time, vec_w)
plt.show()

################ rascunho
'''
começo

Observador:
xbarra = (PHI - GAMMA*K-L*H)*Xbarra
(H*x=y)
#Sensores
## Função sensor retorna y => (y = H*x(xreal) + ruido de medição ==  +randn(1)*sigma)
xbarra = xbarra + L*y
Calcular a acao controle
u = -K*xbarra+Nbarra*ref
### daqui ela iria para o motor
#aqui começa a simulação não linear
# wdot: -a*w-b*w^2+c*u -> motor;
# T: bw^2
zdotdot: 6T - drag(v^2) -mg
# retangular para trás
zdor += zdotdot*Ts
z += zdot*Ts + 0.5*zdotdot*Ts^2

'''
