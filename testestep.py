from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

timestep = 0.001
tf = 5
a = 3
b = 0.5
c = 5
def func_altura(altura, vel_z, ref):
    dz = 1.98125367*altura -0.98176825*vel_z + 0.98175586*ref
    return dz


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
#nossas variaveis
altura = 0
ref = 1


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
    dw = motor_func(omega,volt)
    vec_dw.append(dw)
    #first attempt
    gain = timestep / 12
    aux = 5*dw + 8*dw_k1 - dw_k2
    omega = omega + gain*aux
    vec_w.append(omega)

plt.plot(vec_time, vec_w)
plt.show()



