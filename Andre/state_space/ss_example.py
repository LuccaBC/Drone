################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 11.05.2021
#
#   Info: State space model examples
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np

G = tf([1, 2], [1, 6.5, 8, 2.5])
print(G)
sys_ss = tf2ss(G)
test_ss, T_transform = ctrl.reachable_form(sys_ss)
Ac = np.round(test_ss.A, decimals=2)
Bc = np.round(test_ss.B, decimals=2)
Cc = np.round(test_ss.C, decimals=2)
Dc = np.round(test_ss.D, decimals=2)

print(Ac, Bc, Cc, Dc)

######################

G2 = tf(1, [1, 2, 0])
ss_2 = tf2ss(G2)
ss_2, T_t2 = ctrl.reachable_form(ss_2)
Ac = np.round(ss_2.A, decimals=2)
Bc = np.round(ss_2.B, decimals=2)
Cc = np.round(ss_2.C, decimals=2)
Dc = np.round(ss_2.D, decimals=2)
print(ss_2)
Ts = 0.2

#retangular para frente
Ad = (np.eye(len(Ac))+Ts*Ac)
Bd = Ts*Bc

sys_forward = ss(Ad, Bd, Cc, Dc, Ts)
print(sys_forward)

#retangular para trás
inv_ad = np.linalg.inv(Ad)
Adb = np.round(inv_ad, decimals=2)
Bdb = np.round(np.matmul(inv_ad, Bd), decimals=2)
Cdb = np.round(np.matmul(Cc, inv_ad), decimals=2)
Ddb = np.round(Dc + np.matmul(Cdb, Bd), decimals=2)
sys_backward = ss(Adb, Bdb, Cdb, Ddb, Ts)
print(sys_backward)

#tustin
inv_at = np.linalg.inv(np.eye(len(Ac))-0.5*Ts*Ac)
Adt = np.round(np.matmul(np.eye(len(Ac))+0.5*Ts*Ac, inv_at), decimals=2)
Bdt = np.round(np.matmul(inv_at, Bc*np.sqrt(Ts)), decimals=2)
Cdt = np.round(np.matmul(np.sqrt(Ts)*Cc, inv_at), decimals=2)
Ddt = np.round(Dc + 0.5*Ts*np.matmul(np.matmul(Cc, inv_at), Bc), decimals=2)
sys_tustin = ss(Adt, Bdt, Cdt, Ddt, Ts)
print(sys_tustin)

####################################

sys_ds = c2d(ss_2, Ts)
print(sys_ds)
tf_from_ss = ss2tf(sys_ds)
print(tf_from_ss)
print("compare zoh")
print(c2d(G2, Ts))
