################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 22.01.2021
#
#   Info: What this function should do!
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy

dT = 0.1
T = numpy.linspace(0, 4, 41)

sys = tf([1], [1, 1])
sysd = c2d(sys, dT, method='zoh')
sysd_i = c2d(sys, dT, method='impulse')

y_imp, T = impulse(sys, T, 0)
y_imp_d, T = impulse(sysd_i/dT, T, 0)

yout, T = step(sys, T, 0)
yd, T = step(sysd, T, 0)

plt.figure(1)
plt.plot(T.T, yout.T, 'b:', label="continuous", linewidth=5.0)
plt.step(T.T, yd.T, 'r-', label="discrete", where='post', linewidth=5.0)
plt.show(block=True)

plt.figure(2)
plt.plot(T.T, y_imp.T, 'b:', label="continuous", linewidth=5.0)
plt.step(T.T, y_imp_d.T, 'r-', label="discrete", where='post', linewidth=5.0)
plt.show(block=True)

# rlocus(sys)
# plt.show()
# rlocus(sysd)
# plt.show()