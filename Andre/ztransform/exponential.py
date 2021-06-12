################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 07.03.2021
#
#   Info: Function to show advance
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np


Ts = 1
time_final = 4
n_points = round(time_final / Ts) + 1

time = linspace(0, time_final, n_points)

exp1 = 2 ** time
exp2 = 4*(2 ** time)

fig, axs = plt.subplots()
axs.stem(time, exp1, 'rs')
axs.stem(time+0.1, exp2, 'bo')
plt.show()

