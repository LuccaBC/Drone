################################################################################################################
#
#               DIGITAL CONTROL - EE/UFSCAR
#
#   Author: André Carmona Hernandes
#   Version: 1
#   Last-Update: 06.05.2021
#
#   Info: On hover only, Momentum theory, check motor-blade relationship:
#
#   These codes are used in DIGITAL CONTROL classes. You may use and study by them, however use with caution!
#
################################################################################################################

from control.matlab import *
import matplotlib.pyplot as plt
import numpy as np
#Equations used:
#In Brazil, it is common to use 0.5*pho*Ct*(pi*R^2)*(omega*R)^2, however,
#https://m-selig.ae.illinois.edu/props/propDB.html uses a different convention, thus, to use the measured values,
# a scaling must be done

## Inputs
# Inputs of the blade
#TODO: Read txt from database and que the constants.
c_tu = 0.1059
c_pu = 0.0407
diameter = 11
#https://m-selig.ae.illinois.edu/props/volume-1/plots/apcsf_11x4.7_static_ct.png
#https://m-selig.ae.illinois.edu/props/volume-1/plots/apcsf_11x4.7_static_cp.png
#drone
mass = 1.0
des_control = 0.4 #35% for hover, it can be 40, but as there are some approximations, I am keeping at 35.
v_bat = 12 # 3S
n_motors = 6
#Motor
#https://br.banggood.com/IFlight-XING-E-Pro-2207-1800KV-3-6S-or-2450KV-2750KV-2-4S-Brushless-Motor-for-RC-Drone-FPV-Racing-p-1518196.html?cur_warehouse=CN&ID=47980&rmmds=search
KV = 789 #rpm/V
I_max = 43.6 # max current
#----------------------------
#Enviroment
#input
#São Carlos -----------
height = 856 #meters above seal level
latitude = -22.0104
humidity = 0.39

###São José do Rio Preto --------------
##height = 489
##latitude =  -20.8202
##humidity = 0.3
##
###Ribeirão Preto ---------------
##height = 546.8
##latitude =  -21.1036
##humidity = 0.36
##
###Recife---------------
##height = 4
##latitude =  -8.0403
##humidity = 0.85


#blade
radius = diameter*25.4/2000
ct_br = (8/(np.pi ** 3))*c_tu
cq_br = (8/(np.pi ** 4))*c_pu

#Atmospheric calculation

#constants
vapor_density = humidity*12.83/1000 #kg/m^3
M_dry = 0.0289652
M_water = 0.018016
R_gas = 8.31446
T0 = 288.15
L = 0.0065
p0 = 101325

#equations
temp = T0 - L*height
lat_r = np.deg2rad(latitude)
g = 9.780327*(1+0.0053024*(np.sin(lat_r))**2 - 0.0000058*(np.sin(2*lat_r))**2)
#Buck equation for Vapour pressure
p_vapour = 0.61121*np.exp((18.678 - (temp-273.15)/234.5)*((temp-273.15)/(257.14+(temp-273.15))))
#dry pressure
p_dry = p0*(1 - L*height/T0)**(9.81*M_dry/(R_gas*L))

rho_dry = p_dry*M_dry/(R_gas*temp)
hum_mass = vapor_density/rho_dry
partial_dry = (1 - hum_mass) * p_dry
partial_wet = hum_mass * p_vapour * 100000
air_density = (partial_dry*M_dry+partial_wet*M_water)/(R_gas*temp)
print("Densidade do ar: --------")
print(air_density)

#motor
kw = 30/(np.pi*KV)
kt = 1.5*kw/np.sqrt(3)
r_est = v_bat/I_max
#T_hover = P/n_motors

omega = np.sqrt(2*mass*g/(air_density*n_motors*ct_br*np.pi*(radius**4)))
rpm=omega*30/np.pi
print(omega*30/np.pi)
#ideal torque
Q = 0.5*air_density*cq_br*np.pi*(radius**5)*(omega**2)
Aq = r_est*Q/(kt*v_bat)
#motor BEMF
Q_bemf = kw*omega
Bq = Q_bemf/v_bat

#demanded
Q_used = Aq+Bq

#Desired Q

print((des_control-Q_used)*100)


current = (Q_used*v_bat-kw*omega)/r_est
print(current)

#Calculando contante de trust "b"
const_trust = 2*ct_br*air_density*radius*radius*np.pi*radius*radius*6
print(const_trust)
