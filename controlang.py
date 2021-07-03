from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

###Recife---------------
height = 4
latitude =  -8.0403
humidity = 0.85

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

l = 0.072456 #metade do lado do quadrado que froma o drone - A no witeboard que usamos
Ixx = l**4/12
Iyy = l**4/12

# Sistema

# Constantes do problema

diameter = 11 #inch
radius = diameter*25.4/2000
Ts = 0.02  # Amostragem
w_max = 9380 # Vel. max. rpm
c_tu = 0.1059
A = np.pi*radius**2
A_drone = 0.021
m = 1
Ohm0 = np.pi/30*(0.4*w_max) # convertido para rad/s
ct_br = (8/(np.pi ** 3))*c_tu
cd = 0.8
dz0 = 0.1

# Parametros motor

km = 12.1*(10**(-3))
Cq = 0.8
B = 5.3*(10**(-5))
J = 1*(10**(-4))
Rmotor = 0.275
c_pu = 0.0407
cq_br = (8/(np.pi ** 4))*c_pu


# Velocidade inicial de cada motor 
Ohm1 = np.pi/30*(0.4*w_max) # convertido para rad/s
Ohm2 = np.pi/30*(0.4*w_max) # convertido para rad/s
Ohm3 = np.pi/30*(0.4*w_max) # convertido para rad/s
Ohm4 = np.pi/30*(0.4*w_max) # convertido para rad/s
Ohm5 = np.pi/30*(0.4*w_max) # convertido para rad/s
Ohm6 = np.pi/30*(0.4*w_max) # convertido para rad/s

# Parametros motor
Cq = 0.8
B = 5.3*(10**(-5))
Rmotor = 0.275
c_pu = 0.0407
cq_br = (8/(np.pi ** 4))*c_pu

# Constantes do ss

b = (1/2)*(ct_br*air_density*A*radius**2)/m
d_corpo = (1/2)*(cd*air_density*A_drone)
d_motor = (1/2)*(cq_br*air_density*A*radius**3)

Kfa = 1 # ajustar essa constante 
Ohmr = 1# ajustar essa constante 

# Constantes motor

Am1 = -(km**2/(Rmotor*J) + (2*d_motor*Ohm1)/J)
Am2 = -(km**2/(Rmotor*J) + (2*d_motor*Ohm2)/J)
Am3 = -(km**2/(Rmotor*J) + (2*d_motor*Ohm3)/J)
Am4 = -(km**2/(Rmotor*J) + (2*d_motor*Ohm4)/J)
Am5 = -(km**2/(Rmotor*J) + (2*d_motor*Ohm5)/J)
Am6 = -(km**2/(Rmotor*J) + (2*d_motor*Ohm6)/J)

Bm1 = Bm2 = Bm3 = Bm4 = Bm5 = Bm6 = km/(Rmotor*J)

# Constantes dos estados
const_x1 = (d_corpo*dz0**2)/m
const_x2 = (3*b*Ohm0**2)/m +g
const_x3 = ((d_motor*Ohm0**2)/J)

vet_const = [const_x1, const_x2, const_x1]


# Espaço de estados - Continuo

aux1 = l*b/Ixx
aux2 = l*b*np.sqrt(3)/Ixx

Ac_ang = [[0,1,0,0,0,0,0,0,0,0],[0,(-2*Kfa/Ixx),0,(-J*Ohmr/Ixx),(-aux1*Ohm1),(-2*aux1*Ohm2),(-aux1*Ohm3),(aux1*Ohm4),(2*aux1*Ohm5),(aux1*Ohm6)],[0,0,0,1,0,0,0,0,0,0],[0,(J*Ohmr/Iyy),0,(-2*Kfa/Iyy),(-aux2*Ohm1),0,(aux2*Ohm3),(aux2*Ohm4),0,(-aux2*Ohm6)],[0,0,0,0,Am1,0,0,0,0,0],[0,0,0,0,0,Am2,0,0,0,0],[0,0,0,0,0,0,Am3,0,0,0],[0,0,0,0,0,0,0,Am4,0,0],[0,0,0,0,0,0,0,0,Am5,0],[0,0,0,0,0,0,0,0,0,Am6]]
Bc_ang = [[0],[0],[0],[0],[Bm1],[Bm2],[Bm3],[Bm4],[Bm5],[Bm6]]

Cc_ang = [[1,0,0,0,0,0,0,0,0,0],[0,0,1,0,0,0,0,0,0,0]]
Dc_ang = [[0],[0]]

ang_sis = ss(Ac_ang, Bc_ang, Cc_ang, Dc_ang)
disc_ang_sis = c2d(ang_sis,Ts)
'''
yout_ang,T = step(disc_ang_sis)
plt.step(T, yout_ang, where='post')
plt.show()
'''


#mT_c = np.eye(10)
X0_ang = [[1],[0],[1],[0],[0],[0],[0],[0],[0],[0]]
#init_state_c = np.matmul(mT_c, X0_ang)
final_time_c = 10
time_d_c = linspace(0, int(final_time_c/Ts)*Ts, int(final_time_c/Ts)+1)

yout_cc, Tcc, xout_cc = initial(disc_ang_sis, time_d_c, X0_ang, return_x=True)
#xx = np.matmul(np.linalg.inv(mT_c), xout_cc.T)
plt.step(Tcc, yout_cc, where='post')
#plt.step(Tcc, xx[0,:], where='post')
plt.show()

#controle
#Mp = 0












