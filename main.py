# Projeto de Controle Digital
# Parte 1 - movimento Z

from control.matlab import *
import control as ctrl
import matplotlib.pyplot as plt
import numpy as np
import scipy.signal as sig

#Enviroment input
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

# Constantes do ss

b = (1/2)*(ct_br*air_density*A*radius**2)/m
d_corpo = (1/2)*(cd*air_density*A_drone)
d_motor = (1/2)*(cq_br*air_density*A*radius**2)

# Constantes motor

Am = -(km**2/(Rmotor*J) + (2*d_motor*Ohm0)/J)
Bm = km/(Rmotor*J)

# Constantes dos estados
const_x1 = (d_corpo*dz0**2)/m
const_x2 = (3*b*Ohm0**2)/m +g
const_x3 = ((d_motor*Ohm0**2)/J)

vet_const = [const_x1, const_x2, const_x1]
print("Veotr de constantes:\n", vet_const)
# Espaço de estados - Continuo

Ac = [[0,1,0],[0,(2*d_corpo*dz0)/m,(6*b*Ohm0)/m],[0,0,Am]]
Bc = [[0],[0],[Bm]]
Cc = [1,0,0]
Dc = [0]

#constant_sis = -(6*ct_br*air_density*A*((w0*radius)**2)/m + g) 
#print("Constante:", constant_sis)

cont_ss = ss(Ac,Bc,Cc,Dc)
print("Sistema Continuo:\n", cont_ss)

# Discretização

disc_ss = c2d(cont_ss, Ts)
print("Sistema Discreto:\n", disc_ss)


print(disc_ss.pole())
# Teste

yout, T = step(cont_ss)
plt.step(T, yout, where='post')
plt.show()

# Vel
C_vel = [0,1,0]
cont_ss_vel = ss(Ac,Bc,C_vel,Dc)

yout_vel, T_vel = step(cont_ss_vel)
plt.step(T_vel, yout_vel, where='post')
plt.show()

# Controlador
# Ms = 25%
# Ts de 1% = 10s
# Acc max = 10,2 cm/s^2

zeta = 0.4037
ts_1 = 10

wn = 4.6/(zeta*ts_1)
wd = wn*np.sqrt(1-zeta**2)

# polos continuos
p_c1 = complex(-zeta*wn,wd)
p_c2 = complex(-zeta*wn,-wd)
print("polos Continuos:\n", np.round(p_c1,4), "e" , np.round(p_c2,4))

# polos discretos
p_d1 = np.exp(p_c1*Ts)
p_d2 = np.exp(p_c2*Ts)
print("polos Discretos:\n", np.round(p_d1,4), "e", np.round(p_d2,4))

# Valor inicial - gravidade medida pela IMU ?
X0 = [[-1], [0], [0]]

#termo de "ruido" do resto da expansão de taylor - ficou bem ruim

mT = np.eye(3)
new_ss, mT = ctrl.reachable_form(disc_ss)
init_state = np.matmul(mT, X0)
# colocando os polos
K = place(new_ss.A, new_ss.B, [p_d1, p_d2, 0.6309566]) # Controlabilidade forma canonica
PHI = new_ss.A - np.matmul(new_ss.B,K)
print('K=', np.round(K, 2))

#Parte nova
final_time = 10
offset = 0
print('Matriz nova controle \n', PHI)
cl_ss = ss(PHI, [[0], [0], [0]], new_ss.C, new_ss.D, Ts)
print('A=', np.round(cl_ss.A, 2))
time_d = linspace(0, int(final_time/Ts)*Ts, int(final_time/Ts)+1)
yout, tout, xout = initial(cl_ss, time_d, init_state, return_x=True)
xx = np.matmul(np.linalg.inv(mT), xout.T)
#print('xx=\n', xx)
plt.step(tout, (xx[0, :].T)+offset, 'r', where='post', label='modal')
plt.step(tout, (xx[1, :].T)+offset, 'b', where='post', label='modal')
plt.show()

#Acrescentando ref
all_up = np.concatenate((new_ss.A.A, new_ss.B.A), axis=1)
all_dw = np.concatenate((new_ss.C.A, new_ss.D.A), axis=1)
all_2 = np.concatenate((all_up, all_dw), axis=0)

Nrf = np.matmul(np.linalg.inv(all_2), [[0], [0], [0], [1]])
print("bagulho:\n",np.linalg.inv(all_2))
print('Nrf\n',Nrf,'\nall_2\n',all_2)

Nx = Nrf[0:3]
Nu = -636.38#Nrf[3]
N =  Nu + np.matmul(K, Nx)
print("N:\n",N)
auxB = np.matmul(disc_ss.B, N)
newB = np.concatenate((auxB, auxB), axis=0)
print('matriz A\n', cl_ss.A)
print('matriz auxB\n', auxB)
aug_ss2 = ss(cl_ss.A, auxB, cl_ss.C, [0], Ts)
print('printando a bagaça final\n',aug_ss2)

X0_c = [[0], [0], [0]]
init_state_c = np.matmul(mT, X0_c)
final_time_c = 10
time_d_c = linspace(0, int(final_time_c/Ts)*Ts, int(final_time_c/Ts)+1)

yout_c, tout_c, xout_c = initial(aug_ss2, time_d_c, init_state_c, return_x=True)
xx_c = np.matmul(np.linalg.inv(mT), xout_c.T)
yout_cc, Tcc = step(aug_ss2, time_d_c, init_state_c)
plt.step(Tcc, yout_cc, where='post')
plt.show()


